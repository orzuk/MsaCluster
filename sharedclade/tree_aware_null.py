#!/usr/bin/env python3
"""Tree-aware permutation null for CROSS-METHOD clade agreement.

Why this exists
---------------
``sharedclade/analyse_shared_clade.py`` measures, per fold-switch family,
how many of the 7 fold-preference methods localise their clade-shift signal
to the SAME branch (pairwise Jaccard of ``clade_best_members`` >= 0.5, same
direction; largest such clique = the agreement statistic).  Its null draws
random SUBSETS of the family's cluster set, which is badly anti-conservative
because random subsets are not clades: real trees only offer ~n-2 nested
clades, and nested clades of similar size overlap heavily, so agreement is
far more likely by chance than a subset null admits.

This script computes the proper null: permute each method's per-cluster
values across the leaves of the REAL tree (independently per method), rerun
the clade-shift scan, and recompute the agreement statistic.  It needs the
cluster tree, which lives on the HUJI cluster, so run it there:

    python3 sharedclade/tree_aware_null.py \
        --pairs 2n54B_2hdmA,5jytA_2qkeE,2namA_1uxmK \
        --labels centered \
        --survey-csv docs/fold_diversity_survey.csv \
        --n-perm 10000 --out clade_agreement_null_centered.csv

    # corrected labels
    python3 sharedclade/tree_aware_null.py --pairs ALL --labels corrected \
        --survey-csv docs/fold_diversity_survey_corrected_9method.csv \
        --n-perm 10000 --out clade_agreement_null_corrected.csv

Run it from the repository root (it imports ``utils.phylo_tests`` and
``utils.ancestral_utils``), with the same ``--coarse-rep`` as the scan whose
results you are testing (default ``largest_piece_medoid``, i.e. the
medoid-representative run).

Two nulls are reported per family, both preserving each method's value
distribution and each method's missing-value pattern:

  joint       one leaf permutation applied to ALL methods at once.  Keeps
              the cross-method correlation within a cluster (AF2 and AF3
              agree because they see the same cluster) and destroys only
              the association with the tree.  This is the null that matches
              the question "do the methods agree on a branch beyond what
              their mutual correlation plus tree shape already give?" and is
              the one to quote.
  independent one permutation per method.  Destroys cross-method
              correlation as well, so it is the more liberal of the two;
              report it as a sensitivity check.

Output columns (one row per pair x null):
  pair_id, label, null_type, n_methods, n_perm,
  obs_k, null_mean_k, null_sd_k, p_k              (direction-aware clique)
  obs_k_nodir, null_mean_k_nodir, p_k_nodir       (same branch, any direction)
  obs_edges, null_mean_edges, p_edges             (# agreeing method pairs)
"""
from __future__ import annotations

import argparse
import itertools
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
for p in (REPO, os.path.join(REPO, "..")):
    if p not in sys.path:
        sys.path.insert(0, p)

from utils.phylo_tests import (average_ranks, clade_split_matrix,  # noqa: E402
                               mann_whitney_z)
# ete3-dependent imports are deferred into run_pair()/main() so the scan
# maths below can be unit-tested on a machine without ete3 installed.

METHODS7 = ["AF2", "AF3", "Boltz2", "BioEmu", "ESM", "DDG", "S4PRED"]
MAIN_THR = 0.5
CHUNK = 500


# --------------------------------------------------------------------------
def best_split_batch(M: np.ndarray, x: np.ndarray, order: np.ndarray,
                     min_clade: int = 2):
    """Best split index and its sign for a batch of leaf arrangements.

    M      (S, n) bool clade-membership matrix (from clade_split_matrix)
    x      (n,)   the method's values, aligned to M's columns; NaN = missing
    order  (B, n) leaf-index permutations; row b means "leaf j gets the
                  value x[order[b, j]]"
    Returns (best (B,) int split index, sign (B,) float, ok (B,) bool).
    Mirrors utils.phylo_tests._scan_core / run_clade_shift for one method:
    missing leaves are dropped from the ranking and from the clade sizes,
    and T = max |z| over the valid splits.
    """
    B = order.shape[0]
    finite = np.isfinite(x)
    r = np.zeros_like(x, dtype=float)
    tt = 0.0
    if finite.sum() >= 2:
        rr, tt = average_ranks(x[finite])
        r[finite] = rr
    n_k = float(finite.sum())
    Mt = M.T.astype(float)                       # (n, S)
    R = r[order]                                 # (B, n)
    A = finite[order].astype(float)              # (B, n)
    R1 = R @ Mt                                  # (B, S)
    N1 = A @ Mt                                  # (B, S)
    z = mann_whitney_z(R1, N1, n_k, tt)
    valid = (N1 >= min_clade) & (n_k - N1 >= min_clade) & np.isfinite(z)
    az = np.where(valid, np.abs(z), -np.inf)
    best = az.argmax(axis=1)
    ok = np.isfinite(az[np.arange(B), best]) & (az[np.arange(B), best] > -np.inf)
    sgn = np.sign(z[np.arange(B), best])
    return best, sgn, ok


def agreement_stats(M, best, sgn, ok, masks, thr=MAIN_THR):
    """(clique size with direction, clique size ignoring direction, #edges)
    for each arrangement.  ``masks`` (K, n) marks the leaves each method
    scored: a method's clade is its split intersected with its own leaves,
    exactly as ``clade_best_members`` is reported."""
    K, B = best.shape[0], best.shape[1]
    sets = np.zeros((K, B, M.shape[1]), dtype=bool)
    for k in range(K):
        sets[k] = M[best[k]] & masks[k][None, :]
    pl = list(itertools.combinations(range(K), 2))
    pidx = {p: i for i, p in enumerate(pl)}
    e_dir = np.zeros((B, len(pl)), dtype=bool)
    e_nod = np.zeros((B, len(pl)), dtype=bool)
    for pi, (i, j) in enumerate(pl):
        inter = (sets[i] & sets[j]).sum(1).astype(float)
        union = (sets[i] | sets[j]).sum(1).astype(float)
        with np.errstate(invalid="ignore", divide="ignore"):
            jac = np.where(union > 0, inter / union, 0.0)
        both = ok[i] & ok[j]
        e_nod[:, pi] = both & (jac >= thr - 1e-12)
        e_dir[:, pi] = e_nod[:, pi] & (sgn[i] == sgn[j]) & (sgn[i] != 0)
    out = []
    for E in (e_dir, e_nod):
        st = np.ones(B, dtype=np.int64)
        for size in range(K, 1, -1):
            for sub in itertools.combinations(range(K), size):
                need = [pidx[p] for p in itertools.combinations(sub, 2)]
                st = np.where(E[:, need].all(1) & (st < size), size, st)
        out.append(st)
    return out[0], out[1], e_dir.sum(1)


# --------------------------------------------------------------------------
def run_pair(pair_id, values, rep_strategy, n_perm, rng, min_clade=2):
    from utils.ancestral_utils import build_coarse_tree, _name_internal_nodes
    tree, tag_to_label, label_to_tag = build_coarse_tree(
        pair_id, rep_strategy=rep_strategy)
    _name_internal_nodes(tree)
    tags = [t for t in tag_to_label if t != "DeepMsa"]
    labels = [tag_to_label[t] for t in tags]
    n = len(labels)
    if n < 4:
        return None
    M, _ = clade_split_matrix(tree, labels, min_clade)
    if M.shape[0] == 0:
        return None

    xs, masks, methods = [], [], []
    for m in METHODS7:
        v = np.array([values.get((pair_id, m, t), None) for t in tags],
                     dtype=object)
        x = np.array([np.nan if u is None else float(u) for u in v], float)
        if np.isfinite(x).sum() < 2 * min_clade:
            continue
        xs.append(x)
        masks.append(np.isfinite(x))
        methods.append(m)
    K = len(methods)
    if K < 2:
        return None
    masks = np.vstack(masks)

    ident = np.arange(n)[None, :]
    best0 = np.zeros((K, 1), int)
    sgn0 = np.zeros((K, 1))
    ok0 = np.zeros((K, 1), bool)
    for k in range(K):
        b, s, o = best_split_batch(M, xs[k], ident, min_clade)
        best0[k], sgn0[k], ok0[k] = b, s, o
    obs_k, obs_knd, obs_e = agreement_stats(M, best0, sgn0, ok0, masks)

    res = {}
    for null_type in ("joint", "independent"):
        ks, knds, es = [], [], []
        done = 0
        while done < n_perm:
            b = min(CHUNK, n_perm - done)
            if null_type == "joint":
                order = rng.permuted(np.tile(np.arange(n), (b, 1)), axis=1)
                orders = [order] * K
            else:
                orders = [rng.permuted(np.tile(np.arange(n), (b, 1)), axis=1)
                          for _ in range(K)]
            bb = np.zeros((K, b), int)
            ss = np.zeros((K, b))
            oo = np.zeros((K, b), bool)
            for k in range(K):
                bb[k], ss[k], oo[k] = best_split_batch(M, xs[k], orders[k],
                                                       min_clade)
            a, c, e = agreement_stats(M, bb, ss, oo, masks)
            ks.append(a); knds.append(c); es.append(e)
            done += b
        ks = np.concatenate(ks); knds = np.concatenate(knds)
        es = np.concatenate(es)
        res[null_type] = dict(
            pair_id=pair_id, null_type=null_type, n_methods=K,
            methods=",".join(methods), n_leaves=n, n_splits=int(M.shape[0]),
            n_perm=n_perm,
            obs_k=int(obs_k[0]), null_mean_k=float(ks.mean()),
            null_sd_k=float(ks.std()),
            p_k=(1 + int((ks >= obs_k[0]).sum())) / (1 + n_perm),
            obs_k_nodir=int(obs_knd[0]), null_mean_k_nodir=float(knds.mean()),
            p_k_nodir=(1 + int((knds >= obs_knd[0]).sum())) / (1 + n_perm),
            obs_edges=int(obs_e[0]), null_mean_edges=float(es.mean()),
            p_edges=(1 + int((es >= obs_e[0]).sum())) / (1 + n_perm),
        )
        res[null_type]["_draws_k"] = ks.astype(np.int16)
        res[null_type]["_draws_knd"] = knds.astype(np.int16)
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs", default="ALL",
                    help="comma-separated pair ids, or ALL (every pair in the survey csv)")
    ap.add_argument("--labels", choices=["centered", "corrected"], default="centered")
    ap.add_argument("--survey-csv", required=True)
    ap.add_argument("--coarse-rep", choices=["first", "largest_piece_medoid"],
                    default="largest_piece_medoid")
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--min-clade", type=int, default=2)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    from scripts.phylo_placement import load_fold_value
    value_col = "TMdiff_residual" if args.labels == "corrected" else "TMdiff_centered"
    values = load_fold_value(args.survey_csv, value_col)
    if args.pairs == "ALL":
        pairs = sorted({k[0] for k in values})
    else:
        pairs = [p for p in args.pairs.split(",") if p]

    rng = np.random.default_rng(args.seed)
    rows = []
    draws = {}
    for i, pair in enumerate(pairs, 1):
        print(f"[{i}/{len(pairs)}] {pair}", flush=True)
        try:
            r = run_pair(pair, values, args.coarse_rep, args.n_perm,
                         rng, args.min_clade)
        except Exception as exc:                      # noqa: BLE001
            print(f"  failed: {exc}", flush=True)
            continue
        if r is None:
            print("  skipped (too few leaves / methods)", flush=True)
            continue
        for nt in ("joint", "independent"):
            row = dict(label=args.labels, **r[nt])
            draws[(pair, nt, "k")] = row.pop("_draws_k")
            draws[(pair, nt, "knd")] = row.pop("_draws_knd")
            rows.append(row)
            print(f"  {nt}: obs_k={row['obs_k']} null_mean={row['null_mean_k']:.2f} "
                  f"p={row['p_k']:.4f} | nodir obs={row['obs_k_nodir']} "
                  f"p={row['p_k_nodir']:.4f}", flush=True)

    import csv
    if rows:
        with open(args.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {len(rows)} rows -> {args.out}")
    if draws:
        # Per-permutation clique sizes, so the across-family sum has an
        # empirical null instead of a normal approximation.
        np.savez_compressed(os.path.splitext(args.out)[0] + "_draws.npz",
                            **{"|".join(k): v for k, v in draws.items()})
        print(f"wrote draws -> {os.path.splitext(args.out)[0]}_draws.npz")


if __name__ == "__main__":
    main()
