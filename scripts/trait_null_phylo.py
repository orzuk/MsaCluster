#!/usr/bin/env python3
"""Is fold preference MORE clade-structured than generic sequence traits?

Motivation
----------
The negative-trait baseline (scripts/negative_traits_phylo.py) showed that
generic per-cluster sequence traits (random projections of residue scores,
net charge, hydrophobicity) have significant Moran's I on the coarse cluster
tree about as often as fold preference does. A label-permutation p-value
therefore cannot say whether fold preference is clade-structured for reasons
beyond "every sequence-derived quantity is". This script replaces the
label-permutation null with an empirical null of K random SEQUENCE traits
evaluated on the same tree and the same leaves.

Per pair (coarse tree from ``utils.ancestral_utils.build_coarse_tree`` with
the default ``largest_piece_medoid`` representatives, then
``_name_internal_nodes``; pairs with < --min-clusters leaves are skipped):

1. K random sequence traits (default K = 200). Trait k of a cluster = mean over
   the cluster's member rows (a3m row 0 = query/medoid, skipped) of the mean
   of a fixed random per-residue score over the member's aligned residues
   (upper case = match columns; lower-case insertions and gaps excluded).
   Score table k = 20 amino acids ~ N(0, 1) from
   ``negative_traits_phylo.random_scores(k, --rp-seed)``, i.e.
   ``np.random.default_rng([rp_seed, k])`` -- the same tables as that
   script's ``random_proj_k`` traits (k = 1..3 coincide). Computed with the
   shared vectorised ``negative_traits_phylo.random_proj_values``.
2. For every random trait and every fold-preference method on the same leaf
   set (leaves with a finite fold value AND sequence traits):
     I_norm  Moran's I as a fraction of the achievable maximum, the SAME
             routine as ``utils.phylo_tests.morans_I_test``
             (``morans_I_norm``; tree quantities precomputed once);
     z       permutation-free z of Moran's I from the exact randomization
             moments (Cliff & Ord; ``phylo_tests.morans_analytic_z``). It
             estimates the same quantity as morans_I_test's permutation z
             without the 1000-permutation Monte-Carlo error, and is used for
             BOTH fold preference and random traits, so the comparison is
             self-consistent (values differ slightly from
             docs/phylo_placement*.csv morans_z);
     T       clade-shift statistic max |z| over splits
             (``phylo_tests.clade_T_stat``; no permutations).
3. Fold-preference methods: AF2, AF3, Boltz2, BioEmu, ESM, DDG, S4PRED
   (MSAT / CCMpred excluded: too few clusters). Values:
     raw        TMdiff_centered from docs/fold_diversity_survey.csv
     corrected  TMdiff_residual from docs/fold_diversity_survey_corrected_9method.csv
   (read with ``phylo_placement.load_fold_value``).
   p_trait(m) = (1 + #{random traits with statistic >= observed}) / (K + 1),
   separately for I_norm (primary), z and T.
4. Family level, per label set and statistic: p_fam = min(1, n_m * min_m
   p_trait(m)) (Bonferroni over the n_m methods with a value), then BH over
   families. ``best_method_*`` = argmin p_trait (ties -> larger statistic).
   No "best method vs best random trait" comparison (biased).

Outputs
-------
--out         docs/trait_null_phylo.csv          one row per (pair, labels, method)
--out-family  docs/trait_null_phylo_family.csv   one row per (pair, labels)
and a printed summary (families with q <= 0.05 per statistic and label set).

Usage (repo root)
-----------------
python3 scripts/trait_null_phylo.py --n-jobs 8
python3 scripts/trait_null_phylo.py --pairs 1fzpD_2frhA --n-random 50   # quick
python3 scripts/trait_null_phylo.py --summary-only                       # re-summarise --out
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-cache-" + os.environ.get("USER", "x"))
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)

import argparse
import multiprocessing as mp
import signal
import sys
import time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, "scripts"))

from utils import phylo_tests as pt  # noqa: E402
from utils.ancestral_tests import bh_adjust  # noqa: E402
from negative_traits_phylo import (random_score_matrix, random_proj_values,  # noqa: E402
                                   discover_pairs)
from seq_divergence_correction import load_a3m_seqs  # noqa: E402

METHODS = ("AF2", "AF3", "Boltz2", "BioEmu", "ESM", "DDG", "S4PRED")
LABEL_SOURCES = {
    "raw": ("fold_diversity_survey.csv", "TMdiff_centered"),
    "corrected": ("fold_diversity_survey_corrected_9method.csv", "TMdiff_residual"),
}
STATS = ("inorm", "z", "T")


# ---------------------------------------------------------------------------
# Core: one tree, random-trait matrix, fold-preference values
# ---------------------------------------------------------------------------

# Set by --dump-null: keep every random-trait I_norm draw, not just its
# summary, so the SI figure can show the observed value inside its own null.
DUMP_NULL = False


def tree_precompute(tree, names: Sequence[str], min_clade: int = 2) -> dict:
    W = pt.morans_weights(tree, names)
    M, _, evecs = pt.morans_spectral(W)
    S, _ = pt.clade_split_matrix(tree, names, min_clade)
    return {"W": W, "M": M, "evecs": evecs, "S": S, "min_clade": min_clade}


def trait_stats(pre: dict, x: np.ndarray) -> Tuple[float, float, float, float]:
    """(I, I_norm, z, T) of values ``x`` on the precomputed leaf set."""
    I, I_norm, z = pt.morans_fast_stats(pre["W"], pre["M"], pre["evecs"], x)
    T = pt.clade_T_stat(pre["S"], x, pre["min_clade"])
    return I, I_norm, z, T


def exceed_p(null: np.ndarray, obs: float) -> float:
    """(1 + #{null >= obs}) / (K + 1); NaN null values never exceed."""
    if not np.isfinite(obs):
        return float("nan")
    K = len(null)
    return (1 + int(np.sum(np.isfinite(null) & (null >= obs - 1e-12)))) / (K + 1)


def trait_null_rows(tree, rand_values: Dict[str, np.ndarray],
                    fold_values: Dict[str, Dict[str, Optional[float]]],
                    min_clusters: int = 10, min_clade: int = 2,
                    null_sink: Optional[dict] = None) -> List[dict]:
    """Per-method rows for one tree and one label set.

    rand_values  {leaf_label: (K,) random-trait values}
    fold_values  {method: {leaf_label: fold-preference value or None}}
    The random traits are evaluated on EXACTLY the method's leaf set (finite
    fold value and finite random traits), cached across methods.
    """
    leaf_order = [lf.name for lf in tree.iter_leaves()]
    K = len(next(iter(rand_values.values()))) if rand_values else 0
    cache: Dict[tuple, Tuple[dict, np.ndarray]] = {}
    rows = []
    for m, vals in fold_values.items():
        names = [nm for nm in leaf_order
                 if nm in rand_values and np.all(np.isfinite(rand_values[nm]))
                 and vals.get(nm) is not None and np.isfinite(vals[nm])]
        row = {"method": m, "n_clusters": len(names), "K": K, "I": np.nan,
               "I_norm": np.nan, "z": np.nan, "T": np.nan,
               "p_trait_inorm": np.nan, "p_trait_z": np.nan, "p_trait_T": np.nan,
               "rand_I_norm_median": np.nan, "rand_I_norm_q95": np.nan,
               "rand_z_median": np.nan, "rand_T_median": np.nan, "note": ""}
        if len(names) < min_clusters:
            row["note"] = f"n<{min_clusters}"
            rows.append(row)
            continue
        key = tuple(names)
        if key not in cache:
            pre = tree_precompute(tree, names, min_clade)
            X = np.array([rand_values[nm] for nm in names]).T      # (K, n)
            null = np.array([trait_stats(pre, X[k])[1:] for k in range(K)])
            cache[key] = (pre, null)
        pre, null = cache[key]
        x = np.array([vals[nm] for nm in names], dtype=float)
        I, I_norm, z, T = trait_stats(pre, x)
        if not np.isfinite(I):
            row["note"] = "zero_variance"
        row.update({"I": I, "I_norm": I_norm, "z": z, "T": T,
                    "p_trait_inorm": exceed_p(null[:, 0], I_norm),
                    "p_trait_z": exceed_p(null[:, 1], z),
                    "p_trait_T": exceed_p(null[:, 2], T),
                    "rand_I_norm_median": float(np.nanmedian(null[:, 0])),
                    "rand_I_norm_q95": float(np.nanquantile(null[:, 0], 0.95)),
                    "rand_z_median": float(np.nanmedian(null[:, 1])),
                    "rand_T_median": float(np.nanmedian(null[:, 2]))})
        if null_sink is not None and np.isfinite(I_norm):
            # Keep the K random-trait I_norm draws behind this row so the
            # observed value can be plotted against its own null
            # (SI figure; scripts/plot_trait_null_dist.py).
            null_sink[m] = np.asarray(null[:, 0], dtype=float)
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Per-pair worker (real data)
# ---------------------------------------------------------------------------

def cluster_random_traits(pair_dir: str, tags: Sequence[str], R: np.ndarray
                          ) -> Dict[str, np.ndarray]:
    """{tag: (K,) random-trait values} from output_msa_cluster/<tag>.a3m."""
    out = {}
    for tag in tags:
        f = os.path.join(pair_dir, "output_msa_cluster", f"{tag}.a3m")
        if not os.path.exists(f):
            continue
        seqs = load_a3m_seqs(f)
        if len(seqs) <= 1:
            continue
        v = random_proj_values(seqs[1:], R)
        if np.all(np.isfinite(v)):
            out[tag] = v
    return out


def _alarm(signum, frame):
    raise TimeoutError("timed out")


def analyze_pair(job) -> dict:
    (pair_id, fold_by_label, n_random, rp_seed, min_clusters, min_clade,
     timeout) = job
    from config import DATA_DIR
    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(timeout)
    t0 = time.time()
    try:
        from utils.ancestral_utils import build_coarse_tree, _name_internal_nodes
        tree, tag_to_label, _ = build_coarse_tree(pair_id)
        _name_internal_nodes(tree)
        n_leaves = sum(1 for _ in tree.iter_leaves())
        if n_leaves < min_clusters:
            return {"pair_id": pair_id, "error": f"{n_leaves} leaves < {min_clusters}"}
        R = random_score_matrix(n_random, rp_seed)
        rt = cluster_random_traits(os.path.join(DATA_DIR, pair_id), list(tag_to_label), R)
        rand_values = {tag_to_label[t]: v for t, v in rt.items()}
        rows = []
        nulls = {}
        for labels, by_method in fold_by_label.items():
            fv = {m: {tag_to_label[t]: v for t, v in d.items() if t in tag_to_label}
                  for m, d in by_method.items()}
            sink = {} if DUMP_NULL else None
            for r in trait_null_rows(tree, rand_values, fv, min_clusters,
                                     min_clade, null_sink=sink):
                rows.append({"pair_id": pair_id, "labels": labels,
                             "n_leaves": n_leaves, **r})
            if sink:
                for m, v in sink.items():
                    nulls[f"{pair_id}|{labels}|{m}"] = v
        out = {"pair_id": pair_id, "rows": rows, "secs": time.time() - t0}
        if nulls:
            out["nulls"] = nulls
        return out
    except Exception as e:  # noqa: BLE001
        return {"pair_id": pair_id, "error": f"{type(e).__name__}: {e}"}
    finally:
        signal.alarm(0)


# ---------------------------------------------------------------------------
# Family level + report
# ---------------------------------------------------------------------------

def family_table(df: pd.DataFrame) -> pd.DataFrame:
    """Bonferroni-over-methods p_fam per (pair, labels) and BH q over families."""
    stat_col = {"inorm": "I_norm", "z": "z", "T": "T"}
    fam = []
    for (pair, labels), g in df.groupby(["pair_id", "labels"], sort=True):
        r = {"pair_id": pair, "labels": labels,
             "n_clusters": int(g.n_clusters.max()) if len(g) else 0}
        for s in STATS:
            pc = f"p_trait_{s}"
            ok = g[np.isfinite(g[pc])]
            r[f"n_methods_{s}"] = len(ok)
            if ok.empty:
                r[f"p_fam_{s}"] = np.nan
                r[f"best_method_{s}"] = ""
                continue
            best = ok.sort_values([pc, stat_col[s]], ascending=[True, False]).iloc[0]
            r[f"p_fam_{s}"] = min(1.0, len(ok) * float(best[pc]))
            r[f"best_method_{s}"] = best["method"]
            r[f"best_{stat_col[s]}"] = float(best[stat_col[s]])
        fam.append(r)
    fam = pd.DataFrame(fam)
    if fam.empty:
        return fam
    for s in STATS:
        fam[f"q_{s}"] = np.nan
        for labels, idx in fam.groupby("labels").groups.items():
            fam.loc[idx, f"q_{s}"] = bh_adjust(fam.loc[idx, f"p_fam_{s}"].values)
    return fam


def load_names(path: str) -> Dict[str, str]:
    if not os.path.exists(path):
        return {}
    d = pd.read_csv(path)
    if "display_name" not in d.columns:
        return {}
    return dict(zip(d.pair_id.astype(str), d.display_name.fillna("").astype(str)))


def print_summary(df: pd.DataFrame, fam: pd.DataFrame, names: Dict[str, str],
                  alpha: float = 0.05) -> None:
    print("\n=== Fold preference vs K random sequence traits on the same tree ===")
    ok = df[np.isfinite(df.p_trait_inorm)]
    for labels, g in ok.groupby("labels"):
        print(f"[{labels}] {g.pair_id.nunique()} families, {len(g)} (family, method) "
              f"rows; median p_trait I_norm {g.p_trait_inorm.median():.3f}; "
              f"fraction p_trait_inorm <= {alpha}: {(g.p_trait_inorm <= alpha).mean():.3f}")
        per_m = g.groupby("method").p_trait_inorm.apply(lambda p: (p <= alpha).mean())
        print("   per-method fraction p_trait_inorm <= 0.05: "
              + ", ".join(f"{m} {v:.2f}" for m, v in per_m.items()))
    for labels, g in fam.groupby("labels"):
        for s in STATS:
            sig = g[g[f"q_{s}"] <= alpha].sort_values(f"q_{s}")
            n_t = int(np.isfinite(g[f"p_fam_{s}"]).sum())
            print(f"[{labels}] statistic {s:5s}: {len(sig)}/{n_t} families with "
                  f"BH q <= {alpha}  (p_fam <= {alpha}: "
                  f"{int((g[f'p_fam_{s}'] <= alpha).sum())})")
            for _, r in sig.iterrows():
                nm = names.get(r.pair_id, "")
                print(f"     {r.pair_id:14s} q={r[f'q_{s}']:.3g} best={r[f'best_method_{s}']}"
                      + (f"  {nm}" if nm else ""))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def load_fold_values(docs_dir: str, label_sets: Sequence[str], methods: Sequence[str]
                     ) -> Dict[str, Dict[str, Dict[str, Dict[str, float]]]]:
    """{pair: {labels: {method: {cluster_tag: value}}}} (finite values only)."""
    from phylo_placement import load_fold_value
    out: Dict[str, dict] = {}
    for labels in label_sets:
        fname, col = LABEL_SOURCES[labels]
        path = os.path.join(docs_dir, fname)
        if not os.path.exists(path):
            print(f"[warn] {path} missing; skipping {labels}")
            continue
        for (pair, m, tag), v in load_fold_value(path, col).items():
            if m in methods and v is not None and np.isfinite(v):
                out.setdefault(pair, {}).setdefault(labels, {}).setdefault(m, {})[tag] = v
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--pairs", default="ALL")
    ap.add_argument("--max-pairs", type=int, default=0)
    ap.add_argument("--labels", default="raw,corrected")
    ap.add_argument("--methods", default=",".join(METHODS))
    ap.add_argument("--n-random", type=int, default=200, help="K random sequence traits")
    ap.add_argument("--rp-seed", type=int, default=777)
    ap.add_argument("--min-clusters", type=int, default=10)
    ap.add_argument("--min-clade", type=int, default=2)
    ap.add_argument("--n-jobs", type=int, default=1)
    ap.add_argument("--timeout", type=int, default=1800, help="per-pair seconds")
    ap.add_argument("--docs-dir", default=os.path.join(REPO, "docs"))
    ap.add_argument("--names-csv", default=os.path.join(REPO, "docs", "pair_protein_names.csv"))
    ap.add_argument("--out", default=os.path.join(REPO, "docs", "trait_null_phylo.csv"))
    ap.add_argument("--out-family", default=os.path.join(REPO, "docs",
                                                         "trait_null_phylo_family.csv"))
    ap.add_argument("--summary-only", action="store_true")
    ap.add_argument("--dump-null", default="",
                    help="NPZ path: also save every random-trait I_norm draw "
                         "(key '<pair>|<labels>|<method>'), so the observed "
                         "value can be plotted inside its own null. Use with "
                         "--n-jobs 1 and a short --pairs list.")
    args = ap.parse_args(argv)
    t0 = time.time()
    if args.dump_null:
        globals()["DUMP_NULL"] = True

    if args.summary_only:
        df = pd.read_csv(args.out)
    else:
        labels = [s.strip() for s in args.labels.split(",") if s.strip()]
        methods = [s.strip() for s in args.methods.split(",") if s.strip()]
        fold = load_fold_values(args.docs_dir, labels, methods)
        pairs = (discover_pairs() if args.pairs.upper() == "ALL"
                 else [p.strip() for p in args.pairs.split(",") if p.strip()])
        pairs = [p for p in pairs if p in fold]
        if args.max_pairs:
            pairs = pairs[: args.max_pairs]
        print(f"[pairs] {len(pairs)} pairs with fold values; K={args.n_random}; "
              f"labels={labels}; methods={methods}", flush=True)
        jobs = [(p, fold[p], args.n_random, args.rp_seed, args.min_clusters,
                 args.min_clade, args.timeout) for p in pairs]
        rows = []
        if args.n_jobs > 1:
            with mp.Pool(args.n_jobs) as pool:
                it = pool.imap_unordered(analyze_pair, jobs, chunksize=1)
                res = []
                for k, r in enumerate(it, 1):
                    res.append(r)
                    print(f"  [{k}/{len(jobs)}] {r['pair_id']}: "
                          + (r["error"] if "error" in r else f"{r['secs']:.0f}s"), flush=True)
        else:
            res = []
            for k, j in enumerate(jobs, 1):
                r = analyze_pair(j)
                res.append(r)
                print(f"  [{k}/{len(jobs)}] {r['pair_id']}: "
                      + (r["error"] if "error" in r else f"{r['secs']:.0f}s"), flush=True)
        for r in sorted(res, key=lambda d: d["pair_id"]):
            rows += r.get("rows", [])
        if not rows:
            sys.exit("[error] no pairs analysed")
        if args.dump_null:
            nulls = {}
            for r in res:
                nulls.update(r.get("nulls", {}))
            np.savez_compressed(args.dump_null, **nulls)
            print(f"[dump] {len(nulls)} null vectors -> {args.dump_null}", flush=True)
        df = pd.DataFrame(rows)
        lead = ["pair_id", "method", "labels", "I_norm", "z", "T", "p_trait_inorm",
                "p_trait_z", "p_trait_T", "K", "n_clusters"]
        df = df[lead + [c for c in df.columns if c not in lead]]
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        df.to_csv(args.out, index=False)
        print(f"[out] {args.out}")
    fam = family_table(df)
    names = load_names(args.names_csv)
    if not fam.empty:
        fam.insert(2, "name", fam.pair_id.map(names).fillna(""))
    fam.to_csv(args.out_family, index=False)
    print(f"[out] {args.out_family}")
    print_summary(df, fam, names)
    print(f"[done] {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
