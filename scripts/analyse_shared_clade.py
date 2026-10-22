#!/usr/bin/env python3
"""Do independent fold-preference methods localise the SAME branch in a family?

Reads already-computed clade-shift scan results (medoid-representative,
7-method) and measures cross-method agreement of the best clade per family.

No trees are needed: everything is computed from ``clade_best_members``
(sets of cluster tags).  The null in part (2) is LABEL-FREE (random subsets
of the family's cluster set, matched on size) and therefore ignores tree
structure -- see SHARED_CLADE.md.
"""
from __future__ import annotations

import csv
import itertools
import json
import math
import os
from collections import defaultdict

import numpy as np

S = "/tmp/claude-1000/-mnt-c-Code-Github-MsaCluster/fea976ea-5c76-4e3b-a35d-95918917c7f5/scratchpad"
OUT = os.path.join(S, "sharedclade")
os.makedirs(OUT, exist_ok=True)

METHODS7 = ["AF2", "AF3", "Boltz2", "BioEmu", "ESM", "DDG", "S4PRED"]
STRUCT5 = ["AF2", "AF3", "Boltz2", "BioEmu", "ESM"]
SEQ2 = ["DDG", "S4PRED"]

DATASETS = {
    "raw": dict(
        placement=f"{S}/round3/clade_shift_results/v07/phylo_placement_v07.csv",
        survey=f"{S}/final_docs/docs/fold_diversity_survey.csv",
        value_col="TMdiff_centered",
    ),
    "corrected": dict(
        placement=f"{S}/round3/clade_shift_results/v07/phylo_placement_corrected_v07.csv",
        survey=f"{S}/final_docs/docs/fold_diversity_survey_corrected_9method.csv",
        value_col="TMdiff_residual",
    ),
}

NAMED = {
    "2n54B_2hdmA": "XCL1",
    "5jytA_2qkeE": "KaiB",
    "2namA_1uxmK": "SOD1",
    "5lj3M_3tp2A": "CWC2",
    "5ly6B_5aoeB": "pneumolysin",
    "5wrgC_1wyyB": "SARS-CoV spike",
    "1kctA_3t1pA": "alpha1-antitrypsin",
}

NPERM = 10000
THRESHOLDS = (0.3, 0.5, 0.8)
MAIN_THR = 0.5


# ---------------------------------------------------------------- loading
def fnum(s):
    try:
        v = float(s)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def load_universe(path, value_col):
    """(pair, method) -> set of cluster tags with a finite value."""
    uni = defaultdict(set)
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            if row["method"] not in METHODS7 or row["cluster"] == "DeepMsa":
                continue
            if fnum(row.get(value_col)) is not None:
                uni[(row["pair_id"], row["method"])].add(row["cluster"])
    return uni


def load_placement(path):
    """pair -> {method: record}, plus COMBINED, plus n_clusters."""
    by_pair = defaultdict(dict)
    ncl = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            p, m = row["pair_id"], row["method"]
            if row.get("n_clusters"):
                ncl[p] = int(float(row["n_clusters"]))
            if m not in METHODS7 and m != "COMBINED":
                continue
            mem = [t for t in row["clade_best_members"].split(";") if t]
            by_pair[p][m] = dict(
                members=frozenset(mem),
                size=len(mem),
                direction=int(float(row["clade_best_direction"]))
                if fnum(row["clade_best_direction"]) is not None else 0,
                p=fnum(row["clade_p"]),
                p_bh=fnum(row.get("clade_p_bh")),
                T=fnum(row["clade_T"]),
                n_leaves=int(float(row["clade_n_leaves"]))
                if fnum(row["clade_n_leaves"]) is not None else 0,
                effect=fnum(row["clade_best_effect"]),
                note=row["clade_note"],
            )
    return by_pair, ncl


def jaccard(a, b):
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


# ------------------------------------------------------------ clique code
def max_clique(methods, adj):
    """Largest subset of `methods` that is pairwise connected in dict adj."""
    best = []
    n = len(methods)
    for size in range(n, 1, -1):
        for sub in itertools.combinations(methods, size):
            if all(adj.get((a, b), False) for a, b in itertools.combinations(sub, 2)):
                return list(sub)
    return best


# --------------------------------------------------------------- the null
def null_pvalue(pair_methods, universes, obs_stat, rng, nperm=NPERM):
    """Empirical p for `obs_stat` (= max-clique size at MAIN_THR) under a
    label-free null: each method's best clade is replaced by a uniformly
    random subset of the SAME SIZE drawn from that method's own leaf
    universe; observed directions are kept fixed.  Random subsets are NOT
    clades, so this null is anti-conservative (see report)."""
    ms = list(pair_methods)
    k = len(ms)
    if k < 2 or obs_stat < 2:
        return 1.0, None
    # universe as index space: union of all methods' universes
    allu = sorted(set().union(*[universes[m] for m in ms]))
    idx = {t: i for i, t in enumerate(allu)}
    N = len(allu)
    uni_idx = [np.array(sorted(idx[t] for t in universes[m]), dtype=np.int64) for m in ms]
    sizes = [pair_methods[m]["size"] for m in ms]
    dirs = [pair_methods[m]["direction"] for m in ms]

    B = nperm
    draws = np.zeros((k, B, N), dtype=bool)
    for i, (u, sz) in enumerate(zip(uni_idx, sizes)):
        nu = len(u)
        if sz <= 0 or sz > nu:
            continue
        # random subset of size sz from u, for each of B draws
        r = rng.random((B, nu))
        sel = np.argpartition(r, sz - 1, axis=1)[:, :sz]
        cols = u[sel]                                    # (B, sz)
        rows = np.repeat(np.arange(B)[:, None], sz, axis=1)
        draws[i, rows, cols] = True

    pair_list = list(itertools.combinations(range(k), 2))
    edges = np.zeros((B, len(pair_list)), dtype=bool)
    for pi, (i, j) in enumerate(pair_list):
        if dirs[i] != dirs[j] or dirs[i] == 0:
            continue
        inter = (draws[i] & draws[j]).sum(axis=1).astype(float)
        union = (draws[i] | draws[j]).sum(axis=1).astype(float)
        with np.errstate(invalid="ignore", divide="ignore"):
            jac = np.where(union > 0, inter / union, 0.0)
        edges[:, pi] = jac >= MAIN_THR - 1e-12

    pidx = {pr: n for n, pr in enumerate(pair_list)}
    stat = np.full(B, 1, dtype=np.int64)
    for size in range(k, 1, -1):
        for sub in itertools.combinations(range(k), size):
            need = [pidx[pr] for pr in itertools.combinations(sub, 2)]
            ok = edges[:, need].all(axis=1)
            stat = np.where(ok & (stat < size), size, stat)
    n_ge = int((stat >= obs_stat).sum())
    return (1 + n_ge) / (1 + B), stat


# ------------------------------------------------- secondary null: random trees
def _random_tree_splits(N, rng, min_clade=2):
    """Splits (clade index sets) of a random agglomerative (coalescent-shape)
    tree on N leaves, canonicalised to the smaller side like
    utils.phylo_tests.clade_split_matrix."""
    active = [np.zeros(N, dtype=bool) for _ in range(N)]
    for i in range(N):
        active[i][i] = True
    rows = []
    while len(active) > 1:
        i, j = rng.choice(len(active), size=2, replace=False)
        i, j = (int(i), int(j)) if i < j else (int(j), int(i))
        merged = active[i] | active[j]
        active.pop(j)
        active.pop(i)
        active.append(merged)
        k = int(merged.sum())
        if min_clade <= k <= N - min_clade:
            rows.append(merged if 2 * k <= N else ~merged)
    if not rows:
        return np.zeros((0, N), dtype=bool)
    return np.vstack(rows)


def null_pvalue_randtree(pair_methods, universes, obs_stat, rng,
                         n_trees=100, per_tree=100, ignore_direction=False):
    """Empirical p under a null that keeps clade sizes and directions but
    places each method's clade uniformly among the (size-matched) clades of
    a RANDOM tree on the family's leaf set.  Unlike the label-free null this
    respects the nestedness of real clade sets; it still assumes a random
    (coalescent-shaped) topology rather than the true one."""
    ms = list(pair_methods)
    k = len(ms)
    if k < 2 or obs_stat < 2:
        return 1.0
    allu = sorted(set().union(*[universes[m] for m in ms]))
    idx = {t: i for i, t in enumerate(allu)}
    N = len(allu)
    if N < 4:
        return 1.0
    sizes = [pair_methods[m]["size"] for m in ms]
    dirs = [pair_methods[m]["direction"] for m in ms]
    # leaves each method actually scored (others can never be in its clade)
    masks = [np.zeros(N, dtype=bool) for _ in ms]
    for mi, m in enumerate(ms):
        for t in universes[m]:
            masks[mi][idx[t]] = True

    pair_list = list(itertools.combinations(range(k), 2))
    pidx = {pr: n for n, pr in enumerate(pair_list)}
    stats = []
    for _ in range(n_trees):
        M = _random_tree_splits(N, rng)
        if M.shape[0] == 0:
            continue
        ssz = M.sum(axis=1)
        pick = np.zeros((k, per_tree, N), dtype=bool)
        for mi, sz in enumerate(sizes):
            tol = max(1, int(round(0.25 * sz)))
            elig = np.flatnonzero(np.abs(ssz - sz) <= tol)
            if elig.size == 0:
                d = np.abs(ssz - sz)
                elig = np.flatnonzero(d == d.min())
            ch = rng.choice(elig, size=per_tree, replace=True)
            # restrict to the leaves this method scored
            pick[mi] = M[ch] & masks[mi][None, :]
        edges = np.zeros((per_tree, len(pair_list)), dtype=bool)
        for pi, (i, j) in enumerate(pair_list):
            if not ignore_direction and (dirs[i] != dirs[j] or dirs[i] == 0):
                continue
            inter = (pick[i] & pick[j]).sum(axis=1).astype(float)
            union = (pick[i] | pick[j]).sum(axis=1).astype(float)
            with np.errstate(invalid="ignore", divide="ignore"):
                jac = np.where(union > 0, inter / union, 0.0)
            edges[:, pi] = jac >= MAIN_THR - 1e-12
        st = np.ones(per_tree, dtype=np.int64)
        for size in range(k, 1, -1):
            for sub in itertools.combinations(range(k), size):
                need = [pidx[pr] for pr in itertools.combinations(sub, 2)]
                ok = edges[:, need].all(axis=1)
                st = np.where(ok & (st < size), size, st)
        stats.append(st)
    if not stats:
        return 1.0
    st = np.concatenate(stats)
    return (1 + int((st >= obs_stat).sum())) / (1 + st.size)


# ------------------------------------------------------------------- main
def analyse():
    rng = np.random.default_rng(20261021)
    per_family_rows = []
    pairwise_rows = []
    consensus = {}            # (label, pair) -> dict

    for label, cfg in DATASETS.items():
        placement, ncl = load_placement(cfg["placement"])
        universes = load_universe(cfg["survey"], cfg["value_col"])

        for pair in sorted(placement):
            recs = placement[pair]
            valid = {m: recs[m] for m in METHODS7
                     if m in recs and recs[m]["size"] > 0 and recs[m]["T"] is not None}
            n_clusters = ncl.get(pair, 0)

            # --- pairwise table
            adj = {thr: {} for thr in THRESHOLDS}
            adj_nd = {}
            for a, b in itertools.combinations(METHODS7, 2):
                if a not in valid or b not in valid:
                    continue
                ja = jaccard(valid[a]["members"], valid[b]["members"])
                sd = valid[a]["direction"] == valid[b]["direction"] and valid[a]["direction"] != 0
                for thr in THRESHOLDS:
                    adj[thr][(a, b)] = adj[thr][(b, a)] = bool(sd and ja >= thr - 1e-12)
                adj_nd[(a, b)] = adj_nd[(b, a)] = bool(ja >= MAIN_THR - 1e-12)
                pairwise_rows.append(dict(
                    label=label, pair_id=pair, family=NAMED.get(pair, ""),
                    method_a=a, method_b=b,
                    size_a=valid[a]["size"], size_b=valid[b]["size"],
                    n_intersect=len(valid[a]["members"] & valid[b]["members"]),
                    n_union=len(valid[a]["members"] | valid[b]["members"]),
                    jaccard=round(ja, 4),
                    dir_a=valid[a]["direction"], dir_b=valid[b]["direction"],
                    same_direction=int(sd),
                    agree_0p3=int(adj[0.3][(a, b)]),
                    agree_0p5=int(adj[0.5][(a, b)]),
                    agree_0p8=int(adj[0.8][(a, b)]),
                    agree_0p5_nodir=int(adj_nd[(a, b)]),
                    identical_clade=int(valid[a]["members"] == valid[b]["members"]),
                    p_a=valid[a]["p"], p_b=valid[b]["p"],
                ))

            ms = list(valid)
            cl = {thr: max_clique(ms, adj[thr]) for thr in THRESHOLDS}
            cl_nd = max_clique(ms, adj_nd)
            grp = cl[MAIN_THR]
            k_agree = len(grp)

            # representative clade of the consensus group: the member set of
            # the group method with the smallest clade_p (ties -> largest T)
            if grp:
                rep = min(grp, key=lambda m: (valid[m]["p"] if valid[m]["p"] is not None else 1.0,
                                              -(valid[m]["T"] or 0)))
                rep_members = valid[rep]["members"]
                inter = frozenset.intersection(*[valid[m]["members"] for m in grp])
                union = frozenset.union(*[valid[m]["members"] for m in grp])
                direction = valid[rep]["direction"]
                ps = [valid[m]["p"] for m in grp if valid[m]["p"] is not None]
            else:
                rep, rep_members, inter, union, direction, ps = "", frozenset(), frozenset(), frozenset(), 0, []

            # --- robustness subsets
            def sub_support(subset):
                sm = [m for m in subset if m in valid]
                c = max_clique(sm, adj[MAIN_THR])
                if len(c) < 2 or not rep_members:
                    return 0, c
                # does this independent clique point at the same clade?
                ok = all(jaccard(valid[m]["members"], rep_members) >= MAIN_THR - 1e-12
                         and valid[m]["direction"] == direction for m in c)
                return int(ok), c

            sup_noddg, c_noddg = sub_support([m for m in METHODS7 if m != "DDG"])
            sup_nos4, c_nos4 = sub_support([m for m in METHODS7 if m != "S4PRED"])
            sup_struct, c_struct = sub_support(STRUCT5)
            sup_seq, c_seq = sub_support(SEQ2)

            n_struct_in_grp = sum(1 for m in grp if m in STRUCT5)
            n_seq_in_grp = sum(1 for m in grp if m in SEQ2)
            strong = int(n_struct_in_grp >= 1 and n_seq_in_grp >= 1 and k_agree >= 2)

            # --- null
            uni_pair = {m: universes.get((pair, m), set()) for m in valid}
            if k_agree >= 2:
                pnull, _ = null_pvalue(valid, uni_pair, k_agree, rng)
                pnull_rt = null_pvalue_randtree(valid, uni_pair, k_agree, rng)
            else:
                pnull = pnull_rt = 1.0
            pnull_rt_nd = (null_pvalue_randtree(valid, uni_pair, len(cl_nd), rng,
                                                ignore_direction=True)
                           if len(cl_nd) >= 2 else 1.0)

            comb = recs.get("COMBINED")
            comb_j = jaccard(comb["members"], rep_members) if comb and rep_members else 0.0

            consensus[(label, pair)] = dict(members=rep_members, direction=direction,
                                            group=grp, k=k_agree)

            per_family_rows.append(dict(
                label=label, pair_id=pair, family=NAMED.get(pair, ""),
                n_clusters=n_clusters,
                n_methods_valid=len(valid),
                n_agree=k_agree,
                agree_methods=",".join(grp),
                n_agree_nodir=len(cl_nd), agree_methods_nodir=",".join(cl_nd),
                n_exact_same_clade=max(
                    [sum(1 for m2 in ms if valid[m2]["members"] == valid[m]["members"])
                     for m in ms] or [0]),
                n_agree_j03=len(cl[0.3]), agree_methods_j03=",".join(cl[0.3]),
                n_agree_j08=len(cl[0.8]), agree_methods_j08=",".join(cl[0.8]),
                consensus_rep_method=rep,
                consensus_size=len(rep_members),
                consensus_frac=round(len(rep_members) / n_clusters, 4) if n_clusters else "",
                consensus_direction=direction,
                consensus_fold=("F1" if direction > 0 else "F2" if direction < 0 else ""),
                consensus_intersect_size=len(inter),
                consensus_union_size=len(union),
                p_max_in_group=max(ps) if ps else "",
                p_min_in_group=min(ps) if ps else "",
                p_geomean_in_group=round(
                    math.exp(sum(math.log(max(x, 1e-12)) for x in ps) / len(ps)), 6) if ps else "",
                n_struct_in_group=n_struct_in_grp,
                n_seq_in_group=n_seq_in_grp,
                strong_class=strong,
                support_drop_DDG=sup_noddg, clique_drop_DDG=",".join(c_noddg),
                support_drop_S4PRED=sup_nos4, clique_drop_S4PRED=",".join(c_nos4),
                support_struct5=sup_struct, clique_struct5=",".join(c_struct),
                support_DDG_S4PRED=sup_seq, clique_DDG_S4PRED=",".join(c_seq),
                p_null_labelfree=round(pnull, 6),
                p_null_randtree=round(pnull_rt, 6),
                p_null_randtree_nodir=round(pnull_rt_nd, 6),
                combined_jaccard_vs_consensus=round(comb_j, 4),
                combined_size=comb["size"] if comb else "",
                combined_direction=comb["direction"] if comb else "",
                combined_p=comb["p"] if comb else "",
                consensus_members=";".join(sorted(rep_members)),
            ))

    # raw-vs-corrected consistency
    for row in per_family_rows:
        other = "corrected" if row["label"] == "raw" else "raw"
        a = consensus.get((row["label"], row["pair_id"]))
        b = consensus.get((other, row["pair_id"]))
        if a and b and a["members"] and b["members"]:
            j = jaccard(a["members"], b["members"])
            row["jaccard_vs_other_label"] = round(j, 4)
            row["same_clade_both_labels"] = int(j >= MAIN_THR - 1e-12
                                                and a["direction"] == b["direction"]
                                                and a["k"] >= 2 and b["k"] >= 2)
        else:
            row["jaccard_vs_other_label"] = ""
            row["same_clade_both_labels"] = 0

    # ---- write
    ffields = list(per_family_rows[0].keys())
    with open(f"{OUT}/shared_clade_per_family.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=ffields)
        w.writeheader()
        per_family_rows.sort(key=lambda r: (r["label"], -r["n_agree"],
                                            r["p_max_in_group"] if r["p_max_in_group"] != "" else 1.0))
        w.writerows(per_family_rows)

    with open(f"{OUT}/shared_clade_pairwise.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(pairwise_rows[0].keys()))
        w.writeheader()
        w.writerows(pairwise_rows)

    with open(f"{OUT}/_per_family.json", "w") as fh:
        json.dump(per_family_rows, fh, indent=1)
    print(f"wrote {len(per_family_rows)} family rows, {len(pairwise_rows)} pairwise rows")


if __name__ == "__main__":
    analyse()
