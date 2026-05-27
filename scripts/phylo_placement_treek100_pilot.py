#!/usr/bin/env python3
"""
Phylogenetic-placement p-value for the K=100 + medoid + top-10 pilot.

For each of the 4 pilot modes (AF2/AF3 x full/top10) on a given pair:
  1. Classify each K=100 cluster as F1c / F2c / ambiguous by the sign of
     TMdiff_centered with threshold delta=0.05.
  2. On the K=100 cluster-representative tree (pruned from
     output_phytree/DeepMsa_tree.nwk), compute:
        - 1-nearest-neighbour concordance (NN_1):
            fraction of labelled leaves whose nearest labelled neighbour
            (by tree distance) shares the same label
        - Fitch parsimony minimum changes
  3. Permutation null: shuffle labels across labelled leaves N=10000 times,
     recompute NN_1 and Fitch, count tail-probability vs. observed.

Comparable to the analysis in scripts/phylo_placement.py used in the
existing paper (KaiB at K=22 yielded p_NN1 = 0.005 on DDG).

Usage:
    python3 scripts/phylo_placement_treek100_pilot.py
    python3 scripts/phylo_placement_treek100_pilot.py --pair 5jytA_2qkeE
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy.stats import combine_pvalues

try:
    from Bio import Phylo
except ImportError:
    print("biopython is required.", file=sys.stderr)
    sys.exit(1)

ROOT = Path(__file__).resolve().parents[1]


def read_cluster_members(a3m_path: Path) -> list[str]:
    out = []
    with a3m_path.open() as f:
        for line in f:
            if line.startswith(">"):
                out.append(line[1:].split()[0])
    return out


def pairwise_tree_distance_matrix(tree, leaf_names: list[str]) -> np.ndarray:
    """Compute pairwise tree distances among the given leaves."""
    n = len(leaf_names)
    D = np.full((n, n), np.inf, dtype=float)
    leaf_objs = {l.name: l for l in tree.get_terminals() if l.name in leaf_names}
    for i, ni in enumerate(leaf_names):
        if ni not in leaf_objs:
            continue
        for j, nj in enumerate(leaf_names):
            if i == j:
                D[i, j] = 0.0
                continue
            if nj not in leaf_objs:
                continue
            try:
                D[i, j] = float(tree.distance(leaf_objs[ni], leaf_objs[nj]))
            except Exception:
                pass
    return D


def nn1_concordance(labels: np.ndarray, dist: np.ndarray) -> float:
    """1-NN concordance: fraction of labelled leaves whose nearest labelled
    neighbour (excluding self) shares the same label."""
    mask = labels != ""
    if mask.sum() < 2:
        return float("nan")
    idx = np.where(mask)[0]
    matches = 0
    n = 0
    for i in idx:
        d = dist[i].copy()
        d[i] = np.inf
        # restrict to labelled others
        for j in range(len(labels)):
            if not mask[j] or j == i:
                d[j] = np.inf if j != i else d[j]
        if not np.isfinite(np.nanmin(d)):
            continue
        j_star = int(np.nanargmin(d))
        if labels[j_star] == labels[i]:
            matches += 1
        n += 1
    return matches / n if n > 0 else float("nan")


def fitch_min_changes(tree, labels_by_leaf: dict[str, str]) -> int:
    """Approximate Fitch parsimony: minimum label changes along the tree
    given the labels on leaves. Uses simple postorder Fitch on the
    tree topology.
    """
    # Build a dict of node -> set of possible labels
    def is_leaf(clade):
        return clade.is_terminal()

    label_sets = {}
    changes = [0]

    def post(clade):
        if is_leaf(clade):
            lab = labels_by_leaf.get(clade.name, "")
            if lab:
                label_sets[id(clade)] = {lab}
            else:
                label_sets[id(clade)] = {"F1c", "F2c"}  # unknown
            return
        children = clade.clades
        for c in children:
            post(c)
        # Intersection rule
        sets = [label_sets[id(c)] for c in children if id(c) in label_sets]
        if not sets:
            label_sets[id(clade)] = set()
            return
        intersect = set.intersection(*sets)
        if intersect:
            label_sets[id(clade)] = intersect
        else:
            label_sets[id(clade)] = set.union(*sets)
            changes[0] += 1

    post(tree.root)
    return changes[0]


def run_phylo_test(tree, labels_by_leaf: dict[str, str],
                   leaf_names: list[str], n_perm: int = 10000,
                   seed: int = 42) -> dict:
    """Compute NN_1 and Fitch parsimony observed + permutation p-values."""
    rng = np.random.default_rng(seed)

    # Distance matrix on labelled leaves only (faster)
    D = pairwise_tree_distance_matrix(tree, leaf_names)
    labels = np.array([labels_by_leaf.get(n, "") for n in leaf_names])

    nn1_obs = nn1_concordance(labels, D)
    fitch_obs = fitch_min_changes(tree, labels_by_leaf)

    # Permutation null: shuffle non-empty labels among labelled positions
    non_empty_mask = labels != ""
    non_empty_labels = labels[non_empty_mask].copy()

    nn1_null = []
    fitch_null = []
    for _ in range(n_perm):
        perm = non_empty_labels.copy()
        rng.shuffle(perm)
        labels_perm = labels.copy()
        labels_perm[non_empty_mask] = perm
        nn1_null.append(nn1_concordance(labels_perm, D))
        labels_by_leaf_perm = {n: l for n, l in zip(leaf_names, labels_perm)}
        fitch_null.append(fitch_min_changes(tree, labels_by_leaf_perm))

    nn1_null_arr = np.array(nn1_null)
    fitch_null_arr = np.array(fitch_null)

    # p-values: NN_1 is "higher = better"; Fitch is "lower = better"
    p_nn1 = (np.sum(nn1_null_arr >= nn1_obs) + 1) / (n_perm + 1)
    p_fitch = (np.sum(fitch_null_arr <= fitch_obs) + 1) / (n_perm + 1)

    return {
        "n_labelled": int(non_empty_mask.sum()),
        "n_F1c": int(np.sum(labels == "F1c")),
        "n_F2c": int(np.sum(labels == "F2c")),
        "nn1_obs": float(nn1_obs),
        "nn1_null_mean": float(np.mean(nn1_null_arr)),
        "p_nn1": float(p_nn1),
        "fitch_obs": int(fitch_obs),
        "fitch_null_mean": float(np.mean(fitch_null_arr)),
        "p_fitch": float(p_fitch),
        "n_perm": int(n_perm),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pair", default="5jytA_2qkeE")
    ap.add_argument("--scores-csv",
                    default=str(ROOT / "docs" / "treek100_pilot_tmscores.csv"))
    ap.add_argument("--delta", type=float, default=0.05,
                    help="TMdiff_centered threshold for F1c/F2c classification "
                         "(default 0.05)")
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--output", default=None,
                    help="Output CSV path (default: docs/treek100_pilot_phylo.csv)")
    args = ap.parse_args()

    pair = args.pair
    pair_dir = ROOT / "Pipeline" / "FoldPairs" / pair
    tree_path = pair_dir / "output_phytree" / "DeepMsa_tree.nwk"
    cluster_dir = pair_dir / "output_msa_cluster"
    if not tree_path.is_file():
        sys.exit(f"Tree not found: {tree_path}")
    out_csv = args.output or str(ROOT / "docs" / f"treek100_pilot_phylo_{pair}.csv")

    print(f"[phylo] pair={pair}")
    print(f"[phylo] loading tree: {tree_path}")
    tree = Phylo.read(str(tree_path), "newick")
    tree_leaves = {l.name for l in tree.get_terminals() if l.name}
    print(f"[phylo] tree has {len(tree_leaves)} leaves")

    # Map each cluster -> its representative leaf
    cluster_reps: dict[str, str] = {}
    rep_to_cluster: dict[str, str] = {}
    for a3m in sorted(cluster_dir.glob("ShallowMsa_*.a3m")):
        cluster_tag = a3m.stem
        for m in read_cluster_members(a3m):
            if m in tree_leaves and m not in rep_to_cluster:
                cluster_reps[cluster_tag] = m
                rep_to_cluster[m] = cluster_tag
                break
    print(f"[phylo] resolved {len(cluster_reps)} cluster reps")

    # Prune the tree to the cluster reps
    keep = set(cluster_reps.values())
    drop = [l for l in tree.get_terminals() if l.name not in keep]
    print(f"[phylo] pruning {len(drop)} non-representative leaves...")
    for leaf in drop:
        try:
            tree.prune(leaf)
        except Exception:
            pass
    leaf_names = [l.name for l in tree.get_terminals()]
    print(f"[phylo] post-prune leaf count: {len(leaf_names)}")

    # Load TMdiff_centered, average across chains per (mode, cluster)
    scores = pd.read_csv(args.scores_csv)
    pivot = (scores.dropna(subset=["TMdiff_centered"])
             .groupby(["mode", "cluster"])["TMdiff_centered"]
             .mean().unstack("mode"))

    # For each mode, classify clusters and run the phylo test
    results = []
    for mode in pivot.columns:
        labels_by_leaf: dict[str, str] = {}
        n_f1, n_f2 = 0, 0
        for cluster_tag, rep_leaf in cluster_reps.items():
            if rep_leaf not in leaf_names:
                continue
            if cluster_tag not in pivot.index:
                continue
            v = pivot.loc[cluster_tag, mode]
            if pd.isna(v):
                continue
            if v > args.delta:
                labels_by_leaf[rep_leaf] = "F1c"
                n_f1 += 1
            elif v < -args.delta:
                labels_by_leaf[rep_leaf] = "F2c"
                n_f2 += 1
        if n_f1 < 2 or n_f2 < 2:
            print(f"[phylo] {mode}: n_F1c={n_f1}, n_F2c={n_f2} -- insufficient, skipping")
            continue
        print(f"[phylo] === {mode}: n_F1c={n_f1}, n_F2c={n_f2} ===")
        r = run_phylo_test(tree, labels_by_leaf, leaf_names, n_perm=args.n_perm)
        r["mode"] = mode
        r["delta"] = args.delta
        results.append(r)
        print(f"  NN_1 obs={r['nn1_obs']:.3f}  null_mean={r['nn1_null_mean']:.3f}  "
              f"p={r['p_nn1']:.5f}")
        print(f"  Fitch obs={r['fitch_obs']}  null_mean={r['fitch_null_mean']:.1f}  "
              f"p={r['p_fitch']:.5f}")

    df = pd.DataFrame(results)
    df.to_csv(out_csv, index=False)
    print(f"\n[phylo] wrote {out_csv}")

    # Summary line per mode
    print("\n[phylo] === SUMMARY ===")
    for r in results:
        print(f"  {r['mode']:<28s}  p_NN1={r['p_nn1']:.5f}  "
              f"p_Fitch={r['p_fitch']:.5f}  "
              f"(n_F1c={r['n_F1c']}, n_F2c={r['n_F2c']})")


if __name__ == "__main__":
    main()
