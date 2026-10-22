#!/usr/bin/env python3
"""
How tree-consistent are the MSA clusters?  For every fold-switch pair, compare
the ShallowMsa_* cluster memberships with the pair's DeepMsa_tree.nwk:

  on_tree_frac   fraction of clustered sequences that are leaves of the tree
                 (the tree is built on a <=5000-sequence subsample; the rest are
                 assigned to clusters by sequence identity)
  monophyletic   cluster's on-tree members are exactly the leaves of their MRCA
  purity         |members on tree| / |leaves under their MRCA|  (1 = clade)
  n_pieces       number of maximal single-cluster clades the members split into

Writes one row per cluster and prints a per-pair / overall summary.

  python3 scripts/cluster_monophyly_census.py --pairs ALL --out docs/cluster_monophyly_census.csv
"""
import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import DATA_DIR  # noqa: E402


def _leaf_key(name: str) -> str:
    tok = (name or "").strip().strip("'\"").split()[0]
    return tok.split("/")[0]


def _cluster_members(a3m: str):
    out = []
    with open(a3m) as fh:
        for i, line in enumerate(l for l in fh if l.startswith(">")):
            if i == 0:          # row 0 = the cluster's query/medoid row
                continue
            out.append(_leaf_key(line[1:]))
    return out


def census_pair(pair_dir: str):
    from ete3 import Tree
    nwk = os.path.join(pair_dir, "output_phytree", "DeepMsa_tree.nwk")
    a3ms = sorted(glob.glob(os.path.join(pair_dir, "output_msa_cluster", "ShallowMsa_*.a3m")))
    if not os.path.isfile(nwk) or not a3ms:
        return []
    tree = Tree(nwk, format=1, quoted_node_names=True)
    leaf_of = {}
    for lf in tree.iter_leaves():
        leaf_of[_leaf_key(lf.name)] = lf
    cl_of_leaf = {}
    members = {}
    for f in a3ms:
        tag = re.search(r"ShallowMsa_\d+", f).group(0)
        mem = _cluster_members(f)
        members[tag] = mem
        for m in mem:
            if m in leaf_of:
                cl_of_leaf[m] = tag
    # Tree leaves that belong to no cluster (e.g. dropped by the gap filter
    # after the tree was built) must not count against monophyly: prune them.
    n_unassigned = len(leaf_of) - len(cl_of_leaf)
    if n_unassigned:
        tree.prune([leaf_of[m] for m in cl_of_leaf], preserve_branch_length=True)
    # One postorder pass: number of leaves under each node, and the node's
    # cluster label if ALL its leaves belong to one cluster (else None).
    MIXED = object()
    n_under, label = {}, {}
    for n in tree.traverse("postorder"):
        if n.is_leaf():
            n_under[n] = 1
            label[n] = cl_of_leaf.get(_leaf_key(n.name), MIXED)
        else:
            n_under[n] = sum(n_under[c] for c in n.children)
            labs = {label[c] for c in n.children}
            label[n] = labs.pop() if len(labs) == 1 else MIXED
    pieces_of = {}
    for n in tree.traverse("preorder"):
        lab = label[n]
        if lab is not MIXED and (n.up is None or label[n.up] != lab):
            pieces_of[lab] = pieces_of.get(lab, 0) + 1
    rows = []
    for tag, mem in members.items():
        on = [m for m in mem if m in leaf_of]
        row = {"pair_id": os.path.basename(pair_dir), "cluster": tag, "n_members": len(mem),
               "n_on_tree": len(on), "on_tree_frac": len(on) / len(mem) if mem else np.nan,
               "monophyletic": np.nan, "purity": np.nan, "n_pieces": np.nan}
        if len(on) >= 2:
            mrca = tree.get_common_ancestor([leaf_of[m] for m in on])
            row["purity"] = len(on) / n_under[mrca]
            row["monophyletic"] = bool(n_under[mrca] == len(on))
            row["n_pieces"] = pieces_of.get(tag, 0)
        elif len(on) == 1:
            row.update(monophyletic=True, purity=1.0, n_pieces=1)
        rows.append(row)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs", default="ALL")
    ap.add_argument("--out", default="docs/cluster_monophyly_census.csv")
    args = ap.parse_args()
    if args.pairs == "ALL":
        pairs = sorted(d for d in os.listdir(DATA_DIR)
                       if re.fullmatch(r"[0-9a-zA-Z]{5}_[0-9a-zA-Z]{5}", d))
    else:
        pairs = args.pairs.split(",")
    rows = []
    for p in pairs:
        try:
            r = census_pair(os.path.join(DATA_DIR, p))
        except Exception as e:  # keep going; report the pair
            print(f"[warn] {p}: {e}", file=sys.stderr)
            r = []
        rows += r
        if r:
            d = pd.DataFrame(r)
            print(f"{p}: clusters={len(d)} on_tree={d.n_on_tree.sum()/d.n_members.sum():.2f} "
                  f"monophyletic={np.nanmean(d.monophyletic.astype(float)):.2f} "
                  f"median_purity={d.purity.median():.2f} median_pieces={d.n_pieces.median():.0f}", flush=True)
    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    df.to_csv(args.out, index=False)
    g = df.groupby("pair_id").agg(on_tree=("n_on_tree", "sum"), n=("n_members", "sum"),
                                  mono=("monophyletic", lambda s: np.nanmean(s.astype(float))),
                                  purity=("purity", "median"), pieces=("n_pieces", "median"))
    g["on_tree"] = g.on_tree / g.n
    print("\n=== overall (per-pair medians across", len(g), "pairs) ===")
    print(f"fraction of sequences on the tree: median {g.on_tree.median():.2f} "
          f"(IQR {g.on_tree.quantile(.25):.2f}-{g.on_tree.quantile(.75):.2f})")
    print(f"fraction of clusters monophyletic: median {g.mono.median():.2f} "
          f"(IQR {g.mono.quantile(.25):.2f}-{g.mono.quantile(.75):.2f}); "
          f"pooled {np.nanmean(df.monophyletic.astype(float)):.2f} of {df.monophyletic.notna().sum()} clusters")
    print(f"cluster purity (members / MRCA leaves): median {df.purity.median():.2f}")
    print(f"pieces per cluster: median {df.n_pieces.median():.0f}, 90th pct {df.n_pieces.quantile(.9):.0f}")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
