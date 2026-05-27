#!/usr/bin/env python3
"""
Build a phylogenetic-tree + heatmap figure for the KaiB K=100 pilot.

For 5jytA_2qkeE specifically:
  - Use the per-sequence Deep MSA tree (output_phytree/DeepMsa_tree.nwk)
  - For each K=100 ShallowMsa cluster, pick a representative leaf
    (first sequence in the cluster .a3m that maps to a tree leaf)
  - Prune the tree to those ~99 representative leaves
  - Plot the pruned tree alongside a 4-column heatmap:
        AF2_treek100_full | AF2_treek100_top10 | AF3_treek100_full | AF3_treek100_top10

Result: docs/figs/5jytA_2qkeE_phytree_treek100.png

Usage:
    python3 scripts/plot_treek100_phytree.py
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd

try:
    from Bio import Phylo
except ImportError:
    print("biopython is required (Bio.Phylo). pip install biopython", file=sys.stderr)
    sys.exit(1)

ROOT = Path(__file__).resolve().parents[1]


def read_cluster_members(cluster_a3m: Path) -> list[str]:
    """Return list of headers in a cluster .a3m (without leading '>')."""
    headers = []
    with cluster_a3m.open() as f:
        for line in f:
            if line.startswith(">"):
                headers.append(line[1:].split()[0])
    return headers


def find_cluster_representative(members: list[str], tree_leaves: set[str]
                                 ) -> Optional[str]:
    """Pick the first cluster member that's a leaf of the tree."""
    for m in members:
        if m in tree_leaves:
            return m
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pair", default="5jytA_2qkeE")
    ap.add_argument("--scores-csv", default=str(ROOT / "docs" / "treek100_pilot_tmscores.csv"))
    ap.add_argument("--output", default=str(ROOT / "docs" / "figs" / "5jytA_2qkeE_phytree_treek100.png"))
    args = ap.parse_args()

    pair_dir = ROOT / "Pipeline" / "FoldPairs" / args.pair
    tree_path = pair_dir / "output_phytree" / "DeepMsa_tree.nwk"
    cluster_dir = pair_dir / "output_msa_cluster"

    if not tree_path.is_file():
        sys.exit(f"Tree not found: {tree_path}")
    if not cluster_dir.is_dir():
        sys.exit(f"Cluster dir not found: {cluster_dir}")

    # Load tree
    print(f"[plot] loading tree: {tree_path}")
    tree = Phylo.read(str(tree_path), "newick")
    tree_leaves = {l.name for l in tree.get_terminals() if l.name}
    print(f"[plot] tree has {len(tree_leaves)} leaves")

    # Map each cluster -> its representative leaf
    cluster_reps: dict[str, str] = {}  # cluster_tag -> rep leaf name
    cluster_a3ms = sorted(cluster_dir.glob("ShallowMsa_*.a3m"))
    for cluster_a3m in cluster_a3ms:
        cluster_tag = cluster_a3m.stem  # ShallowMsa_NNN
        members = read_cluster_members(cluster_a3m)
        rep = find_cluster_representative(members, tree_leaves)
        if rep is None:
            print(f"[plot] WARN {cluster_tag}: no member in tree; skipping")
            continue
        cluster_reps[cluster_tag] = rep
    print(f"[plot] {len(cluster_reps)} clusters have tree-leaf representatives")

    # Prune the tree to just the representative leaves
    keep_leaves = set(cluster_reps.values())
    leaves_to_drop = [l for l in tree.get_terminals() if l.name not in keep_leaves]
    print(f"[plot] pruning {len(leaves_to_drop)} non-representative leaves...")
    for leaf in leaves_to_drop:
        try:
            tree.prune(leaf)
        except Exception:
            pass

    # Reverse map: rep leaf name -> cluster_tag (for labelling)
    rep_to_cluster = {v: k for k, v in cluster_reps.items()}
    # Short labels: S000, S001, ...
    def short_label(cluster_tag: str) -> str:
        n = cluster_tag.replace("ShallowMsa_", "")
        return f"S{int(n):03d}"

    # Load TMdiff_centered scores
    scores = pd.read_csv(args.scores_csv)
    # Aggregate per (mode, cluster): mean over the 2 chains
    pivot = (scores.dropna(subset=["TMdiff_centered"])
                  .groupby(["mode", "cluster"])["TMdiff_centered"]
                  .mean()
                  .unstack("mode"))
    # Order columns: full first, top10 second, AF2 above AF3
    col_order = ["AF2_treek100_full", "AF2_treek100_top10",
                 "AF3_treek100_full", "AF3_treek100_top10"]
    col_order = [c for c in col_order if c in pivot.columns]
    pivot = pivot[col_order]
    print(f"[plot] heatmap shape: {pivot.shape}")
    print(f"[plot] columns: {list(pivot.columns)}")

    # ----- Plot tree + heatmap -----
    fig = plt.figure(figsize=(11, 14))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.3, 1.0], wspace=0.05)
    ax_tree = fig.add_subplot(gs[0, 0])
    ax_hm = fig.add_subplot(gs[0, 1])

    # Plot the tree using Bio.Phylo with custom labels
    # We need leaf order to match heatmap row order
    def get_leaf_order(tree) -> list[str]:
        return [l.name for l in tree.get_terminals()]
    leaf_order = get_leaf_order(tree)
    print(f"[plot] post-prune leaf count: {len(leaf_order)}")

    # Phylo.draw works in the standard ladderized order
    def relabel(node):
        if node.is_terminal() and node.name in rep_to_cluster:
            return short_label(rep_to_cluster[node.name])
        return ""
    Phylo.draw(tree, axes=ax_tree, do_show=False, label_func=relabel,
               show_confidence=False)
    ax_tree.set_xlabel("")
    ax_tree.set_ylabel("")
    ax_tree.set_title("K=100 tree-cut clusters\n(rep leaves)", fontsize=11)
    # Hide axis ticks for cleaner look
    ax_tree.set_yticks([])
    # Build heatmap data in leaf order
    cluster_order = []
    for leaf_name in leaf_order:
        if leaf_name in rep_to_cluster:
            cluster_order.append(rep_to_cluster[leaf_name])
    # Bio.Phylo draws leaves bottom-to-top by default; we need to align
    # the heatmap accordingly
    cluster_order_for_hm = list(reversed(cluster_order))
    hm_data = pivot.reindex(cluster_order_for_hm).to_numpy(dtype=float)

    # Diverging colormap red <-> white <-> blue
    cmap = LinearSegmentedColormap.from_list(
        "fold_pref",
        ["#1f4abf", "#85a4d3", "#ffffff", "#d99097", "#bf1f1f"],
        N=256,
    )
    vmax = max(0.1, np.nanmax(np.abs(hm_data)))
    im = ax_hm.imshow(hm_data, aspect="auto", cmap=cmap, vmin=-vmax, vmax=vmax,
                      interpolation="nearest")
    ax_hm.set_yticks(range(len(cluster_order_for_hm)))
    ax_hm.set_yticklabels([short_label(c) for c in cluster_order_for_hm], fontsize=6)
    ax_hm.yaxis.tick_right()
    ax_hm.set_xticks(range(len(col_order)))
    col_short = {
        "AF2_treek100_full": "AF2 full",
        "AF2_treek100_top10": "AF2 top10",
        "AF3_treek100_full": "AF3 full",
        "AF3_treek100_top10": "AF3 top10",
    }
    ax_hm.set_xticklabels([col_short.get(c, c) for c in col_order],
                          rotation=30, ha="right", fontsize=10)
    ax_hm.set_title("$\\Delta$TM$_{F_1 - F_2}$ centered\n(per cluster, mean over chains)",
                    fontsize=11)

    # Colorbar
    cbar = plt.colorbar(im, ax=ax_hm, fraction=0.04, pad=0.10)
    cbar.set_label("$\\Delta$TM centered\n(red = F1, blue = F2)", fontsize=9)

    fig.suptitle("KaiB (5jytA_2qkeE) K=100 + medoid-query pilot:\n"
                 "AF2/AF3 per-cluster fold preference, full vs top-10 MSA",
                 fontsize=12, y=0.995)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] wrote {out_path}")

    # Also save a CSV of cluster-ordered scores so caller can verify
    csv_out = out_path.with_suffix("").with_name(out_path.stem + "_data.csv")
    df_out = pd.DataFrame({"cluster": cluster_order_for_hm,
                           "leaf_label": [short_label(c) for c in cluster_order_for_hm]})
    for c in col_order:
        df_out[c] = [pivot.loc[cl, c] if cl in pivot.index else np.nan
                     for cl in cluster_order_for_hm]
    df_out.to_csv(csv_out, index=False)
    print(f"[plot] wrote {csv_out}")


if __name__ == "__main__":
    main()
