#!/usr/bin/env python3
"""PCA of the per-cluster fold-preference matrix for one family.

Each cluster is a point in the space of the seven analysed methods; the clade
picked out by the clade-shift scan is coloured separately. Methods are
z-scored first, so the decomposition is of the correlation matrix and no
method dominates through its units (TM-score, dddG, secondary-structure
fraction).

Usage:
    python3 scripts/plot_xcl1_pca.py [PAIR_ID] [OUT_PNG]
"""
import csv
import os
import sys

import numpy as np

SURVEY = "docs/fold_diversity_survey_corrected_9method.csv"
PLACEMENT = "docs/phylo_placement_medoidrep.csv"
METHODS = ["AF2", "AF3", "ESM", "Boltz2", "BioEmu", "DDG", "S4PRED"]
VALUE_COL = "TMdiff_centered"


def load_matrix(pair_id, value_col=VALUE_COL):
    """Return (clusters, methods, X) with X[i, j] finite for every entry."""
    by_cluster = {}
    with open(SURVEY) as fh:
        for row in csv.DictReader(fh):
            if row["pair_id"] != pair_id or row["method"] not in METHODS:
                continue
            raw = row.get(value_col, "")
            if raw in ("", "nan", "NaN"):
                continue
            by_cluster.setdefault(row["cluster"], {})[row["method"]] = float(raw)
    clusters = sorted(c for c, v in by_cluster.items() if len(v) == len(METHODS))
    X = np.array([[by_cluster[c][m] for m in METHODS] for c in clusters], float)
    return clusters, METHODS, X


def load_clade(pair_id, method="AF2"):
    with open(PLACEMENT) as fh:
        for row in csv.DictReader(fh):
            if row["pair_id"] == pair_id and row["method"] == method:
                members = (row.get("clade_best_members") or "").split(";")
                return [m for m in members if m], row.get("clade_p", "")
    return [], ""


def pca(X):
    """Correlation PCA: z-score the columns, then SVD."""
    Z = (X - X.mean(0)) / np.where(X.std(0) > 0, X.std(0), 1.0)
    U, S, Vt = np.linalg.svd(Z, full_matrices=False)
    scores = U * S
    var = S ** 2 / (S ** 2).sum()
    return scores, Vt, var


def main(pair_id="2n54B_2hdmA", out_png="poster/figs/xcl1_pca.png"):
    clusters, methods, X = load_matrix(pair_id)
    clade, clade_p = load_clade(pair_id)
    scores, loadings, var = pca(X)
    in_clade = np.array([c in clade for c in clusters])

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(3.6, 3.0), dpi=300)
    ax.axhline(0, color="0.85", lw=0.6, zorder=0)
    ax.axvline(0, color="0.85", lw=0.6, zorder=0)
    ax.scatter(scores[~in_clade, 0], scores[~in_clade, 1], s=16, c="#1F77B4",
               alpha=0.65, linewidths=0, label=f"Other clusters ({(~in_clade).sum()})")
    ax.scatter(scores[in_clade, 0], scores[in_clade, 1], s=34, c="#D62728",
               edgecolors="white", linewidths=0.5, zorder=3,
               label=f"Shifted clade ({in_clade.sum()})")
    ax.set_xlabel(f"PC1 ({100 * var[0]:.0f}% of variance)", fontsize=7)
    ax.set_ylabel(f"PC2 ({100 * var[1]:.0f}%)", fontsize=7)
    ax.tick_params(labelsize=6)
    ax.legend(fontsize=6, frameon=False, loc="best")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.tight_layout(pad=0.3)
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    fig.savefig(out_png)
    print(f"{pair_id}: {len(clusters)} clusters x {len(methods)} methods; "
          f"clade {in_clade.sum()} clusters (p = {clade_p}); "
          f"PC1 {100 * var[0]:.1f}%, PC2 {100 * var[1]:.1f}%")
    print("wrote", out_png)


if __name__ == "__main__":
    main(*sys.argv[1:])
