#!/usr/bin/env python3
"""7x7 Spearman correlation of per-cluster centered fold preference across methods.

Reads docs/fold_diversity_survey.csv, pools all (pair_id, cluster) rows
EXCLUDING the DeepMsa rows, builds a (pair_id, cluster) x method matrix of
TMdiff_centered values, computes pairwise Spearman rank-correlation, and
saves a labeled 7x7 heatmap to docs/figs/cross_method_correlation.png.

Why this exists
---------------
The per-pair tree-heatmaps (utils/protein_plot_utils.py) show each method's
per-cluster signal. Cross-pair, we want to know: which methods systematically
agree (likely shared bias) and which are independent (orthogonal evidence).
AF2 and AF3 are expected to track tightly (shared training-data bias); MSAT
and CCMpred should be more orthogonal; ESM stands alone (single-sequence
transformer); DDG is energy-based and likely orthogonal to the rest.

This figure quantifies those expectations across all 93 pairs.
"""

from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.method_labels import METHOD_ORDER, disp  # noqa: E402

DOCS = "docs"
OUT = os.path.join(DOCS, "figs", "cross_method_correlation.png")


def load_long() -> pd.DataFrame:
    csv = os.path.join(DOCS, "fold_diversity_survey.csv")
    df = pd.read_csv(csv)
    # Drop the un-clustered Deep MSA rows; correlations are per-cluster.
    df = df[~df["cluster"].astype(str).str.lower().str.startswith("deep")]
    df["TMdiff_centered"] = pd.to_numeric(df["TMdiff_centered"], errors="coerce")
    return df


def pivot_wide(df: pd.DataFrame) -> pd.DataFrame:
    """(pair_id, cluster) -> method matrix of TMdiff_centered."""
    wide = df.pivot_table(
        index=["pair_id", "cluster"], columns="method",
        values="TMdiff_centered", aggfunc="first",
    )
    cols = [m for m in METHOD_ORDER if m in wide.columns]
    return wide[cols]


def spearman_matrix(wide: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (corr 7x7, pairwise n 7x7) using all available (pair,cluster) rows."""
    from scipy.stats import spearmanr
    cols = list(wide.columns)
    corr = pd.DataFrame(index=cols, columns=cols, dtype=float)
    counts = pd.DataFrame(index=cols, columns=cols, dtype=int)
    arr = wide.to_numpy(dtype=float)
    for i, a in enumerate(cols):
        for j, b in enumerate(cols):
            x = arr[:, i]
            y = arr[:, j]
            mask = np.isfinite(x) & np.isfinite(y)
            n = int(mask.sum())
            counts.loc[a, b] = n
            if n < 5:
                corr.loc[a, b] = float("nan")
                continue
            if i == j:
                corr.loc[a, b] = 1.0
                continue
            rho, _ = spearmanr(x[mask], y[mask])
            corr.loc[a, b] = float(rho)
    return corr, counts


def render(corr: pd.DataFrame, counts: pd.DataFrame, out_path: str) -> None:
    cols = list(corr.columns)
    labels = [disp(c) for c in cols]
    M = corr.to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(7.0, 6.0))
    im = ax.imshow(M, cmap="RdBu_r", vmin=-1.0, vmax=1.0,
                   interpolation="nearest", aspect="equal")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_yticks(range(len(cols)))
    ax.set_yticklabels(labels)
    # cell labels
    for i in range(len(cols)):
        for j in range(len(cols)):
            v = M[i, j]
            n = int(counts.iloc[i, j])
            if not np.isfinite(v):
                txt = "—"
                color = "#888"
            else:
                txt = f"{v:.2f}"
                color = "white" if abs(v) > 0.5 else "black"
            ax.text(j, i, txt, ha="center", va="center",
                    fontsize=8, color=color)
    cbar = plt.colorbar(im, ax=ax, shrink=0.85)
    cbar.set_label("Spearman ρ of TMdiff_centered\n(across all per-pair × cluster rows)",
                   fontsize=9)
    ax.set_title("Cross-method correlation of per-cluster centered fold preference\n"
                 f"({M.shape[0]} methods, pooled across 85 pairs)",
                 fontsize=10, pad=8)
    plt.tight_layout()
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    plt.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"[ok] wrote {out_path}")


def main() -> None:
    df = load_long()
    print(f"Loaded {len(df)} (pair, cluster, method) rows (Deep MSA excluded)")
    wide = pivot_wide(df)
    print(f"Wide matrix: {wide.shape[0]} (pair, cluster) rows x {wide.shape[1]} methods")
    print(f"Non-NaN per method:\n{wide.notna().sum()}")
    corr, counts = spearman_matrix(wide)
    print("\nSpearman correlation matrix:")
    print(corr.round(2).to_string())
    # Persist the numbers (not just the PNG) so the matrix is readable/diffable.
    csv_out = OUT.rsplit(".", 1)[0] + ".csv"
    n_out = OUT.rsplit(".", 1)[0] + "_n.csv"
    os.makedirs(os.path.dirname(csv_out), exist_ok=True)
    corr.round(4).to_csv(csv_out)
    counts.to_csv(n_out)
    print(f"[ok] wrote {csv_out} (correlation matrix) and {n_out} (pairwise n)")
    render(corr, counts, OUT)


if __name__ == "__main__":
    main()
