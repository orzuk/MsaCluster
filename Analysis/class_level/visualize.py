"""Class-level analysis figures.

Per-method scatter of (seq_identity, statistic) coloured by trigger_class,
plus a stacked-panel summary of Fisher / permutation p-values per analysis
mode (raw vs residualised vs identity-filtered). One figure per statistic.

Usage:
    python3 -m Analysis.class_level.visualize
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Stable per-class colours so figures across statistics are visually consistent.
CLASS_COLOURS = {
    "equilibrium_or_unknown": "#808080",  # gray
    "ligand":                 "#1f77b4",  # blue
    "mutation":               "#d62728",  # red
    "oligomerization":        "#2ca02c",  # green
    "protein_binding":        "#9467bd",  # purple
}


def scatter_seq_identity_by_class(
    stats: pd.DataFrame,
    statistic_col: str,
    method_col: str = "method",
    class_col: str = "trigger_class",
    identity_col: str = "seq_identity",
    output_path: str = "docs/class_level/scatter_seqid_by_class.png",
) -> None:
    """One subplot per method, x=seq_identity, y=statistic, coloured by class.
    Tests visually whether class colours segregate vertically at matched x,
    or whether classes just trace out the same x-y curve in different
    x-regions (the latter = the class effect was a seq_identity confound).
    """
    sub = stats.dropna(subset=[statistic_col, identity_col, class_col]).copy()
    methods = sorted(sub[method_col].unique())
    n = len(methods)
    ncol = 4
    nrow = (n + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.5 * ncol, 3.2 * nrow),
                             squeeze=False)
    for ax, m in zip(axes.flat, methods):
        gm = sub[sub[method_col] == m]
        for cls, group in gm.groupby(class_col):
            ax.scatter(group[identity_col], group[statistic_col],
                       c=CLASS_COLOURS.get(cls, "black"),
                       label=cls, alpha=0.65, s=22, edgecolor="none")
        # OLS reference line to visualise the seq_identity confound
        mask = gm[identity_col].notna() & gm[statistic_col].notna()
        if mask.sum() >= 5:
            x = gm.loc[mask, identity_col].astype(float).values
            y = gm.loc[mask, statistic_col].astype(float).values
            beta, alpha_ = np.polyfit(x, y, 1)
            xfit = np.linspace(x.min(), x.max(), 100)
            ax.plot(xfit, alpha_ + beta * xfit, ls="--", color="black",
                    alpha=0.5, lw=1.0)
        ax.set_title(m, fontsize=11)
        ax.set_xlabel("seq_identity")
        ax.set_ylabel(statistic_col, fontsize=9)
    # Hide unused panels
    for j in range(n, nrow * ncol):
        axes.flat[j].axis("off")
    # One legend on the figure
    handles = [plt.Line2D([0], [0], marker='o', color='w',
                          markerfacecolor=c, label=k, markersize=8)
               for k, c in CLASS_COLOURS.items() if k in sub[class_col].unique()]
    fig.legend(handles=handles, loc="lower center", ncol=5,
               frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle(f"Per-method scatter of {statistic_col} vs seq_identity "
                 f"(colour = trigger_class)", fontsize=13, y=1.01)
    fig.tight_layout()
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"[viz] wrote {output_path}")


def summary_pvalue_plot(
    summary_csv: str = "docs/class_level/summary.csv",
    output_path: str = "docs/class_level/summary_pvalues.png",
) -> None:
    """Bar plot summarising Fisher and permutation p-values across the four
    analysis modes (raw / residualised / identity-filtered / lmm-covariate),
    grouped by statistic. Visualises whether the class effect survives the
    seq_identity controls.
    """
    df = pd.read_csv(summary_csv)
    if "mode_prefix" not in df.columns:
        print(f"[viz] {summary_csv} doesn't have mode_prefix column; skipping")
        return
    stats_to_plot = ["stddev_centered", "stddev_residual",
                     "max_abs_centered", "max_abs_residual", "abs_beta"]
    modes_to_plot = ["raw", "resid", "highid"]
    mode_labels = {
        "raw": "Raw",
        "resid": "Residualised\non seq_identity",
        "highid": "Filtered to\nidentity≥0.95",
    }
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    width = 0.25
    x_positions = np.arange(len(stats_to_plot))
    for ax, p_col, title in zip(
        axes,
        ["combined_kruskal_fisher_p", "permutation_p"],
        ["Fisher combined p (per-method Kruskal-Wallis)",
         "Multi-method permutation p (sum of per-method H)"],
    ):
        for i, mode in enumerate(modes_to_plot):
            heights = []
            for s in stats_to_plot:
                row = df[(df.mode_prefix == mode) & (df.statistic == s)]
                if row.empty or pd.isna(row[p_col].values[0]):
                    heights.append(np.nan)
                else:
                    p = float(row[p_col].values[0])
                    # -log10 transform so small p shows tall bar
                    heights.append(-np.log10(max(p, 1e-7)))
            offset = (i - 1) * width
            ax.bar(x_positions + offset, heights, width=width,
                   label=mode_labels[mode],
                   color=plt.cm.viridis(0.2 + 0.3 * i),
                   edgecolor="black", linewidth=0.5)
        ax.axhline(-np.log10(0.05), ls="--", color="red", lw=1.0,
                   label="$\\alpha = 0.05$")
        ax.set_xticks(x_positions)
        ax.set_xticklabels(stats_to_plot, rotation=15)
        ax.set_ylabel("$-\\log_{10} p$")
        ax.set_title(title, fontsize=11)
        ax.legend(loc="upper right", fontsize=8, framealpha=0.85)
    fig.suptitle("Class-level effects: do they survive seq_identity controls?",
                 fontsize=13)
    fig.tight_layout()
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"[viz] wrote {output_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stats", default="docs/class_level/per_pair_per_method_stats.csv")
    ap.add_argument("--summary", default="docs/class_level/summary.csv")
    ap.add_argument("--outdir", default="docs/class_level")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    stats = pd.read_csv(args.stats)
    for stat_col in ["stddev_centered", "stddev_residual",
                     "max_abs_centered", "max_abs_residual", "abs_beta"]:
        if stat_col not in stats.columns:
            continue
        scatter_seq_identity_by_class(
            stats, statistic_col=stat_col,
            output_path=str(outdir / f"scatter_seqid_{stat_col}.png"),
        )

    summary_pvalue_plot(
        summary_csv=args.summary,
        output_path=str(outdir / "summary_pvalues.png"),
    )


if __name__ == "__main__":
    main()
