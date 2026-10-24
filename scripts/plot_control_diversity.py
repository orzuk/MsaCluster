#!/usr/bin/env python3
"""SI figure: structural spread of per-cluster predictions, switchers vs controls.

Every target's MSA clusters are folded one representative sequence at a time by
ESMFold2, and TM-align is run over all pairs of the resulting structures. A
single-fold protein should give one structure in every cluster; a family whose
clusters prefer different folds should not. The figure shows the distribution of
the smallest pairwise TM-score per target, which is where a second fold appears,
and then shows that the separation is not a by-product of how far apart the
folded sequences are.

Inputs
  docs/esmfold2_cluster_diversity.csv   one row per target: group, n_clusters,
                                        mean/median/min/p10 pairwise TM and the
                                        mean/min identity of its representatives
                                        (scripts/esmfold2_diversity.py)

Usage (repo root)
  python3 scripts/plot_control_diversity.py
  python3 scripts/plot_control_diversity.py --out "/mnt/g/.../figs" --stem control_diversity
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ORDER = ["FoldSwitch", "MatchedSingleFold", "PorterIDP", "RandomSingleFold"]
LABEL = {"FoldSwitch": "fold-switching pairs",
         "MatchedSingleFold": "Pfam-matched single-fold",
         "PorterIDP": "disordered (DisProt)",
         "RandomSingleFold": "random single-fold"}
COLOR = {"FoldSwitch": "#c1272d", "MatchedSingleFold": "#4878a8",
         "PorterIDP": "#6a9a52", "RandomSingleFold": "#8c6bb1"}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--csv", default=os.path.join(REPO, "docs",
                                                  "esmfold2_cluster_diversity.csv"))
    ap.add_argument("--out", default=os.path.join(REPO, "docs", "figures", "si"))
    ap.add_argument("--stem", default="control_diversity")
    args = ap.parse_args(argv)

    if not os.path.isfile(args.csv):
        sys.exit(f"[error] missing {args.csv}")
    df = pd.read_csv(args.csv)
    df = df[df["group"].isin(ORDER)]

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.7))

    # Panel A: distribution of the per-target minimum pairwise TM-score.
    ax = axes[0]
    for i, g in enumerate(ORDER):
        v = df.loc[df["group"] == g, "min_tm"].dropna().values
        if not len(v):
            continue
        x = np.sort(v)
        ax.step(x, np.arange(1, len(x) + 1) / len(x), where="post",
                color=COLOR[g], lw=2.0 if g == "FoldSwitch" else 1.4,
                label=f"{LABEL[g]} ($n$ = {len(v)})")
        ax.plot([np.median(v)], [0.5], marker="o", ms=4.5, color=COLOR[g])
    ax.set_xlabel("smallest pairwise TM-score between a target's clusters", fontsize=8.5)
    ax.set_ylabel("cumulative fraction of targets", fontsize=8.5)
    ax.set_title("A   Clusters of single-fold proteins stay together",
                 fontsize=9, loc="left")
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    ax.tick_params(labelsize=7.5)

    # Panel B: the spread does not follow how diverged the folded sequences are.
    ax = axes[1]
    for g in ORDER:
        s = df[df["group"] == g].dropna(subset=["mean_rep_id", "min_tm"])
        ax.scatter(s["mean_rep_id"], s["min_tm"], s=16,
                   facecolor=COLOR[g], edgecolor="none",
                   alpha=0.85 if g == "FoldSwitch" else 0.55)
    sw = df[df["group"] == "FoldSwitch"].dropna(subset=["mean_rep_id", "min_tm"])
    ct = df[df["group"] != "FoldSwitch"].dropna(subset=["mean_rep_id", "min_tm"])
    for s, col in ((sw, COLOR["FoldSwitch"]), (ct, "0.35")):
        if len(s) > 2:
            b = np.polyfit(s["mean_rep_id"], s["min_tm"], 1)
            xs = np.linspace(s["mean_rep_id"].min(), s["mean_rep_id"].max(), 20)
            ax.plot(xs, np.polyval(b, xs), color=col, lw=1.3, ls="--")
    ax.set_xlabel("mean identity between the folded cluster representatives", fontsize=8.5)
    ax.set_ylabel("smallest pairwise TM-score", fontsize=8.5)
    ax.set_title("B   Spread does not follow sequence divergence",
                 fontsize=9, loc="left")
    ax.tick_params(labelsize=7.5)

    for ax in axes:
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    fig.tight_layout()

    os.makedirs(args.out, exist_ok=True)
    for ext in ("png", "pdf"):
        path = os.path.join(args.out, f"{args.stem}.{ext}")
        fig.savefig(path, dpi=300 if ext == "png" else None, bbox_inches="tight")
        print(f"[plot] wrote {path}", flush=True)
    plt.close(fig)
    return 0


if __name__ == "__main__":
    sys.exit(main())
