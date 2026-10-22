#!/usr/bin/env python3
"""SI figure: observed fold-preference I_norm inside its own random-trait null.

The per-family random-trait test (scripts/trait_null_phylo.py) ranks each
method's Moran's I_norm against K random projections of the clusters'
amino-acid composition on the same tree and the same leaves. Seven families
pass it after multiplicity correction, and all seven sit at the smallest
attainable family p (K+2 over K+1 draws), so the p-values cannot rank them.
This figure shows what the p-values cannot: where the observed value falls in
the null it is tested against, family by family.

One panel per family. Grey: the K null I_norm values, as a histogram. Coloured
line: the observed I_norm of the family's best method. The panel title gives
the family, the method and the number of clusters; the annotation gives the
per-method p and the fraction of the null the observed value exceeds.

Inputs
  docs/trait_null_draws_7fam.npz   null draws, key "<pair>|<labels>|<method>"
                                   (scripts/trait_null_phylo.py --dump-null)
  docs/trait_null_phylo.csv        observed I_norm and p_trait_inorm per method
  docs/trait_null_phylo_family.csv family-level p, q and best method
  docs/pair_protein_names.csv      display names

Usage (repo root)
  python3 scripts/plot_trait_null_dist.py
  python3 scripts/plot_trait_null_dist.py --labels corrected --out docs/figures/si
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
sys.path.insert(0, REPO)

try:
    from utils.method_labels import method_color, method_label
except Exception:  # keep the script runnable without the registry
    def method_label(m):
        return m

    def method_color(m):
        return "#1f77b4"


def load_names(path: str) -> dict:
    if not os.path.isfile(path):
        return {}
    df = pd.read_csv(path)
    col = "name" if "name" in df.columns else df.columns[1]
    return dict(zip(df["pair_id"], df[col]))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--draws", default=os.path.join(REPO, "docs",
                                                    "trait_null_draws_7fam.npz"))
    ap.add_argument("--per-method", default=os.path.join(REPO, "docs",
                                                         "trait_null_phylo.csv"))
    ap.add_argument("--family", default=os.path.join(REPO, "docs",
                                                     "trait_null_phylo_family.csv"))
    ap.add_argument("--names-csv", default=os.path.join(REPO, "docs",
                                                        "pair_protein_names.csv"))
    ap.add_argument("--labels", default="corrected")
    ap.add_argument("--out", default=os.path.join(REPO, "docs", "figures", "si"))
    ap.add_argument("--stem", default="trait_null_distributions")
    args = ap.parse_args(argv)

    if not os.path.isfile(args.draws):
        sys.exit(f"[error] missing {args.draws}; run trait_null_phylo.py --dump-null")
    draws = np.load(args.draws)
    per = pd.read_csv(args.per_method)
    fam = pd.read_csv(args.family)
    names = load_names(args.names_csv)

    fam = fam[fam["labels"] == args.labels].copy()
    fam = fam[pd.to_numeric(fam["q_inorm"], errors="coerce") <= 0.05]
    fam = fam.sort_values("best_I_norm", ascending=False)
    if fam.empty:
        sys.exit("[error] no families pass the null at q <= 0.05")

    n = len(fam)
    ncol = 3 if n > 4 else 2
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.6 * ncol, 2.9 * nrow))
    axes = np.atleast_1d(axes).ravel()

    for ax, (_, row) in zip(axes, fam.iterrows()):
        pair, method = row["pair_id"], row["best_method_inorm"]
        key = f"{pair}|{args.labels}|{method}"
        if key not in draws:
            ax.set_visible(False)
            print(f"[warn] no draws for {key}", flush=True)
            continue
        null = np.asarray(draws[key], dtype=float)
        null = null[np.isfinite(null)]
        obs = float(row["best_I_norm"])

        pm = per[(per["pair_id"] == pair) & (per["labels"] == args.labels)
                 & (per["method"] == method)]
        p_trait = float(pm["p_trait_inorm"].iloc[0]) if len(pm) else np.nan
        n_clusters = int(pm["n_clusters"].iloc[0]) if len(pm) else int(row["n_clusters"])

        col = method_color(method)
        q95 = float(np.quantile(null, 0.95))
        ax.hist(null, bins=45, color="0.72", edgecolor="none")
        ax.axvline(obs, color=col, lw=2.2, zorder=5)
        ax.axvline(q95, color="0.35", lw=1.0, ls=":")

        # Leave room on both sides so the observed line never sits on the
        # frame, and put the label on whichever side the line is not.
        lo, hi = min(float(null.min()), obs), max(float(null.max()), obs)
        pad = 0.10 * (hi - lo if hi > lo else 1.0)
        ax.set_xlim(lo - pad, hi + pad)
        right = obs > 0.5 * (lo + hi)
        ax.annotate(f"observed {obs:.2f}\n$p$ = {p_trait:.4f}\n"
                    f"null 95th pct {q95:.2f}, max {float(null.max()):.2f}",
                    xy=(0.03 if right else 0.97, 0.95), xycoords="axes fraction",
                    ha="left" if right else "right", va="top",
                    fontsize=7.5, color=col)
        title = names.get(pair, pair)
        ax.set_title(f"{title}\n{method_label(method)}, {n_clusters} clusters",
                     fontsize=8.5)
        ax.set_xlabel(r"$I_{\mathrm{norm}}$", fontsize=8)
        ax.set_ylabel("random traits", fontsize=8)
        ax.tick_params(labelsize=7)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)

    for ax in axes[n:]:
        ax.set_visible(False)

    k = int(np.median([len(np.asarray(draws[k_])) for k_ in draws.files]))
    fig.suptitle(f"Observed fold-preference $I_{{\\mathrm{{norm}}}}$ against the "
                 f"family's own null of {k:,} random composition projections "
                 f"({args.labels} preferences)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))

    os.makedirs(args.out, exist_ok=True)
    for ext in ("png", "pdf"):
        path = os.path.join(args.out, f"{args.stem}.{ext}")
        fig.savefig(path, dpi=300 if ext == "png" else None, bbox_inches="tight")
        print(f"[plot] wrote {path}", flush=True)
    plt.close(fig)
    return 0


if __name__ == "__main__":
    sys.exit(main())
