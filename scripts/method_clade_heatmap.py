#!/usr/bin/env python3
"""Methods x pairs heatmap of clade-signal strength, to see which method gives
the strongest phylogenetic (clade-specific) signal per family.

Reads docs/phylo_placement_corrected.csv (falls back to raw). Uses morans_I_norm
(ceiling-normalized effect size; the fair cross-method/cross-pair metric) if
present, else morans_z. Renders a heatmap (pairs x methods), prints the per-pair
argmax method and the tally ("DDG strongest in N/85 pairs"), and writes a CSV.

Usage:
    python3 scripts/method_clade_heatmap.py            # corrected, I_norm if present
    python3 scripts/method_clade_heatmap.py --metric morans_z --labels raw
"""
import argparse
import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCS = os.path.join(ROOT, "docs")
METHOD_ORDER = ["AF2", "AF3", "ESM", "Boltz2", "DDG", "MSAT", "CCMpred", "S4PRED"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", choices=["raw", "corrected"], default="corrected")
    ap.add_argument("--metric", default=None,
                    help="morans_I_norm (default if present) or morans_z")
    args = ap.parse_args()

    fn = ("phylo_placement_corrected.csv" if args.labels == "corrected"
          else "phylo_placement.csv")
    d = pd.read_csv(os.path.join(DOCS, fn))
    metric = args.metric or ("morans_I_norm" if "morans_I_norm" in d.columns
                             else "morans_z")
    if metric not in d.columns:
        raise SystemExit(f"{metric} not in {fn} — re-run phylo_placement first")
    print(f"[heatmap] metric={metric}  labels={args.labels}")

    wide = d.pivot_table(index="pair_id", columns="method", values=metric, aggfunc="first")
    cols = [m for m in METHOD_ORDER if m in wide.columns]
    wide = wide[cols]

    # Per-pair winner (argmax over methods) + tally
    win = wide.idxmax(axis=1)
    tally = win.value_counts()
    print(f"\n=== Strongest method per pair ({metric}) — tally over {win.notna().sum()} pairs ===")
    for m in cols:
        print(f"  {m:9} strongest in {int(tally.get(m, 0)):>3} pairs")
    print(f"\n  median {metric} per method:")
    for m in cols:
        print(f"  {m:9} {wide[m].median():+.3f}")

    out_csv = os.path.join(DOCS, f"method_clade_{metric}_{args.labels}.csv")
    wide.to_csv(out_csv)
    print(f"\n[ok] wrote {out_csv}")

    # Heatmap (sorted by DDG strength so the pattern is visible)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        sort_col = "DDG" if "DDG" in wide.columns else cols[0]
        w = wide.sort_values(sort_col, ascending=False)
        M = w.to_numpy(dtype=float)
        vlim = np.nanpercentile(np.abs(M), 98) or 1.0
        fig, ax = plt.subplots(figsize=(6, max(4, 0.16 * len(w))))
        im = ax.imshow(M, cmap="RdBu_r", vmin=-vlim, vmax=vlim, aspect="auto")
        ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=45, ha="right")
        ax.set_yticks(range(len(w))); ax.set_yticklabels(w.index, fontsize=4)
        ax.set_title(f"Clade signal per method x pair ({metric}, {args.labels})\n"
                     f"rows sorted by DDG", fontsize=10)
        plt.colorbar(im, ax=ax, shrink=0.6, label=metric)
        plt.tight_layout()
        out_png = os.path.join(DOCS, "figs", f"method_clade_{metric}_{args.labels}.png")
        os.makedirs(os.path.dirname(out_png), exist_ok=True)
        plt.savefig(out_png, dpi=200); plt.close(fig)
        print(f"[ok] wrote {out_png}")
    except Exception as e:
        print(f"[warn] heatmap render failed: {e}")


if __name__ == "__main__":
    main()
