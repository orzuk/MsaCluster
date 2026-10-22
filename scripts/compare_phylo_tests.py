#!/usr/bin/env python3
"""Compare Moran's I vs the clade-shift scan in a phylo_placement CSV.

Input: a docs/phylo_placement[_corrected].csv produced with
    scripts/phylo_placement.py --tests morans,clade_shift,clade_shift_combined
(columns morans_p, morans_p_bh, clade_p, clade_p_bh, clade_best_* ...).

Per method (and for COMBINED) it reports how many families (pairs) are
significant at raw p<0.05 and BH<0.05 by Moran's I vs by the clade-shift scan,
their overlap (both / only Moran / only clade) and the Spearman correlation of
-log10 p. For the COMBINED row (one per pair) the Moran comparator is the
pair's best Moran's I over methods (min raw p, min BH) -- "significant by
Moran's in any method".

Outputs
  docs/phylo_tests_comparison<suffix>.csv   per-method summary table
  docs/phylo_tests_comparison<suffix>.txt   short text summary
  docs/figures/global/phylo_tests_morans_vs_clade<suffix>.png
      small multiples: -log10 p Moran (x) vs clade-shift (y), per method

Usage
  python3 scripts/compare_phylo_tests.py                       # raw labels
  python3 scripts/compare_phylo_tests.py --labels corrected    # suffix _corrected
  python3 scripts/compare_phylo_tests.py --pair 2n54B_2hdmA    # detail for one pair
"""

from __future__ import annotations

import argparse
import os
import sys

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from config import TABLES_RES  # noqa: E402
from utils.method_labels import order_methods  # noqa: E402

ALPHA = 0.05
P_FLOOR = 1e-6          # for -log10 of p == 0 (never happens with +1 perms)


def _nlp(p) -> np.ndarray:
    p = pd.to_numeric(p, errors="coerce").to_numpy(dtype=float)
    return -np.log10(np.clip(p, P_FLOOR, 1.0))


def _spearman(x: np.ndarray, y: np.ndarray):
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return np.nan, np.nan
    try:
        from scipy.stats import spearmanr
        r = spearmanr(x[ok], y[ok])
        return float(r.statistic if hasattr(r, "statistic") else r[0]), float(r.pvalue)
    except ImportError:
        return float(pd.Series(x[ok]).corr(pd.Series(y[ok]), method="spearman")), np.nan


def paired_table(df: pd.DataFrame) -> pd.DataFrame:
    """One row per (pair, method) with Moran and clade p's side by side,
    including COMBINED rows (Moran = best over the pair's methods)."""
    per = df[df["method"] != "COMBINED"][
        ["pair_id", "method", "morans_p", "morans_p_bh", "clade_p", "clade_p_bh"]
    ].copy()
    comb = df[df["method"] == "COMBINED"][
        ["pair_id", "method", "clade_p", "clade_p_bh"]].copy()
    if not comb.empty:
        best = (per.groupby("pair_id")[["morans_p", "morans_p_bh"]].min()
                .reset_index())
        comb = comb.merge(best, on="pair_id", how="left")
    out = pd.concat([per, comb], ignore_index=True)
    for c in ("morans_p", "morans_p_bh", "clade_p", "clade_p_bh"):
        out[c] = pd.to_numeric(out[c], errors="coerce")
    return out


def summarize(pt: pd.DataFrame) -> pd.DataFrame:
    rows = []
    methods = order_methods([m for m in pt["method"].unique() if m != "COMBINED"])
    if (pt["method"] == "COMBINED").any():
        methods.append("COMBINED")
    for m in methods:
        sub = pt[(pt["method"] == m) & pt["morans_p"].notna() & pt["clade_p"].notna()]
        r = {"method": m, "n_tested": len(sub)}
        for tag, mc, cc in (("raw", "morans_p", "clade_p"),
                            ("bh", "morans_p_bh", "clade_p_bh")):
            sm = sub[mc] < ALPHA
            sc = sub[cc] < ALPHA
            r.update({f"n_morans_{tag}": int(sm.sum()),
                      f"n_clade_{tag}": int(sc.sum()),
                      f"n_both_{tag}": int((sm & sc).sum()),
                      f"n_only_morans_{tag}": int((sm & ~sc).sum()),
                      f"n_only_clade_{tag}": int((~sm & sc).sum())})
        rho, rp = _spearman(_nlp(sub["morans_p"]), _nlp(sub["clade_p"]))
        r["spearman_rho_nlog10p"] = round(rho, 3) if np.isfinite(rho) else np.nan
        r["spearman_p"] = rp
        if len(sub):                      # skip methods absent from this CSV
            rows.append(r)
    return pd.DataFrame(rows)


def summary_text(summ: pd.DataFrame, src: str) -> str:
    lines = [f"Moran's I vs clade-shift scan ({src}); alpha={ALPHA}",
             "COMBINED: Moran comparator = best Moran's I over the pair's methods.",
             "",
             f"{'method':10} {'n':>4} | {'raw: Mor clade both onlyM onlyC':>32} | "
             f"{'BH: Mor clade both onlyM onlyC':>31} | rho(-log p)"]
    for _, r in summ.iterrows():
        raw = (f"{r.n_morans_raw:>8} {r.n_clade_raw:>5} {r.n_both_raw:>4} "
               f"{r.n_only_morans_raw:>5} {r.n_only_clade_raw:>5}")
        bh = (f"{r.n_morans_bh:>7} {r.n_clade_bh:>5} {r.n_both_bh:>4} "
              f"{r.n_only_morans_bh:>5} {r.n_only_clade_bh:>5}")
        rho = f"{r.spearman_rho_nlog10p:.2f}" if np.isfinite(r.spearman_rho_nlog10p) else "  -"
        lines.append(f"{r.method:10} {r.n_tested:>4} | {raw:>32} | {bh:>31} | {rho}")
    return "\n".join(lines)


def scatter(pt: pd.DataFrame, summ: pd.DataFrame, out_png: str, title: str) -> None:
    import matplotlib.pyplot as plt
    methods = list(summ["method"])
    if not methods:
        return
    ncol = min(5, len(methods))
    nrow = int(np.ceil(len(methods) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.9 * ncol, 2.9 * nrow),
                             squeeze=False, sharex=True, sharey=True)
    fig.patch.set_facecolor("#fcfcfb")
    x_all, y_all = _nlp(pt["morans_p"]), _nlp(pt["clade_p"])
    lim = max(1.5, np.nanmax(np.r_[x_all, y_all, 1.0]) * 1.05)
    thr = -np.log10(ALPHA)
    for ax, m in zip(axes.ravel(), methods):
        sub = pt[pt["method"] == m]
        x, y = _nlp(sub["morans_p"]), _nlp(sub["clade_p"])
        color = "#eb6834" if m == "COMBINED" else "#2a78d6"
        ax.set_facecolor("#fcfcfb")
        ax.plot([0, lim], [0, lim], color="#c3c2b7", lw=1, zorder=1)
        ax.axvline(thr, color="#c3c2b7", lw=1, ls="--", zorder=1)
        ax.axhline(thr, color="#c3c2b7", lw=1, ls="--", zorder=1)
        ax.scatter(x, y, s=22, color=color, alpha=0.8, edgecolor="#fcfcfb",
                   linewidth=0.8, zorder=2)
        r = summ[summ["method"] == m].iloc[0]
        lab = "COMBINED (x: best Moran)" if m == "COMBINED" else m
        ax.set_title(lab, fontsize=9, color="#0b0b0b")
        ax.text(0.97, 0.03,
                f"raw p<.05  both {r.n_both_raw} | M {r.n_only_morans_raw} | "
                f"C {r.n_only_clade_raw}\nrho={r.spearman_rho_nlog10p:.2f}",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=7,
                color="#52514e")
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color("#c3c2b7")
        ax.tick_params(colors="#52514e", labelsize=7)
    for ax in axes.ravel()[len(methods):]:
        ax.axis("off")
    for ax in axes[-1]:
        ax.set_xlabel("-log10 p  Moran's I", fontsize=8, color="#52514e")
    for ax in axes[:, 0]:
        ax.set_ylabel("-log10 p  clade shift", fontsize=8, color="#52514e")
    fig.suptitle(title, fontsize=10, color="#0b0b0b")
    fig.tight_layout()
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    fig.savefig(out_png, dpi=150, facecolor=fig.get_facecolor())
    plt.close(fig)


def show_pair(df: pd.DataFrame, pair_id: str) -> None:
    sub = df[df["pair_id"] == pair_id].copy()
    if sub.empty:
        print(f"[error] pair {pair_id} not in CSV")
        sys.exit(2)
    order = order_methods([m for m in sub["method"].unique() if m != "COMBINED"])
    order += ["COMBINED"] if (sub["method"] == "COMBINED").any() else []
    sub = sub.set_index("method").loc[order].reset_index()
    cols = ["method", "morans_n", "morans_I", "morans_I_norm", "morans_p",
            "morans_p_bh", "clade_n_leaves", "clade_T", "clade_p", "clade_p_bh",
            "clade_best_size", "clade_best_direction", "clade_best_effect",
            "clade_best_rbc", "clade_note"]
    cols = [c for c in cols if c in sub.columns]
    with pd.option_context("display.width", 250, "display.max_columns", 40):
        print(f"=== {pair_id} ===")
        print(sub[cols].to_string(index=False))
    print("\nBest clade per method (members as cluster tags):")
    for _, r in sub.iterrows():
        mem = r.get("clade_best_members")
        if not isinstance(mem, str) or not mem:
            continue
        extra = r.get("clade_best_method_z")
        extra = f"  [{extra}]" if isinstance(extra, str) and extra else ""
        print(f"  {r['method']:9} size={int(r['clade_best_size'])} "
              f"dir={int(r['clade_best_direction']):+d}: {mem}{extra}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--labels", choices=["centered", "corrected"], default="centered",
                    help="picks the default input and output suffix")
    ap.add_argument("--input", default=None,
                    help="phylo_placement CSV (default docs/phylo_placement"
                         "[_corrected].csv by --labels)")
    ap.add_argument("--suffix", default=None,
                    help="output suffix (default '' for centered, '_corrected')")
    ap.add_argument("--out-dir", default=TABLES_RES,
                    help="where CSV/TXT go (figure under <out-dir>/figures/global)")
    ap.add_argument("--pair", default=None,
                    help="print detailed per-method results for one pair and exit")
    args = ap.parse_args()

    corrected = args.labels == "corrected"
    in_csv = args.input or os.path.join(
        TABLES_RES, "phylo_placement_corrected.csv" if corrected else "phylo_placement.csv")
    suffix = args.suffix if args.suffix is not None else ("_corrected" if corrected else "")
    df = pd.read_csv(in_csv)
    missing = [c for c in ("morans_p", "clade_p") if c not in df.columns]
    if missing:
        sys.exit(f"[error] {in_csv} lacks {missing}; rerun scripts/phylo_placement.py "
                 "with --tests morans,clade_shift,clade_shift_combined")
    if "clade_p_bh" not in df.columns:
        df["clade_p_bh"] = np.nan

    if args.pair:
        show_pair(df, args.pair)
        return

    pt = paired_table(df)
    summ = summarize(pt)
    out_csv = os.path.join(args.out_dir, f"phylo_tests_comparison{suffix}.csv")
    out_txt = os.path.join(args.out_dir, f"phylo_tests_comparison{suffix}.txt")
    os.makedirs(args.out_dir, exist_ok=True)
    out_png = os.path.join(args.out_dir, "figures", "global",
                           f"phylo_tests_morans_vs_clade{suffix}.png")
    summ.to_csv(out_csv, index=False)
    txt = summary_text(summ, os.path.basename(in_csv))
    with open(out_txt, "w") as fh:
        fh.write(txt + "\n")
    print(txt)
    scatter(pt, summ, out_png,
            f"Moran's I vs clade-shift scan ({os.path.basename(in_csv)})")
    print(f"\nWrote {out_csv}\n      {out_txt}\n      {out_png}")


if __name__ == "__main__":
    main()
