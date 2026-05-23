#!/usr/bin/env python3
"""Compare RAW vs CORRECTED stratified-by-trigger BH FDR analysis.

Why this exists
---------------
The headline stratified-by-trigger analysis in
scripts/stratified_concordance_by_trigger.py uses the POST-correction
empirical permutation p-values. That tells us which pairs survive
sequence-divergence correction within their trigger class.

A natural complementary question: which pairs would have been called
significant by the RAW (pre-correction) analysis with the same
within-class BH? The difference between the two lists is the set of
pairs the correction filtered out — a direct demonstration of why the
correction matters.

This script computes both and prints the comparison.

Reads:
   docs/fold_diversity_concordance_perm.csv              (raw)
   docs/fold_diversity_concordance_corrected_perm.csv    (corrected)
   docs/triggers_from_pdb.csv

Output: stdout summary (no CSV — readers can `tee` if they want a log).
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

_THIS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _THIS)
sys.path.insert(0, os.path.dirname(_THIS))

from permutation_null_concordance import bh_adjust  # noqa: E402

DOCS = "docs"
ALPHA = 0.05


def load_or_die(path):
    if not os.path.isfile(path):
        print(f"[error] missing CSV: {path}", file=sys.stderr)
        sys.exit(2)
    return pd.read_csv(path)


def stratified_bh(df: pd.DataFrame, pcol: str, label: str) -> dict:
    """Per-class BH FDR. Returns dict {class -> list of (pair_id, raw_p, bh_p)} for hits."""
    print(f"\n=== {label}: per-class BH FDR <= {ALPHA} ===")
    print(f"{'class':<25} {'n':>4} {'rank-1 hit':<22} {'p_raw':>8} {'p_BH':>8}")
    print("-" * 75)
    hits_by_class: dict = {}
    for cls, sub in df.groupby("trigger_class"):
        sub = sub.copy()
        sub["p_bh"] = bh_adjust(sub[pcol].to_numpy(dtype=float))
        hits = sub[sub["p_bh"] <= ALPHA].sort_values("p_bh")
        hits_by_class[cls] = []
        if hits.empty:
            print(f"{cls:<25} {len(sub):>4d}  (no within-class hits)")
            continue
        for _, r in hits.iterrows():
            print(f"{cls:<25} {len(sub):>4d}  {r['pair_id']:<22} "
                  f"{r[pcol]:>8.4f} {r['p_bh']:>8.4f}")
            hits_by_class[cls].append(
                (r["pair_id"], float(r[pcol]), float(r["p_bh"]))
            )
    return hits_by_class


def main() -> None:
    raw = load_or_die(os.path.join(DOCS, "fold_diversity_concordance_perm.csv"))
    cor = load_or_die(os.path.join(DOCS, "fold_diversity_concordance_corrected_perm.csv"))
    trig = load_or_die(os.path.join(DOCS, "triggers_from_pdb.csv"))

    raw_p_col = "perm_mean_concordance_p"
    # Corrected CSV uses the same column name when produced by
    # permutation_null_concordance.py --mode corrected.
    cor_p_col = "perm_mean_concordance_p" if raw_p_col in cor.columns \
                else "perm_mean_concordance_corrected_p"

    df_raw = raw.merge(trig[["pair_id", "trigger_class"]], on="pair_id", how="inner")
    df_cor = cor.merge(trig[["pair_id", "trigger_class"]], on="pair_id", how="inner")

    print(f"RAW   pairs in BH pool: {len(df_raw)}")
    print(f"COR   pairs in BH pool: {len(df_cor)}")

    raw_hits = stratified_bh(df_raw, raw_p_col,
                             label="RAW (pre-correction) within-class BH")
    cor_hits = stratified_bh(df_cor, cor_p_col,
                             label="CORRECTED within-class BH")

    # Reconcile
    print("\n=== Comparison: pairs gained / lost between RAW and CORRECTED ===")
    raw_set = {p for v in raw_hits.values() for (p, _, _) in v}
    cor_set = {p for v in cor_hits.values() for (p, _, _) in v}
    only_raw = raw_set - cor_set
    only_cor = cor_set - raw_set
    both = raw_set & cor_set

    def _line(label, pairs):
        if not pairs:
            print(f"  {label}: (none)")
        else:
            for p in sorted(pairs):
                cls = trig[trig["pair_id"] == p]["trigger_class"].iloc[0]
                raw_p = df_raw[df_raw["pair_id"] == p][raw_p_col].iloc[0] if p in df_raw["pair_id"].values else float("nan")
                cor_p = df_cor[df_cor["pair_id"] == p][cor_p_col].iloc[0] if p in df_cor["pair_id"].values else float("nan")
                print(f"  {label}: {p:<22} class={cls:<22} "
                      f"raw_p={raw_p:.4f} cor_p={cor_p:.4f}")

    print()
    _line("BOTH (robust)        ", both)
    _line("ONLY RAW (filtered)  ", only_raw)
    _line("ONLY CORRECTED (new) ", only_cor)

    print()
    print("Interpretation:")
    print("  BOTH    : signal robust to sequence-divergence correction.")
    print("  ONLY RAW: pre-correction signal driven (in part) by per-cluster")
    print("            PID-to-query gradient; demoted by the regression correction.")
    print("            These are the cases that motivate the correction step.")
    print("  ONLY COR: rare; corrected signal stronger than raw (concordance was")
    print("            masked pre-correction by method-specific PID dependencies).")


if __name__ == "__main__":
    main()
