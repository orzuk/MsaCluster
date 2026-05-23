#!/usr/bin/env python3
"""Stratified BH FDR on per-pair concordance permutation p-values, by trigger class.

This script tests a pre-specified hypothesis from the fold-switching
literature (Chakravarty et al. 2023, Distinguishing fold-switchers from
non-fold-switchers): TRIGGERED / stratified fold-switchers (where the
switch depends on ligand, mutation, oligomerization, or binding partner)
should show stronger clade-level fold-preference structure than
EQUILIBRIUM switchers (where both folds interconvert reversibly within
each clade).

Statistical anchor
------------------
The headline analysis applies BH FDR over the full evaluable pool (~43
pairs after filtering), which is dramatically over-conservative if the
biological signal is concentrated in a structural subclass. By
pre-specifying a trigger-class stratification (5 classes from
`docs/triggers_from_pdb.csv`, auto-generated from PDB metadata by
`scripts/classify_triggers_from_pdb.py`), we accomplish two things:

  1. Per-class BH FDR (over n=13-29 pairs depending on class). A smaller
     pool is a less penalizing multiple-testing burden, so candidates
     diluted by the full-screen BH may surface within their class.

  2. Pooled 1-d.o.f. test of the Chakravarty hypothesis: pool the four
     triggered classes (78 pairs) and compare the BH-significant rate
     against the equilibrium pool (15 pairs) via Fisher's exact test.

To keep this defensible (NOT post-hoc p-hacking), the trigger
classification MUST be locked in before any candidates are inspected
class-by-class. The classification in `docs/triggers_from_pdb.csv` was
generated from PDB metadata (ligands, mutation, chain count) by
`classify_triggers_from_pdb.py` on 2026-05-19 and is committed to the
repo. This script reads it as-is.

Inputs (defaults are the post-correction CSVs; pass --perm-csv for raw):
  docs/fold_diversity_concordance_corrected_perm.csv
                     (per-pair empirical perm p-value from
                      permutation_null_concordance.py --mode corrected)
  docs/triggers_from_pdb.csv
                     (per-pair trigger_class label)

Outputs:
  docs/stratified_concordance.csv
                     (per trigger class: n_pairs, n_BH_sig at alpha,
                      hit_rate, comma-separated hit list)

Usage:
  python3 scripts/stratified_concordance_by_trigger.py
  python3 scripts/stratified_concordance_by_trigger.py \\
          --perm-csv docs/fold_diversity_concordance_perm.csv  # raw
  python3 scripts/stratified_concordance_by_trigger.py --alpha 0.10
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Reuse the BH function from the consolidated permutation_null script
_THIS = Path(__file__).resolve().parent
sys.path.insert(0, str(_THIS))
from permutation_null_concordance import bh_adjust  # noqa: E402


TRIGGERED_CLASSES = {"ligand", "oligomerization", "protein_binding", "mutation"}
EQUILIBRIUM_CLASSES = {"equilibrium_or_unknown"}

DEFAULT_PERM_CSV = "docs/fold_diversity_concordance_corrected_perm.csv"
DEFAULT_TRIGGERS_CSV = "docs/triggers_from_pdb.csv"
DEFAULT_OUTPUT = "docs/stratified_concordance.csv"


def _detect_p_column(df: pd.DataFrame) -> str:
    """Find the raw permutation p-value column.

    permutation_null_concordance.py uses 'perm_mean_concordance_p' in both
    raw and corrected modes. Older variants may also expose
    'perm_mean_concordance_corrected_p'.
    """
    for cand in ("perm_mean_concordance_p", "perm_mean_concordance_corrected_p"):
        if cand in df.columns:
            return cand
    raise SystemExit(
        f"no recognised p-value column in input (looked for "
        f"perm_mean_concordance_p). columns: {list(df.columns)}"
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--perm-csv", default=DEFAULT_PERM_CSV,
                    help=f"Per-pair permutation-p CSV "
                         f"(default: {DEFAULT_PERM_CSV})")
    ap.add_argument("--triggers-csv", default=DEFAULT_TRIGGERS_CSV,
                    help=f"Trigger classification CSV "
                         f"(default: {DEFAULT_TRIGGERS_CSV})")
    ap.add_argument("--output", default=DEFAULT_OUTPUT,
                    help=f"Output CSV (default: {DEFAULT_OUTPUT})")
    ap.add_argument("--alpha", type=float, default=0.05,
                    help="BH FDR threshold (default: 0.05)")
    ap.add_argument("--min-class-size", type=int, default=10,
                    help="Skip classes with fewer than N pairs (default: 10)")
    args = ap.parse_args()

    if not os.path.isfile(args.perm_csv):
        raise SystemExit(f"perm CSV not found: {args.perm_csv}")
    if not os.path.isfile(args.triggers_csv):
        raise SystemExit(f"triggers CSV not found: {args.triggers_csv}")

    perm = pd.read_csv(args.perm_csv)
    triggers = pd.read_csv(args.triggers_csv)

    if "trigger_class" not in triggers.columns:
        raise SystemExit(f"triggers CSV missing 'trigger_class' column: {args.triggers_csv}")
    if "pair_id" not in triggers.columns:
        raise SystemExit(f"triggers CSV missing 'pair_id' column: {args.triggers_csv}")

    p_col = _detect_p_column(perm)

    # Merge — INNER JOIN ensures only pairs present in BOTH CSVs are tested.
    df = perm.merge(triggers[["pair_id", "trigger_class"]],
                    on="pair_id", how="inner")
    if df.empty:
        raise SystemExit("no overlap between perm CSV and triggers CSV "
                         "after merge on pair_id")
    print(f"Loaded {len(perm)} pairs from {args.perm_csv}")
    print(f"Loaded {len(triggers)} pairs from {args.triggers_csv}")
    print(f"After inner merge on pair_id: {len(df)} pairs evaluable")
    print(f"Using p-value column: {p_col}")
    print(f"Trigger-class distribution in merged set:")
    for cls, n in df["trigger_class"].value_counts().items():
        marker = "TRIGGERED" if cls in TRIGGERED_CLASSES else \
                 "EQUILIBRIUM" if cls in EQUILIBRIUM_CLASSES else "OTHER"
        print(f"    {cls:<25s} n={n:>3d}  ({marker})")
    print()

    # ----- (1) Per-class BH -----
    print(f"=== Per-class Benjamini-Hochberg FDR <= {args.alpha} ===")
    rows = []
    for cls, sub in df.groupby("trigger_class"):
        n = len(sub)
        if n < args.min_class_size:
            print(f"  {cls:<25s} n={n} < {args.min_class_size}; skipped")
            continue
        sub = sub.copy()
        sub["p_bh_within_class"] = bh_adjust(sub[p_col].to_numpy(dtype=float))
        hits = sub[sub["p_bh_within_class"] <= args.alpha].sort_values("p_bh_within_class")
        rows.append({
            "trigger_class": cls,
            "n_pairs": n,
            "n_BH_sig_at_alpha": len(hits),
            "hit_rate": round(len(hits) / n, 3),
            "hits": ",".join(hits["pair_id"].tolist()),
            "alpha": args.alpha,
        })
        hit_str = ", ".join(hits["pair_id"].tolist()) or "(none)"
        print(f"  {cls:<25s} n={n:>3d}  BH-sig={len(hits):>2d}  hits: {hit_str}")

    out_df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    out_df.to_csv(args.output, index=False)
    print(f"\nWrote {args.output}")

    # ----- (2) Chakravarty 2023 pooled test -----
    triggered = df[df["trigger_class"].isin(TRIGGERED_CLASSES)].copy()
    equilibrium = df[df["trigger_class"].isin(EQUILIBRIUM_CLASSES)].copy()

    if len(triggered) and len(equilibrium):
        triggered["p_bh"] = bh_adjust(triggered[p_col].to_numpy(dtype=float))
        equilibrium["p_bh"] = bh_adjust(equilibrium[p_col].to_numpy(dtype=float))

        t_hits = int((triggered["p_bh"] <= args.alpha).sum())
        e_hits = int((equilibrium["p_bh"] <= args.alpha).sum())
        t_n, e_n = len(triggered), len(equilibrium)

        print(f"\n=== Chakravarty 2023 prediction test (pre-specified) ===")
        print(f"  Triggered pool (ligand+oligomerization+protein_binding+mutation):")
        print(f"    n={t_n}  BH-sig={t_hits}  hit_rate={100*t_hits/t_n:.1f}%")
        print(f"  Equilibrium pool (equilibrium_or_unknown):")
        print(f"    n={e_n}  BH-sig={e_hits}  hit_rate={100*e_hits/e_n:.1f}%")

        try:
            from scipy.stats import fisher_exact
            table = [[t_hits, t_n - t_hits],
                     [e_hits, e_n - e_hits]]
            odds, p_two = fisher_exact(table, alternative="two-sided")
            _, p_greater = fisher_exact(table, alternative="greater")
            print(f"  Fisher's exact test:")
            print(f"    odds ratio (triggered / equilibrium): {odds:.3f}")
            print(f"    p (two-sided):                         {p_two:.4f}")
            print(f"    p (one-sided, triggered > equilibrium): {p_greater:.4f}")
        except ImportError:
            print("    [scipy not available; skipping Fisher's exact]")
    else:
        print(f"\n[warn] cannot run pooled Chakravarty test: "
              f"triggered n={len(triggered)}, equilibrium n={len(equilibrium)}")


if __name__ == "__main__":
    main()
