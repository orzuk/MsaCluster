#!/usr/bin/env python3
"""Summarize the v2 downstream results in one shot, from the docs/ CSVs.

Prints, with no arguments:
  1. Per-method Moran's phylogenetic-signal counts (RAW and CORRECTED):
     tested pairs, raw p<.05, BH<.05, median I.
  2. The corrected BH-significant (pair, method) hits — the candidate list.
  3. ESMFold2 scoring coverage (valid vs failed pairs).
  4. The cross-method Spearman correlation matrix of per-cluster TMdiff_centered
     (sanity check: AF2/AF3/ESM(=ESMFold2)/Boltz2 should be mutually positive).

Usage (from anywhere):
    python3 scripts/show_v2_phylo_results.py
"""
from __future__ import annotations
import glob
import os
import sys

import numpy as np
import pandas as pd

DOCS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "docs")
METHOD_ORDER = ["AF2", "AF3", "ESM", "Boltz2", "DDG", "MSAT", "CCMpred", "S4PRED"]


def morans_summary():
    for tag, fn in [("RAW", "phylo_placement.csv"),
                    ("CORRECTED", "phylo_placement_corrected.csv")]:
        f = os.path.join(DOCS, fn)
        print(f"\n=== Moran's {tag}  ({fn}) ===")
        if not os.path.isfile(f):
            print(f"  MISSING: {f}")
            continue
        d = pd.read_csv(f)
        if "morans_p" not in d.columns:
            print("  NO morans_p column (stale / pre-Moran's run)")
            continue
        print(f"  {'method':10}{'tested':>7}{'raw<.05':>9}{'BH<.05':>8}"
              f"{'medI':>8}{'p90|I|':>8}{'max|I|':>8}")
        for m in sorted(d.method.unique()):
            s = d[(d.method == m) & d.morans_p.notna()]
            if not len(s):
                continue
            bh = int((s["morans_p_bh"] < 0.05).sum()) if "morans_p_bh" in s.columns else 0
            ai = s.morans_I.abs()
            print(f"  {m:10}{len(s):>7}{int((s.morans_p < 0.05).sum()):>9}"
                  f"{bh:>8}{s.morans_I.median():>8.3f}{ai.quantile(0.9):>8.3f}{ai.max():>8.3f}")
        if tag == "CORRECTED" and "morans_p_bh" in d.columns:
            # sort by p, breaking the 1/(n_perm+1) floor ties by SIGNED I
            # (most positive = strongest clade clustering first)
            sig = d[d.morans_p_bh < 0.05].sort_values(
                ["morans_p", "morans_I"], ascending=[True, False])
            print(f"\n  Corrected BH<0.05: {len(sig)} (pair, method, I, p, p_bh)")
            for _, r in sig.head(40).iterrows():
                print(f"    {r.pair_id:14} {r.method:8} I={r.morans_I:6.3f} "
                      f"p={r.morans_p:.4f} bh={r.morans_p_bh:.4f}")
            if not len(sig):
                print("    (none survive correction)")


def esmfold2_coverage():
    ok = bad = 0
    for f in glob.glob(os.path.join(DOCS, "esmfold2_tmscores_*.csv")):
        d = pd.read_csv(f)
        v = pd.to_numeric(d.get("TM_F1"), errors="coerce").notna().sum() \
            if "TM_F1" in d.columns else 0
        ok += v > 0
        bad += v == 0
    print(f"\n=== ESMFold2 scoring: {ok} valid / {bad} failed (of {ok + bad}) ===")


def correlation_matrix():
    f = os.path.join(DOCS, "fold_diversity_survey.csv")
    print("\n=== Cross-method Spearman correlation (per-cluster TMdiff_centered) ===")
    if not os.path.isfile(f):
        print(f"  MISSING: {f}")
        return
    from scipy.stats import spearmanr
    d = pd.read_csv(f)
    d = d[d["cluster"].astype(str) != "DeepMsa"]
    wide = d.pivot_table(index=["pair_id", "cluster"], columns="method",
                         values="TMdiff_centered", aggfunc="first")
    cols = [m for m in METHOD_ORDER if m in wide.columns]
    arr = wide[cols].to_numpy(dtype=float)
    M = pd.DataFrame(index=cols, columns=cols, dtype=float)
    for i in range(len(cols)):
        for j in range(len(cols)):
            x, y = arr[:, i], arr[:, j]
            mask = np.isfinite(x) & np.isfinite(y)
            if mask.sum() < 5:
                M.iloc[i, j] = np.nan
            elif i == j:
                M.iloc[i, j] = 1.0
            else:
                M.iloc[i, j] = spearmanr(x[mask], y[mask]).correlation
    print(M.round(2).to_string())


if __name__ == "__main__":
    morans_summary()
    esmfold2_coverage()
    correlation_matrix()
