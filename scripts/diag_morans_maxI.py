#!/usr/bin/env python3
"""Diagnose Moran's results by SIGNED effect size (not p, which floors at
1/(n_perm+1)=0.001 and ties everything).

Sign matters: I>0 = clade clustering, I<0 (below the -1/(n-1) null) = dispersion
/ sister-leaf flipping. We check whether the tail is positively biased (real
signal in a subset) or symmetric (noise), flag small-morans_n pairs (large |I|
can be a small-n artifact), and test whether the corrected max-|I| pair has
method-IDENTICAL residuals (the suspicious identical 0.456 -> a correction bug).
"""
import os
import numpy as np
import pandas as pd

DOCS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "docs")
CORR = "phylo_placement_corrected.csv"
RAW = "phylo_placement.csv"


def load(fn):
    f = os.path.join(DOCS, fn)
    if not os.path.isfile(f):
        print(f"  MISSING {f}"); return None
    d = pd.read_csv(f)
    return d if "morans_I" in d.columns else None


def signed_distribution(d, tag):
    print(f"\n=== {tag}: SIGNED Moran's I distribution per method ===")
    print(f"  (null E[I] = -1/(n-1); positive-tail = clade clustering, "
          f"negative = dispersion)")
    print(f"  {'method':9}{'n':>5}{'meanI':>8}{'medI':>8}{'I>0.2':>7}"
          f"{'0.1..0.2':>9}{'<-0.1':>7}{'maxI':>8}{'minI':>8}")
    for m in sorted(d.method.unique()):
        s = d[(d.method == m) & d.morans_I.notna()]
        if not len(s):
            continue
        I = s.morans_I
        print(f"  {m:9}{len(s):>5}{I.mean():>8.3f}{I.median():>8.3f}"
              f"{int((I > 0.2).sum()):>7}{int(((I > 0.1) & (I <= 0.2)).sum()):>9}"
              f"{int((I < -0.1).sum()):>7}{I.max():>8.3f}{I.min():>8.3f}")


def top_signed(d, tag, k=20):
    print(f"\n=== {tag}: top {k} pairs by SIGNED I (most positive = strongest "
          f"clade clustering) ===")
    cols = [c for c in ["pair_id", "method", "morans_n", "morans_I",
                        "morans_p_bh", "n_clusters"] if c in d.columns]
    print(d.sort_values("morans_I", ascending=False)[cols].head(k).to_string(index=False))
    print(f"\n--- {tag}: most negative 8 (dispersion) ---")
    print(d.sort_values("morans_I")[cols].head(8).to_string(index=False))


def cross_method_signed(d, k=12):
    """Pair x method matrix of signed I for the pairs with the largest
    mean positive I across methods -- shows whether methods AGREE (all +)."""
    piv = d.pivot_table(index="pair_id", columns="method", values="morans_I", aggfunc="first")
    # rank pairs by how positive they are on average across methods
    piv["_mean"] = piv.mean(axis=1, numeric_only=True)
    top = piv.sort_values("_mean", ascending=False).head(k).drop(columns="_mean")
    print(f"\n=== Cross-method SIGNED I for top {k} pairs (by mean across methods) ===")
    print(top.round(3).to_string())


def identical_residual_check(d):
    d = d.copy(); d["absI"] = d.morans_I.abs()
    pid = d.sort_values("absI", ascending=False)["pair_id"].iloc[0]
    nmax = d.sort_values("absI", ascending=False)["morans_n"].iloc[0] if "morans_n" in d else "?"
    print(f"\n=== max-|I| corrected pair = {pid}  (morans_n={nmax}) ===")
    sc = os.path.join(DOCS, "fold_diversity_survey_corrected.csv")
    if not os.path.isfile(sc):
        print(f"  MISSING {sc}"); return
    s = pd.read_csv(sc); s = s[s.pair_id == pid]
    ccol = "cluster_norm" if "cluster_norm" in s.columns else "cluster"
    piv = s.pivot_table(index=ccol, columns="method", values="TMdiff_residual", aggfunc="first")
    print(f"  per-method residuals (first 6 of {len(piv)} clusters):")
    print(piv.head(6).round(4).to_string())
    methods = list(piv.columns)
    if len(methods) >= 2:
        ref = piv[methods[0]].fillna(-999).round(6)
        identical = all(piv[m].fillna(-999).round(6).equals(ref) for m in methods[1:])
        print(f"\n  All methods' residuals identical for {pid}? {identical}")
        if identical:
            print("  >>> BUG: seq_divergence_correction emits method-agnostic "
                  "TMdiff_residual -> identical Moran's I across methods.")


if __name__ == "__main__":
    dc = load(CORR)
    dr = load(RAW)
    if dr is not None:
        signed_distribution(dr, "RAW")
    if dc is not None:
        signed_distribution(dc, "CORRECTED")
        top_signed(dc, "CORRECTED")
        cross_method_signed(dc)
        identical_residual_check(dc)
