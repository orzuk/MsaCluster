#!/usr/bin/env python3
"""Diagnose Moran's results by SIGNED standardized effect size.

z = (I_obs - mean(I_perm)) / sd(I_perm) references each pair to its OWN
permutation null, so it is comparable across pairs of different n (raw I's
null E[I]=-1/(n-1) shifts with n) and does not floor like the p-value.

Key checks:
  1. Per-method signed z distribution (positive = clade clustering).
  2. Top pairs by z (the proper candidate ranking).
  3. Confirm large-negative RAW I are just small-n nulls: list I < -0.1 with
     their morans_n and z (expect tiny n and z ~ 0).
  4. Cross-method signed-z matrix for the top pairs (who AGREES).
"""
import os
import pandas as pd

DOCS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "docs")


def load(fn):
    f = os.path.join(DOCS, fn)
    if not os.path.isfile(f):
        print(f"  MISSING {f}"); return None
    d = pd.read_csv(f)
    return d if "morans_I" in d.columns else None


def has_z(d):
    return "morans_z" in d.columns


def signed_distribution(d, tag):
    print(f"\n=== {tag}: signed distribution per method ===")
    z = has_z(d)
    hdr = f"  {'method':9}{'n':>5}{'medI':>8}{'maxI':>8}{'minI':>8}"
    if z:
        hdr += f"{'medZ':>8}{'z>2':>6}{'z<-2':>7}{'maxZ':>8}"
    print(hdr)
    for m in sorted(d.method.unique()):
        s = d[(d.method == m) & d.morans_I.notna()]
        if not len(s):
            continue
        I = s.morans_I
        line = f"  {m:9}{len(s):>5}{I.median():>8.3f}{I.max():>8.3f}{I.min():>8.3f}"
        if z:
            Z = s.morans_z.dropna()
            line += (f"{Z.median():>8.2f}{int((Z > 2).sum()):>6}"
                     f"{int((Z < -2).sum()):>7}{Z.max():>8.2f}")
        print(line)


def top_by_z(d, tag, k=20):
    sort_col = "morans_z" if has_z(d) else "morans_I"
    cols = [c for c in ["pair_id", "method", "morans_n", "morans_I",
                        "morans_z", "morans_p_bh", "n_clusters"] if c in d.columns]
    print(f"\n=== {tag}: top {k} by {sort_col} (strongest clade clustering) ===")
    print(d.sort_values(sort_col, ascending=False)[cols].head(k).to_string(index=False))


def confirm_neg_small_n(d):
    cols = [c for c in ["pair_id", "method", "morans_n", "morans_I",
                        "morans_z", "morans_p_bh"] if c in d.columns]
    neg = d[d.morans_I < -0.1].sort_values("morans_I")
    print(f"\n=== pairs with RAW I < -0.1 ({len(neg)}) — are they all small-n / z~0? ===")
    print(neg[cols].head(20).to_string(index=False))
    if "morans_n" in neg.columns and len(neg):
        print(f"  morans_n among these: max={neg.morans_n.max()}, "
              f"median={neg.morans_n.median()}  (small n -> null is very negative)")


def cross_method_z(d, k=12):
    val = "morans_z" if has_z(d) else "morans_I"
    piv = d.pivot_table(index="pair_id", columns="method", values=val, aggfunc="first")
    piv["_mean"] = piv.mean(axis=1, numeric_only=True)
    top = piv.sort_values("_mean", ascending=False).head(k).drop(columns="_mean")
    print(f"\n=== cross-method {val} for top {k} pairs (by mean across methods) ===")
    print(top.round(2).to_string())


if __name__ == "__main__":
    dc = load("phylo_placement_corrected.csv")
    dr = load("phylo_placement.csv")
    if dr is not None:
        signed_distribution(dr, "RAW")
    if dc is not None:
        signed_distribution(dc, "CORRECTED")
        confirm_neg_small_n(dc)
        top_by_z(dc, "CORRECTED")
        cross_method_z(dc)
