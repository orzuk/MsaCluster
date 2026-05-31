#!/usr/bin/env python3
"""Inspect ONE pair's per-method per-cluster fold preference + Moran's stats.

Tells you whether a pair (e.g. KaiB 5jytA_2qkeE) shows flat per-cluster
values (methods didn't capture per-cluster variation -> no signal) or has
variance that simply isn't clade-aligned (Moran's I ~ 0 despite spread).

Usage:
    python3 scripts/inspect_pair.py 5jytA_2qkeE
"""
import os
import sys
import pandas as pd

DOCS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "docs")


def main(pair):
    # 1. per-cluster fold-preference spread per method (is it flat?)
    for fn, val, label in [("fold_diversity_survey.csv", "TMdiff_centered", "RAW"),
                           ("fold_diversity_survey_corrected.csv", "TMdiff_residual", "CORRECTED")]:
        f = os.path.join(DOCS, fn)
        if not os.path.isfile(f):
            continue
        d = pd.read_csv(f)
        d = d[d.pair_id == pair]
        if "cluster" in d.columns:
            d = d[d["cluster"].astype(str) != "DeepMsa"]
        if not len(d) or val not in d.columns:
            print(f"\n=== {label}: no data / no {val} for {pair} ===")
            continue
        print(f"\n=== {label}: per-cluster {val} spread for {pair} ===")
        print(f"  {'method':9}{'n_clu':>6}{'std':>10}{'min':>9}{'max':>9}{'frac|v|>0.1':>12}")
        for m, g in d.groupby("method"):
            v = pd.to_numeric(g[val], errors="coerce").dropna()
            if not len(v):
                continue
            print(f"  {m:9}{len(v):>6}{v.std():>10.4f}{v.min():>9.3f}"
                  f"{v.max():>9.3f}{(v.abs() > 0.1).mean():>12.2f}")

    # 2. Moran's stats for this pair (raw + corrected)
    for fn, label in [("phylo_placement.csv", "RAW"),
                      ("phylo_placement_corrected.csv", "CORRECTED")]:
        f = os.path.join(DOCS, fn)
        if not os.path.isfile(f):
            continue
        d = pd.read_csv(f)
        d = d[d.pair_id == pair]
        if not len(d):
            continue
        cols = [c for c in ["method", "morans_n", "morans_I", "morans_z",
                            "morans_p", "morans_p_bh", "n_F1c_leaves",
                            "n_F2c_leaves"] if c in d.columns]
        print(f"\n=== {label}: Moran's for {pair} ===")
        print(d[cols].to_string(index=False))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("usage: python3 scripts/inspect_pair.py <pair_id>")
    main(sys.argv[1])
