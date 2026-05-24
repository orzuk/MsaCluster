#!/usr/bin/env python3
"""Per-method, per-pair within-pair variance of TMdiff_centered.

For each of the 8 methods, computes (across the 93 pairs) the distribution of
per-pair standard deviations of TMdiff_centered across the pair's shallow
clusters (DeepMsa rows excluded). A method whose per-pair stddev is consistently
small across pairs is producing flat per-cluster predictions (saturating to one
fold; potentially noise). A method with larger per-pair stddev is producing
genuinely cluster-discriminating signal.

Compares all 8 methods on the same scale so we can see how the 8th method
(S4PRED) ranks relative to the other 7 in terms of within-pair informativeness.

Reads
-----
  docs/fold_diversity_survey.csv

Output
------
  stdout: a table per method showing median, mean, IQR, min, max of per-pair
          stddev, plus number of pairs covered.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

SURVEY = os.path.join("docs", "fold_diversity_survey.csv")
METHOD_ORDER = ["AF2", "AF3", "ESM", "MSAT", "CCMpred", "DDG", "Boltz2", "S4PRED"]


def main() -> None:
    if not os.path.isfile(SURVEY):
        print(f"[error] missing {SURVEY}", file=sys.stderr)
        sys.exit(2)
    df = pd.read_csv(SURVEY)
    df = df[~df["cluster"].astype(str).str.lower().str.startswith("deep")]
    df["TMdiff_centered"] = pd.to_numeric(df["TMdiff_centered"], errors="coerce")

    print(f"Survey: {len(df)} (pair, cluster, method) rows (DeepMsa excluded)\n")
    print(f"{'method':<10} {'n_pairs':>7} {'median_sd':>10} {'mean_sd':>10} "
          f"{'q25_sd':>9} {'q75_sd':>9} {'min_sd':>9} {'max_sd':>9}")
    print("-" * 80)

    for method in METHOD_ORDER:
        sub = df[df["method"] == method]
        if sub.empty:
            print(f"{method:<10} {'(no rows)':>7}")
            continue
        # Per-pair stddev across shallow clusters
        sd = sub.groupby("pair_id")["TMdiff_centered"].std(ddof=0)
        # Only include pairs with >= 2 clusters (so stddev is defined)
        sd = sd.dropna()
        if sd.empty:
            print(f"{method:<10} {0:>7} (all NaN)")
            continue
        print(f"{method:<10} {len(sd):>7d} "
              f"{sd.median():>10.4f} {sd.mean():>10.4f} "
              f"{sd.quantile(0.25):>9.4f} {sd.quantile(0.75):>9.4f} "
              f"{sd.min():>9.4f} {sd.max():>9.4f}")

    # Also dump per-pair counts of non-NaN clusters per method (signal coverage)
    print("\nCluster coverage per method (median across pairs):")
    print(f"{'method':<10} {'median_clusters':>15} {'pairs_with_any':>16}")
    for method in METHOD_ORDER:
        sub = df[df["method"] == method]
        if sub.empty:
            print(f"{method:<10} {'-':>15} {'-':>16}")
            continue
        per_pair_n = sub.groupby("pair_id")["TMdiff_centered"].apply(
            lambda x: x.notna().sum())
        n_with = (per_pair_n > 0).sum()
        print(f"{method:<10} {per_pair_n.median():>15.1f} {n_with:>16d}")

    # Per-method, what fraction of clusters cross a threshold (visually call F1c/F2c)?
    threshold = 0.05
    print(f"\nFraction of (pair, cluster) rows with |TMdiff_centered| > {threshold}:")
    print(f"{'method':<10} {'frac_nonzero':>13}")
    for method in METHOD_ORDER:
        sub = df[df["method"] == method]["TMdiff_centered"]
        sub = sub.dropna()
        if sub.empty:
            print(f"{method:<10} {'-':>13}")
            continue
        frac = (sub.abs() > threshold).mean()
        print(f"{method:<10} {frac:>13.3f}")


if __name__ == "__main__":
    main()
