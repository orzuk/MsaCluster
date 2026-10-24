#!/usr/bin/env python3
"""Cluster-count-matched comparison of predicted structural dispersion.

Fold switchers have many more clusters than controls (median 100 against 17),
and the smallest pairwise TM-score among a target's clusters is an extreme
value: its expectation falls as the number of compared pairs grows, so part of
the switcher-control gap in that statistic is arithmetic rather than
biological. Resampling a fixed number of clusters from every target puts both
groups on the same number of pairwise comparisons and removes that artifact.

For each cluster budget k, each target with at least k clusters contributes the
mean over draws of the minimum, mean and 10th-percentile TM-score among k
randomly chosen clusters. Switchers and controls are then compared with a
Mann-Whitney U test.

Usage: tm_subsample_matched.py MATRIX_DIR OUT_CSV [--draws 200]
"""
from __future__ import annotations

import argparse
import glob
import os

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


def target_stats(mat: np.ndarray, k: int, draws: int, rng) -> dict | None:
    """Mean over draws of min / mean / p10 TM among k sampled clusters."""
    n = mat.shape[0]
    if n < k:
        return None
    iu = np.triu_indices(k, 1)
    mins, means, p10s = [], [], []
    for _ in range(draws):
        idx = rng.choice(n, size=k, replace=False)
        sub = mat[np.ix_(idx, idx)][iu]
        sub = sub[np.isfinite(sub)]
        if sub.size == 0:
            continue
        mins.append(sub.min())
        means.append(sub.mean())
        p10s.append(np.quantile(sub, 0.1))
    if not mins:
        return None
    return dict(min_tm=float(np.mean(mins)), mean_tm=float(np.mean(means)),
                p10_tm=float(np.mean(p10s)))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("matrix_dir")
    ap.add_argument("out_csv")
    ap.add_argument("--draws", type=int, default=200)
    ap.add_argument("--seed", type=int, default=20261024)
    a = ap.parse_args()

    files = sorted(glob.glob(os.path.join(a.matrix_dir, "*.npz")))
    targets = []
    for f in files:
        d = np.load(f, allow_pickle=True)
        mat = np.asarray(d["tm"], dtype=float)
        np.fill_diagonal(mat, np.nan)
        grp = str(d["group"])
        targets.append((str(d["target"]),
                        "FoldSwitch" if grp == "FoldSwitch" else "Control",
                        mat))
    print(f"{len(targets)} targets "
          f"({sum(t[1] == 'FoldSwitch' for t in targets)} switchers, "
          f"{sum(t[1] == 'Control' for t in targets)} controls)")

    rows = []
    for k in (3, 5, 8, 10, 15, 17, 20, 30, 50):
        rng = np.random.default_rng([a.seed, k])
        per = {"FoldSwitch": [], "Control": []}
        for name, grp, mat in targets:
            st = target_stats(mat, k, a.draws, rng)
            if st:
                per[grp].append(st)
        if len(per["FoldSwitch"]) < 5 or len(per["Control"]) < 5:
            continue
        row = dict(k=k, n_switch=len(per["FoldSwitch"]),
                   n_control=len(per["Control"]))
        for stat in ("min_tm", "mean_tm", "p10_tm"):
            sw = np.array([d[stat] for d in per["FoldSwitch"]])
            ct = np.array([d[stat] for d in per["Control"]])
            u, p = mannwhitneyu(sw, ct, alternative="two-sided")
            row[f"{stat}_switch"] = float(np.median(sw))
            row[f"{stat}_control"] = float(np.median(ct))
            row[f"{stat}_p"] = float(p)
        rows.append(row)

    out = pd.DataFrame(rows)
    out.to_csv(a.out_csv, index=False)
    pd.set_option("display.width", 200)
    print(out.to_string(index=False, float_format=lambda v: f"{v:.4g}"))


if __name__ == "__main__":
    main()
