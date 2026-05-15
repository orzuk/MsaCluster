"""Benjamini-Hochberg multiple-testing correction on per-pair concordance p-values.

Reads docs/fold_diversity_concordance_perm.csv (empirical permutation p-values)
and applies Benjamini-Hochberg FDR correction across the per-pair
mean-concordance p-values. Reports adjusted p-values for the 4 candidate
pairs and the corrected high-confidence list.

Usage:
    python3 scripts/bh_correct_concordance_p.py
    python3 scripts/bh_correct_concordance_p.py --alpha 0.10
"""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd

INPUT = Path("docs/fold_diversity_concordance_perm.csv")
OUTPUT = Path("docs/fold_diversity_concordance_perm_bh.csv")
PAIRS = ["2n54B_2hdmA", "2namA_1uxmK", "1zk9A_3jv6A", "3kdsG_2ce7C",
         "4qhhA_4qhfA"]


def bh_adjust(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (one-step procedure)."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    # ensure monotone non-decreasing from large to small p
    adj_sorted = np.minimum.accumulate(ranked[::-1])[::-1]
    adj_sorted = np.minimum(adj_sorted, 1.0)
    adj = np.empty_like(adj_sorted)
    adj[order] = adj_sorted
    return adj


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, default=INPUT)
    ap.add_argument("--output", type=Path, default=OUTPUT)
    ap.add_argument("--alpha", type=float, default=0.05,
                    help="FDR threshold for highlighting (default 0.05)")
    args = ap.parse_args()

    df = pd.read_csv(args.input)
    if "perm_mean_concordance_p" not in df.columns:
        raise SystemExit("input missing perm_mean_concordance_p column")
    sub = df.dropna(subset=["perm_mean_concordance_p"]).copy()
    sub["p_bh"] = bh_adjust(sub["perm_mean_concordance_p"].to_numpy())
    sub = sub.sort_values("perm_mean_concordance_p")
    sub.to_csv(args.output, index=False)
    print(f"wrote {args.output}: {len(sub)} pairs tested, BH-adjusted")

    print(f"\n=== Pairs passing BH FDR <= {args.alpha} ===")
    passing = sub[sub["p_bh"] <= args.alpha]
    for _, r in passing.iterrows():
        print(f"  {r['pair_id']:<18}  raw_p={r['perm_mean_concordance_p']:.4f}  "
              f"BH_p={r['p_bh']:.4f}  n_clusters={int(r['n_clusters'])}  "
              f"observed_mean={r['observed_mean_concordance']:.3f}")
    if not len(passing):
        print("  (none)")

    print("\n=== Candidate pairs (whether or not they pass BH) ===")
    print(f"{'pair_id':<18} {'raw_p':>8} {'BH_p':>8} {'n_clu':>6} {'obs_mean':>9}")
    for p in PAIRS:
        r = sub[sub["pair_id"] == p]
        if not len(r):
            print(f"  {p}: NOT in perm CSV"); continue
        r = r.iloc[0]
        flag = "  PASS BH" if r["p_bh"] <= args.alpha else ""
        print(f"  {p:<18} {r['perm_mean_concordance_p']:>8.4f} "
              f"{r['p_bh']:>8.4f} {int(r['n_clusters']):>6d} "
              f"{r['observed_mean_concordance']:>9.3f}{flag}")


if __name__ == "__main__":
    main()
