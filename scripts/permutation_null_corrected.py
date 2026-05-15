"""Empirical permutation null on POST-correction concordance + BH FDR.

The original `permutation_null_concordance.py` operates on TMdiff_centered
(pre-correction). After applying the sequence-divergence regression
correction (`scripts/seq_divergence_correction.py`), the per-cluster
quantity that drives cross-method concordance is TMdiff_residual, not
TMdiff_centered, so the original p-values do not apply.

This script computes empirical permutation null p-values for the
post-correction mean concordance per pair, then applies BH FDR across
all evaluated pairs.

For each pair p:
  1. For each method m, gather the residual vector r_{m,c} from
     docs/fold_diversity_survey_corrected.csv.
  2. Compute observed pairwise mean sign-concordance bar(A_p).
  3. N_PERM times, independently shuffle each method's residual vector
     within the pair (preserving per-method marginals but breaking
     per-cluster sign alignment). Recompute bar(A_p)^shuf.
  4. Empirical one-sided p = Pr[bar(A)^shuf >= bar(A)^obs] / N_PERM.
  5. BH-correct over pairs (default alpha = 0.05).

Reads:
    docs/fold_diversity_survey_corrected.csv
Writes:
    docs/fold_diversity_concordance_corrected_perm.csv  (empirical p)
    docs/fold_diversity_concordance_corrected_perm_bh.csv  (+ BH adj p)

Usage:
    python3 scripts/permutation_null_corrected.py
    python3 scripts/permutation_null_corrected.py --n-perm 5000 --alpha 0.10
"""

from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

METHODS = ["AF2", "ESM", "MSAT", "DDG"]
SURVEY = Path("docs/fold_diversity_survey_corrected.csv")
OUT_PERM = Path("docs/fold_diversity_concordance_corrected_perm.csv")
OUT_BH = Path("docs/fold_diversity_concordance_corrected_perm_bh.csv")
CANDIDATES = ["2n54B_2hdmA", "2namA_1uxmK", "1zk9A_3jv6A",
              "3kdsG_2ce7C", "4qhhA_4qhfA"]


def mean_concordance_from_residuals(per_method: dict[str, np.ndarray]) -> tuple[float, int]:
    """Mean pairwise sign agreement across method-pairs. Returns (mean_A, mean_n)."""
    present = [m for m in METHODS if m in per_method and len(per_method[m]) > 0]
    if len(present) < 2:
        return float("nan"), 0
    agreements: list[float] = []
    ns: list[int] = []
    for a, b in combinations(present, 2):
        va, vb = per_method[a], per_method[b]
        n = min(len(va), len(vb))
        if n < 1:
            continue
        sa = np.sign(va[:n])
        sb = np.sign(vb[:n])
        mask = (sa != 0) & (sb != 0)
        if mask.sum() < 1:
            continue
        agreements.append((sa[mask] == sb[mask]).mean())
        ns.append(int(mask.sum()))
    if not agreements:
        return float("nan"), 0
    return float(np.mean(agreements)), int(np.mean(ns))


def bh_adjust(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    adj_sorted = np.minimum.accumulate(ranked[::-1])[::-1]
    adj_sorted = np.minimum(adj_sorted, 1.0)
    adj = np.empty_like(adj_sorted)
    adj[order] = adj_sorted
    return adj


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, default=SURVEY)
    ap.add_argument("--out-perm", type=Path, default=OUT_PERM)
    ap.add_argument("--out-bh", type=Path, default=OUT_BH)
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--alpha", type=float, default=0.05,
                    help="BH FDR threshold for highlighting (default 0.05)")
    args = ap.parse_args()

    df = pd.read_csv(args.input)
    needed = {"pair_id", "method", "TMdiff_residual"}
    if not needed.issubset(df.columns):
        raise SystemExit(f"missing columns in {args.input}: need {needed - set(df.columns)}")
    if "cluster_norm" not in df.columns:
        df["cluster_norm"] = df["cluster"]
    df = df[df["cluster_norm"] != "DeepMsa"].copy()

    rng = np.random.default_rng(args.seed)
    rows = []
    for pid, sub in df.groupby("pair_id"):
        per_method_obs: dict[str, np.ndarray] = {}
        for m in METHODS:
            ms = sub[sub["method"] == m]
            vals = ms["TMdiff_residual"].dropna().to_numpy(dtype=float)
            if len(vals):
                per_method_obs[m] = vals
        if len(per_method_obs) < 2:
            continue
        obs_mean, obs_n = mean_concordance_from_residuals(per_method_obs)
        if np.isnan(obs_mean):
            continue
        n_clusters = max(len(v) for v in per_method_obs.values())

        n_ge = 0
        for _ in range(args.n_perm):
            permuted = {m: rng.permutation(v) for m, v in per_method_obs.items()}
            perm_mean, _ = mean_concordance_from_residuals(permuted)
            if not np.isnan(perm_mean) and perm_mean >= obs_mean:
                n_ge += 1
        p_emp = (n_ge + 1) / (args.n_perm + 1)

        rows.append({
            "pair_id": pid,
            "n_clusters": int(n_clusters),
            "n_methods_present": len(per_method_obs),
            "observed_mean_concordance_corrected": obs_mean,
            "observed_mean_n_clusters": obs_n,
            "n_perm": args.n_perm,
            "perm_mean_concordance_p": p_emp,
        })

    out = pd.DataFrame(rows).sort_values("perm_mean_concordance_p")
    out["p_bh"] = bh_adjust(out["perm_mean_concordance_p"].to_numpy())
    args.out_perm.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_perm, index=False)
    out.to_csv(args.out_bh, index=False)
    print(f"wrote {args.out_perm} and {args.out_bh}: {len(out)} pairs, N_PERM={args.n_perm}")

    passing = out[out["p_bh"] <= args.alpha]
    print(f"\n=== Pairs passing post-correction BH FDR <= {args.alpha} ===")
    for _, r in passing.iterrows():
        print(f"  {r['pair_id']:<18}  obs={r['observed_mean_concordance_corrected']:.3f}  "
              f"raw_p={r['perm_mean_concordance_p']:.4f}  BH_p={r['p_bh']:.4f}  "
              f"n_clu={int(r['n_clusters'])}")
    if not len(passing):
        print("  (none)")

    print("\n=== Candidate pairs (post-correction, with empirical perm + BH) ===")
    print(f"{'pair_id':<18} {'n_clu':>5} {'obs':>6} {'raw_p':>8} {'BH_p':>8}")
    out_idx = out.set_index("pair_id")
    for p in CANDIDATES:
        if p not in out_idx.index:
            print(f"  {p}: not in output"); continue
        r = out_idx.loc[p]
        flag = "  PASS BH" if r["p_bh"] <= args.alpha else ""
        print(f"  {p:<18} {int(r['n_clusters']):>5d} "
              f"{r['observed_mean_concordance_corrected']:>6.3f} "
              f"{r['perm_mean_concordance_p']:>8.4f} {r['p_bh']:>8.4f}{flag}")


if __name__ == "__main__":
    main()
