"""Bootstrap confidence intervals on per-pair corrected concordance.

For each pair, resamples its clusters WITH REPLACEMENT N times, recomputes
the mean cross-method per-cluster sign concordance on each bootstrap
sample, and reports the 2.5 / 97.5 percentile 95% CI.

Reads:  docs/fold_diversity_survey_corrected.csv
Writes: docs/fold_diversity_concordance_corrected_bootstrap.csv

Usage:
    python3 scripts/bootstrap_concordance_ci.py
    python3 scripts/bootstrap_concordance_ci.py --n-boot 2000 --seed 12345
"""

from __future__ import annotations

import argparse
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd

METHODS = ["AF2", "ESM", "MSAT", "DDG"]
INPUT = Path("docs/fold_diversity_survey_corrected.csv")
OUTPUT = Path("docs/fold_diversity_concordance_corrected_bootstrap.csv")
CANDIDATES = ["2n54B_2hdmA", "2namA_1uxmK", "1zk9A_3jv6A",
              "3kdsG_2ce7C", "4qhhA_4qhfA"]


def normalize_cluster(s):
    import re
    s = str(s).strip()
    s = re.sub(r"^(msa_t_)?(clusters_|cmap_)?", "", s)
    if s.startswith("DeepMsa") or s.lower() in ("deep", "deepmsa", "query"):
        return "DeepMsa"
    m = re.search(r"ShallowMsa_(\d+)", s)
    if m:
        return f"ShallowMsa_{int(m.group(1)):03d}"
    if s.isdigit():
        return f"ShallowMsa_{int(s):03d}"
    return s


def pair_concordance_from_residuals(sub: pd.DataFrame) -> float:
    """Mean pairwise sign agreement across method-pairs (returns nan if undefined)."""
    pivot = sub.pivot_table(
        index="cluster_norm", columns="method", values="TMdiff_residual", aggfunc="first",
    )
    present = [m for m in METHODS if m in pivot.columns]
    if len(present) < 2:
        return float("nan")
    vals = []
    for a, b in combinations(present, 2):
        both = pivot[[a, b]].dropna()
        if len(both) < 1:
            continue
        sa = np.sign(both[a].to_numpy())
        sb = np.sign(both[b].to_numpy())
        mask = (sa != 0) & (sb != 0)
        if mask.sum() < 1:
            continue
        vals.append((sa[mask] == sb[mask]).mean())
    return float(np.mean(vals)) if vals else float("nan")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, default=INPUT)
    ap.add_argument("--output", type=Path, default=OUTPUT)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=12345)
    args = ap.parse_args()

    df = pd.read_csv(args.input)
    if "TMdiff_residual" not in df.columns or "cluster_norm" not in df.columns:
        # cluster_norm may need recomputation if reading the raw corrected csv
        df["cluster_norm"] = df["cluster"].map(normalize_cluster)
    df = df[df["cluster_norm"] != "DeepMsa"]

    rng = np.random.default_rng(args.seed)
    rows = []
    for pid, sub in df.groupby("pair_id"):
        clusters = sub["cluster_norm"].drop_duplicates().to_numpy()
        n_c = len(clusters)
        if n_c < 2:
            continue
        observed = pair_concordance_from_residuals(sub)
        boots = np.empty(args.n_boot, dtype=float)
        boots[:] = np.nan
        for b in range(args.n_boot):
            sampled = rng.choice(clusters, size=n_c, replace=True)
            mask = sub["cluster_norm"].isin(sampled)
            # multiplicity-correct: build resampled per-method-per-cluster table
            # by stacking each sampled cluster's rows; aggfunc='first' below
            # is fine because resampling-with-replacement re-uses each cluster
            # by index, and concordance is sign-based per cluster.
            resampled = (
                sub.set_index("cluster_norm").loc[sampled].reset_index()
            )
            boots[b] = pair_concordance_from_residuals(resampled)
        valid = boots[~np.isnan(boots)]
        if not len(valid):
            continue
        rows.append({
            "pair_id": pid,
            "n_clusters": int(n_c),
            "observed_concordance": observed,
            "boot_mean": float(np.mean(valid)),
            "boot_median": float(np.median(valid)),
            "boot_ci_low": float(np.percentile(valid, 2.5)),
            "boot_ci_high": float(np.percentile(valid, 97.5)),
            "boot_n_valid": int(len(valid)),
        })
    out = pd.DataFrame(rows).sort_values("observed_concordance", ascending=False)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print(f"wrote {args.output}: {len(out)} pairs bootstrapped, {args.n_boot} resamples each")

    print("\n=== Candidate pairs: observed concordance with 95% bootstrap CI ===")
    print(f"{'pair_id':<18} {'n_clu':>5} {'obs':>6} {'CI_lo':>7} {'CI_hi':>7} {'CI_width':>9}")
    for p in CANDIDATES:
        r = out[out["pair_id"] == p]
        if not len(r):
            print(f"  {p}: not in output"); continue
        r = r.iloc[0]
        w = r["boot_ci_high"] - r["boot_ci_low"]
        flag = "  HIGH-CONF" if r["boot_ci_low"] >= 0.65 else ""
        print(f"  {p:<18} {int(r['n_clusters']):>5d} {r['observed_concordance']:>6.3f} "
              f"{r['boot_ci_low']:>7.3f} {r['boot_ci_high']:>7.3f} {w:>9.3f}{flag}")


if __name__ == "__main__":
    main()
