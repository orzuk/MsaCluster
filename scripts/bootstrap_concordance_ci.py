"""Bootstrap confidence intervals on per-pair corrected concordance.

For each pair, resamples its clusters WITH REPLACEMENT N times, recomputes
the mean cross-method per-cluster sign concordance on each bootstrap
sample, and reports the 2.5 / 97.5 percentile 95% CI.

Reads:  docs/fold_diversity_survey_corrected.csv
Writes: docs/fold_diversity_concordance_corrected_bootstrap.csv

Usage:
    python3 scripts/bootstrap_concordance_ci.py
    python3 scripts/bootstrap_concordance_ci.py --n-boot 2000 --seed 12345
    python3 scripts/bootstrap_concordance_ci.py --n-boot 200    # quick

Prints a per-pair progress line so you can see it advancing.
"""

from __future__ import annotations

import argparse
import re
import time
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd

METHODS = ["AF2", "ESM", "MSAT", "DDG"]
INPUT = Path("docs/fold_diversity_survey_corrected.csv")
OUTPUT = Path("docs/fold_diversity_concordance_corrected_bootstrap.csv")
CANDIDATES = ["2n54B_2hdmA", "2namA_1uxmK", "1zk9A_3jv6A",
              "3kdsG_2ce7C", "4qhhA_4qhfA"]


def normalize_cluster(s: str) -> str:
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


def concordance_from_sign_matrix(sign_mat: np.ndarray) -> float:
    """sign_mat: (n_clusters, n_methods) of {-1, 0, +1}. Returns mean pairwise sign agreement."""
    n_methods = sign_mat.shape[1]
    if n_methods < 2:
        return float("nan")
    vals: list[float] = []
    for a, b in combinations(range(n_methods), 2):
        sa = sign_mat[:, a]
        sb = sign_mat[:, b]
        mask = (sa != 0) & (sb != 0) & ~np.isnan(sa) & ~np.isnan(sb)
        if mask.sum() < 1:
            continue
        vals.append((sa[mask] == sb[mask]).mean())
    if not vals:
        return float("nan")
    return float(np.mean(vals))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, default=INPUT)
    ap.add_argument("--output", type=Path, default=OUTPUT)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--progress-every", type=int, default=1,
                    help="print progress every N pairs (default 1)")
    args = ap.parse_args()

    print(f"Loading {args.input} ...")
    df = pd.read_csv(args.input)
    if "cluster_norm" not in df.columns:
        df["cluster_norm"] = df["cluster"].map(normalize_cluster)
    df = df[df["cluster_norm"] != "DeepMsa"]
    print(f"  {len(df)} rows across {df['pair_id'].nunique()} pairs after dropping DeepMsa")

    rng = np.random.default_rng(args.seed)
    rows = []
    pairs = sorted(df["pair_id"].unique())
    t0 = time.time()
    for i, pid in enumerate(pairs):
        sub = df[df["pair_id"] == pid]
        # Build per-cluster sign-matrix: rows = clusters, cols = METHODS
        clusters = sorted(sub["cluster_norm"].unique())
        n_c = len(clusters)
        if n_c < 2:
            if i % args.progress_every == 0:
                print(f"  [{i+1:>3}/{len(pairs)}] {pid:<18}  skipped (n_clusters={n_c})")
            continue
        cluster_to_idx = {c: k for k, c in enumerate(clusters)}
        sign_mat = np.full((n_c, len(METHODS)), np.nan, dtype=float)
        for j, m in enumerate(METHODS):
            ms = sub[sub["method"] == m]
            for _, r in ms.iterrows():
                idx = cluster_to_idx.get(r["cluster_norm"])
                v = r.get("TMdiff_residual", np.nan)
                if idx is not None and not pd.isna(v):
                    sign_mat[idx, j] = np.sign(float(v))
        observed = concordance_from_sign_matrix(sign_mat)
        if np.isnan(observed):
            if i % args.progress_every == 0:
                print(f"  [{i+1:>3}/{len(pairs)}] {pid:<18}  no valid concordance")
            continue

        # Bootstrap: resample cluster indices with replacement, slice sign_mat
        idxs = np.arange(n_c)
        boots = np.empty(args.n_boot, dtype=float)
        for b in range(args.n_boot):
            sampled = rng.choice(idxs, size=n_c, replace=True)
            boots[b] = concordance_from_sign_matrix(sign_mat[sampled, :])
        valid = boots[~np.isnan(boots)]
        if not len(valid):
            continue
        ci_lo = float(np.percentile(valid, 2.5))
        ci_hi = float(np.percentile(valid, 97.5))
        rows.append({
            "pair_id": pid,
            "n_clusters": int(n_c),
            "observed_concordance": observed,
            "boot_mean": float(np.mean(valid)),
            "boot_median": float(np.median(valid)),
            "boot_ci_low": ci_lo,
            "boot_ci_high": ci_hi,
            "boot_n_valid": int(len(valid)),
        })
        if i % args.progress_every == 0:
            elapsed = time.time() - t0
            eta = elapsed / (i + 1) * (len(pairs) - i - 1)
            print(f"  [{i+1:>3}/{len(pairs)}] {pid:<18}  obs={observed:.3f}  "
                  f"CI=[{ci_lo:.3f}, {ci_hi:.3f}]  n_clu={n_c}  "
                  f"elapsed={elapsed:.0f}s  eta={eta:.0f}s")

    out = pd.DataFrame(rows).sort_values("observed_concordance", ascending=False)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print(f"\nwrote {args.output}: {len(out)} pairs bootstrapped, {args.n_boot} resamples each")
    print(f"total wall time: {time.time() - t0:.0f}s")

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
