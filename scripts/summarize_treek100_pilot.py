#!/usr/bin/env python3
"""
Aggregate all per-pair K=100 + medoid + top-10 pilot CSVs into a single
side-by-side comparison table.

Reads:
  docs/treek100_pilot_tmscores.csv                   (KaiB, default name)
  docs/treek100_pilot_tmscores_<pair>.csv            (other pairs)

For each pair x (AF2/AF3) computes:
  - n_clusters scored
  - mean/stddev/min/max of TMdiff_centered (chain-averaged)
  - fraction of clusters with |TMdiff_centered| > 0.10  ("signal" fraction)

Prints a single table to stdout and writes docs/treek100_pilot_summary.csv.

Usage:
    python3 scripts/summarize_treek100_pilot.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

PAIR_LABELS = {
    "5jytA_2qkeE": "KaiB",
    "2n54B_2hdmA": "XCL1",
    "1zk9A_3jv6A": "RelB",
    "2namA_1uxmK": "SOD1",
    "6c6sD_2ougC": "RfaH",
    "4qhhA_4qhfA": "Selecase",
}


def load_pair(csv_path: Path, pair_id: str) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    df["pair"] = pair_id
    return df


def find_csvs(docs_dir: Path) -> dict[str, Path]:
    """Return {pair_id: csv_path}. The KaiB run uses the default
    filename (no pair suffix); others have _<pair>.csv."""
    out: dict[str, Path] = {}
    default = docs_dir / "treek100_pilot_tmscores.csv"
    if default.is_file():
        # Heuristic: identify by which chain tags appear in the file.
        try:
            head = pd.read_csv(default, nrows=5)
            chains = set(head["chain"].astype(str).tolist())
            for pid in PAIR_LABELS:
                p1, p2 = pid.split("_")
                if p1 in chains or p2 in chains:
                    out[pid] = default
                    break
        except Exception:
            out["5jytA_2qkeE"] = default  # fall back to KaiB
    for f in sorted(docs_dir.glob("treek100_pilot_tmscores_*.csv")):
        stem = f.stem  # treek100_pilot_tmscores_<pair>
        pid = stem.replace("treek100_pilot_tmscores_", "")
        out[pid] = f
    return out


def per_pair_stats(df: pd.DataFrame) -> pd.DataFrame:
    """Reduce to one row per (pair, mode) with summary stats.

    Aggregation: average across the 2 chains per (pair, mode, cluster)
    *before* computing stddev — this matches the figure convention.
    """
    if "TMdiff_centered" not in df.columns:
        return pd.DataFrame()
    g = (df.dropna(subset=["TMdiff_centered"])
           .groupby(["pair", "mode", "cluster"])["TMdiff_centered"]
           .mean()
           .reset_index())
    out = (g.groupby(["pair", "mode"])["TMdiff_centered"]
             .agg(n="count", mean="mean", std="std",
                  vmin="min", vmax="max")
             .reset_index())
    # Fraction with |TMdiff_centered| > 0.10 -- "non-trivial" signal
    sig = (g.assign(big=lambda d: d["TMdiff_centered"].abs() > 0.10)
             .groupby(["pair", "mode"])["big"].mean()
             .reset_index(name="frac_abs_gt_0.10"))
    out = out.merge(sig, on=["pair", "mode"], how="left")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--docs-dir", default=str(ROOT / "docs"))
    ap.add_argument("--output", default=str(ROOT / "docs" / "treek100_pilot_summary.csv"))
    args = ap.parse_args()

    docs_dir = Path(args.docs_dir)
    csvs = find_csvs(docs_dir)
    if not csvs:
        sys.exit(f"No pilot CSVs found under {docs_dir}")
    print(f"[summary] found {len(csvs)} CSV(s):")
    for pid, f in csvs.items():
        print(f"  {PAIR_LABELS.get(pid, '?'):<10s} {pid:<15s} {f.name}")

    frames = []
    for pid, csv in csvs.items():
        try:
            frames.append(load_pair(csv, pid))
        except Exception as e:
            print(f"[summary] WARN failed to load {csv}: {e}")
    if not frames:
        sys.exit("No CSVs loaded.")
    full = pd.concat(frames, ignore_index=True)
    summary = per_pair_stats(full)
    if summary.empty:
        sys.exit("No TMdiff_centered column found in any CSV.")

    # Add human-readable pair label
    summary["protein"] = summary["pair"].map(PAIR_LABELS).fillna(summary["pair"])
    # Reorder columns and rows
    summary = summary[["protein", "pair", "mode", "n", "mean", "std",
                       "vmin", "vmax", "frac_abs_gt_0.10"]]
    summary = summary.sort_values(["protein", "mode"]).reset_index(drop=True)
    # Round for display
    disp = summary.copy()
    for c in ("mean", "std", "vmin", "vmax", "frac_abs_gt_0.10"):
        disp[c] = disp[c].round(3)

    print("\n[summary] === Per-(pair, mode) TMdiff_centered stats "
          "(chain-averaged per cluster) ===\n")
    print(disp.to_string(index=False))

    # Compact comparison view: pivot std + frac side-by-side
    pivot_std = summary.pivot(index="protein", columns="mode", values="std").round(3)
    pivot_frac = summary.pivot(index="protein", columns="mode", values="frac_abs_gt_0.10").round(3)
    print("\n[summary] === stddev pivot (rows=protein, cols=mode) ===\n")
    print(pivot_std.to_string())
    print("\n[summary] === fraction |TMdiff_centered|>0.10 ===\n")
    print(pivot_frac.to_string())

    out_csv = Path(args.output)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out_csv, index=False)
    print(f"\n[summary] wrote {out_csv}")


if __name__ == "__main__":
    main()
