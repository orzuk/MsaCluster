#!/usr/bin/env python3
"""Score AF2/AF3 per-cluster predictions vs the two truth folds.

Thin wrapper over utils.score_against_folds.score_predictions_vs_truth.
Handles enumerating AF2/AF3 prediction PDBs (different directory layouts
per AF version + per-chain subdirs) and delegates the TM-align loop +
DataFrame construction to the shared util.

Usage:
    python3 scripts/score_treek100_pilot.py --pair 5jytA_2qkeE
    python3 scripts/score_treek100_pilot.py --pair 5jytA_2qkeE \\
        --modes AF2_newmethod_top10,AF3_newmethod_top10
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Optional

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from utils.score_against_folds import (
    find_truth_pdbs, score_predictions_vs_truth,
)


def find_af2_rank1(cluster_dir: Path, chain: str) -> Optional[Path]:
    chain_dir = cluster_dir / chain
    if not chain_dir.is_dir():
        return None
    hits = list(chain_dir.glob("*unrelaxed_rank_001*.pdb"))
    return hits[0] if hits else None


def find_af3_rank1(cluster_dir: Path, chain: str) -> Optional[Path]:
    chain_dir = cluster_dir / chain
    if not chain_dir.is_dir():
        return None
    hits = list(chain_dir.glob("tmp_*_first/tmp_*_first_rank1_model.pdb"))
    return hits[0] if hits else None


def collect_predictions(mode_dir: Path, is_af3: bool, chain_F1: str,
                         chain_F2: str) -> list[tuple[str, str, Path]]:
    """Return (cluster_id, chain, pred_pdb) tuples for every prediction
    found under mode_dir, across both chains."""
    finder = find_af3_rank1 if is_af3 else find_af2_rank1
    out = []
    for cdir in sorted(d for d in mode_dir.iterdir()
                       if d.is_dir() and d.name.startswith("ShallowMsa_")):
        cluster_tag = cdir.name
        for chain in (chain_F1, chain_F2):
            pred = finder(cdir, chain)
            if pred is not None:
                out.append((cluster_tag, chain, pred))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pair", default="5jytA_2qkeE")
    ap.add_argument("--output", default=None,
                    help="Output CSV path (default depends on --modes)")
    ap.add_argument("--modes",
                    default="AF2_treek100_full,AF2_treek100_top10,"
                            "AF3_treek100_full,AF3_treek100_top10",
                    help="Comma-separated AF output-dir suffixes to score.")
    args = ap.parse_args()

    pair_id = args.pair
    pair_dir = ROOT / "Pipeline" / "FoldPairs" / pair_id
    output_af_dir = pair_dir / "output_AF"

    # Locate truth PDBs (so we can pre-print + reuse downstream)
    (truth_F1, truth_F2), (ch1, ch2), (pdb1, pdb2) = find_truth_pdbs(pair_id, pair_dir)
    chain_F1 = f"{pdb1}{ch1}"
    chain_F2 = f"{pdb2}{ch2}"
    print(f"[score] pair={pair_id}")
    print(f"[score] truth_F1={truth_F1.name}  chain_tag={chain_F1}")
    print(f"[score] truth_F2={truth_F2.name}  chain_tag={chain_F2}")

    modes = [m.strip() for m in args.modes.split(",") if m.strip()]
    out_csv = Path(args.output) if args.output else (
        ROOT / "docs" / f"newmethod_tmscores_{pair_id}.csv"
    )

    all_dfs = []
    for mode in modes:
        mode_dir = output_af_dir / mode
        if not mode_dir.is_dir():
            print(f"[score] SKIP {mode}: dir not found ({mode_dir})")
            continue
        is_af3 = mode.startswith("AF3")
        version = "AF3" if is_af3 else "AF2"
        preds = collect_predictions(mode_dir, is_af3, chain_F1, chain_F2)
        if not preds:
            print(f"[score] {mode}: no predictions found")
            continue
        print(f"[score] === {mode}: {len(preds)} predictions ===")

        # Split per chain so TMdiff_centered is per (mode, chain) group
        # (the util groups by (mode, chain) for centering)
        cluster_ids = [c for c, _, _ in preds]
        chains = [ch for _, ch, _ in preds]
        pred_paths = [p for _, _, p in preds]

        df = score_predictions_vs_truth(
            pair_id=pair_id,
            pred_paths=pred_paths,
            cluster_ids=cluster_ids,
            mode_label=mode,
            version_label=version,
            pair_dir=pair_dir,
            chain_label=None,  # let the util default, BUT we pass per-row chain below
            repo_root=ROOT,
        )
        # Override the chain column with per-row chain tags
        df["chain"] = chains
        # Recompute TMdiff_centered now that chain is correct
        if "TMdiff" in df.columns:
            df["TMdiff_centered"] = df.groupby(["mode", "chain"])["TMdiff"].transform(
                lambda g: g - g.median()
            )
        all_dfs.append(df)

    if not all_dfs:
        print("[score] no predictions scored; nothing to write")
        return
    full = pd.concat(all_dfs, ignore_index=True)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    full.to_csv(out_csv, index=False)
    print(f"\n[score] wrote {out_csv}  ({len(full)} rows)")

    # Summary
    if "TMdiff_centered" in full.columns:
        summary = (full.dropna(subset=["TMdiff_centered"])
                       .groupby(["mode", "chain"])["TMdiff_centered"]
                       .agg(["count", "mean", "std", "min", "max"])
                       .round(4))
        print("\n[score] Per-(mode, chain) within-pair stddev of TMdiff_centered:")
        print(summary.to_string())


if __name__ == "__main__":
    main()
