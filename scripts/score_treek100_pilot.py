#!/usr/bin/env python3
"""
Score the K=100 + medoid-query AF2/AF3 pilot for the KaiB pair
(5jytA_2qkeE).

For each of the 4 result sets:
  output_AF/AF2_treek100_full/
  output_AF/AF2_treek100_top10/
  output_AF/AF3_treek100_full/
  output_AF/AF3_treek100_top10/

and for each (cluster, chain), find the rank-1 PDB and TM-align it against
the two ground-truth structures (5jyt and 2qke), record TM_F1, TM_F2,
TMdiff = TM_F1 - TM_F2, and TMdiff_centered = TMdiff - median(TMdiff per
mode-version-chain).

Output: docs/treek100_pilot_tmscores.csv

Usage:
    python3 scripts/score_treek100_pilot.py
    python3 scripts/score_treek100_pilot.py --pair 5jytA_2qkeE
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from utils.align_utils import _run_tmalign_binary


def find_af2_rank1(cluster_dir: Path, chain: str) -> Optional[Path]:
    """Find AF2 rank-1 PDB. Pattern:
       <cluster_dir>/<chain>/*unrelaxed_rank_001*.pdb
    """
    chain_dir = cluster_dir / chain
    if not chain_dir.is_dir():
        return None
    hits = list(chain_dir.glob("*unrelaxed_rank_001*.pdb"))
    return hits[0] if hits else None


def find_af3_rank1(cluster_dir: Path, chain: str) -> Optional[Path]:
    """Find AF3 rank-1 PDB. Pattern:
       <cluster_dir>/<chain>/tmp_<cluster_tag>__<chain>_first/
                              tmp_<cluster_tag>__<chain>_first_rank1_model.pdb
    """
    chain_dir = cluster_dir / chain
    if not chain_dir.is_dir():
        return None
    # The "first" subdir name uses the cluster's full tag
    candidates = list(chain_dir.glob("tmp_*_first/tmp_*_first_rank1_model.pdb"))
    return candidates[0] if candidates else None


def find_truth_pdbs(pair_dir: Path, pair_id: str) -> tuple[Path, Path]:
    """Locate the two truth PDB files for the pair.
    Tries: <pair_dir>/<pdb1>.pdb, <pair_dir>/chain_pdb_files/<pdb1>.pdb
    """
    m = re.match(r"^([0-9a-zA-Z]{4})([A-Za-z])_([0-9a-zA-Z]{4})([A-Za-z])$", pair_id)
    if not m:
        raise ValueError(f"Cannot parse pair_id {pair_id}")
    pdb1, ch1, pdb2, ch2 = m.groups()
    candidates_per_pdb = []
    for pdb in (pdb1, pdb2):
        found = None
        for c in [
            pair_dir / f"{pdb}.pdb",
            pair_dir / f"{pdb}_cif.pdb",
            pair_dir / "chain_pdb_files" / f"{pdb}.pdb",
        ]:
            if c.is_file():
                found = c
                break
        if found is None:
            raise FileNotFoundError(f"Truth PDB for {pdb} not found in {pair_dir}")
        candidates_per_pdb.append(found)
    return tuple(candidates_per_pdb), (ch1, ch2), (pdb1, pdb2)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pair", default="5jytA_2qkeE")
    ap.add_argument("--output", default=str(ROOT / "docs" / "treek100_pilot_tmscores.csv"))
    ap.add_argument("--modes", default="AF2_treek100_full,AF2_treek100_top10,AF3_treek100_full,AF3_treek100_top10",
                    help="Comma-separated AF output-dir suffixes to score.")
    ap.add_argument("--checkpoint-every", type=int, default=20,
                    help="Save partial CSV every N rows (default 20)")
    args = ap.parse_args()

    pair = args.pair
    pair_dir = ROOT / "Pipeline" / "FoldPairs" / pair
    output_af_dir = pair_dir / "output_AF"

    (truth_F1, truth_F2), (ch1, ch2), (pdb1, pdb2) = find_truth_pdbs(pair_dir, pair)
    chain_F1 = f"{pdb1}{ch1}"  # e.g. 5jytA
    chain_F2 = f"{pdb2}{ch2}"  # e.g. 2qkeE
    print(f"[score] pair={pair}")
    print(f"[score] truth_F1={truth_F1}  (chain_tag={chain_F1})")
    print(f"[score] truth_F2={truth_F2}  (chain_tag={chain_F2})")

    modes = [m.strip() for m in args.modes.split(",") if m.strip()]
    rows = []
    out_csv = Path(args.output)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    for mode in modes:
        mode_dir = output_af_dir / mode
        if not mode_dir.is_dir():
            print(f"[score] SKIP {mode}: dir not found ({mode_dir})")
            continue
        is_af3 = mode.startswith("AF3")
        finder = find_af3_rank1 if is_af3 else find_af2_rank1
        cluster_dirs = sorted([d for d in mode_dir.iterdir()
                               if d.is_dir() and d.name.startswith("ShallowMsa_")])
        print(f"[score] === {mode}: {len(cluster_dirs)} clusters ===")
        for ci, cdir in enumerate(cluster_dirs):
            cluster_tag = cdir.name  # ShallowMsa_NNN
            for chain_tag in (chain_F1, chain_F2):
                pred_pdb = finder(cdir, chain_tag)
                if pred_pdb is None:
                    print(f"  [{mode}/{cluster_tag}/{chain_tag}] no rank-1 PDB; skipping")
                    rows.append({
                        "mode": mode, "version": "AF3" if is_af3 else "AF2",
                        "cluster": cluster_tag, "chain": chain_tag,
                        "TM_F1": None, "TM_F2": None, "TMdiff": None,
                        "pred_pdb": "", "error": "no_rank1_pdb",
                    })
                    continue
                try:
                    # TM-align uses the first chain in the file if no chain specified
                    res_F1 = _run_tmalign_binary(str(pred_pdb), str(truth_F1), None, ch1)
                    res_F2 = _run_tmalign_binary(str(pred_pdb), str(truth_F2), None, ch2)
                    # Use max of by-1 and by-2 normalization, matching the
                    # rest of the pipeline convention
                    tm_F1 = max(res_F1.get("tm_by_1") or 0.0,
                                res_F1.get("tm_by_2") or 0.0)
                    tm_F2 = max(res_F2.get("tm_by_1") or 0.0,
                                res_F2.get("tm_by_2") or 0.0)
                    rows.append({
                        "mode": mode, "version": "AF3" if is_af3 else "AF2",
                        "cluster": cluster_tag, "chain": chain_tag,
                        "TM_F1": tm_F1, "TM_F2": tm_F2,
                        "TMdiff": tm_F1 - tm_F2,
                        "pred_pdb": str(pred_pdb.relative_to(ROOT)),
                        "error": "",
                    })
                except Exception as e:
                    print(f"  [{mode}/{cluster_tag}/{chain_tag}] TM-align failed: {e}")
                    rows.append({
                        "mode": mode, "version": "AF3" if is_af3 else "AF2",
                        "cluster": cluster_tag, "chain": chain_tag,
                        "TM_F1": None, "TM_F2": None, "TMdiff": None,
                        "pred_pdb": str(pred_pdb.relative_to(ROOT)),
                        "error": str(e)[:100],
                    })
            if (ci + 1) % args.checkpoint_every == 0:
                pd.DataFrame(rows).to_csv(out_csv, index=False)
                print(f"  [{mode}] checkpoint {ci+1}/{len(cluster_dirs)} -> {out_csv}")

    df = pd.DataFrame(rows)

    # Compute TMdiff_centered per (mode, chain): subtract within-group median
    if not df.empty and "TMdiff" in df.columns:
        df["TMdiff_centered"] = df.groupby(["mode", "chain"])["TMdiff"].transform(
            lambda g: g - g.median()
        )

    df.to_csv(out_csv, index=False)
    print(f"\n[score] wrote {out_csv}  ({len(df)} rows)")

    # Summary
    print("\n[score] Per-(mode, chain) within-pair stddev of TMdiff_centered:")
    if "TMdiff_centered" in df.columns:
        summary = (df.dropna(subset=["TMdiff_centered"])
                   .groupby(["mode", "chain"])["TMdiff_centered"]
                   .agg(["count", "mean", "std", "min", "max"])
                   .round(4))
        print(summary.to_string())


if __name__ == "__main__":
    main()
