#!/usr/bin/env python3
"""Score ESMFold2 cluster-medoid predictions against fold-1 + fold-2 truth.

ESMFold2 produces one PDB per cluster medoid:
    Pipeline/<pair>/output_esmfold2/ShallowMsa_NNN.pdb

For each cluster, TM-align the prediction against the two truth PDBs and
record TM_F1, TM_F2, TMdiff = TM_F1 - TM_F2.

Output: docs/esmfold2_tmscores_<pair>.csv with columns matching the
existing pilot scoring CSVs:
    mode, version, cluster, chain, TM_F1, TM_F2, TMdiff, TMdiff_centered,
    pred_pdb, error

Since ESMFold2 is single-sequence (only one prediction per cluster, not
per chain), `chain` column is set to a single fold label (the F1 chain
tag) for compatibility with the pilot summarizer. The TMdiff_centered
is computed per-pair (no chain stratification needed).

Usage:
    python3 scripts/score_esmfold2.py --pair 5jytA_2qkeE
    python3 scripts/score_esmfold2.py --foldpair_ids ALL
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

try:
    from utils.align_utils import tmalign_unified as _tmalign
except ImportError:
    from utils.align_utils import _run_tmalign_binary as _tmalign  # type: ignore


def find_truth_pdbs(pair_dir: Path, pair_id: str):
    """Locate the two truth PDB files. Same logic as score_treek100_pilot.py."""
    m = re.match(r"^([0-9a-zA-Z]{4})([A-Za-z])_([0-9a-zA-Z]{4})([A-Za-z])$", pair_id)
    if not m:
        raise ValueError(f"Cannot parse pair_id {pair_id}")
    pdb1, ch1, pdb2, ch2 = m.groups()
    found = []
    for pdb in (pdb1, pdb2):
        for c in [pair_dir / f"{pdb}.pdb",
                  pair_dir / f"{pdb}_cif.pdb",
                  pair_dir / "chain_pdb_files" / f"{pdb}.pdb"]:
            if c.is_file():
                found.append(c)
                break
        else:
            raise FileNotFoundError(f"Truth PDB for {pdb} not found in {pair_dir}")
    return tuple(found), (ch1, ch2), (pdb1, pdb2)


def score_one_pair(pair_id: str, mode_label: str,
                    output_dir: Path, pred_subdir: str = "output_esmfold2",
                    ) -> pd.DataFrame:
    pair_dir = ROOT / "Pipeline" / "FoldPairs" / pair_id
    preds_dir = pair_dir / pred_subdir
    if not preds_dir.is_dir():
        print(f"[score-esmfold2] {pair_id}: no {pred_subdir}/ -- skipping")
        return pd.DataFrame()

    (truth_F1, truth_F2), (ch1, ch2), (pdb1, pdb2) = find_truth_pdbs(pair_dir, pair_id)
    chain_F1 = f"{pdb1}{ch1}"
    chain_F2 = f"{pdb2}{ch2}"
    print(f"[score-esmfold2] {pair_id}: truth_F1={truth_F1.name} ({chain_F1}); "
          f"truth_F2={truth_F2.name} ({chain_F2})")

    rows = []
    pred_files = sorted(preds_dir.glob("ShallowMsa_*.pdb"))
    print(f"[score-esmfold2] {pair_id}: {len(pred_files)} predictions")
    for pred_pdb in pred_files:
        cluster_tag = pred_pdb.stem  # e.g. ShallowMsa_007
        try:
            res_F1 = _tmalign(str(pred_pdb), str(truth_F1), chain1=None, chain2=ch1)
            res_F2 = _tmalign(str(pred_pdb), str(truth_F2), chain1=None, chain2=ch2)
            tm_F1 = max(res_F1.get("tm_by_1") or 0.0, res_F1.get("tm_by_2") or 0.0)
            tm_F2 = max(res_F2.get("tm_by_1") or 0.0, res_F2.get("tm_by_2") or 0.0)
            rows.append({
                "mode":     mode_label,
                "version":  "ESMFold2",
                "cluster":  cluster_tag,
                "chain":    chain_F1,  # nominal; ESMFold2 is single-seq, no chain split
                "TM_F1":    tm_F1,
                "TM_F2":    tm_F2,
                "TMdiff":   tm_F1 - tm_F2,
                "pred_pdb": str(pred_pdb.relative_to(ROOT)),
                "error":    "",
            })
        except Exception as e:
            print(f"  [{cluster_tag}] TM-align failed: {e}")
            rows.append({
                "mode": mode_label, "version": "ESMFold2",
                "cluster": cluster_tag, "chain": chain_F1,
                "TM_F1": None, "TM_F2": None, "TMdiff": None,
                "pred_pdb": str(pred_pdb.relative_to(ROOT)),
                "error": str(e)[:120],
            })

    df = pd.DataFrame(rows)
    if not df.empty and "TMdiff" in df.columns:
        df["TMdiff_centered"] = df["TMdiff"] - df["TMdiff"].median()

    out_csv = output_dir / f"esmfold2_tmscores_{pair_id}.csv"
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"[score-esmfold2] {pair_id}: wrote {out_csv}  ({len(df)} rows)")
    if not df.empty and "TMdiff_centered" in df.columns:
        clean = df.dropna(subset=["TMdiff_centered"])
        if not clean.empty:
            print(f"  n={len(clean)}  mean={clean.TMdiff_centered.mean():.4f}  "
                  f"std={clean.TMdiff_centered.std():.4f}  "
                  f"min={clean.TMdiff_centered.min():.4f}  "
                  f"max={clean.TMdiff_centered.max():.4f}")
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pair", help="Single pair id")
    ap.add_argument("--foldpair_ids", default=None,
                    help="Comma-separated pair ids, or 'ALL'")
    ap.add_argument("--mode_label", default="ESMFold2_newmethod",
                    help="Label written to the 'mode' column (default ESMFold2_newmethod)")
    ap.add_argument("--output_dir", default=str(ROOT / "docs"),
                    help="Where to write per-pair CSVs (default: docs/)")
    ap.add_argument("--pred_subdir", default="output_esmfold2",
                    help="Sub-directory under Pipeline/<pair>/ containing the PDB outputs")
    args = ap.parse_args()

    if args.pair:
        pair_list = [args.pair]
    elif args.foldpair_ids:
        if args.foldpair_ids == "ALL":
            pair_list = sorted(p.name for p in (ROOT / "Pipeline" / "FoldPairs").glob("*/")
                               if p.name not in ("_s4pred_work", "jobs"))
        else:
            pair_list = [s.strip() for s in args.foldpair_ids.split(",") if s.strip()]
    else:
        sys.exit("Must give --pair or --foldpair_ids")

    output_dir = Path(args.output_dir)
    for pid in pair_list:
        try:
            score_one_pair(pid, args.mode_label, output_dir,
                           pred_subdir=args.pred_subdir)
        except Exception as e:
            print(f"[score-esmfold2] {pid}: ERROR {e}")


if __name__ == "__main__":
    main()
