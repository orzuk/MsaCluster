#!/usr/bin/env python3
"""Score ESMFold2 per-cluster medoid predictions vs the two truth folds.

Thin wrapper over utils.score_against_folds.score_predictions_vs_truth.
ESMFold2 emits ONE PDB per cluster (not per-chain-per-rank like AF), so
enumeration is straightforward.

Usage:
    python3 scripts/score_esmfold2.py --pair 5jytA_2qkeE
    python3 scripts/score_esmfold2.py --foldpair_ids ALL
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from utils.score_against_folds import score_predictions_vs_truth


def score_one_pair(pair_id: str, mode_label: str, output_dir: Path,
                    pred_subdir: str = "output_esmfold2") -> None:
    pair_dir = ROOT / "Pipeline" / "FoldPairs" / pair_id
    preds_dir = pair_dir / pred_subdir
    if not preds_dir.is_dir():
        print(f"[score-esmfold2] {pair_id}: no {pred_subdir}/ -- skip")
        return

    # ESMFold2 writes .cif (mmCIF) when the result has multi-char chain IDs.
    # TMalign's auto-chain-detection misfires on mmCIF, so convert .cif -> .pdb
    # next to the original (cached; idempotent) and score the .pdb version.
    cif_files = sorted(preds_dir.glob("ShallowMsa_*.cif"))
    if cif_files:
        try:
            from Bio.PDB.MMCIFParser import MMCIFParser
            from Bio.PDB.PDBIO import PDBIO
            parser = MMCIFParser(QUIET=True)
            io = PDBIO()
            n_conv = 0
            for cif in cif_files:
                pdb = cif.with_suffix(".pdb")
                if pdb.is_file() and pdb.stat().st_size > 0:
                    continue
                try:
                    s = parser.get_structure("s", str(cif))
                    io.set_structure(s)
                    io.save(str(pdb))
                    n_conv += 1
                except Exception as e:
                    print(f"[score-esmfold2] WARN: .cif->.pdb failed for {cif.name}: {e}")
            if n_conv > 0:
                print(f"[score-esmfold2] {pair_id}: converted {n_conv} .cif -> .pdb")
        except ImportError:
            print("[score-esmfold2] WARN: biopython not available; falling back to .cif (may NaN)")

    pred_files = sorted(preds_dir.glob("ShallowMsa_*.pdb"))
    if not pred_files:
        # Fallback: if .pdb didn't materialize (e.g., biopython missing),
        # use .cif files directly (will likely produce NaN TM-scores)
        pred_files = sorted(preds_dir.glob("ShallowMsa_*.cif"))
    if not pred_files:
        print(f"[score-esmfold2] {pair_id}: no predictions in {preds_dir}")
        return
    cluster_ids = [p.stem for p in pred_files]
    print(f"[score-esmfold2] {pair_id}: {len(pred_files)} predictions")

    df = score_predictions_vs_truth(
        pair_id=pair_id,
        pred_paths=pred_files,
        cluster_ids=cluster_ids,
        mode_label=mode_label,
        version_label="ESMFold2",
        pair_dir=pair_dir,
        repo_root=ROOT,
    )

    out_csv = output_dir / f"esmfold2_tmscores_{pair_id}.csv"
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"[score-esmfold2] {pair_id}: wrote {out_csv}  ({len(df)} rows)")
    if "TMdiff_centered" in df.columns:
        clean = df.dropna(subset=["TMdiff_centered"])
        if not clean.empty:
            print(f"  n={len(clean)}  mean={clean.TMdiff_centered.mean():.4f}  "
                  f"std={clean.TMdiff_centered.std():.4f}  "
                  f"min={clean.TMdiff_centered.min():.4f}  "
                  f"max={clean.TMdiff_centered.max():.4f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pair", help="Single pair id")
    ap.add_argument("--foldpair_ids", default=None,
                    help="Comma-separated pair ids, or 'ALL'")
    ap.add_argument("--mode_label", default="ESMFold2_newmethod")
    ap.add_argument("--output_dir", default=str(ROOT / "docs"))
    ap.add_argument("--pred_subdir", default="output_esmfold2")
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
