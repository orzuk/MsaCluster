#!/usr/bin/env python3
"""
Run ESMFold2 (Biohub, May 2026) on cluster medoid sequences.

For each ShallowMsa_NNN.a3m in --input_dir, extract the first sequence
(which is the cluster medoid after the --query_type=medoid patch in
run_foldswitch_pipeline.py), run ESMFold2 in single-sequence mode, and
write the predicted structure to <output_dir>/ShallowMsa_NNN.pdb.

This is the v2 paper's single-sequence baseline alongside the
AF2/AF3 + top-10 MSA pilot, replacing the older ESM-Fold (Meta 2022)
baseline.

Usage (single pair):
    python3 run_ESMFold2.py \
        --input_dir Pipeline/FoldPairs/5jytA_2qkeE/output_msa_cluster \
        --output_dir Pipeline/FoldPairs/5jytA_2qkeE/output_esmfold2

Usage (driven by run_foldswitch_pipeline.py-style ids):
    python3 run_ESMFold2.py --foldpair_ids 5jytA_2qkeE
"""
from __future__ import annotations

import argparse
import os
import sys
import time
import traceback
from pathlib import Path
from typing import Optional


def _read_first_seq_from_a3m(path: Path) -> tuple[Optional[str], str]:
    """Return (header, sequence) for the first record in an a3m.

    a3m files use uppercase for match-state columns and lowercase for
    insertions; ESMFold2 wants a clean uppercase sequence with gaps
    stripped (matching what the cluster medoid would look like as a
    single sequence input).
    """
    header: Optional[str] = None
    chunks: list[str] = []
    with path.open() as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    break  # finished the first record
                header = line[1:].split()[0]
                continue
            if header is None:
                continue
            chunks.append(line)
    seq = "".join(chunks).upper().replace("-", "").replace(".", "")
    return header, seq


def _load_model(device: str):
    """Load Biohub/ESMFold2 from HuggingFace. Wrapped in a function so
    the heavy import doesn't run at module import time (lets --help work
    without GPU)."""
    print(f"[esmfold2] loading biohub/ESMFold2 onto {device}...", flush=True)
    from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model
    model = ESMFold2Model.from_pretrained("biohub/ESMFold2")
    model = model.to(device).eval()
    return model


def _fold(model, name: str, sequence: str,
          num_loops: int, num_sampling_steps: int) -> str:
    """Run one fold; return PDB text."""
    from esm.models.esmfold2 import (
        ESMFold2InputBuilder, ProteinInput, StructurePredictionInput,
    )
    spi = StructurePredictionInput(
        sequences=[ProteinInput(id=name, sequence=sequence)]
    )
    result = ESMFold2InputBuilder().fold(
        model, spi,
        num_loops=num_loops,
        num_sampling_steps=num_sampling_steps,
    )
    # The result object exposes PDB output in one of a few shapes
    # depending on the package version; try the common ones.
    for attr in ("to_pdb", "pdb"):
        v = getattr(result, attr, None)
        if v is None:
            continue
        return v() if callable(v) else v
    preds = getattr(result, "predictions", None)
    if preds:
        first = preds[0]
        for attr in ("to_pdb", "pdb"):
            v = getattr(first, attr, None)
            if v is None:
                continue
            return v() if callable(v) else v
    raise RuntimeError(
        f"Cannot extract PDB from result of type {type(result).__name__}; "
        f"attributes: {dir(result)[:40]}..."
    )


def _resolve_input_dirs(args) -> list[tuple[str, Path, Path]]:
    """Return list of (label, input_dir, output_dir) to process."""
    if args.input_dir and args.output_dir:
        return [("(custom)", Path(args.input_dir), Path(args.output_dir))]
    if args.foldpair_ids:
        repo = Path(__file__).resolve().parent
        out = []
        for pid in args.foldpair_ids.split(","):
            pid = pid.strip()
            if not pid:
                continue
            pair_dir = repo / "Pipeline" / "FoldPairs" / pid
            inp = pair_dir / "output_msa_cluster"
            outp = pair_dir / "output_esmfold2"
            out.append((pid, inp, outp))
        return out
    raise SystemExit("Must give --input_dir+--output_dir or --foldpair_ids")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--input_dir",
                    help="Dir of ShallowMsa_*.a3m files")
    ap.add_argument("--output_dir",
                    help="Dir to write ShallowMsa_*.pdb files")
    ap.add_argument("--foldpair_ids",
                    help="Comma-separated foldpair ids; replaces "
                         "--input_dir/--output_dir")
    ap.add_argument("--pattern", default="ShallowMsa_*.a3m")
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--num_loops", type=int, default=3)
    ap.add_argument("--num_sampling_steps", type=int, default=50)
    ap.add_argument("--skip_existing", action="store_true",
                    help="Skip clusters whose .pdb already exists")
    ap.add_argument("--max_clusters", type=int, default=None,
                    help="Process only the first N clusters (debug)")
    args = ap.parse_args()

    targets = _resolve_input_dirs(args)
    print(f"[esmfold2] {len(targets)} target dir(s)")

    model = _load_model(args.device)

    total_ok, total_skip, total_fail = 0, 0, 0
    for label, in_dir, out_dir in targets:
        if not in_dir.is_dir():
            print(f"[esmfold2] SKIP {label}: input dir not found ({in_dir})")
            continue
        out_dir.mkdir(parents=True, exist_ok=True)
        fastas = sorted(in_dir.glob(args.pattern))
        if args.max_clusters is not None:
            fastas = fastas[:args.max_clusters]
        print(f"[esmfold2] === {label}: {len(fastas)} input(s) "
              f"-> {out_dir} ===")
        for i, fa in enumerate(fastas):
            name = fa.stem  # ShallowMsa_NNN
            out_pdb = out_dir / f"{name}.pdb"
            if args.skip_existing and out_pdb.is_file() and out_pdb.stat().st_size > 0:
                total_skip += 1
                continue
            header, seq = _read_first_seq_from_a3m(fa)
            if not seq:
                print(f"  [{name}] empty sequence; skip")
                total_fail += 1
                continue
            t0 = time.time()
            try:
                pdb_text = _fold(
                    model, name=header or name, sequence=seq,
                    num_loops=args.num_loops,
                    num_sampling_steps=args.num_sampling_steps,
                )
                out_pdb.write_text(pdb_text)
                dt = time.time() - t0
                total_ok += 1
                if (i + 1) % 10 == 0 or i == 0:
                    print(f"  [{i+1}/{len(fastas)}] {name} L={len(seq)} "
                          f"OK in {dt:.1f}s -> {out_pdb.name}", flush=True)
            except Exception as e:
                total_fail += 1
                print(f"  [{name}] FAIL: {e}")
                traceback.print_exc(limit=2)

    print(f"\n[esmfold2] DONE  ok={total_ok}  skip={total_skip}  fail={total_fail}")


if __name__ == "__main__":
    main()
