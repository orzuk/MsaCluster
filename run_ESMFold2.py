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


def _read_representative_from_a3m(
    path: Path,
    method: str = "medoid",
) -> tuple[Optional[str], str]:
    """Extract the cluster's representative sequence (medoid by default).

    Delegates to utils.cluster_representative.get_cluster_representative so
    AF2/AF3/Boltz-2/ESMFold2 all share the same logic. Returns
    (header, sequence). Sequence is ungapped, uppercased.
    """
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from utils.cluster_representative import get_cluster_representative
    header, seq = get_cluster_representative(path, method=method)
    return header, seq or ""


def _load_model(device: str):
    """Load Biohub/ESMFold2 from HuggingFace. Wrapped in a function so
    the heavy import doesn't run at module import time (lets --help work
    without GPU)."""
    print(f"[esmfold2] loading biohub/ESMFold2 onto {device}...", flush=True)
    from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model
    model = ESMFold2Model.from_pretrained("biohub/ESMFold2")
    model = model.to(device).eval()
    return model


def _extract_pdb_text(obj) -> Optional[str]:
    """Try common PDB-extraction APIs on a structure-bearing object.
    Returns PDB text on success, None otherwise."""
    for attr in (
        "to_pdb_string", "pdb_string", "write_pdb_string",
        "to_pdb", "pdb",
    ):
        v = getattr(obj, attr, None)
        if v is None:
            continue
        out = v() if callable(v) else v
        if isinstance(out, str) and out.strip().startswith(("ATOM", "HEADER", "MODEL", "REMARK")):
            return out
        if isinstance(out, (bytes, bytearray)):
            return out.decode("utf-8", errors="replace")
    return None


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
    # Try direct extraction from result (older API shapes).
    pdb = _extract_pdb_text(result)
    if pdb:
        return pdb
    # MolecularComplexResult: structure lives on .complex (a MolecularComplex).
    # ESMFold2's MolecularComplex outputs mmCIF natively. The legacy ESM
    # ProteinComplex API (with PDB methods) is reachable via to_protein_complex().
    cplx = getattr(result, "complex", None)
    if cplx is not None:
        pdb = _extract_pdb_text(cplx)
        if pdb:
            return pdb
        # Convert MolecularComplex -> ProteinComplex (legacy API with PDB methods).
        to_pc = getattr(cplx, "to_protein_complex", None)
        if callable(to_pc):
            try:
                pc = to_pc()
                pdb = _extract_pdb_text(pc)
                if pdb:
                    return pdb
            except Exception as e:
                print(f"  [warn] to_protein_complex() failed: {e}", flush=True)
        # Last-resort: emit mmCIF (caller writes .cif instead of .pdb).
        to_mmcif = getattr(cplx, "to_mmcif", None)
        if callable(to_mmcif):
            mmcif = to_mmcif()
            if isinstance(mmcif, (bytes, bytearray)):
                mmcif = mmcif.decode("utf-8", errors="replace")
            if isinstance(mmcif, str) and mmcif.strip():
                return ("__MMCIF__\n" + mmcif)
        attrs = [a for a in dir(cplx) if not a.startswith("_")]
        raise RuntimeError(
            f"Could not find PDB or mmCIF extractor on result.complex "
            f"(type={type(cplx).__name__}). Public attrs: {attrs}"
        )
    # Predictions list fallback (kept from earlier API shapes).
    preds = getattr(result, "predictions", None)
    if preds:
        pdb = _extract_pdb_text(preds[0])
        if pdb:
            return pdb
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
    ap.add_argument("--representative_method", default="medoid",
                    choices=["medoid", "consensus", "first"],
                    help="How to pick one sequence per cluster (default: medoid).")
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
            out_cif = out_dir / f"{name}.cif"
            if args.skip_existing and (
                (out_pdb.is_file() and out_pdb.stat().st_size > 0)
                or (out_cif.is_file() and out_cif.stat().st_size > 0)
            ):
                total_skip += 1
                continue
            header, seq = _read_representative_from_a3m(fa, method=args.representative_method)
            if not seq:
                print(f"  [{name}] empty sequence; skip")
                total_fail += 1
                continue
            t0 = time.time()
            try:
                text = _fold(
                    model, name=header or name, sequence=seq,
                    num_loops=args.num_loops,
                    num_sampling_steps=args.num_sampling_steps,
                )
                if text.startswith("__MMCIF__\n"):
                    out_path = out_cif
                    out_path.write_text(text[len("__MMCIF__\n"):])
                else:
                    out_path = out_pdb
                    out_path.write_text(text)
                dt = time.time() - t0
                total_ok += 1
                if (i + 1) % 10 == 0 or i == 0:
                    print(f"  [{i+1}/{len(fastas)}] {name} L={len(seq)} "
                          f"OK in {dt:.1f}s -> {out_path.name}", flush=True)
            except Exception as e:
                total_fail += 1
                print(f"  [{name}] FAIL: {e}")
                traceback.print_exc(limit=2)

    print(f"\n[esmfold2] DONE  ok={total_ok}  skip={total_skip}  fail={total_fail}")


if __name__ == "__main__":
    main()
