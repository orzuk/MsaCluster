#!/usr/bin/env python3
"""bioemu2_sample.py -- sample one BioEmu2 equilibrium ensemble from an A3M.

Runs INSIDE bioemu2-venv (python 3.12), called by scripts/shell/RunBioEmu2.sh.
BioEmu2 (microsoft/bioemu2, checkpoint bioemu2-preview) exposes no CLI for a
local MSA, so this is the thin Python entry point: build a System from the
query sequence, attach the cluster A3M as the MSA, sample, and write the frames
as PDBs where run_bioemu.py expects them.

Unlike BioEmu1 this model is all-atom native and writes topology.cif, so we
split the trajectory here (mdtraj ships with bioemu2) instead of leaving a cif
for the caller to interpret.

Usage
  python3 bioemu2_sample.py --a3m cluster.a3m --out_dir DIR --num_samples 50
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def read_query(a3m: Path) -> str:
    """First sequence of the a3m, with insertions and gaps removed.

    get_msa validates the query row against the System sequence after the same
    stripping, so the System must be built from exactly this string.
    """
    seq = []
    seen_header = False
    with open(a3m) as fh:
        for line in fh:
            if line.startswith(">"):
                if seen_header:
                    break
                seen_header = True
                continue
            seq.append(line.strip())
    q = "".join(seq)
    return "".join(c for c in q if c.isupper() and c != "-")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--a3m", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--num_samples", type=int, default=50)
    ap.add_argument("--msa_rows", type=int, default=0,
                    help="Keep at most this many MSA rows including the query "
                         "(0 = all). BioEmu2's own docs subsample with msa[:N].")
    ap.add_argument("--num_steps", type=int, default=0,
                    help="Denoiser steps (0 = model default, 200). Lowering it "
                         "trades sample quality for throughput; the per-cluster "
                         "signal we read off is ordinal.")
    args = ap.parse_args()

    a3m = Path(args.a3m)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    from bioemu2 import EmuSession, System, get_msa, load_model, save_frames

    query = read_query(a3m)
    if not query:
        print(f"[bioemu2] empty query in {a3m}", file=sys.stderr)
        return 2
    print(f"[bioemu2] {a3m.name}: L={len(query)}, n={args.num_samples}", flush=True)

    system = System.from_description(proteins=[query])
    msa = get_msa(system, source=a3m)
    if args.msa_rows and args.msa_rows > 0:
        msa = msa[: args.msa_rows]

    emu = EmuSession(system, model=load_model(), mode="equilibrium", msa=msa)
    kw = {"num_samples": args.num_samples}
    if args.num_steps:
        kw["num_steps"] = args.num_steps
    samples = emu.sample(**kw)
    save_frames(samples, out_dir)

    # Per-frame PDBs, the format utils.score_against_folds TM-scores.
    frames_dir = out_dir / "samples"
    frames_dir.mkdir(exist_ok=True)
    for i in range(samples.n_frames):
        samples[i].save_pdb(str(frames_dir / f"frame_{i:04d}.pdb"))
    print(f"[bioemu2] wrote {samples.n_frames} frames to {frames_dir}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
