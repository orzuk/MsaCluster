#!/usr/bin/env bash
# RunBioEmu2.sh -- wrapper around BioEmu2 equilibrium-ensemble sampling.
# Called from run_bioemu.py --engine bioemu2, mirroring RunBioEmu.sh.
#
# Usage:
#   RunBioEmu2.sh <cluster.a3m> <out_dir> <num_samples> [num_steps]
#
# BioEmu2 is trained from scratch (MSA module -> pairformer -> all-atom EDM
# denoiser) rather than reading frozen AlphaFold2 Evoformer embeddings the way
# BioEmu1 does, so none of BioEmu1's jax/protobuf ordering workarounds apply
# and no colabfold AF2 parameters are needed. Output: samples.xtc +
# topology.cif + samples/frame_*.pdb (written by scripts/bioemu2_sample.py).

set -euo pipefail

BIOEMU2_VENV="${BIOEMU2_VENV:-/sci/labs/orzuk/orzuk/venvs/bioemu2-venv}"
# Same lab-storage HF cache BioEmu1 uses; home has a ~2 GB quota. The
# bioemu2-preview checkpoint is pre-staged there by scripts/install_bioemu2.sbatch
# on a glacier node, because GPU nodes have no outbound internet.
export HF_HOME="${HF_HOME:-/sci/labs/orzuk/orzuk/hf_cache}"
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <cluster.a3m> <out_dir> <num_samples> [num_steps]" >&2
  exit 2
fi
A3M="$1"; OUT="$2"; N="$3"; STEPS="${4:-0}"

mkdir -p "$OUT"
. "$BIOEMU2_VENV/bin/activate"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CMD="python3 $HERE/../bioemu2_sample.py --a3m $A3M --out_dir $OUT --num_samples $N"
[[ "$STEPS" != "0" ]] && CMD="$CMD --num_steps $STEPS"

echo "[bioemu2] $CMD"
exec $CMD
