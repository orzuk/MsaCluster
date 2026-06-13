#!/usr/bin/env bash
# RunBioEmu.sh -- thin wrapper around BioEmu equilibrium-ensemble sampling.
# Called from run_bioemu.py. Edit the variables below after installing bioemu
# in a DEDICATED venv (see INSTALL_bioemu.md -- its torch/AF-evoformer stack
# conflicts with my-python-venv / af-venv / esmfold2-venv / torch2-venv).
#
# Usage:
#   RunBioEmu.sh <cluster.a3m> <out_dir> <num_samples>
#
# BioEmu samples the approximate equilibrium ensemble of ONE chain. We pass the
# cluster's ShallowMsa A3M as the input so BioEmu builds its AF2-evoformer
# embedding from the SAME per-cluster MSA that AF2/AF3/Boltz-2 use (rather than
# fetching a generic colabfold MSA) -- that is the whole point of a per-cluster
# prediction. Output: topology.pdb + samples.xtc (split into per-frame PDBs by
# run_bioemu.py via mdtraj).

set -euo pipefail

# --- CONFIG (edit after install) ----------------------------------------
BIOEMU_VENV="${BIOEMU_VENV:-/sci/labs/orzuk/orzuk/venvs/bioemu-venv}"
BIOEMU_CACHE="${BIOEMU_CACHE:-/sci/labs/orzuk/orzuk/bioemu_cache}"
# ------------------------------------------------------------------------

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <cluster.a3m> <out_dir> <num_samples>" >&2
  exit 2
fi
A3M="$1"; OUT="$2"; N="$3"

mkdir -p "$OUT" "$BIOEMU_CACHE"

. "$BIOEMU_VENV/bin/activate"

# jaxlib + tensorflow both register xla_data.proto into protobuf's C++ pool;
# jaxlib aborts when it is the SECOND registrant. sitecustomize.py in jaxfirst/
# forces jax to init first (see that file). Put it at the FRONT of PYTHONPATH.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="$HERE/jaxfirst${PYTHONPATH:+:$PYTHONPATH}"
export TF_CPP_MIN_LOG_LEVEL=3

# Home has ~2 GB quota; keep ALL caches in lab storage. bioemu downloads AF2
# params (~3.5 GB) to ~/.cache/colabfold (a symlink to lab) -- make sure the
# symlink target exists so the download doesn't fail on a dangling link.
BIOEMU_SO3_CACHE="${BIOEMU_SO3_CACHE:-$BIOEMU_CACHE/so3}"
mkdir -p "$BIOEMU_SO3_CACHE"
CF_DIR="$(readlink -f "$HOME/.cache/colabfold" 2>/dev/null || true)"
[ -n "$CF_DIR" ] && mkdir -p "$CF_DIR/params"

# BioEmu CLI: `python -m bioemu.sample`. The exact flag names can vary across
# bioemu releases -- verify with `python -m bioemu.sample --help` after install
# and adjust here. As of the public release it accepts a sequence OR an a3m via
# --sequence, the sample count via --num_samples, and an output dir via
# --output_dir. We pass the A3M so the embedding uses our cluster MSA.
CMD="python -m bioemu.sample \
    --sequence $A3M \
    --num_samples $N \
    --output_dir $OUT \
    --cache_embeds_dir $BIOEMU_CACHE \
    --cache_so3_dir $BIOEMU_SO3_CACHE"

echo "[bioemu] $CMD"
exec $CMD
