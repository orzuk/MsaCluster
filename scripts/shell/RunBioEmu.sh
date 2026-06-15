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
# HF_HOME must point at lab storage (home has ~2 GB quota) AND match where the
# bioemu checkpoint was first downloaded, so per-cluster runs hit the local
# cache instead of re-downloading from the Hub (compute nodes can be offline).
# Hard-set (do NOT defer to an inherited HF_HOME): the bioemu checkpoint was
# cached under hf_cache, and a stray HF_HOME from the login shell/another venv
# (e.g. .../tmp/huggingface) points the cache lookup at the wrong dir -> the
# --ckpt_path bypass comes up empty -> bioemu tries to download -> offline fail.
export HF_HOME="/sci/labs/orzuk/orzuk/hf_cache"
# GPU compute nodes have no outbound internet, so hf_hub_download's HEAD call
# fails and raises LocalEntryNotFoundError even though the checkpoint IS cached.
# Force offline mode: use the cached snapshot, never hit the network. All
# required files are pre-cached by the smoke test (HF checkpoint+config; AF2
# params live in the colabfold dir; the MSA is the local cluster a3m).
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
# ------------------------------------------------------------------------

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <cluster.a3m> <out_dir> <num_samples> [batch_size_100]" >&2
  exit 2
fi
A3M="$1"; OUT="$2"; N="$3"; BS100="${4:-}"   # BS100 optional (see run_bioemu.py)

mkdir -p "$OUT" "$BIOEMU_CACHE"

. "$BIOEMU_VENV/bin/activate"

# jaxlib + tensorflow both register xla_data.proto into protobuf's C++ pool;
# jaxlib aborts when it is the SECOND registrant. sitecustomize.py in jaxfirst/
# forces jax to init first (see that file). Put it at the FRONT of PYTHONPATH.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="$HERE/jaxfirst${PYTHONPATH:+:$PYTHONPATH}"
export TF_CPP_MIN_LOG_LEVEL=3

# By default jax pre-allocates ~75% of the GPU on init, which starves co-scheduled
# bioemu jobs (mass 5s OOM failures in the sweep). PREALLOCATE=false makes jax
# grab memory ON DEMAND, so small jobs use little and coexist on a shared card --
# that alone fixes the contention. Do NOT hard-cap MEM_FRACTION low: long chains
# need quadratically more memory, and a tight cap made them OOM. Leave the ceiling
# at the jax default so big proteins can use the card when they have it.
export XLA_PYTHON_CLIENT_PREALLOCATE=false

# Home has ~2 GB quota; keep ALL caches in lab storage. bioemu downloads AF2
# params (~3.5 GB) to ~/.cache/colabfold (a symlink to lab) -- make sure the
# symlink target exists so the download doesn't fail on a dangling link.
BIOEMU_SO3_CACHE="${BIOEMU_SO3_CACHE:-$BIOEMU_CACHE/so3}"
mkdir -p "$BIOEMU_SO3_CACHE"
CF_DIR="$(readlink -f "$HOME/.cache/colabfold" 2>/dev/null || true)"
[ -n "$CF_DIR" ] && mkdir -p "$CF_DIR/params"

# Pass the cached checkpoint + config EXPLICITLY so bioemu's
# maybe_download_checkpoint() short-circuits (returns the given paths) instead
# of calling hf_hub_download -- which fails on the offline GPU nodes even though
# the files are cached. Resolve them from the HF cache snapshot dir.
BIOEMU_MODEL="${BIOEMU_MODEL:-bioemu-v1.1}"
HF_BIOEMU="$HF_HOME/hub/models--microsoft--bioemu"
# Resolve via nullglob arrays (NOT `ls | head`, which returns non-zero under
# `set -euo pipefail` and would abort the script with exit 2).
shopt -s nullglob
_ckpts=( "$HF_BIOEMU"/snapshots/*/checkpoints/"$BIOEMU_MODEL"/checkpoint.ckpt )
_cfgs=(  "$HF_BIOEMU"/snapshots/*/checkpoints/"$BIOEMU_MODEL"/config.yaml )
shopt -u nullglob
CKPT="${_ckpts[0]:-}"
CFG="${_cfgs[0]:-}"
CKPT_ARGS=""
if [[ -f "$CKPT" && -f "$CFG" ]]; then
    CKPT_ARGS="--ckpt_path $CKPT --model_config_path $CFG"
    echo "[bioemu] using cached checkpoint: $CKPT"
else
    echo "[bioemu] WARNING: cached checkpoint/config not found under $HF_BIOEMU;" \
         "bioemu will try to download (will fail on an offline node)." >&2
fi

# BioEmu CLI: `python -m bioemu.sample`. The exact flag names can vary across
# bioemu releases -- verify with `python -m bioemu.sample --help` after install
# and adjust here. As of the public release it accepts a sequence OR an a3m via
# --sequence, the sample count via --num_samples, and an output dir via
# --output_dir. We pass the A3M so the embedding uses our cluster MSA.
# batch_size_100 (optional 4th arg): run_bioemu.py scales it up for long chains
# so BioEmu's int(batch_size_100*(100/L)^2) batch never underflows to 0.
BS_ARG=""
[[ -n "$BS100" ]] && BS_ARG="--batch_size_100 $BS100"
# --filter_samples False: BioEmu's physicality filter crashes the whole cluster
# with `ValueError: Invalid suffix '_unphysical.xtc'` (buggy with_suffix() call)
# whenever it has unphysical samples to write -- the #1 cause of missing clusters
# (1408 of them). Disabling it keeps the full ensemble; clashy samples just
# TM-score low and don't bias the population fraction.
CMD="python -m bioemu.sample \
    --sequence $A3M \
    --num_samples $N \
    --output_dir $OUT \
    --cache_embeds_dir $BIOEMU_CACHE \
    --cache_so3_dir $BIOEMU_SO3_CACHE \
    --filter_samples=False \
    $CKPT_ARGS \
    $BS_ARG"

echo "[bioemu] $CMD"
exec $CMD
