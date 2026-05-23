#!/usr/bin/env bash
# RunBoltz2.sh -- thin wrapper around Boltz-2 inference.
# Designed to be called from run_Boltz2.py. User must edit the variables
# below after installing boltz in torch2-venv.
#
# Usage:
#   RunBoltz2.sh <input.yaml_or_fasta> <out_dir> [extra flags...]
#
# Boltz-2 accepts a YAML input that can include protein chains, ligands,
# nucleic acid partners, etc. For our use case (per-cluster fold
# prediction) we pass a simple single-chain YAML. To include binding
# partners (e.g. RfaH-RNAP fold-switching), edit the YAML accordingly --
# Boltz-2 supports protein + protein, protein + ligand, protein + RNA.

set -euo pipefail

# --- CONFIG (edit after install) ----------------------------------------
BOLTZ_VENV="${BOLTZ_VENV:-/sci/labs/orzuk/orzuk/venvs/torch2-venv}"
BOLTZ_CACHE="${BOLTZ_CACHE:-/sci/labs/orzuk/orzuk/boltz_cache}"
# ------------------------------------------------------------------------

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <input.yaml_or_fasta> <out_dir> [extra flags...]" >&2
  exit 2
fi
INP="$1"; OUT="$2"; shift 2

mkdir -p "$OUT" "$BOLTZ_CACHE"

. "$BOLTZ_VENV/bin/activate"

# NOTE: --use_msa_server is intentionally NOT passed. The caller
# (run_Boltz2.py) embeds the cluster-specific ShallowMsa A3M into the
# YAML (`msa: <path>` per chain), which is the whole point of per-cluster
# Boltz-2 prediction. With --use_msa_server, Boltz would refetch a
# generic colabfold MSA and wash out the cluster signal. To force a
# server fetch (e.g. for an apo single-sequence baseline), re-add
# --use_msa_server below or pass it via the positional `$*` extras.
CMD="boltz predict \
    $INP \
    --out_dir $OUT \
    --cache $BOLTZ_CACHE \
    --output_format pdb \
    $*"

echo "[boltz2] $CMD"
exec $CMD
