#!/usr/bin/env bash
# RunAlphaFlow.sh -- thin wrapper around AlphaFlow inference.
# Designed to be called from run_AlphaFlow.py. User must edit the variables
# below after installing AlphaFlow in torch2-venv to point at the correct
# repo / weights / command.
#
# Usage:
#   RunAlphaFlow.sh <input.csv> <out_dir> [samples] [mode]
# Args:
#   input.csv   -- two-column CSV (name, sequence) following AlphaFlow's
#                   predict.py convention
#   out_dir     -- destination directory for predicted PDBs
#   samples     -- number of ensemble samples to draw (default 8)
#   mode        -- AlphaFlow inference mode (default esmflow_md_distilled)

set -euo pipefail

# --- CONFIG (edit after install) ----------------------------------------
ALPHAFLOW_VENV="${ALPHAFLOW_VENV:-/sci/labs/orzuk/orzuk/venvs/torch2-venv}"
ALPHAFLOW_REPO="${ALPHAFLOW_REPO:-/sci/labs/orzuk/orzuk/github/alphaflow}"
ALPHAFLOW_WEIGHTS="${ALPHAFLOW_WEIGHTS:-$ALPHAFLOW_REPO/weights/esmflow_md_distilled_202402.pt}"
ALPHAFLOW_MODE_DEFAULT="${ALPHAFLOW_MODE_DEFAULT:-esmflow_md_distilled}"
ALPHAFLOW_NO_DIFFUSION="${ALPHAFLOW_NO_DIFFUSION:-1}"     # 1 = pass --no_diffusion
# ------------------------------------------------------------------------

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <input.csv> <out_dir> [samples] [mode]" >&2
  exit 2
fi
INP="$1"; OUT="$2"; shift 2
SAMPLES="${1:-8}"; shift || true
MODE="${1:-$ALPHAFLOW_MODE_DEFAULT}"

mkdir -p "$OUT"

# Use POSIX `.` (not `source`) for /bin/sh-compatible sbatch wraps
. "$ALPHAFLOW_VENV/bin/activate"

# Build the command. The exact predict.py flags will differ slightly across
# AlphaFlow versions; this matches the recipe in the README at time of
# writing (2024). User can override via $ALPHAFLOW_CMD if needed.
EXTRA=""
if [[ "$ALPHAFLOW_NO_DIFFUSION" == "1" ]]; then
  EXTRA="--no_diffusion"
fi

CMD="python3 $ALPHAFLOW_REPO/predict.py \
    --mode $MODE \
    --input_csv $INP \
    --weights $ALPHAFLOW_WEIGHTS \
    --samples $SAMPLES \
    --outpdb $OUT \
    $EXTRA"

echo "[alphaflow] $CMD"
exec $CMD
