#!/usr/bin/env bash
# Pipeline/colabfold_env.sh
set -euo pipefail

# --- Optional: point to a ColabFold weights cache directory (leave empty to auto-download/cache) ---
: "${COLABFOLD_WEIGHTS:=}"   # e.g., COLABFOLD_WEIGHTS="/sci/labs/orzuk/orzuk/colabfold_weights"

# Resolve the colabfold_batch binary:
# Order:
#  1) $COLABFOLD_BIN if explicitly set and executable
#  2) Whatever is first on PATH
#  3) Common local & cluster locations
resolve_colabfold() {
  if [ -n "${COLABFOLD_BIN:-}" ] && [ -x "$COLABFOLD_BIN" ]; then return 0; fi
  if command -v colabfold_batch >/dev/null 2>&1; then
    COLABFOLD_BIN="$(command -v colabfold_batch)"; return 0; fi
  local candidates=(
    "/sci/labs/orzuk/orzuk/my-python-venv/bin/colabfold_batch"
    "$HOME/localcolabfold/colabfold-conda/bin/colabfold_batch"
    "$HOME/.local/bin/colabfold_batch"
    "/sci/labs/orzuk/orzuk/localcolabfold/colabfold-conda/bin/colabfold_batch"
  )
  for c in "${candidates[@]}"; do
    if [ -x "$c" ]; then COLABFOLD_BIN="$c"; return 0; fi
  done
  echo "ERROR: colabfold_batch not found. Install via 'pip install colabfold' or LocalColabFold," \
       "or export COLABFOLD_BIN=/full/path/to/colabfold_batch" >&2
  exit 127
}
resolve_colabfold

# --- JAX runtime knobs (safe defaults; can be overridden by caller) ---
export XLA_PYTHON_CLIENT_ALLOCATOR=platform
: "${XLA_PYTHON_CLIENT_MEM_FRACTION:=4.0}"

# --- Optional Hugging Face cache for ESM and friends ---
: "${HF_HOME:=$HOME/.cache/huggingface}"

# Build optional --data flag only if the weights dir exists
CF_DATA_FLAGS=()
if [ -n "${COLABFOLD_WEIGHTS}" ] && [ -d "${COLABFOLD_WEIGHTS}" ]; then
  CF_DATA_FLAGS=(--data="${COLABFOLD_WEIGHTS}")
fi

# Helpful log so scripts show what they're using (comment out if too chatty)
echo "[colabfold_env] using: ${COLABFOLD_BIN}"
