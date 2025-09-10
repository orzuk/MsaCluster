#!/usr/bin/env bash
set -euo pipefail

# Activate your Python venv
source /sci/labs/orzuk/orzuk/af2-venv/bin/activate

# CUDA toolchain & libdevice location (from system CUDA)
CUDA_DIR="$(dirname "$(dirname "$(readlink -f "$(which ptxas)")")")"
export PATH="$CUDA_DIR/bin:$PATH"
export XLA_FLAGS="--xla_gpu_cuda_data_dir=$CUDA_DIR"


# Make directory for colabfold params
export COLABFOLD_PARAMS_DIR="${COLABFOLD_PARAMS_DIR:-$HOME/.cache/colabfold}"
mkdir -p "$COLABFOLD_PARAMS_DIR"

# JAX memory knobs + wheel libs (safe defaults)
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_MEM_FRACTION=.80
export LD_LIBRARY_PATH="/sci/labs/orzuk/orzuk/af2-venv/lib/python3.11/site-packages/nvidia/cudnn/lib:\
/sci/labs/orzuk/orzuk/af2-venv/lib/python3.11/site-packages/nvidia/cuda_nvrtc/lib:\
/sci/labs/orzuk/orzuk/af2-venv/lib/python3.11/site-packages/nvidia/cuda_runtime/lib:\
/sci/labs/orzuk/orzuk/af2-venv/lib/python3.11/site-packages/nvidia/cublas/lib:${LD_LIBRARY_PATH:-}"


# If asked to run Python (for get_msa.py), do that; otherwise default to colabfold_batch.
if [[ "${1:-}" == "--python" ]]; then
  shift
  python3 "$@"
elif [[ "${1:-}" == "python3" || "${1:-}" == "python" || "${1:-}" =~ \.py$ ]]; then
  python3 "$@"
else
  colabfold_batch "$@"
fi