#!/usr/bin/env bash
set -euo pipefail

# Activate your Python venv
source /sci/labs/orzuk/orzuk/af2-venv/bin/activate


# --- ensure OUTDIR exists if user passed: <input> <outdir> [extra args]
INP=""
OUT=""
if [[ $# -ge 2 && -e "$1" ]]; then
  INP="$1"
  OUT="$2"
  # If second arg looks like a directory, create it
  if [[ -n "$OUT" ]]; then
    mkdir -p "$OUT"
  fi
fi


# --- CUDA toolchain discovery (robust) ---
if [[ -n "${CUDA_DIR_OVERRIDE:-}" && -x "${CUDA_DIR_OVERRIDE}/bin/ptxas" ]]; then
  CUDA_DIR="${CUDA_DIR_OVERRIDE}"
elif command -v ptxas >/dev/null 2>&1; then
  CUDA_DIR="$(dirname "$(dirname "$(readlink -f "$(command -v ptxas)")")")"
elif [[ -x "/usr/local/APP/nvidia/cuda/12.4.1/bin/ptxas" ]]; then
  CUDA_DIR="/usr/local/APP/nvidia/cuda/12.4.1"
elif [[ -n "${CUDA_HOME:-}" && -x "${CUDA_HOME}/bin/ptxas" ]]; then
  CUDA_DIR="${CUDA_HOME}"
else
  echo "[fatal] Couldn't locate ptxas. Set CUDA_DIR_OVERRIDE=/path/to/cuda" >&2
  exit 2
fi
export PATH="${CUDA_DIR}/bin:${PATH}"
export XLA_FLAGS="--xla_gpu_cuda_data_dir=${CUDA_DIR}"



# --- if no GPU is visible and we're not already inside a SLURM job, relaunch under srun
if [[ "${USE_SRUN:-1}" == "1" && -z "${CUDA_VISIBLE_DEVICES:-}" && -z "${SLURM_JOB_ID:-}" ]]; then
  exec srun --gres=gpu:a100:1 --time=2:00:00 --mem=16G --cpus-per-task=4 --pty bash -lc \
    "source /sci/labs/orzuk/orzuk/af2-venv/bin/activate; CUDA_DIR=\"$(dirname "$(dirname "$(readlink -f "$(which ptxas)")")")\"; export PATH=\"\$CUDA_DIR/bin:\$PATH\"; export XLA_FLAGS=\"--xla_gpu_cuda_data_dir=\$CUDA_DIR\"; $(printf '%q ' "$0" "$@")"
fi

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


echo "[diag] nvidia-smi:"
nvidia-smi || true
echo "[diag] which ptxas: $(command -v ptxas || echo 'not found')"
echo "[diag] ptxas --version:"; (ptxas --version || true)
echo "[diag] CUDA_DIR=${CUDA_DIR}"

# If asked to run Python (for get_msa.py), do that; otherwise default to colabfold_batch.
if [[ "${1:-}" == "--python" ]]; then
  shift
  python3 "$@"
elif [[ "${1:-}" == "python3" || "${1:-}" == "python" || "${1:-}" =~ \.py$ ]]; then
  python3 "$@"
else
  # If user passed INP+OUT, feed them explicitly; otherwise fall back to raw args
  if [[ -n "$INP" && -n "$OUT" ]]; then
    echo "[run] colabfold_batch $INP $OUT ${@:3}"
    colabfold_batch "$INP" "$OUT" "${@:3}"
  else
    colabfold_batch "$@"
  fi
fi


# --- cleanup: remove per-run temp A3M if it lives under tmp_af2_pairs
# Only runs on success because 'set -euo pipefail' is enabled.
if [[ -n "${INP:-}" ]]; then
  case "$INP" in
    */tmp_af2_pairs/*)
      rm -f -- "$INP" || true
      # Try to remove the tmp directory if now empty (safe with concurrent jobs)
      rmdir -p "$(dirname "$INP")" 2>/dev/null || true
      ;;
  esac
fi