#!/usr/bin/env bash
set -euo pipefail

AF2_VENV="${AF2_VENV:-/sci/labs/orzuk/orzuk/af2-venv}"
CF_BATCH="${CF_BATCH:-$AF2_VENV/bin/colabfold_batch}"

# Activate AF2 venv (same as your current flow)
source "$AF2_VENV/bin/activate"

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <input.{fasta|fa|a3m}> <outdir> [extra colabfold flags...]"
  exit 2
fi
INP="$1"; OUT="$2"; shift 2
mkdir -p "$OUT"

# CUDA env (kept from your runner)
if [[ -n "${CUDA_DIR_OVERRIDE:-}" && -x "${CUDA_DIR_OVERRIDE}/bin/ptxas" ]]; then
  CUDA_DIR="${CUDA_DIR_OVERRIDE}"
elif command -v ptxas >/dev/null 2>&1; then
  CUDA_DIR="$(dirname "$(dirname "$(readlink -f "$(command -v ptxas)")")")"
elif [[ -x "/usr/local/APP/nvidia/cuda/12.4.1/bin/ptxas" ]]; then
  CUDA_DIR="/usr/local/APP/nvidia/cuda/12.4.1"
else
  echo "[fatal] Couldn't locate ptxas. Set CUDA_DIR_OVERRIDE=/path/to/cuda" >&2; exit 2
fi
export PATH="${CUDA_DIR}/bin:${PATH}"
export XLA_FLAGS="--xla_gpu_cuda_data_dir=${CUDA_DIR}"
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_MEM_FRACTION=.80

PYSP=$(python - <<'PY'
import site; print([p for p in site.getsitepackages() if 'site-packages' in p][0])
PY
)
export LD_LIBRARY_PATH="$PYSP/nvidia/cudnn/lib:$PYSP/nvidia/cuda_nvrtc/lib:$PYSP/nvidia/cuda_runtime/lib:$PYSP/nvidia/cublas/lib:${CUDA_DIR}/lib64:${LD_LIBRARY_PATH:-}"

echo "[diag] nvidia-smi:"; nvidia-smi || true
echo "[diag] which ptxas: $(command -v ptxas || echo not-found)"
echo "[diag] CUDA_DIR=$CUDA_DIR"

case "$INP" in
  *.a3m)
    echo "[mode] A3M given → run AF2 with provided MSA"
    echo "[run] $CF_BATCH $INP $OUT $*"
    "$CF_BATCH" "$INP" "$OUT" "$@"
    ;;
  *.fa|*.fasta)
    echo "[mode] FASTA given → first call fasta2MSA_ColabFold.sh → then run AF2 on produced A3M"
    MSADIR="$OUT/msas"; mkdir -p "$MSADIR"
    bash "$(dirname "$0")/fasta2MSA_ColabFold.sh" "$INP" "$MSADIR" --jobname-prefix "$(basename "${INP%.*}")"_ 
    A3M="$(ls -1 "$MSADIR"/$(basename "${INP%.*}")_*.a3m 2>/dev/null | head -n1 || true)"
    [[ -s "${A3M:-}" ]] || { echo "[fatal] No A3M in $MSADIR"; exit 3; }
    echo "[run] $CF_BATCH $A3M $OUT $*"
    "$CF_BATCH" "$A3M" "$OUT" "$@"
    ;;
  *)
    echo "[fatal] Input must be .fasta/.fa or .a3m"; exit 2;;
esac
