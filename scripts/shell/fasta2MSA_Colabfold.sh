#!/usr/bin/env bash
set -euo pipefail

# One job: FASTA -> A3M (using the *working* AF2 venv colabfold_batch --msa-only)
AF2_VENV="${AF2_VENV:-/sci/labs/orzuk/orzuk/af2-venv}"
CF_BATCH="${CF_BATCH:-$AF2_VENV/bin/colabfold_batch}"

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <input.{fasta|fa}> <out_msa_dir> [--jobname-prefix NAME_ | other colabfold flags...]"
  exit 2
fi
INP="$1"; OUT="$2"; shift 2
mkdir -p "$OUT"

# sanity
: "${CF_BATCH:?colabfold_batch not found (set CF_BATCH=/path/to/colabfold_batch)}"

# CUDA discovery (same style as your runners)
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

echo "[diag] nvidia-smi:"; nvidia-smi || true
echo "[diag] which ptxas: $(command -v ptxas || echo not-found)"
echo "[diag] CUDA_DIR=$CUDA_DIR"

case "$INP" in
  *.fa|*.fasta)
    NAME="$(basename "${INP%.*}")"
    # Let caller pass flags like --jobname-prefix, --host-url, etc.
    echo "[run] $CF_BATCH --msa-only $* $INP $OUT"
    "$CF_BATCH" --msa-only "$@" "$INP" "$OUT"

    # Find output A3M (colabfold names: ${NAME}_*.a3m)
    A3M="$(ls -1 "$OUT"/${NAME}_*.a3m 2>/dev/null | head -n1 || true)"
    [[ -s "${A3M:-}" ]] || { echo "[fatal] No A3M produced in $OUT"; exit 3; }
    echo "[ok] A3M: $A3M"
    ;;
  *)
    echo "[fatal] Input must be .fasta/.fa"; exit 2;;
esac
