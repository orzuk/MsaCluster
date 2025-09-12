#!/usr/bin/env bash
set -euo pipefail

# --- CONFIG: paths (edit if your layout differs)
AF2_VENV="/sci/labs/orzuk/orzuk/af2-venv"

# Prefer whatever is on PATH; otherwise fall back to AF2 venv tools
CF_SEARCH="${CF_SEARCH:-$(command -v colabfold_search || true)}"
CF_BATCH="${CF_BATCH:-$(command -v colabfold_batch || true)}"
if [[ -z "$CF_SEARCH" && -x "$AF2_VENV/bin/colabfold_search" ]]; then CF_SEARCH="$AF2_VENV/bin/colabfold_search"; fi
if [[ -z "$CF_BATCH"  && -x "$AF2_VENV/bin/colabfold_batch"  ]]; then CF_BATCH="$AF2_VENV/bin/colabfold_batch";  fi
: "${CF_SEARCH:?colabfold_search not found (set CF_SEARCH=/path/to/colabfold_search)}"
: "${CF_BATCH:?colabfold_batch not found (set CF_BATCH=/path/to/colabfold_batch)}"

# --- Activate AF2 venv (so XLA + libs are right for colabfold)
source "$AF2_VENV/bin/activate"

# --- Args
if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <input.{fasta|fa|a3m}> <outdir> [extra colabfold flags...]"
  exit 2
fi
INP="$1"; OUT="$2"; shift 2
mkdir -p "$OUT"

# --- CUDA discovery
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

# --- Add pip CUDA/CUDNN/CUBLAS/NVRTC wheels to runtime path (+ system lib64)
PYSP=$(python - <<'PY'
import site; print([p for p in site.getsitepackages() if 'site-packages' in p][0])
PY
)
export LD_LIBRARY_PATH="$PYSP/nvidia/cudnn/lib:$PYSP/nvidia/cuda_nvrtc/lib:$PYSP/nvidia/cuda_runtime/lib:$PYSP/nvidia/cublas/lib:${CUDA_DIR}/lib64:${LD_LIBRARY_PATH:-}"

echo "[diag] nvidia-smi:"; nvidia-smi || true
echo "[diag] which ptxas: $(command -v ptxas || echo not-found)"
echo "[diag] CUDA_DIR=$CUDA_DIR"

# --- Helper: extract first sequence from A3M to a tmp FASTA
a3m_to_first_fa() {
  local a3m="$1" outfa="$2"
  python - "$a3m" "$outfa" <<'PY'
import sys, re, pathlib
a3m = pathlib.Path(sys.argv[1]).read_text().splitlines()
seq=[]; started=False
for ln in a3m:
    if ln.startswith('>'):
        if started: break
        started=True; continue
    if not started: continue
    seq.append(re.sub(r'[a-z.]','',ln.strip()))
s=''.join(seq)
with open(sys.argv[2],'w') as f:
    f.write('>query\n'); f.write(s+'\n')
print("wrote", sys.argv[2], "len", len(s))
PY
}


# Force binaries
CF_SEARCH="${CF_SEARCH:-}"
CF_BATCH="${CF_BATCH:-}"
if [[ -x "$AF2_VENV/bin/colabfold_search" ]]; then CF_SEARCH="$AF2_VENV/bin/colabfold_search"; fi
if [[ -x "$AF2_VENV/bin/colabfold_batch"  ]]; then CF_BATCH="$AF2_VENV/bin/colabfold_batch";  fi
# If still empty, try PATH:
: "${CF_BATCH:=$(command -v colabfold_batch || true)}"
: "${CF_SEARCH:=$(command -v colabfold_search || true)}"
: "${CF_BATCH:?colabfold_batch not found}"

# --- Decide mode
case "$INP" in
  *.a3m)
    echo "[mode] A3M → AF2 using given MSA"
    echo "[run] $CF_BATCH $INP $OUT $*"
    "$CF_BATCH" "$INP" "$OUT" "$@"
    ;;
 *.fa|*.fasta)
    echo "[mode] FASTA → AF2 via colabfold_batch (does its own MSA)"
    echo "[run] $CF_BATCH $INP $OUT $*"
    "$CF_BATCH" "$INP" "$OUT" "$@"
    ;;
  *)
    echo "[fatal] Input must be .fasta/.fa or .a3m"; exit 2;;
esac
