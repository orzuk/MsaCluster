#!/usr/bin/env bash
set -euo pipefail

# --- Activate AF3 venv
source /sci/labs/orzuk/orzuk/af3-venv/bin/activate

# --- Parse minimal args
if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <fasta|af3_json> <outdir> [--af3-json <path_from_colabfold_search>] [--inference-only] [extra run_alphafold.py flags...]"
  exit 2
fi
INP="$1"        # FASTA or AF3-JSON
OUT="$2"; shift 2
mkdir -p "$OUT"

# --- CUDA toolchain discovery (like AF2 wrapper)
if [[ -n "${CUDA_DIR_OVERRIDE:-}" && -x "${CUDA_DIR_OVERRIDE}/bin/ptxas" ]]; then
  CUDA_DIR="${CUDA_DIR_OVERRIDE}"
elif command -v ptxas >/dev/null 2>&1; then
  CUDA_DIR="$(dirname "$(dirname "$(readlink -f "$(command -v ptxas)")")")"
elif [[ -x "/usr/local/APP/nvidia/cuda/12.4.1/bin/ptxas" ]]; then
  CUDA_DIR="/usr/local/APP/nvidia/cuda/12.4.1"
else
  echo "[fatal] Couldn't locate ptxas. Set CUDA_DIR_OVERRIDE=/path/to/cuda" >&2
  exit 2
fi
export PATH="${CUDA_DIR}/bin:${PATH}"
export XLA_FLAGS="--xla_gpu_cuda_data_dir=${CUDA_DIR}"

# --- Big tmp/cache to avoid HOME filling up
export TMPDIR=/sci/labs/orzuk/orzuk/alphafold3/tmp
export TMP=$TMPDIR; export TEMP=$TMPDIR
export XDG_CACHE_HOME=/sci/labs/orzuk/orzuk/alphafold3/cache
export JAX_CACHE_DIR=$XDG_CACHE_HOME/jax
export HF_HOME=/sci/labs/orzuk/orzuk/alphafold3/hf_cache
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_MEM_FRACTION=.80

# --- Auto-relaunch under srun if we’re not on a GPU shell
if [[ "${USE_SRUN:-1}" == "1" && -z "${CUDA_VISIBLE_DEVICES:-}" && -z "${SLURM_JOB_ID:-}" ]]; then
  exec srun --gres=gpu:1 --time=02:00:00 --mem=32G --cpus-per-task=8 --pty bash -lc \
    "source /sci/labs/orzuk/orzuk/af3-venv/bin/activate; CUDA_DIR=\"$(dirname "$(dirname "$(readlink -f "$(which ptxas)")")")\"; export PATH=\"\$CUDA_DIR/bin:\$PATH\"; export XLA_FLAGS=\"--xla_gpu_cuda_data_dir=\$CUDA_DIR\"; $(printf '%q ' "$0" "$INP" "$OUT" "$@")"
fi

echo "[diag] nvidia-smi:"; nvidia-smi || true
echo "[diag] which ptxas: $(command -v ptxas || echo 'not found')"
echo "[diag] CUDA_DIR=${CUDA_DIR}"
echo "[diag] TMPDIR=$TMPDIR  XDG_CACHE_HOME=$XDG_CACHE_HOME"

# --- Ensure the CCD pickle is present (symlink from repo if missing)
PYSP=$(python -c "import site,sys; print([p for p in site.getsitepackages() if 'site-packages' in p][0])")
PICKLE="$PYSP/alphafold3/constants/converters/chemical_component_sets.pickle"
if [[ ! -s "$PICKLE" && -s "/sci/labs/orzuk/orzuk/github/alphafold3/alphafold3/constants/converters/chemical_component_sets.pickle" ]]; then
  mkdir -p "$(dirname "$PICKLE")"
  ln -sf /sci/labs/orzuk/orzuk/github/alphafold3/alphafold3/constants/converters/chemical_component_sets.pickle "$PICKLE"
fi

# --- Build or take AF3 JSON
JSON=""
AF3_JSON_FROM_CF=""
INFERENCE_ONLY=0

# parse optional flags we support here
while [[ $# -gt 0 ]]; do
  case "$1" in
    --af3-json) AF3_JSON_FROM_CF="$2"; shift 2 ;;
    --inference-only) INFERENCE_ONLY=1; shift ;;
    *) break ;; # leave remaining flags for run_alphafold.py
  esac
done

if [[ -n "$AF3_JSON_FROM_CF" ]]; then
  # Use ColabFold-produced AF3 JSON with MSAs
  JSON="$OUT/input_af3.json"
  python - <<'PY'
import json,sys,shutil,os
src=sys.argv[1]; dst=sys.argv[2]
os.makedirs(os.path.dirname(dst), exist_ok=True)
with open(src) as f: data=json.load(f)
with open(dst,'w') as f: json.dump(data,f,indent=2)
print("Wrote",dst)
PY "$AF3_JSON_FROM_CF" "$JSON"
  INFERENCE_ONLY=1
else
  # If INP already ends with .json, use it directly
  if [[ "$INP" == *.json ]]; then
    JSON="$INP"
  else
    # Create minimal JSON from FASTA (single-chain A)
    NAME="$(basename "${INP%.*}")"
    JSON="$OUT/${NAME}.json"
    python - <<'PY'
import json,sys,Pathlib
from pathlib import Path
fa=Path(sys.argv[1]).read_text().splitlines()
seq="".join(l.strip() for l in fa if not l.startswith(">"))
name=sys.argv[2]; out=Path(sys.argv[3]); out.parent.mkdir(parents=True, exist_ok=True)
payload={"name":name,"sequences":[{"protein":{"id":["A"],"sequence":seq}}],"modelSeeds":[1],"dialect":"alphafold3","version":1}
out.write_text(json.dumps(payload,indent=2)); print("Wrote",out)
PY "$INP" "$NAME" "$JSON"
  fi
fi

# --- Run AF3
RUN=/sci/labs/orzuk/orzuk/github/alphafold3/run_alphafold.py
if [[ $INFERENCE_ONLY -eq 1 ]]; then
  echo "[run] AF3 inference-only"
  python "$RUN" --json_path "$JSON" --model_dir /sci/labs/orzuk/orzuk/trained_models/alphafold3 --output_dir "$OUT" --run_data_pipeline=false --run_inference=true "$@"
else
  echo "[run] AF3 full pipeline"
  python "$RUN" --json_path "$JSON" --model_dir /sci/labs/orzuk/orzuk/trained_models/alphafold3 --output_dir "$OUT" "$@"
fi
