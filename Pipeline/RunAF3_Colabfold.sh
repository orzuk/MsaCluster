#!/usr/bin/env bash
set -euo pipefail

# --- CONFIG
AF3_VENV="/sci/labs/orzuk/orzuk/af3-venv"
AF3_REPO="/sci/labs/orzuk/orzuk/github/alphafold3"
AF3_MODELDIR="/sci/labs/orzuk/orzuk/trained_models/alphafold3"
SCR_BASE="/sci/labs/orzuk/orzuk/alphafold3"               # big caches
FAKE_HOME="/sci/labs/orzuk/orzuk/software/conda_tmp/home_fake"

# Path to your A3M→AF3-JSON converter
CONVERTER="${CONVERTER:-/sci/labs/orzuk/orzuk/github/MsaCluster/a3m_toaf3json.py}"


# ColabFold tools (prefer PATH, else AF2 venv fallbacks)
CF_SEARCH="${CF_SEARCH:-$(command -v colabfold_search || true)}"
if [[ -z "$CF_SEARCH" && -x "/sci/labs/orzuk/orzuk/af2-venv/bin/colabfold_search" ]]; then
  CF_SEARCH="/sci/labs/orzuk/orzuk/af2-venv/bin/colabfold_search"
fi
: "${CF_SEARCH:?colabfold_search not found (set CF_SEARCH=/path/to/colabfold_search)}"

# --- Activate AF3 venv
source "$AF3_VENV/bin/activate"

# --- Args
if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <input.{fasta|fa|a3m|json}> <outdir> [extra run_alphafold.py flags...]"
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

# --- Big caches away from $HOME (which is small on your cluster)
export TMPDIR="$SCR_BASE/tmp"
export TMP="$TMPDIR"; export TEMP="$TMPDIR"
export XDG_CACHE_HOME="$SCR_BASE/cache"
export JAX_CACHE_DIR="$XDG_CACHE_HOME/jax"
export HF_HOME="$SCR_BASE/hf_cache"
export HOME="$FAKE_HOME"; mkdir -p "$HOME" "$TMPDIR" "$XDG_CACHE_HOME" "$HF_HOME"

echo "[diag] nvidia-smi:"; nvidia-smi || true
echo "[diag] which ptxas: $(command -v ptxas || echo not-found)"
echo "[diag] CUDA_DIR=$CUDA_DIR"
echo "[diag] TMPDIR=$TMPDIR  XDG_CACHE_HOME=$XDG_CACHE_HOME"

# --- Ensure CCD pickle exists; build once if missing (uses wwPDB CCD)
PYSP=$(python -c "import site; print([p for p in site.getsitepackages() if 'site-packages' in p][0])")
PICKLE="$PYSP/alphafold3/constants/converters/chemical_component_sets.pickle"
if [[ ! -s "$PICKLE" ]]; then
  echo "[fix] CCD pickle missing; building"
  CCD="/sci/labs/orzuk/orzuk/software/conda_tmp/ccd/components.cif"
  if [[ ! -s "$CCD" ]]; then
    echo "[fix] Downloading CCD"
    mkdir -p "$(dirname "$CCD")"
    curl -L -o "${CCD}.gz" https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz || \
    wget -O "${CCD}.gz"    https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz
    gunzip -f "${CCD}.gz"
  fi
  # try console script then module forms
  (command -v build_data >/dev/null && build_data --components-cif "$CCD") || \
  python -m alphafold3.build_data --components-cif "$CCD" || \
  python -m alphafold3.scripts.build_data --components-cif "$CCD"
fi

# --- Helper: extract first sequence from A3M
a3m_to_first_seq() {
  python - "$1" <<'PY'
import sys, re, pathlib
a3m = pathlib.Path(sys.argv[1]).read_text().splitlines()
seq=[]; started=False
for ln in a3m:
    if ln.startswith('>'):
        if started: break
        started=True; continue
    if not started: continue
    seq.append(re.sub(r'[a-z.]','',ln.strip()))
print(''.join(seq))
PY
}

# --- Helper: A3M -> AF3 JSON (unpaired MSA for chain A)
a3m_to_af3json() {
  local fa="$1" a3m="$2" outjson="$3"
  python - "$fa" "$a3m" "$outjson" <<'PY'
import sys, json, pathlib, re
fa = pathlib.Path(sys.argv[1]).read_text().splitlines()
seq = "".join(l.strip() for l in fa if not l.startswith(">"))
def read_a3m(p):
    aln=[]; seq=None
    with open(p) as f:
        for line in f:
            line=line.rstrip("\n")
            if not line: continue
            if line.startswith(">"):
                if seq is not None: aln.append(seq)
                seq=""
            else:
                seq += re.sub(r"[a-z.]", "", line)
        if seq is not None: aln.append(seq)
    return aln
msa = read_a3m(sys.argv[2])
if not msa or msa[0] != seq:
    msa = [seq] + [s for s in msa if s != seq]
delmat = [[0]*len(s) for s in msa]
payload = {
  "name": pathlib.Path(sys.argv[1]).stem,
  "sequences": [{"protein": {"id": ["A"], "sequence": seq}}],
  "unpairedMSAs": [[{
    "descriptions": ["colabfold_a3m"]*len(msa),
    "sequences": msa,
    "deletionMatrix": delmat
  }]],
  "modelSeeds": [1],
  "dialect": "alphafold3",
  "version": 1
}
out = pathlib.Path(sys.argv[3]); out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(payload, indent=2))
print("Wrote", out)
PY
}

# --- Decide mode and run AF3 (inference-only; no local DBs)
RUN="$AF3_REPO/run_alphafold.py"
JSON=""

case "$INP" in
  *.json)
    echo "[mode] AF3-JSON provided → inference-only"
    JSON="$INP"
    ;;
  *.a3m)
    echo "[mode] A3M provided → convert to AF3-JSON → AF3 inference-only"
    NAME="$(basename "${INP%.*}")"
    TMPFA="$OUT/${NAME}_first.fa"
    # helper to grab the first sequence of the A3M and write a small FASTA
    FIRSTSEQ="$(a3m_to_first_seq "$INP")"
    printf ">query\n%s\n" "$FIRSTSEQ" > "$TMPFA"
    JSON="$OUT/${NAME}_af3.json"
    echo "[conv] python $CONVERTER $TMPFA $INP $JSON"
    python "$CONVERTER" "$TMPFA" "$INP" "$JSON"
    ;;
  *.fa|*.fasta)
    echo "[mode] FASTA → build MSA with colabfold_search (3-arg) → convert to AF3-JSON → AF3 inference-only"
    MSADIR="$OUT/msas"; mkdir -p "$MSADIR"
    NAME="$(basename "${INP%.*}")"
    BASE="$MSADIR/$NAME"
    echo "[run] $CF_SEARCH $INP $MSADIR $BASE"
    "$CF_SEARCH" "$INP" "$MSADIR" "$BASE"
    A3M="${BASE}.a3m"
    [[ -s "$A3M" ]] || { echo "[fatal] A3M not found at $A3M"; exit 2; }
    JSON="$OUT/${NAME}_af3.json"
    echo "[conv] python $CONVERTER $INP $A3M $JSON"
    python "$CONVERTER" "$INP" "$A3M" "$JSON"
    ;;
  *)
    echo "[fatal] Input must be .json/.a3m/.fasta/.fa"; exit 2;;
esac

echo "[run] AF3 inference-only"
python "$RUN" \
  --json_path "$JSON" \
  --model_dir "$AF3_MODELDIR" \
  --output_dir "$OUT" \
  --run_data_pipeline=false \
  --run_inference=true \
  "$@"
