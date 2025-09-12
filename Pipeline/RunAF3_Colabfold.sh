#!/usr/bin/env bash
set -euo pipefail

# --- CONFIG
AF3_VENV="${AF3_VENV:-/sci/labs/orzuk/orzuk/af3-venv}"
AF3_REPO="${AF3_REPO:-/sci/labs/orzuk/orzuk/github/alphafold3}"
AF3_MODELDIR="${AF3_MODELDIR:-/sci/labs/orzuk/orzuk/trained_models/alphafold3}"
SCR_BASE="${SCR_BASE:-/sci/labs/orzuk/orzuk/alphafold3}"    # big caches
FAKE_HOME="${FAKE_HOME:-/sci/labs/orzuk/orzuk/software/conda_tmp/home_fake}"

CONVERTER="${CONVERTER:-/sci/labs/orzuk/orzuk/github/MsaCluster/a3m_toaf3json.py}"  # your converter

# Activate AF3 env
source "$AF3_VENV/bin/activate"

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <input.{fasta|fa|a3m|json}> <outdir> [extra run_alphafold.py flags...]"
  exit 2
fi
INP="$1"; OUT="$2"; shift 2
mkdir -p "$OUT"


# --- Parse optional --pdb=all|rank1|none (default: all)
PDB_MODE="all"
AF3_ARGS=()
for arg in "$@"; do
  case "$arg" in
    --pdb=all)   PDB_MODE="all" ;;
    --pdb=rank1) PDB_MODE="rank1" ;;
    --pdb=none|--no-pdb) PDB_MODE="none" ;;
    *) AF3_ARGS+=("$arg") ;;
  esac
done
set -- "${AF3_ARGS[@]}"


# CUDA (same style as your AF3 runner)
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

# Caches off $HOME
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

# Ensure CCD pickle exists (kept from your AF3 script)
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
  (command -v build_data >/dev/null && build_data --components-cif "$CCD") || \
  python -m alphafold3.build_data --components-cif "$CCD" || \
  python -m alphafold3.scripts.build_data --components-cif "$CCD"
fi

RUN="$AF3_REPO/run_alphafold.py"
JSON=""

# Helper: first sequence from A3M -> tmp FASTA (for converter sanity)
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

case "$INP" in
  *.json)
    echo "[mode] AF3-JSON provided → inference-only"
    JSON="$INP"
    ;;
  *.a3m)
    echo "[mode] A3M provided → convert to AF3-JSON → AF3 inference-only"
    NAME="$(basename "${INP%.*}")"
    TMPFA="$OUT/${NAME}_first.fa"
    a3m_to_first_fa "$INP" "$TMPFA"
    JSON="$OUT/${NAME}_af3.json"
    echo "[conv] python $CONVERTER $TMPFA $INP $JSON"
    python "$CONVERTER" "$TMPFA" "$INP" "$JSON"
    ;;
  *.fa|*.fasta)
    echo "[mode] FASTA provided → call shared fasta2MSA_Colabfold.sh → convert → AF3 inference-only"
    MSADIR="$OUT/msas"; mkdir -p "$MSADIR"
    bash "$(dirname "$0")/fasta2MSA_Colabfold.sh" "$INP" "$MSADIR" --jobname-prefix "$(basename "${INP%.*}")"_
    A3M="$(ls -1 "$MSADIR"/$(basename "${INP%.*}")_*.a3m 2>/dev/null | head -n1 || true)"
    [[ -s "${A3M:-}" ]] || { echo "[fatal] No A3M in $MSADIR"; exit 3; }
    NAME="$(basename "${INP%.*}")"
    TMPFA="$OUT/${NAME}_first.fa"
    a3m_to_first_fa "$A3M" "$TMPFA"
    JSON="$OUT/${NAME}_af3.json"
    echo "[conv] python $CONVERTER $TMPFA $A3M $JSON"
    python "$CONVERTER" "$TMPFA" "$A3M" "$JSON"
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


# --- Optional: export PDBs from mmCIFs (default: all)
if [[ "${PDB_MODE}" != "none" ]]; then
  # job name from the json we just ran
  JOB_NAME="$(python - <<'PY'
import json; print(json.load(open(r"""$JSON""")).get("name","job"))
PY
)"
  # find newest AF3 output dir for this job
  OUTDIR="$(ls -1dt "$OUT/${JOB_NAME}_"* 2>/dev/null | head -n1)"
  if [[ ! -d "$OUTDIR" ]]; then
    echo "[warn] cannot locate AF3 output dir for $JOB_NAME"; exit 0
  fi

  convert_cif_to_pdb() {
    local cif="$1" pdb="$2"
    if command -v gemmi >/dev/null 2>&1; then
      gemmi convert "$cif" "$pdb"
    else
      python - <<'PY' || exit 1
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
import sys, os
cif, pdb = r"""$1""", r"""$2"""
parser=MMCIFParser(QUIET=True)
s=parser.get_structure("model", cif)
io=PDBIO(); io.set_structure(s); io.save(pdb)
PY
    fi
  }

  if [[ "$PDB_MODE" == "rank1" ]]; then
    CSV="$OUTDIR/${JOB_NAME}_ranking_scores.csv"
    if [[ -s "$CSV" ]]; then
      BEST_MODEL="$(python - <<'PY'
import csv, sys
p=r"""$CSV"""
rows=list(csv.DictReader(open(p)))
keys=("ranking_score","mean_plddt","plddt","confidence")
k=next((x for x in keys if x in rows[0]), list(rows[0].keys())[-1])
best=max(rows, key=lambda r: float(r[k]))
print(best.get("model","seed-1_sample-0"))
PY
)"
      CIF="$OUTDIR/$BEST_MODEL/${JOB_NAME}_${BEST_MODEL}_model.cif"
      PDB="$OUTDIR/$BEST_MODEL/${JOB_NAME}_${BEST_MODEL}_model.pdb"
      if [[ -f "$CIF" ]]; then
        convert_cif_to_pdb "$CIF" "$PDB" && echo "[ok] wrote $PDB"
        # also drop a convenience copy at root
        cp -f "$PDB" "$OUTDIR/${JOB_NAME}_rank1_model.pdb"
      else
        echo "[warn] missing CIF for best model: $CIF"
      fi
    else
      echo "[warn] ranking csv not found, skipping rank1 export"
    fi
  else
    # PDB_MODE=all: convert every *_model.cif we find
    while IFS= read -r -d '' CIF; do
      PDB="${CIF%.cif}.pdb"
      convert_cif_to_pdb "$CIF" "$PDB" && echo "[ok] wrote $PDB"
    done < <(find "$OUTDIR" -type f -name "*_model.cif" -print0)
  fi
fi


