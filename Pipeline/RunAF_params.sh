#!/bin/bash
#SBATCH --job-name=af2_pair
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --gres=gpu:a100:1
# (Optional, site-specific — uncomment if your cluster requires it)
# #SBATCH --partition=a100
# #SBATCH --constraint=a100

set -euo pipefail

# ------------------ args & paths ------------------
PAIR_ID="${1:?Usage: RunAF_params.sh <PAIR_ID> [BASE_DIR] }"
BASE_DIR="${2:-Pipeline}"
PAIR_DIR="${BASE_DIR}/${PAIR_ID}"

MSA_DEEP="${PAIR_DIR}/output_get_msa"
MSA_CLUST="${PAIR_DIR}/output_msa_cluster"
OUT_DIR="${PAIR_DIR}/output_AF/AF2"
mkdir -p "${OUT_DIR}"

echo "[AF2] Pair:           ${PAIR_ID}"
echo "[AF2] Deep MSA dir:   ${MSA_DEEP}"
echo "[AF2] Cluster MSA dir:${MSA_CLUST}"
echo "[AF2] Output dir:     ${OUT_DIR}"

# ------------------ environment -------------------
# If you have a venv/conda to activate, set this env var before sbatch:
#   export COLABFOLD_ENV_ACTIVATE=/path/to/your/env/bin/activate
if [[ -n "${COLABFOLD_ENV_ACTIVATE:-}" && -f "${COLABFOLD_ENV_ACTIVATE}" ]]; then
  # shellcheck disable=SC1090
  source "${COLABFOLD_ENV_ACTIVATE}"
fi

# Or use site modules instead of a venv (example):
# module load cuda/12.1
# module load python/3.11

COLABFOLD_BIN="${COLABFOLD_BIN:-colabfold_batch}"
CF_FLAGS=${CF_FLAGS:-}

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-4}}"

# ------------------ run AF2 -----------------------
# Run once on the Deep MSA (full alignment) and once on the clustered MSAs.
# Adjust these two lines if your AF frontend expects different inputs.
"${COLABFOLD_BIN}" "${MSA_DEEP}"  "${OUT_DIR}" ${CF_FLAGS}
"${COLABFOLD_BIN}" "${MSA_CLUST}" "${OUT_DIR}" ${CF_FLAGS}

# ------------------ post-process (optional) -------
shopt -s nullglob
if compgen -G "${OUT_DIR}/*.json" > /dev/null; then gzip -f "${OUT_DIR}"/*.json; fi
if compgen -G "${OUT_DIR}/*.pdb"  > /dev/null; then gzip -f "${OUT_DIR}"/*.pdb;  fi

# Convert PNG to JPG if ImageMagick is available (optional)
if command -v convert >/dev/null 2>&1; then
  for p in "${OUT_DIR}"/*.png; do
    [[ -e "$p" ]] || continue
    convert "$p" "${p%.*}.jpg" && rm -f "$p"
  done
fi

echo "[AF2] Done for ${PAIR_ID}"
