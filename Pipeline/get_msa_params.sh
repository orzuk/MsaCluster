#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
set -euo pipefail

SEED_A3M="${1:?missing seed a3m path}"
PAIR_ID="${2:?missing pair id}"

REPO_DIR="/sci/labs/orzuk/orzuk/github/MsaCluster"
cd "$REPO_DIR"

mkdir -p "Pipeline/${PAIR_ID}/output_get_msa"

# Run get_msa.py inside the AF2 venv via the unified wrapper
bash ./Pipeline/RunAF2_Colabfold.sh --python ./get_msa.py \
  "$SEED_A3M" "Pipeline/${PAIR_ID}/output_get_msa" --name DeepMsa
