#!/usr/bin/env bash
# Pipeline/RunAF.sh
# SLURM options (ignored on local runs)
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --gres=gpu:1            # change to your cluster's syntax if needed (e.g., gpu:a100:1)

set -euo pipefail

# Resolve ColabFold binary + optional --data flag and JAX envs
. "$(dirname "$0")/colabfold_env.sh"

# If the cluster uses Environment Modules, try loading CUDA/CuDNN; harmless locally
if command -v module >/dev/null 2>&1; then
  module purge || true
  module load cuda/12.1 2>/dev/null || module load cuda/11.8 2>/dev/null || module load cuda/11.1 2>/dev/null || true
  module load cudnn/8    2>/dev/null || true
fi

# Inputs (override by passing args): ./Pipeline/RunAF.sh <INPUT_DIR> <OUTPUT_DIR>
INPUT_DIR="${1:-./output/output_msa_cluster}"
OUTPUT_DIR="${2:-./output/AF_preds}"
mkdir -p "$OUTPUT_DIR"

# Run ColabFold (add flags you use regularly below)
"$COLABFOLD_BIN" "$INPUT_DIR" "$OUTPUT_DIR" "${CF_DATA_FLAGS[@]}" --amber --templates
