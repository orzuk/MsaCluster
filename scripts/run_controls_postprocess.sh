#!/bin/bash
# ============================================================================
# run_controls_postprocess.sh — sbatch wrapper for controls postprocess.
#
# Computes per-protein cluster structural diversity (df_diversity.csv) for
# the requested control set. CPU-only (TM-align via tmtools) — runs on the
# `glacier` partition, no GPU.
#
# Usage (from repo root):
#   sbatch -A course-52017-25 -p glacier \
#          -J pp_PorterIDP \
#          --output Pipeline/Controls/jobs/postprocess_PorterIDP.out \
#          scripts/run_controls_postprocess.sh PorterIDP
#
#   # All three control types in parallel:
#   for ct in RandomSingleFold PorterIDP MatchedSingleFold; do
#     sbatch -A course-52017-25 -p glacier \
#            -J pp_$ct \
#            --output "Pipeline/Controls/jobs/postprocess_${ct}.out" \
#            scripts/run_controls_postprocess.sh "$ct"
#   done
# ============================================================================
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 06:00:00
set -euo pipefail

source /sci/labs/orzuk/orzuk/venvs/my-python-venv/bin/activate
cd /sci/labs/orzuk/orzuk/github/MsaCluster

CT="${1:?usage: $0 <ControlType: PorterIDP | MatchedSingleFold | RandomSingleFold>}"

python3 run_controls_pipeline.py \
    --run_mode postprocess \
    --protein_ids ALL \
    --control_type "$CT"
