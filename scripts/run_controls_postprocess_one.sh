#!/bin/bash
# ============================================================================
# run_controls_postprocess_one.sh — sbatch wrapper for ONE protein's postprocess
#
# Usage (from repo root):
#   sbatch -J pp_PorterIDP_Q9NZC2 \
#          -c 2 --mem=8G -t 02:00:00 \
#          -A course-52017-25 -p glacier \
#          --output Pipeline/Controls/PorterIDP/jobs/postprocess_Q9NZC2.out \
#          scripts/run_controls_postprocess_one.sh PorterIDP Q9NZC2
#
# Designed for granular per-protein submission so 100+ proteins can run
# in parallel on glacier instead of serial in one job. Avoids `--wrap`
# multi-line shell-quoting issues.
# ============================================================================
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 02:00:00
set -euo pipefail

source /sci/labs/orzuk/orzuk/venvs/my-python-venv/bin/activate
cd /sci/labs/orzuk/orzuk/github/MsaCluster

CT="${1:?usage: $0 <ControlType> <ProteinID>}"
PID="${2:?usage: $0 <ControlType> <ProteinID>}"

python3 run_controls_pipeline.py \
    --run_mode postprocess \
    --protein_ids "$PID" \
    --control_type "$CT"
