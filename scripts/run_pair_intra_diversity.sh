#!/bin/bash
# ============================================================================
# run_pair_intra_diversity.sh — sbatch wrapper for fold-pair intra-cluster
#                               pairwise TM-score computation.
#
# Computes the same statistic as run_controls_pipeline.task_postprocess
# (mean/median/min/max/std of pairwise TM-score between cluster prediction
# PDBs) but for the 93 fold-switching pairs instead of single-protein
# controls. Output: docs/foldpair_intra_diversity.csv (one row per pair).
#
# Usage (from repo root):
#   sbatch -A course-52017-25 -p glacier \
#          -J pair_diversity \
#          --output Pipeline/FoldPairs/jobs/pair_diversity.out \
#          scripts/run_pair_intra_diversity.sh
# ============================================================================
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 06:00:00
set -euo pipefail

source /sci/labs/orzuk/orzuk/venvs/my-python-venv/bin/activate
cd /sci/labs/orzuk/orzuk/github/MsaCluster

mkdir -p Pipeline/FoldPairs/jobs

python3 scripts/compute_pair_intra_diversity.py --pairs ALL
