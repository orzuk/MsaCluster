#!/bin/bash
# ============================================================================
# run_permutation_null.sh — sbatch wrapper for the empirical permutation-null
#                           test on per-pair cross-method concordance and
#                           centered-diversity statistics.
#
# Reads the per-cluster CSV (docs/fold_diversity_survey.csv) produced by
# scripts/fold_diversity_survey.py and writes
# docs/fold_diversity_concordance_perm.csv with empirical p-values that
# replace the binomial-approximation p-values in the manuscript.
#
# This is a pure re-statistical computation on existing CSV data — it does
# NOT re-run AF2/ESM/MSAT/DDG. Wall time on 93 pairs * 1000 permutations is
# ~1-10 minutes single-threaded; you can also run inline:
#
#   python3 scripts/permutation_null_concordance.py
#
# Usage (cluster, from repo root):
#   sbatch -A course-52017-25 -p glacier \
#          -J perm_null \
#          --output Pipeline/FoldPairs/jobs/perm_null.out \
#          scripts/run_permutation_null.sh
# ============================================================================
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 01:00:00
set -euo pipefail

source /sci/labs/orzuk/orzuk/venvs/my-python-venv/bin/activate
cd /sci/labs/orzuk/orzuk/github/MsaCluster

mkdir -p Pipeline/FoldPairs/jobs

python3 scripts/permutation_null_concordance.py \
    --n-perm 1000 \
    --input  docs/fold_diversity_survey.csv \
    --output docs/fold_diversity_concordance_perm.csv \
    --seed 12345
