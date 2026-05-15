#!/bin/bash
# ============================================================================
# run_ancestral.sh — sbatch wrapper for ancestral fold-preference
#                    reconstruction across the 93 fold-switching pairs.
#
# Runs scripts/run_ancestral_analysis.py on all eligible pairs (those with
# a per-pair UPGMA tree and df_af.csv available). For each pair the script:
#   - builds a coarse cluster-level tree (one leaf per ShallowMsa cluster),
#   - assigns each cluster a fold-preference label (F1/F2/Amb) from the
#     centered AF2 TM-score with delta = 0.05,
#   - runs ancestral state reconstruction (pastml MPPA, Fitch fallback),
#     counts transitions, classifies the pattern (single-origin / convergent
#     / mixed / uniform), and runs an empirical permutation null test,
#   - renders an ancestral-tree PNG colored by fold preference.
#
# Reads:
#   - Pipeline/FoldPairs/<pair>/output_phytree/DeepMsa_tree.nwk
#   - Pipeline/FoldPairs/<pair>/output_msa_cluster/ShallowMsa_*.a3m
#   - Pipeline/FoldPairs/<pair>/Analysis/df_af.csv
#
# Writes:
#   - docs/ancestral_summary.csv                       (one row per pair)
#   - Pipeline/FoldPairs/<pair>/output_phytree/ancestral_tree.png (per pair)
#
# Usage (from repo root):
#   sbatch -A course-52017-25 -p glacier \
#          -J ancestral \
#          --output Pipeline/FoldPairs/jobs/ancestral.out \
#          scripts/run_ancestral.sh
# ============================================================================
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 06:00:00
set -euo pipefail

source /sci/labs/orzuk/orzuk/venvs/my-python-venv/bin/activate
cd /sci/labs/orzuk/orzuk/github/MsaCluster

mkdir -p Pipeline/FoldPairs/jobs

python3 scripts/run_ancestral_analysis.py \
    --pairs ALL \
    --delta 0.05
