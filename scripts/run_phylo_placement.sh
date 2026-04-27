#!/bin/bash
# ============================================================================
# run_phylo_placement.sh — sbatch wrapper for the phylogenetic-structure test
#                          of cluster fold-preference (F1c/F2c/Amc).
#
# Reads:
#   - Pipeline/FoldPairs/<pair>/output_phytree/DeepMsa_tree.nwk   (per pair)
#   - Pipeline/FoldPairs/<pair>/output_msa_cluster/ShallowMsa_*.a3m
#   - docs/fold_diversity_survey.csv
#
# Writes:
#   - docs/phylo_placement.csv               (one row per pair x method)
#   - docs/phylo_placement_xcl1_clusters.csv (XCL1 cross-check, if pair present)
#   - docs/figs/phylo_placement_<pair>_<method>.png  (top-N hits, per method)
#
# Usage (from repo root):
#   sbatch -A course-52017-25 -p glacier \
#          -J phylo_placement \
#          --output Pipeline/FoldPairs/jobs/phylo_placement.out \
#          scripts/run_phylo_placement.sh
# ============================================================================
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 04:00:00
set -euo pipefail

source /sci/labs/orzuk/orzuk/venvs/my-python-venv/bin/activate
cd /sci/labs/orzuk/orzuk/github/MsaCluster

mkdir -p Pipeline/FoldPairs/jobs docs/figs

python3 scripts/phylo_placement.py \
    --pairs ALL \
    --make-figures 10 \
    --seed 12345
