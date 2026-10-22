#!/bin/bash
#SBATCH --job-name=trait_null_anc
#SBATCH --account=orzuk
#SBATCH --partition=glacier
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --chdir=/sci/labs/orzuk/orzuk/github/MsaCluster
#SBATCH --output=trait_null_ancestral.out
# (A) Fold preference vs K=200 random sequence traits on the same coarse tree
#     -> docs/trait_null_phylo.csv, docs/trait_null_phylo_family.csv
# (B) 4-state ancestral gain/loss test + sequence-trait control
#     -> docs/ancestral_summary.csv (overwritten), docs/ancestral_events.csv
# Override the working tree with: sbatch --chdir=<worktree> scripts/run_trait_null_ancestral.sbatch
# (needs Pipeline/ (symlink ok) and docs/fold_diversity_survey*.csv).
set -e
source /sci/labs/orzuk/orzuk/venvs/my-python-venv/bin/activate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1   # 8 procs, 1 thread each
export MPLBACKEND=Agg QT_QPA_PLATFORM=offscreen PYTHONUNBUFFERED=1
echo "host $(hostname)  cwd $(pwd)  start $(date)"

echo "=== unit tests ==="
python3 tests/test_trait_null_ancestral.py
python3 tests/test_phylo_tests.py
python3 tests/test_power_negtraits.py

echo "=== (A) trait null: raw + corrected, K=200 ==="
python3 scripts/trait_null_phylo.py --pairs ALL --labels raw,corrected \
    --n-random 200 --min-clusters 10 --n-jobs 8

echo "=== (B) ancestral 4-state gain/loss + trait control (K_d=50) ==="
python3 scripts/run_ancestral_analysis.py --pairs ALL --n-perm 1000 \
    --n-trait 50 --timeout 300 --n-jobs 8

echo "=== DONE $(date) ==="
