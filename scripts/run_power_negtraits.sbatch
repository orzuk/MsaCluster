#!/bin/bash
#SBATCH --job-name=power_negtraits
#SBATCH --account=orzuk
#SBATCH --partition=glacier
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --chdir=/sci/labs/orzuk/orzuk/github/MsaCluster
#SBATCH --output=power_negtraits.out
# Referee analyses #5 (power simulation) and #6 (negative-trait baseline).
# Outputs: docs/power_simulation*.csv, docs/negative_traits_*.csv,
#          docs/figures/global/{power_simulation,negative_traits_morans}.png
set -e
source /sci/labs/orzuk/orzuk/venvs/my-python-venv/bin/activate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1   # 8 procs, 1 thread each
export MPLBACKEND=Agg QT_QPA_PLATFORM=offscreen

echo "=== unit tests ==="
python3 tests/test_power_negtraits.py

echo "=== #6 negative traits (all pairs, Moran 1000 perms as in the paper) ==="
python3 scripts/negative_traits_phylo.py --pairs ALL --n-jobs 8

echo "=== #5 power simulation: plan + timing ==="
# --n-trees 0 = every pair with >= 20 clusters (~74 trees); still well inside 6h.
# Use --n-trees 20 for the referee-spec subset.
python3 scripts/power_simulation_phylo.py --pairs ALL --n-trees 0 --n-jobs 8 --dry-run

echo "=== #5 power simulation: run ==="
python3 scripts/power_simulation_phylo.py --pairs ALL --n-trees 0 --n-jobs 8 \
    --n-reps 100 --morans-perm 199 --clade-perm 999
echo "=== DONE ==="
