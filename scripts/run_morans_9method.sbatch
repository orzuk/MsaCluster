#!/bin/bash
#SBATCH --job-name=morans_9m
#SBATCH --account=course-52017-25
#SBATCH --partition=glacier
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:30:00
#SBATCH --chdir=/sci/labs/orzuk/orzuk/github/MsaCluster
#SBATCH --output=morans_9method.out
set -e
source /sci/labs/orzuk/orzuk/venvs/my-python-venv/bin/activate
ALL9="AF2,AF3,Boltz2,BioEmu,ESM,MSAT,CCMpred,DDG,S4PRED"

echo "=== backup current phylo CSVs ==="
cp -f docs/phylo_placement.csv           docs/phylo_placement.csv.bak8 2>/dev/null || true
cp -f docs/phylo_placement_corrected.csv docs/phylo_placement_corrected.csv.bak8 2>/dev/null || true

echo "=== STEP 1: 9-method corrected survey (separate file, canonical untouched) ==="
python3 scripts/seq_divergence_correction.py --methods "$ALL9" --out-suffix _9method

echo "=== STEP 2: RAW Moran's, 9 methods (BioEmu now in SCORED_METHODS) ==="
python3 scripts/phylo_placement.py --pairs ALL --no-d-statistic

echo "=== STEP 3: CORRECTED Moran's, 9 methods (from the 9-method corrected survey) ==="
python3 scripts/phylo_placement.py --pairs ALL --labels corrected --no-d-statistic \
    --survey-csv docs/fold_diversity_survey_corrected_9method.csv

echo "=== STEP 4: regenerate global plots ==="
python3 run_foldswitch_pipeline.py --run_mode plot --foldpair_ids ALL \
    --run_job_mode inline --plot_scope global

echo "=== DONE; methods now in phylo CSVs: ==="
python3 -c "import pandas as pd; print('raw      :', sorted(pd.read_csv('docs/phylo_placement.csv').method.unique())); print('corrected:', sorted(pd.read_csv('docs/phylo_placement_corrected.csv').method.unique()))"
