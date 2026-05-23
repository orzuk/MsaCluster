#!/bin/bash
#SBATCH --job-name=postproc_all
#SBATCH --partition=glacier
#SBATCH --account=course-52017-25
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=Pipeline/FoldPairs/jobs/postproc_all.%j.out
#SBATCH --error=Pipeline/FoldPairs/jobs/postproc_all.%j.err
#
# run_postprocess_all.sh
# ----------------------
# Run Analysis/postprocess_unified.py on all 93 fold-switching pairs.
# Recomputes per-pair Analysis CSVs (df_af.csv with AF2+AF3 rows,
# df_esm.csv, df_cmap.csv with MSA-Transformer, df_cmap_ccmpred.csv
# with CCMpred, df_ddg.csv) then rebuilds the unified docs/ tables.
#
# After this finishes, the next step is to run scripts/fold_diversity_survey.py
# which produces docs/fold_diversity_survey.csv from these per-pair CSVs.
#
# Usage:
#   sbatch scripts/run_postprocess_all.sh
#
# Monitor:
#   squeue -u $USER
#   tail -f Pipeline/FoldPairs/jobs/postproc_all.<jobid>.out

set -euo pipefail
cd /sci/labs/orzuk/orzuk/github/MsaCluster
. /sci/labs/orzuk/orzuk/venvs/my-python-venv/bin/activate

echo "== run_postprocess_all start: $(date) =="
python3 Analysis/postprocess_unified.py --pairs ALL --force_rerun
echo "== run_postprocess_all done: $(date) =="

# Quick coverage report after the run
echo
echo "== per-pair CSV coverage =="
for fn in df_af.csv df_esm.csv df_cmap.csv df_cmap_ccmpred.csv df_ddg.csv; do
  n=$(ls Pipeline/FoldPairs/*/Analysis/$fn 2>/dev/null | wc -l)
  printf "  %-25s : %3d / 93\n" "$fn" "$n"
done
