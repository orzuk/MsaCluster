#!/bin/bash
#SBATCH --job-name=plot_global
#SBATCH --account=course-52017-25
#SBATCH --partition=glacier
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --chdir=/sci/labs/orzuk/orzuk/github/MsaCluster
#SBATCH --output=plot_global.out
# Regenerate global analysis figures + page from the current 9-method survey.
source /sci/labs/orzuk/orzuk/venvs/my-python-venv/bin/activate
python3 -u run_foldswitch_pipeline.py --run_mode plot --foldpair_ids ALL \
    --run_job_mode inline --plot_scope global
