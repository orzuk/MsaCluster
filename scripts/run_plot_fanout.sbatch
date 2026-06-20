#!/bin/bash
#SBATCH --job-name=plot_fanout
#SBATCH --account=course-52017-25
#SBATCH --partition=glacier
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --chdir=/sci/labs/orzuk/orzuk/github/MsaCluster
#SBATCH --output=plot_fanout.out
# Submitter job: fans out one per-pair plot job per pair (scope=pair skips the
# global block). Runs as a job because the login node OOM-kills the pipeline.
# Pass partition/account through so the child jobs inherit them too.
source /sci/labs/orzuk/orzuk/venvs/my-python-venv/bin/activate
python3 -u run_foldswitch_pipeline.py --run_mode plot --foldpair_ids ALL \
    --run_job_mode sbatch --plot_scope pair --plot_null_permutation TRUE \
    --sbatch_partition glacier --sbatch_account course-52017-25
