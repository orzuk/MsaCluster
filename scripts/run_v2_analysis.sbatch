#!/bin/bash
#SBATCH --job-name=v2_analysis
#SBATCH --account=course-52017-25
#SBATCH --partition=glacier
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --output=Pipeline/FoldPairs/jobs/v2_analysis_%j.out
# ===========================================================================
# Batch wrapper for the v2 downstream analysis
# (ESMFold2 scoring -> survey -> S4PRED -> concordance perm/BH -> stratified
#  -> seq-divergence correction -> Moran's phylo placement raw+corrected
#  -> cross-method correlation).
#
# Use this instead of `sbatch --wrap "bash -c '...'"`: a pasted newline
# between `bash` and the script path there caused ExitCode 126
# (Permission denied). A real job script has no such quoting hazard.
#
# Submit from the repo root:
#   sbatch scripts/run_v2_analysis.sbatch                  # defaults: --no-html --tag v2_82
#   sbatch scripts/run_v2_analysis.sbatch --no-html --tag v2_90   # override args
#
# Log: Pipeline/FoldPairs/jobs/v2_analysis_<jobid>.out
# ===========================================================================
set -uo pipefail

REPO=/sci/labs/orzuk/orzuk/github/MsaCluster
VENV=/sci/labs/orzuk/orzuk/venvs/my-python-venv

cd "$REPO"
# shellcheck disable=SC1091
source "$VENV/bin/activate"

# Pass-through args; default to the standard stats-only run if none given.
ARGS=("$@")
if [[ ${#ARGS[@]} -eq 0 ]]; then
    ARGS=(--no-html --tag v2_82)
fi

echo "[v2-analysis-job] host=$(hostname)  start=$(date)"
echo "[v2-analysis-job] commit=$(git rev-parse --short HEAD)"
echo "[v2-analysis-job] args: ${ARGS[*]}"
echo

bash scripts/run_v2_subset_analysis.sh "${ARGS[@]}"
rc=$?

echo
echo "[v2-analysis-job] run_v2_subset_analysis.sh exit=$rc  done=$(date)"
exit "$rc"
