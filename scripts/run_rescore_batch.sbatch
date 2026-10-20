#!/bin/bash
#SBATCH --job-name=rescore_batch
#SBATCH --account=orzuk
#SBATCH --partition=glacier
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --chdir=/sci/labs/orzuk/orzuk/github/MsaCluster
#SBATCH --output=Pipeline/FoldPairs/jobs/rescore_batch_%j.out
# ===========================================================================
# Re-score existing Boltz-2 or BioEmu predictions (no prediction / sampling)
# for many pairs inside ONE job, running pairs in parallel. Use instead of the
# per-pair --run_job_mode sbatch fan-out when the cluster's jobs-per-hour
# submission limit would be hit (~60/h on Moriah).
#
#   sbatch scripts/run_rescore_batch.sbatch boltz2 <pair_id> [pair_id ...]
#   sbatch scripts/run_rescore_batch.sbatch bioemu ALL
#
# Each pair runs `run_foldswitch_pipeline.py --run_mode run_<method>
# --run_job_mode inline --rescore_only TRUE`; per-pair logs go to
# Pipeline/FoldPairs/<pair>/jobs/rescore_<method>_<pair>.out.
# ===========================================================================
set -euo pipefail
source /sci/labs/orzuk/orzuk/venvs/my-python-venv/bin/activate

METHOD="${1:?usage: run_rescore_batch.sbatch {boltz2|bioemu} <pairs...|ALL>}"
shift
case "$METHOD" in boltz2|bioemu) ;; *) echo "unknown method $METHOD"; exit 2 ;; esac

if [[ "${1:-}" == "ALL" ]]; then
    mapfile -t PAIRS < <(ls -d Pipeline/FoldPairs/*/ | xargs -n1 basename \
                         | grep -E '^[0-9a-zA-Z]{5}_[0-9a-zA-Z]{5}$')
else
    PAIRS=("$@")
fi
NPAR="${SLURM_CPUS_PER_TASK:-8}"
echo "[rescore-batch] method=$METHOD pairs=${#PAIRS[@]} parallel=$NPAR host=$(hostname) start=$(date)"

rescore_one() {
    local method="$1" pair="$2"
    mkdir -p "Pipeline/FoldPairs/$pair/jobs"
    if python3 run_foldswitch_pipeline.py --run_mode "run_$method" --foldpair_ids "$pair" \
            --run_job_mode inline --rescore_only TRUE \
            > "Pipeline/FoldPairs/$pair/jobs/rescore_${method}_${pair}.out" 2>&1; then
        echo "OK   $pair"
    else
        echo "FAIL $pair (see Pipeline/FoldPairs/$pair/jobs/rescore_${method}_${pair}.out)"
    fi
}
export -f rescore_one

printf '%s\n' "${PAIRS[@]}" | xargs -P "$NPAR" -I{} bash -c 'rescore_one "$0" "$1"' "$METHOD" {}
echo "[rescore-batch] done=$(date)"
