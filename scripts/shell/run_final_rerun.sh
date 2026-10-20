#!/bin/bash
# Final CPU-only re-run: re-score Boltz-2 + BioEmu outputs (no prediction), then
# the downstream analysis chain (scripts/run_final_chain.sbatch) once all the
# re-score jobs have finished (afterany).
# Run from a compute node (srun), not the login node:
#   srun -A orzuk -p glacier --mem=8G -t 00:30:00 bash scripts/shell/run_final_rerun.sh
set -euo pipefail
source /sci/labs/orzuk/orzuk/venvs/my-python-venv/bin/activate
cd /sci/labs/orzuk/orzuk/github/MsaCluster

ACCOUNT=orzuk
PARTITION=glacier
mkdir -p Pipeline/FoldPairs/jobs
JOBIDS=Pipeline/FoldPairs/jobs/final_rerun_jobids_$(date +%Y%m%d_%H%M%S).txt
: > "$JOBIDS"

for mode in run_boltz2 run_bioemu; do
  echo "=== submitting $mode --rescore_only TRUE fan-out ==="
  python3 run_foldswitch_pipeline.py --run_mode "$mode" --foldpair_ids ALL \
      --rescore_only TRUE --run_job_mode sbatch \
      --sbatch_partition "$PARTITION" --sbatch_account "$ACCOUNT" \
      --sbatch_jobid_file "$JOBIDS"
done

N_JOBS=$(grep -c . "$JOBIDS" || true)
if [[ "$N_JOBS" -eq 0 ]]; then
  echo "[final-rerun] ERROR: no re-score jobs were submitted; not submitting the chain." >&2
  exit 1
fi
DEP="afterany:$(paste -sd: "$JOBIDS")"
CHAIN_ID=$(sbatch --parsable --dependency="$DEP" scripts/run_final_chain.sbatch)

echo
echo "=== submitted ==="
echo "  re-score jobs : $N_JOBS (ids in $JOBIDS)"
echo "  chain job     : $CHAIN_ID  (scripts/run_final_chain.sbatch, afterany on all re-score jobs)"
echo "  chain log     : Pipeline/FoldPairs/jobs/final_chain_${CHAIN_ID}.out"
echo "  re-score logs : Pipeline/FoldPairs/<pair>/jobs/run_boltz2_<pair>.out, run_bioemu_<pair>.out"
