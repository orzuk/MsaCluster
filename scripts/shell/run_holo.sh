#!/bin/bash
# Submit condition-aware (holo) AF3 + Boltz-2 runs for the given fold pairs.
# Run from a compute node (srun), not the login node:
#   srun -A orzuk -p glacier --mem=8G -t 00:30:00 bash scripts/shell/run_holo.sh 2lelA_2k0qA [more pairs...]
# Methods to run: HOLO_METHODS="af3 boltz2" (default both).
# Outputs: output_AF/AF3_holo + Analysis/df_af_holo.csv; output_boltz2_holo + Analysis/df_boltz2_holo.csv.
# Then: python3 scripts/fold_diversity_survey.py --condition-switch
set -euo pipefail
source /sci/labs/orzuk/orzuk/venvs/my-python-venv/bin/activate
cd /sci/labs/orzuk/orzuk/github/MsaCluster

PAIRS="$*"
[[ -z "$PAIRS" ]] && { echo "usage: $0 <pair_id> [pair_id ...]"; exit 2; }
METHODS="${HOLO_METHODS:-af3 boltz2}"
SB=(--run_job_mode sbatch --sbatch_partition "${HOLO_PARTITION:-salmon}"
    --sbatch_gres "${HOLO_GRES:-gpu:l40s:1}" --sbatch_account "${HOLO_ACCOUNT:-orzuk}")

for p in $PAIRS; do
  python3 -c "from utils.trigger_utils import condition_spec; print('$p', condition_spec('$p'))"
done
if [[ " $METHODS " == *" af3 "* ]]; then
  python3 run_foldswitch_pipeline.py --run_mode run_AF --foldpair_ids $PAIRS \
      --af_ver 3 --af3_mode holo --af_output_suffix holo "${SB[@]}"
fi
if [[ " $METHODS " == *" boltz2 "* ]]; then
  python3 run_foldswitch_pipeline.py --run_mode run_boltz2 --foldpair_ids $PAIRS \
      --boltz2_mode holo "${SB[@]}"
fi
