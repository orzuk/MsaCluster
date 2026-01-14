#!/usr/bin/env bash
set -euo pipefail

# Usage: ./scripts/smoke_test.sh [PAIR_ID]
# Requires: ColabFold/MMseqs setup for get_msa; optional PyRosetta for deltaG (not used here).

PAIR_ID="${1:-1dzlA_5keqF}"

echo "[smoke] Pair: ${PAIR_ID}"

python3 run_foldswitch_pipeline.py --run_mode load --foldpair_ids "${PAIR_ID}"
python3 run_foldswitch_pipeline.py --run_mode get_msa --foldpair_ids "${PAIR_ID}" --run_job_mode inline
python3 run_foldswitch_pipeline.py --run_mode tree --foldpair_ids "${PAIR_ID}"
python3 run_foldswitch_pipeline.py --run_mode cluster_msa --foldpair_ids "${PAIR_ID}" --cluster_alg tree
python3 run_foldswitch_pipeline.py --run_mode postprocess --foldpair_ids "${PAIR_ID}" --cached_only FALSE --reports none

echo "[smoke] Done."
