#!/bin/bash
# ============================================================================
# esmfold2_smoke.sh -- ESMFold2 smoke test on one KaiB cluster medoid.
#
# Submits an sbatch GPU job that:
#   1. Activates the esmfold2-venv
#   2. cd's into the repo
#   3. Runs run_ESMFold2.py on cluster 0 of KaiB
#   4. Should produce Pipeline/FoldPairs/5jytA_2qkeE/output_esmfold2/ShallowMsa_000.pdb
#
# First run downloads ~12 GB of biohub/ESMFold2 weights from HuggingFace;
# expect ~5-10 minutes before the first fold starts. Subsequent jobs reuse
# the cached weights.
#
# Usage:
#   bash scripts/esmfold2_smoke.sh
# ============================================================================
set -euo pipefail

ACCOUNT="course-52017-25"
PARTITION="salmon"
GRES="gpu:l40s:1"
VENV="/sci/labs/orzuk/orzuk/venvs/esmfold2-venv"
REPO="/sci/labs/orzuk/orzuk/github/MsaCluster"
LOG="${REPO}/Pipeline/FoldPairs/jobs/esmfold2_smoke.out"

mkdir -p "$(dirname "$LOG")"

# Build the inner command as a single line to avoid newline-in-quoted-string bugs
INNER="source ${VENV}/bin/activate && cd ${REPO} && python3 run_ESMFold2.py --foldpair_ids 5jytA_2qkeE --max_clusters 1 --device cuda"

echo "Submitting sbatch ESMFold2 smoke test..."
echo "  log: $LOG"
echo

sbatch -J esmfold2_smoke \
    --account="$ACCOUNT" --partition="$PARTITION" --gres="$GRES" \
    --cpus-per-task=4 --mem=24G --time=0:30:00 \
    --output="$LOG" \
    --wrap="bash -c \"$INNER\""

echo
echo "Watch progress:"
echo "  squeue -u \$USER"
echo "  tail -f $LOG"
echo
echo "When the job finishes:"
echo "  cat $LOG"
echo "  ls -lh Pipeline/FoldPairs/5jytA_2qkeE/output_esmfold2/"
