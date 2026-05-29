#!/bin/bash
# ============================================================================
# run_pdb_mine.sh -- launch the PDB fold-switch candidate mining as sbatch.
#
# Wraps scripts/pdb_mine_foldswitch_candidates.py for one-command submission.
# CPU-only (glacier partition); 8 workers; resume-safe (idempotent on rerun).
#
# Usage:
#   bash scripts/run_pdb_mine.sh                # full PDB scan, ~6-24h
#   bash scripts/run_pdb_mine.sh --max 100      # test run, ~30 min
#   bash scripts/run_pdb_mine.sh --workers 16   # more parallelism
#   bash scripts/run_pdb_mine.sh --dry-run      # print sbatch command, don't submit
#
# Outputs go to:
#   log:     Pipeline/FoldPairs/jobs/pdb_mine.out
#   results: data/pdb_mining_candidates.csv
#            data/pdb_mining_candidates.new_only.csv  (deduped vs existing 93)
#
# Monitor:
#   squeue -u $USER
#   tail -f Pipeline/FoldPairs/jobs/pdb_mine.out
# ============================================================================
set -eo pipefail

WORKERS=8
MAX_CLUSTERS=""
DRY_RUN=0
EXTRA_ARGS=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --workers)  WORKERS="$2"; shift 2 ;;
        --max)      MAX_CLUSTERS="--max-clusters $2"; shift 2 ;;
        --dry-run)  DRY_RUN=1; shift ;;
        *)          EXTRA_ARGS+=" $1"; shift ;;
    esac
done

REPO_DIR="/sci/labs/orzuk/orzuk/github/MsaCluster"
VENV_DIR="/sci/labs/orzuk/orzuk/venvs/my-python-venv"
LOG="$REPO_DIR/Pipeline/FoldPairs/jobs/pdb_mine.out"
mkdir -p "$(dirname "$LOG")"

# Inner command runs in the sbatch job.
INNER="source $VENV_DIR/bin/activate \
    && cd $REPO_DIR \
    && python3 scripts/pdb_mine_foldswitch_candidates.py \
       --workers $WORKERS --resume $MAX_CLUSTERS $EXTRA_ARGS"

SBATCH_CMD="sbatch -J pdb_mine \
    --partition=glacier --cpus-per-task=$WORKERS --mem=16G --time=24:00:00 \
    --output=$LOG \
    --wrap=\"bash -c '$INNER'\""

echo "[pdb-mine] log:    $LOG"
echo "[pdb-mine] inner:  $INNER"
echo

if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "[DRY-RUN] $SBATCH_CMD"
    exit 0
fi

eval "$SBATCH_CMD"
echo
echo "Monitor:"
echo "  squeue -u \$USER"
echo "  tail -f $LOG"
echo
echo "When it finishes, view candidates:"
echo "  wc -l data/pdb_mining_candidates.new_only.csv"
echo "  head -20 data/pdb_mining_candidates.new_only.csv"
