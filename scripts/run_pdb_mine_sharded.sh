#!/bin/bash
# ============================================================================
# run_pdb_mine_sharded.sh -- launch N parallel sbatch jobs for PDB mining.
#
# Each shard processes 1/N of the clusters (cluster_idx % N == shard-1).
# Outputs go to data/pdb_mining_candidates.shard_XofN.csv per shard.
# Merge with scripts/merge_pdb_mining_shards.sh after all shards finish.
#
# Resume-safe: re-running this script picks up where each shard left off
# (because each shard's --resume reads its own .shard_XofN.csv).
# Per-shard wallclock ~ (total / N), so 8 shards finish ~8x faster
# than a single sequential job.
#
# Usage:
#   bash scripts/run_pdb_mine_sharded.sh              # default N=8
#   bash scripts/run_pdb_mine_sharded.sh 16           # N=16
#   bash scripts/run_pdb_mine_sharded.sh 8 --max 1000 # 8 shards, 1000-cluster test
# ============================================================================
set -eo pipefail

N=8
EXTRA_ARGS=""

# First positional arg is N if it's a number; otherwise treat as extra
if [[ "$1" =~ ^[0-9]+$ ]]; then
    N="$1"; shift
fi
# Remaining args pass through to python
EXTRA_ARGS="$@"

REPO_DIR="/sci/labs/orzuk/orzuk/github/MsaCluster"
VENV_DIR="/sci/labs/orzuk/orzuk/venvs/my-python-venv"
JOBS_DIR="$REPO_DIR/Pipeline/FoldPairs/jobs"
mkdir -p "$JOBS_DIR"

echo "[sharded] Launching $N parallel pdb_mine shards"
echo "[sharded] Extra args: $EXTRA_ARGS"
echo

SUBMITTED=()
for X in $(seq 1 "$N"); do
    LOG="$JOBS_DIR/pdb_mine.shard_${X}of${N}.out"
    INNER="source $VENV_DIR/bin/activate \
        && cd $REPO_DIR \
        && python3 scripts/pdb_mine_foldswitch_candidates.py \
           --workers 8 --resume --shard ${X}/${N} $EXTRA_ARGS"
    JOB=$(sbatch -J "pdb_mine_${X}of${N}" \
        --partition=glacier --cpus-per-task=8 --mem=16G --time=24:00:00 \
        --output="$LOG" \
        --wrap="bash -c '$INNER'" | awk '{print $NF}')
    SUBMITTED+=("$JOB")
    echo "[sharded] shard $X/$N -> job $JOB  log=$LOG"
done

echo
echo "Submitted ${#SUBMITTED[@]} shards: ${SUBMITTED[*]}"
echo
echo "Monitor:"
echo "  squeue -u \$USER -h --format='%i %T %M' | grep pdb_mine"
echo "  for X in \$(seq 1 $N); do echo === shard \$X ===; tail -3 $JOBS_DIR/pdb_mine.shard_\${X}of${N}.out; done"
echo
echo "When all shards finish, merge with:"
echo "  bash scripts/merge_pdb_mining_shards.sh $N"
