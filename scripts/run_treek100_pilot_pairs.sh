#!/bin/bash
# ============================================================================
# run_treek100_pilot_pairs.sh -- launch K=100 + medoid + top-10 AF pilots
# on multiple pairs at once.
#
# For each pair in PAIRS:
#   1. (Optional) re-cluster at K=100 via tree-cut (skipped if dir exists)
#   2. Submit 1 sbatch job: AF2+AF3, --query_type medoid --af_msa_top_n 10
#   3. (Optional, for comparison) submit a second sbatch job with full MSA
#
# The clustering is fast (~1 min CPU). The AF jobs each run ~5-7h on
# salmon. If the cluster has GPU availability, all PAIRS x MODES jobs
# run in parallel.
#
# Usage:
#   bash scripts/run_treek100_pilot_pairs.sh                     # default 5 pairs, top10 only
#   bash scripts/run_treek100_pilot_pairs.sh --include-full      # also full-MSA runs
#   bash scripts/run_treek100_pilot_pairs.sh --pairs 2n54B_2hdmA # custom list
#   bash scripts/run_treek100_pilot_pairs.sh --dry-run           # print commands, don't submit
#
# After all jobs finish, score with:
#   for P in <pairs>; do
#       python3 scripts/score_treek100_pilot.py --pair $P
#   done
# ============================================================================
set -euo pipefail

# ---- defaults ---------------------------------------------------------------
PAIRS_DEFAULT="2n54B_2hdmA 1zk9A_3jv6A 2namA_1uxmK 6c6sD_2ougC 4qhhA_4qhfA"
INCLUDE_FULL=0
DRY_RUN=0
PAIRS=""
K=100

# ---- parse args -------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --include-full) INCLUDE_FULL=1; shift ;;
        --dry-run)      DRY_RUN=1; shift ;;
        --pairs)        PAIRS="$2"; shift 2 ;;
        --k)            K="$2"; shift 2 ;;
        -h|--help)
            echo "Usage: $0 [--include-full] [--dry-run] [--pairs 'P1 P2 ...'] [--k K]"
            exit 0
            ;;
        *)              echo "Unknown arg: $1"; exit 1 ;;
    esac
done
PAIRS="${PAIRS:-$PAIRS_DEFAULT}"

# ---- paths + slurm ---------------------------------------------------------
REPO_DIR="/sci/labs/orzuk/orzuk/github/MsaCluster"
VENV_DIR="/sci/labs/orzuk/orzuk/venvs/my-python-venv"
ACCOUNT="course-52017-25"
PARTITION="salmon"
GRES="gpu:l40s:1"
CPUS=8
MEM="40G"
TIME="12:00:00"

cd "$REPO_DIR"
# shellcheck disable=SC1091
source "$VENV_DIR/bin/activate"
mkdir -p Pipeline/FoldPairs/jobs

echo "[pilot-batch] pairs: $PAIRS"
echo "[pilot-batch] K=$K  include_full=$INCLUDE_FULL  dry_run=$DRY_RUN"
echo

# ---- Step 1: cluster each pair if needed -----------------------------------
echo "[pilot-batch] === Step 1: tree-cut clustering at K=$K ==="
for P in $PAIRS; do
    CL_DIR="Pipeline/FoldPairs/$P/output_msa_cluster"
    # check if K=$K dir already exists (heuristic: cluster count near K)
    n_existing=$(ls -d "$CL_DIR"/ShallowMsa_*/ 2>/dev/null | wc -l || echo 0)
    if [[ "$n_existing" -ge $((K - 10)) ]] && [[ "$n_existing" -le $((K + 10)) ]]; then
        echo "[pilot-batch]   $P: already has $n_existing clusters (~K=$K); SKIP re-cluster"
        continue
    fi
    echo "[pilot-batch]   $P: re-cluster (was $n_existing, target K=$K)"
    CMD="python3 run_foldswitch_pipeline.py \
        --run_mode cluster_msa --foldpair_ids $P \
        --cluster_alg tree \
        --cluster_tree_min_num_clusters $K \
        --cluster_max_clusters $K \
        --cluster_min_output_size 30 \
        --cluster_min_neff 20 \
        --run_job_mode inline"
    echo "    $CMD"
    if [[ "$DRY_RUN" -eq 0 ]]; then
        eval "$CMD" 2>&1 | tail -3
    fi
done
echo

# ---- Step 2: submit AF runs ------------------------------------------------
echo "[pilot-batch] === Step 2: submit AF2+AF3 sbatch jobs ==="
N_JOBS=0
for P in $PAIRS; do
    # Top-10 mode (the headline methodology)
    CMD_INNER_TOP10="python3 run_foldswitch_pipeline.py \
        --run_mode run_AF --foldpair_ids $P \
        --af_ver both --query_type medoid --af_msa_top_n 10 \
        --af_output_suffix treek${K}_top10 \
        --run_job_mode inline --allow_inline_af"
    LOG_TOP10="Pipeline/FoldPairs/jobs/${P}_treek${K}_top10.out"
    SBATCH_TOP10="sbatch -J ${P}_top10 \
        --account=$ACCOUNT --partition=$PARTITION --gres=$GRES \
        --cpus-per-task=$CPUS --mem=$MEM --time=$TIME \
        --output=$LOG_TOP10 \
        --wrap=\"$CMD_INNER_TOP10\""
    echo "[pilot-batch]   $P top10:"
    echo "    $SBATCH_TOP10"
    if [[ "$DRY_RUN" -eq 0 ]]; then
        eval "$SBATCH_TOP10"
        N_JOBS=$((N_JOBS + 1))
    fi

    # Optional: full-MSA mode (for full-vs-top10 comparison)
    if [[ "$INCLUDE_FULL" -eq 1 ]]; then
        CMD_INNER_FULL="python3 run_foldswitch_pipeline.py \
            --run_mode run_AF --foldpair_ids $P \
            --af_ver both --query_type medoid \
            --af_output_suffix treek${K}_full \
            --run_job_mode inline --allow_inline_af"
        LOG_FULL="Pipeline/FoldPairs/jobs/${P}_treek${K}_full.out"
        SBATCH_FULL="sbatch -J ${P}_full \
            --account=$ACCOUNT --partition=$PARTITION --gres=$GRES \
            --cpus-per-task=$CPUS --mem=$MEM --time=$TIME \
            --output=$LOG_FULL \
            --wrap=\"$CMD_INNER_FULL\""
        echo "[pilot-batch]   $P full:"
        echo "    $SBATCH_FULL"
        if [[ "$DRY_RUN" -eq 0 ]]; then
            eval "$SBATCH_FULL"
            N_JOBS=$((N_JOBS + 1))
        fi
    fi
done

echo
echo "[pilot-batch] Submitted $N_JOBS sbatch jobs"
echo
cat <<EOF
Monitor:
  squeue -u \$USER -h

After completion, score each pair:
  for P in $PAIRS; do
      python3 scripts/score_treek100_pilot.py --pair \$P
  done

Then build per-pair phytree figures:
  for P in $PAIRS; do
      python3 scripts/plot_treek100_phytree.py --pair \$P
  done
EOF
