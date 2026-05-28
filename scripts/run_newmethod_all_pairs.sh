#!/bin/bash
# ============================================================================
# run_newmethod_all_pairs.sh -- v2 methodology on all 85 fold-switch pairs.
#
# Methodology (v2 paper):
#   1. Tree-cut clustering with adaptive K = clip(N_seqs/30, 20, 100)
#   2. Variant-as-query: cluster MEDOID replaces the original chain
#   3. Top-10 nearest-by-Hamming sub-sample as AF MSA input
#
# Stages (each is a separate sbatch dispatch):
#   --stage cluster   : cluster_msa for all pairs (CPU, ~1 min/pair)
#   --stage af        : run_AF (AF2+AF3, both chains) for all pairs (GPU)
#   --stage esmfold2  : run_ESMFold2 on cluster medoids (GPU, ~10 min/pair)
#   --stage all       : all three stages in dependency order
#
# Each pair's AF outputs land in:
#   Pipeline/FoldPairs/<pair>/output_AF/AF2_newmethod_top10/
#   Pipeline/FoldPairs/<pair>/output_AF/AF3_newmethod_top10/
# ESMFold2 outputs:
#   Pipeline/FoldPairs/<pair>/output_esmfold2/
#
# Usage:
#   bash scripts/run_newmethod_all_pairs.sh --stage cluster
#   bash scripts/run_newmethod_all_pairs.sh --stage af
#   bash scripts/run_newmethod_all_pairs.sh --stage esmfold2
#   bash scripts/run_newmethod_all_pairs.sh --stage all
#   bash scripts/run_newmethod_all_pairs.sh --stage cluster --dry-run
#   bash scripts/run_newmethod_all_pairs.sh --stage af --pairs '5jytA_2qkeE 2n54B_2hdmA'
# ============================================================================
set -euo pipefail

# ---- defaults ---------------------------------------------------------------
STAGE=""
DRY_RUN=0
PAIRS=""
SEQS_PER_CLUSTER=30
K_MIN=20
K_MAX=100
SUFFIX="newmethod_top10"
# v2 methodology: top-10 sampling means clusters with >=10 seqs are usable,
# and Neff requirements relax because AF input is a shallow MSA.
MIN_OUTPUT_SIZE=10
MIN_NEFF=5
# Coarse-resolution params (for CCMpred / MSATransformer).
# K_min=2 / min_size=30 / min_neff=10 chosen so every fine-viable pair
# also has a coarse clustering. Depth (~200 seqs/cluster) comes from the
# K_max=30 cap not from filters.
SEQS_PER_CLUSTER_COARSE=200
K_MIN_COARSE=2
K_MAX_COARSE=30
MIN_OUTPUT_SIZE_COARSE=30
MIN_NEFF_COARSE=10

# ---- parse args -------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --stage)               STAGE="$2"; shift 2 ;;
        --dry-run)             DRY_RUN=1; shift ;;
        --pairs)               PAIRS="$2"; shift 2 ;;
        --seqs-per-cluster)    SEQS_PER_CLUSTER="$2"; shift 2 ;;
        --k-min)               K_MIN="$2"; shift 2 ;;
        --k-max)               K_MAX="$2"; shift 2 ;;
        --min-output-size)     MIN_OUTPUT_SIZE="$2"; shift 2 ;;
        --min-neff)            MIN_NEFF="$2"; shift 2 ;;
        --suffix)              SUFFIX="$2"; shift 2 ;;
        -h|--help)
            sed -n '2,30p' "$0"
            exit 0
            ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done
[[ -z "$STAGE" ]] && { echo "Must give --stage cluster|af|esmfold2|all"; exit 1; }

# ---- paths + slurm ----------------------------------------------------------
REPO_DIR="/sci/labs/orzuk/orzuk/github/MsaCluster"
VENV_DIR="/sci/labs/orzuk/orzuk/venvs/my-python-venv"
ESMFOLD2_VENV="/sci/labs/orzuk/orzuk/venvs/esmfold2-venv"
ACCOUNT="course-52017-25"
PARTITION="salmon"
GRES_GPU="gpu:l40s:1"

cd "$REPO_DIR"
# shellcheck disable=SC1091
source "$VENV_DIR/bin/activate"
mkdir -p Pipeline/FoldPairs/jobs

# ---- resolve which pairs to process ----------------------------------------
if [[ -n "$PAIRS" ]]; then
    PAIRS_ARG="$PAIRS"
    echo "[v2] Restricted to pairs: $PAIRS"
else
    PAIRS_ARG="ALL"
    echo "[v2] Using ALL pairs (auto-discovered from Pipeline/FoldPairs)"
fi
echo "[v2] Adaptive-K policy: clip(N_seqs/$SEQS_PER_CLUSTER, $K_MIN, $K_MAX)"
echo "[v2] Output suffix: $SUFFIX"
echo "[v2] Dry-run: $DRY_RUN"
echo

run_or_echo() {
    if [[ "$DRY_RUN" -eq 1 ]]; then
        echo "[DRY-RUN] $*"
    else
        echo "[run]     $*"
        eval "$@"
    fi
}

# ---- Stage 1: fine-resolution cluster_msa ----------------------------------
if [[ "$STAGE" == "cluster" ]] || [[ "$STAGE" == "cluster_fine" ]] \
   || [[ "$STAGE" == "all" ]]; then
    echo "==================================================================="
    echo "[v2] STAGE 1a: cluster_msa FINE resolution (K=$K_MIN-$K_MAX)"
    echo "==================================================================="
    CMD="python3 run_foldswitch_pipeline.py \
        --run_mode cluster_msa --foldpair_ids $PAIRS_ARG \
        --cluster_alg tree \
        --cluster_resolution fine \
        --cluster_adaptive_k_seqs_per_cluster $SEQS_PER_CLUSTER \
        --cluster_adaptive_k_min $K_MIN \
        --cluster_adaptive_k_max $K_MAX \
        --cluster_min_output_size $MIN_OUTPUT_SIZE \
        --cluster_min_neff $MIN_NEFF \
        --run_job_mode sbatch"
    run_or_echo "$CMD"
    echo
fi

# ---- Stage 1b: coarse-resolution cluster_msa -------------------------------
if [[ "$STAGE" == "cluster_coarse" ]] || [[ "$STAGE" == "cluster_both" ]] \
   || [[ "$STAGE" == "all" ]]; then
    echo "==================================================================="
    echo "[v2] STAGE 1b: cluster_msa COARSE resolution (K=$K_MIN_COARSE-$K_MAX_COARSE)"
    echo "==================================================================="
    CMD="python3 run_foldswitch_pipeline.py \
        --run_mode cluster_msa --foldpair_ids $PAIRS_ARG \
        --cluster_alg tree \
        --cluster_resolution coarse \
        --cluster_adaptive_k_seqs_per_cluster_coarse $SEQS_PER_CLUSTER_COARSE \
        --cluster_adaptive_k_min_coarse $K_MIN_COARSE \
        --cluster_adaptive_k_max_coarse $K_MAX_COARSE \
        --cluster_min_output_size_coarse $MIN_OUTPUT_SIZE_COARSE \
        --cluster_min_neff_coarse $MIN_NEFF_COARSE \
        --run_job_mode sbatch"
    run_or_echo "$CMD"
    echo
fi

# ---- Stage 1c: fine->coarse mapping (after both finish) --------------------
if [[ "$STAGE" == "cluster_mapping" ]] || [[ "$STAGE" == "cluster_both" ]] \
   || [[ "$STAGE" == "all" ]]; then
    echo "==================================================================="
    echo "[v2] STAGE 1c: build fine->coarse mapping per pair"
    echo "==================================================================="
    if [[ "$PAIRS_ARG" == "ALL" ]]; then
        PAIR_LIST=$(ls -d Pipeline/FoldPairs/*/ 2>/dev/null \
                    | sed -E 's#.*FoldPairs/##; s#/##' \
                    | grep -Ev '^(_s4pred_work|jobs)$')
    else
        PAIR_LIST="$PAIRS_ARG"
    fi
    for P in $PAIR_LIST; do
        FINE="Pipeline/FoldPairs/$P/output_msa_cluster"
        COARSE="Pipeline/FoldPairs/$P/output_msa_cluster_coarse"
        if [[ ! -d "$FINE" ]] || [[ ! -d "$COARSE" ]]; then continue; fi
        CMD="python3 scripts/build_fine_to_coarse_mapping.py \
            --fine-dir $FINE --coarse-dir $COARSE \
            --output $FINE/fine_to_coarse.csv"
        run_or_echo "$CMD"
    done
    echo
fi

# ---- Stage 2: CCMpred at COARSE resolution ---------------------------------
if [[ "$STAGE" == "ccmpred_coarse" ]] || [[ "$STAGE" == "ccmpred" ]]; then
    echo "==================================================================="
    echo "[v2] STAGE: CCMpred at COARSE resolution"
    echo "==================================================================="
    CMD="python3 run_foldswitch_pipeline.py \
        --run_mode run_cmap_ccmpred --foldpair_ids $PAIRS_ARG \
        --cluster_resolution coarse \
        --run_job_mode sbatch"
    run_or_echo "$CMD"
    echo
fi

# ---- Stage 2: MSA-Transformer at COARSE resolution -------------------------
if [[ "$STAGE" == "msat_coarse" ]] || [[ "$STAGE" == "msat" ]]; then
    echo "==================================================================="
    echo "[v2] STAGE: MSA-Transformer at COARSE resolution"
    echo "==================================================================="
    CMD="python3 run_foldswitch_pipeline.py \
        --run_mode run_cmap_msa_transformer --foldpair_ids $PAIRS_ARG \
        --cluster_resolution coarse \
        --run_job_mode sbatch"
    run_or_echo "$CMD"
    echo
fi

# ---- Stage 2: AF2+AF3 with medoid + top-10 ---------------------------------
if [[ "$STAGE" == "af" ]] || [[ "$STAGE" == "all" ]]; then
    echo "==================================================================="
    echo "[v2] STAGE 2: AF2+AF3 with medoid-as-query + top-10 sampling"
    echo "==================================================================="
    CMD="python3 run_foldswitch_pipeline.py \
        --run_mode run_AF --foldpair_ids $PAIRS_ARG \
        --af_ver both \
        --query_type medoid \
        --af_msa_top_n 10 \
        --af_output_suffix $SUFFIX \
        --run_job_mode sbatch"
    run_or_echo "$CMD"
    echo
fi

# ---- Stage 3: ESMFold2 on cluster medoids ----------------------------------
if [[ "$STAGE" == "esmfold2" ]] || [[ "$STAGE" == "all" ]]; then
    echo "==================================================================="
    echo "[v2] STAGE 3: ESMFold2 on cluster medoids (single-sequence baseline)"
    echo "==================================================================="
    if [[ "$PAIRS_ARG" == "ALL" ]]; then
        # Enumerate pairs from FoldPairs/ (those with output_msa_cluster/)
        PAIR_LIST=$(ls -d Pipeline/FoldPairs/*/ 2>/dev/null \
                    | sed -E 's#.*FoldPairs/##; s#/##' \
                    | while read p; do
                        [[ -d "Pipeline/FoldPairs/$p/output_msa_cluster" ]] && echo "$p"
                      done)
    else
        PAIR_LIST="$PAIRS_ARG"
    fi
    n_jobs=0
    for P in $PAIR_LIST; do
        LOG="Pipeline/FoldPairs/jobs/${P}_esmfold2.out"
        CMD="sbatch -J ${P}_esmfold2 \
            --account=$ACCOUNT --partition=$PARTITION --gres=$GRES_GPU \
            --cpus-per-task=4 --mem=24G --time=2:00:00 \
            --output=$LOG \
            --wrap=\"bash -c 'source $ESMFOLD2_VENV/bin/activate && cd $REPO_DIR && python3 run_ESMFold2.py --foldpair_ids $P --skip_existing'\""
        run_or_echo "$CMD"
        n_jobs=$((n_jobs+1))
    done
    echo "[v2] Submitted $n_jobs ESMFold2 jobs"
    echo
fi

echo "==================================================================="
cat <<EOF
[v2] Done with stage(s): $STAGE

Monitor:
  squeue -u \$USER -h | wc -l
  squeue -u \$USER -h --format='%i %T %M %j' | head

After Stage 1 finishes, run Stage 2:
  bash scripts/run_newmethod_all_pairs.sh --stage af

After Stage 2 finishes, score:
  for P in \$(ls -d Pipeline/FoldPairs/*/ | sed -E 's#.*FoldPairs/##; s#/##'); do
      python3 scripts/score_treek100_pilot.py --pair \$P \\
          --modes AF2_${SUFFIX},AF3_${SUFFIX} \\
          --output docs/newmethod_tmscores_\${P}.csv
  done

Then summarize:
  python3 scripts/summarize_treek100_pilot.py
EOF
