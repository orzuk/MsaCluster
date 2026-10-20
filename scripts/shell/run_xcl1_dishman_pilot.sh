#!/bin/bash
# ===========================================================================
# XCL1 / Dishman curated-MSA pilot (pair 2n54B_2hdmA: F1 = 2n54B alternative
# fold, F2 = 2hdmA chemokine fold).
#
# Re-runs the per-cluster fold-preference analysis on a CURATED vertebrate
# chemokine set (data/pilots/xcl1_dishman/, built by
# scripts/build_curated_chemokine_msa.py; cf. Dishman et al., Science 2021)
# instead of the ~8000-sequence MMseqs2/UniRef MSA, with TRUE-CLADE
# (monophyletic) clustering, to test whether the clade-shift scan localizes
# the gain of metamorphism on the XCL lineage.
#
# Production outputs are untouched: every step runs with
# --data_root <pilot_root> (env MSACLUSTER_DATA_DIR), so all pipeline paths
# resolve to <pilot_root>/<pair>/ instead of Pipeline/FoldPairs/<pair>/.
#
# Run from a COMPUTE node (never the login node, never ~):
#   srun -A orzuk -p glacier --mem=8G -t 01:00:00 \
#        bash scripts/shell/run_xcl1_dishman_pilot.sh
#
# Options:
#   --dry-run              print every command (setup + sbatch), change nothing
#   --setup-only           prepare the pilot dir + tree + clusters, submit nothing
#   --skip-setup           assume the pilot dir / clusters exist, only submit
#   --methods "a b c"      subset of: af esmfold2 boltz2 bioemu ddg s4pred
#                          (default: all; 'af' = AF2 + AF3 in one job)
#   --clusters K           target number of monophyletic clusters (default 30)
#   --min-cluster-size N   min sequences per cluster (default 5)
#   --pilot-root DIR       default Pipeline/Pilots/xcl1_dishman under the repo
#   --no-analysis          skip the final dependent analysis job
# ===========================================================================
set -euo pipefail

REPO="${XCL1_PILOT_REPO:-/sci/labs/orzuk/orzuk/github/MsaCluster}"
if [[ ! -d "$REPO" ]]; then
    # Not on the cluster (e.g. a local --dry-run): use this checkout.
    REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
    echo "[pilot] NOTE: cluster repo path not found; using $REPO"
fi
VENV=/sci/labs/orzuk/orzuk/venvs/my-python-venv
ACCOUNT=orzuk
GPU_PARTITION=salmon
GPU_GRES=gpu:l40s:1
CPU_PARTITION=glacier
PAIR=2n54B_2hdmA
PILOT_ROOT="$REPO/Pipeline/Pilots/xcl1_dishman"
CURATED_DIR="$REPO/data/pilots/xcl1_dishman"
METHODS="af esmfold2 boltz2 bioemu ddg s4pred"
K_CLUSTERS=30
MIN_CLUSTER_SIZE=5
DRY_RUN=0
SETUP_ONLY=0
SKIP_SETUP=0
RUN_ANALYSIS=1

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)           DRY_RUN=1; shift ;;
        --setup-only)        SETUP_ONLY=1; shift ;;
        --skip-setup)        SKIP_SETUP=1; shift ;;
        --no-analysis)       RUN_ANALYSIS=0; shift ;;
        --methods)           METHODS="$2"; shift 2 ;;
        --clusters)          K_CLUSTERS="$2"; shift 2 ;;
        --min-cluster-size)  MIN_CLUSTER_SIZE="$2"; shift 2 ;;
        --pilot-root)        PILOT_ROOT="$2"; shift 2 ;;
        -h|--help)           sed -n '2,36p' "$0"; exit 0 ;;
        *) echo "Unknown arg: $1" >&2; exit 2 ;;
    esac
done

PROD_PAIR_DIR="$REPO/Pipeline/FoldPairs/$PAIR"
PILOT_PAIR_DIR="$PILOT_ROOT/$PAIR"
JOBS_DIR="$PILOT_ROOT/jobs"

run() {   # run (or, with --dry-run, just print) a command
    if [[ "$DRY_RUN" -eq 1 ]]; then
        echo "[dry-run] $*"
    else
        echo "[run]     $*"
        "$@"
    fi
}

echo "=== XCL1 Dishman pilot ==="
echo "  pair          : $PAIR"
echo "  pilot root    : $PILOT_ROOT   (production $PROD_PAIR_DIR is read-only here)"
echo "  curated MSA   : $CURATED_DIR"
echo "  clustering    : tree_mono, K=$K_CLUSTERS, min_size=$MIN_CLUSTER_SIZE"
echo "  methods       : $METHODS"
echo "  dry-run=$DRY_RUN setup_only=$SETUP_ONLY skip_setup=$SKIP_SETUP analysis=$RUN_ANALYSIS"
echo

cd "$REPO"
if [[ "$DRY_RUN" -eq 0 ]]; then
    # shellcheck disable=SC1091
    source "$VENV/bin/activate"
fi
export MSACLUSTER_DATA_DIR="$PILOT_ROOT"
export PYTHONUNBUFFERED=1
run mkdir -p "$JOBS_DIR" "$PILOT_PAIR_DIR/output_get_msa" "$PILOT_ROOT/results"

# ---------------------------------------------------------------------------
# Setup: mirror the pair layout (truth PDBs + chain files), install the
# curated DeepMsa.a3m, build the tree, cut it into monophyletic clades.
# ---------------------------------------------------------------------------
if [[ "$SKIP_SETUP" -eq 0 ]]; then
    if [[ ! -d "$PROD_PAIR_DIR" ]]; then
        if [[ "$DRY_RUN" -eq 1 ]]; then
            echo "[dry-run] NOTE: production pair dir $PROD_PAIR_DIR not present here"
        else
            echo "Missing production pair dir $PROD_PAIR_DIR" >&2
            exit 1
        fi
    fi

    echo "--- setup 1: copy truth structures + chain FASTAs from the production pair ---"
    for sub in chain_pdb_files fasta_chain_files; do
        [[ -d "$PROD_PAIR_DIR/$sub" ]] && run cp -r "$PROD_PAIR_DIR/$sub" "$PILOT_PAIR_DIR/"
    done
    if [[ "$DRY_RUN" -eq 1 ]]; then
        echo "[dry-run] cp $PROD_PAIR_DIR/*.pdb $PILOT_PAIR_DIR/"
    else
        shopt -s nullglob
        for f in "$PROD_PAIR_DIR"/*.pdb "$PROD_PAIR_DIR"/*.fasta "$PROD_PAIR_DIR"/_seed_both.a3m; do
            cp "$f" "$PILOT_PAIR_DIR/"
        done
        shopt -u nullglob
    fi

    echo "--- setup 2: build DeepMsa.a3m from the curated chemokine set ---"
    # Align against the PRODUCTION chain FASTAs (what every scorer uses as the
    # truth sequences), not the RCSB entity sequences.
    run python3 scripts/build_curated_chemokine_msa.py align \
        --fasta "$CURATED_DIR/chemokines_curated.fasta" \
        --query_fasta "$PILOT_PAIR_DIR/fasta_chain_files/2n54B.fasta" \
        --partner_fasta "$PILOT_PAIR_DIR/fasta_chain_files/2hdmA.fasta" \
        --query_tag 2n54B --partner_tag 2hdmA \
        --out "$PILOT_PAIR_DIR/output_get_msa/DeepMsa.a3m"

    echo "--- setup 3: phylogenetic tree over the curated MSA (all sequences) ---"
    run python3 run_foldswitch_pipeline.py --run_mode tree --foldpair_ids "$PAIR" \
        --data_root "$PILOT_ROOT" --run_job_mode inline --tree_max_seqs 0

    echo "--- setup 4: monophyletic (true-clade) clustering, K=$K_CLUSTERS ---"
    run python3 run_foldswitch_pipeline.py --run_mode cluster_msa --foldpair_ids "$PAIR" \
        --data_root "$PILOT_ROOT" --run_job_mode inline \
        --cluster_alg tree_mono \
        --cluster_tree_min_num_clusters "$K_CLUSTERS" \
        --cluster_max_clusters "$K_CLUSTERS" \
        --cluster_min_output_size "$MIN_CLUSTER_SIZE" \
        --cluster_min_neff 1 \
        --cluster_frac_gaps_cutoff 0.9

    if [[ "$DRY_RUN" -eq 0 ]]; then
        N_CL=$(ls "$PILOT_PAIR_DIR/output_msa_cluster"/ShallowMsa_*.a3m 2>/dev/null | wc -l)
        echo "[setup] clusters produced: $N_CL"
        [[ "$N_CL" -ge 2 ]] || { echo "[setup] ERROR: fewer than 2 clusters; stopping." >&2; exit 1; }
    fi
fi

if [[ "$SETUP_ONLY" -eq 1 ]]; then
    echo "[pilot] --setup-only: nothing submitted."
    exit 0
fi

# ---------------------------------------------------------------------------
# Submit ONE job per method (per-cluster predictions run inline inside it).
# ---------------------------------------------------------------------------
submit() {  # submit <jobname> <sbatch opts...> -- <script> <args...>
    local name="$1"; shift
    local cmd=(sbatch --parsable -J "$name" -A "$ACCOUNT" \
               --output="$JOBS_DIR/${name}_%j.out" "$@")
    if [[ "$DRY_RUN" -eq 1 ]]; then
        echo "[dry-run] ${cmd[*]}" >&2
        echo "DRYRUN"
    else
        echo "[sbatch]  ${cmd[*]}" >&2
        "${cmd[@]}"
    fi
}

JOBIDS=()
for M in $METHODS; do
    case "$M" in
        af)       OPTS=(-p "$GPU_PARTITION" --gres="$GPU_GRES" -c 8 --mem=48G -t 16:00:00) ;;
        esmfold2) OPTS=(-p "$GPU_PARTITION" --gres="$GPU_GRES" -c 4 --mem=24G -t 04:00:00) ;;
        boltz2)   OPTS=(-p "$GPU_PARTITION" --gres="$GPU_GRES" -c 8 --mem=32G -t 08:00:00) ;;
        bioemu)   OPTS=(-p "$GPU_PARTITION" --gres="$GPU_GRES" -c 8 --mem=32G -t 10:00:00) ;;
        ddg)      OPTS=(-p "$GPU_PARTITION" --gres="$GPU_GRES" -c 4 --mem=16G -t 02:00:00) ;;
        s4pred)   OPTS=(-p "$CPU_PARTITION" -c 4 --mem=8G -t 02:00:00) ;;
        *) echo "Unknown method: $M" >&2; exit 2 ;;
    esac
    JID=$(submit "xcl1_$M" "${OPTS[@]}" \
          "$REPO/scripts/xcl1_dishman_method.sbatch" "$M" "$PILOT_ROOT" "$PAIR")
    JOBIDS+=("$JID")
    echo "  submitted $M -> $JID"
done

if [[ "$RUN_ANALYSIS" -eq 1 ]]; then
    DEP="afterok:$(IFS=:; echo "${JOBIDS[*]}")"
    AID=$(submit xcl1_analysis --dependency="$DEP" \
          "$REPO/scripts/xcl1_dishman_analysis.sbatch" "$PILOT_ROOT" "$PAIR")
    echo "  submitted analysis -> $AID ($DEP)"
fi

cat <<EOF

=== submitted ${#JOBIDS[@]} method job(s) ===
  logs    : $JOBS_DIR/
  results : $PILOT_ROOT/results/ (survey, phylo_placement, cluster_family_composition)
  monitor : squeue -u \$USER
EOF
