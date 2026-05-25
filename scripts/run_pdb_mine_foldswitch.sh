#!/bin/bash
# ============================================================================
# run_pdb_mine_foldswitch.sh -- sbatch wrapper for the PDB mining script that
#                               searches for new fold-switching pair
#                               candidates beyond the existing 93.
#
# Implements Porter & Looger 2018's methodology on the current PDB:
#   1. Download RCSB sequence-cluster file at --min-seq-identity (default 95%)
#   2. Filter clusters by metadata (X-ray + cryo-EM, resolution <= 3.0 A,
#      length 50-500 AA)
#   3. Pairwise TM-align within each cluster (capped at 20 representatives)
#   4. Keep pairs with TM-score <= 0.5 AND DSSP secondary-structure H<->E
#      switch fraction >= 0.10 (excludes domain swaps and hinge motions)
#   5. Deduplicate against data/foldswitch_PDB_IDs_full.txt
#
# Outputs:
#   data/pdb_mining_candidates.csv           (all TM<=0.5+SS-pass pairs)
#   data/pdb_mining_candidates.new_only.csv  (subset not in existing 93)
#
# Cache:
#   data/cache/pdb_mining/bc-95.out          (RCSB cluster file)
#   data/cache/pdb_mining/metadata/*.json    (per-PDB metadata, cached)
#   data/cache/pdb_mining/pdb_files/*.pdb    (per-PDB structures, cached)
#
# Usage (from repo root):
#   # 100-cluster test run (5-15 min, validates pipeline + caches PDBs):
#   sbatch scripts/run_pdb_mine_foldswitch.sh test
#
#   # Full run (all clusters at 95% identity; 1-3 hours on 16 workers):
#   sbatch scripts/run_pdb_mine_foldswitch.sh full
#
#   # Strict 99% identity, no SS filter (for comparison):
#   sbatch scripts/run_pdb_mine_foldswitch.sh strict
#
#   # Resume an interrupted full run:
#   sbatch scripts/run_pdb_mine_foldswitch.sh resume
#
# Custom invocation -- pass arbitrary args to the Python script:
#   sbatch scripts/run_pdb_mine_foldswitch.sh custom --min-seq-identity 0.90 \
#          --tm-threshold 0.45 --workers 32
#
# Monitor:
#   tail -f Pipeline/FoldPairs/jobs/pdb_mine.out
#   squeue -u $USER -j <jobid>
#
# After completion, check the candidate list:
#   wc -l data/pdb_mining_candidates.new_only.csv
#   head data/pdb_mining_candidates.new_only.csv
# ============================================================================
#SBATCH -J pdb_mine_foldswitch
#SBATCH --account=course-52017-25
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 06:00:00
#SBATCH --output=Pipeline/FoldPairs/jobs/pdb_mine.out
#SBATCH --error=Pipeline/FoldPairs/jobs/pdb_mine.err

set -euo pipefail

# ---- paths ------------------------------------------------------------------
REPO_DIR="/sci/labs/orzuk/orzuk/github/MsaCluster"
VENV_DIR="/sci/labs/orzuk/orzuk/venvs/my-python-venv"

cd "$REPO_DIR"
# shellcheck disable=SC1091
source "$VENV_DIR/bin/activate"

mkdir -p Pipeline/FoldPairs/jobs data/cache/pdb_mining

# ---- pick mode --------------------------------------------------------------
MODE="${1:-full}"

echo "[run_pdb_mine] mode=$MODE"
echo "[run_pdb_mine] node=$(hostname)  cpus=${SLURM_CPUS_PER_TASK:-16}"
echo "[run_pdb_mine] start: $(date -Iseconds)"
echo

case "$MODE" in
    test)
        # 100 largest clusters -- 5-15 min, validates pipeline end-to-end
        python3 scripts/pdb_mine_foldswitch_candidates.py \
            --max-clusters 100 \
            --workers "${SLURM_CPUS_PER_TASK:-16}" \
            --output data/pdb_mining_candidates_test.csv
        ;;
    full)
        # Full PDB scan at default 95% identity, default thresholds
        python3 scripts/pdb_mine_foldswitch_candidates.py \
            --workers "${SLURM_CPUS_PER_TASK:-16}"
        ;;
    strict)
        # Strict: 99% identity, no SS filter (compare TM-only vs Porter SS)
        python3 scripts/pdb_mine_foldswitch_candidates.py \
            --min-seq-identity 0.99 \
            --min-ss-switch-fraction 0.0 \
            --workers "${SLURM_CPUS_PER_TASK:-16}" \
            --output data/pdb_mining_candidates_strict99.csv
        ;;
    resume)
        # Pick up after an interrupted full run
        python3 scripts/pdb_mine_foldswitch_candidates.py \
            --workers "${SLURM_CPUS_PER_TASK:-16}" \
            --resume
        ;;
    custom)
        shift  # drop the "custom" arg, pass the rest through
        python3 scripts/pdb_mine_foldswitch_candidates.py \
            --workers "${SLURM_CPUS_PER_TASK:-16}" \
            "$@"
        ;;
    *)
        echo "Unknown mode: $MODE"
        echo "Usage: sbatch $0 {test|full|strict|resume|custom [args...]}"
        exit 1
        ;;
esac

echo
echo "[run_pdb_mine] end: $(date -Iseconds)"

cat <<'EOF'

Inspect results:
  wc -l data/pdb_mining_candidates.csv
  wc -l data/pdb_mining_candidates.new_only.csv
  head -20 data/pdb_mining_candidates.new_only.csv

Next steps after first run:
  1. Review data/pdb_mining_candidates.new_only.csv for biologically
     plausible candidates (sort by tm_max ascending + ss_switch_fraction
     descending to put the strongest fold-switchers at the top).
  2. To add accepted candidates to the dataset, append their pair IDs to
     data/foldswitch_PDB_IDs_full.txt and re-run the pipeline on the
     extended set.
EOF
