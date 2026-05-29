#!/bin/bash
# ============================================================================
# clean_v1_for_v2.sh -- delete v1 standard-output method dirs in
# Pipeline/FoldPairs/*/ so v2 can write into clean dirs.
#
# Safety:
#   * Refuses to run unless a backups/results_v1_*/ dir exists
#     (the rsync'd 61 GB copy is the only safety net).
#   * --dry-run prints what would be deleted, no actual rm.
#   * Per-pair, per-method confirmation prints (set -v style).
#
# Method dirs that are deleted (the v2 8-method set + the orphan
# v1 ESMFold dir, which v2 replaces with output_esmfold2/):
#   output_AF, output_esmfold2, output_esm_fold, output_ESMFold,
#   output_esmfold, output_boltz2, output_Boltz2, output_ddg, output_DDG,
#   output_s4pred, output_cmaps
#
# NOT deleted (these are *inputs* to v2):
#   output_msa_cluster, output_msa_cluster_coarse, output_phytree,
#   output_get_msa, anything else.
#
# Usage:
#   bash scripts/clean_v1_for_v2.sh --dry-run
#   bash scripts/clean_v1_for_v2.sh
#   bash scripts/clean_v1_for_v2.sh --pairs '5jytA_2qkeE 2n54B_2hdmA'
# ============================================================================
set -euo pipefail

DRY_RUN=0
PAIRS=""
ANALYSIS_ONLY=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)       DRY_RUN=1; shift ;;
        --pairs)         PAIRS="$2"; shift 2 ;;
        --analysis-only) ANALYSIS_ONLY=1; shift ;;
        -h|--help) sed -n '2,28p' "$0"; exit 0 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

# Safety: require a backup to exist before deleting.
BACKUP=$(ls -dt backups/results_v1_* 2>/dev/null | head -1)
if [[ -z "$BACKUP" ]]; then
    echo "[clean] ERROR: no backups/results_v1_* dir found."
    echo "[clean] Refusing to delete without a backup."
    exit 1
fi
BACKUP_SIZE=$(du -sh "$BACKUP" 2>/dev/null | awk '{print $1}')
echo "[clean] Backup present: $BACKUP ($BACKUP_SIZE)"
echo "[clean] Dry-run: $DRY_RUN"
echo

METHOD_DIRS=(
    output_AF
    output_esmfold2
    output_esm_fold
    output_ESMFold
    output_esmfold
    output_boltz2
    output_Boltz2
    output_ddg
    output_DDG
    output_s4pred
    output_cmaps
)
# Per-method aggregator CSVs / artifacts in Analysis/. Several method
# wrappers self-skip if these exist, so v2 must clear them too.
ANALYSIS_FILES=(
    Analysis/df_ddg.csv
    Analysis/df_s4pred.csv
    Analysis/df_boltz2.csv
    Analysis/df_esmfold2.csv
    Analysis/df_af.csv
)
# Per-cluster S4PRED .ss2 files (one per cluster, sit in Analysis/).
ANALYSIS_GLOBS=(
    "Analysis/s4pred_*.ss2"
)

if [[ -n "$PAIRS" ]]; then
    PAIR_LIST="$PAIRS"
else
    PAIR_LIST=$(ls -d Pipeline/FoldPairs/*/ 2>/dev/null \
                | sed -E 's#.*FoldPairs/##; s#/##' \
                | grep -vE '^(_s4pred_work|jobs)$')
fi

n_pairs=0
n_dirs_deleted=0
total_size_freed=0
for P in $PAIR_LIST; do
    PAIR_DIR="Pipeline/FoldPairs/$P"
    [[ ! -d "$PAIR_DIR" ]] && continue
    n_pairs=$((n_pairs+1))
    for M in "${METHOD_DIRS[@]}"; do
        D="$PAIR_DIR/$M"
        # --analysis-only mode: skip method-dir deletion entirely (only clear
        # Analysis/df_*.csv below). Use this when v2 outputs are already
        # present in output_* dirs and only the aggregator CSVs need clearing.
        [[ "$ANALYSIS_ONLY" -eq 1 ]] && continue
        if [[ -d "$D" ]]; then
            SZ=$(du -sb "$D" 2>/dev/null | awk '{print $1}')
            [[ -z "$SZ" ]] && SZ=0
            SZ_HUMAN=$(du -sh "$D" 2>/dev/null | awk '{print $1}')
            if [[ "$DRY_RUN" -eq 1 ]]; then
                echo "[DRY-RUN] rm -rf $D  (${SZ_HUMAN})"
            else
                rm -rf "$D"
                echo "[clean] deleted $D  (${SZ_HUMAN})"
            fi
            n_dirs_deleted=$((n_dirs_deleted+1))
            # Use awk for the running total -- bash arithmetic comparisons can
            # misbehave on >2 GB values on some systems; awk uses double-prec
            # floats so the sum is always faithful up to 2^53 bytes (~8 PB).
            total_size_freed=$(awk "BEGIN{print $total_size_freed + $SZ}")
        fi
    done
    # Per-method aggregator CSVs (the pipeline skips methods if these exist).
    for F in "${ANALYSIS_FILES[@]}"; do
        FP="$PAIR_DIR/$F"
        if [[ -f "$FP" ]]; then
            SZ=$(stat -c %s "$FP" 2>/dev/null || echo 0)
            if [[ "$DRY_RUN" -eq 1 ]]; then
                echo "[DRY-RUN] rm -f $FP  (${SZ} bytes)"
            else
                rm -f "$FP"
                echo "[clean] deleted $FP  (${SZ} bytes)"
            fi
            n_dirs_deleted=$((n_dirs_deleted+1))
            total_size_freed=$(awk "BEGIN{print $total_size_freed + $SZ}")
        fi
    done
    for G in "${ANALYSIS_GLOBS[@]}"; do
        shopt -s nullglob
        for FP in $PAIR_DIR/$G; do
            SZ=$(stat -c %s "$FP" 2>/dev/null || echo 0)
            if [[ "$DRY_RUN" -eq 1 ]]; then
                echo "[DRY-RUN] rm -f $FP"
            else
                rm -f "$FP"
            fi
            n_dirs_deleted=$((n_dirs_deleted+1))
            total_size_freed=$(awk "BEGIN{print $total_size_freed + $SZ}")
        done
        shopt -u nullglob
    done
done

# Convert bytes -> human-readable. numfmt is part of coreutils (always present
# on cluster Linux); awk fallback handles weirder environments.
SIZE_HUMAN=$(numfmt --to=iec --suffix=B --format='%.1f' "$total_size_freed" 2>/dev/null \
             || awk -v b="$total_size_freed" 'BEGIN{
                 if (b>=1073741824) printf "%.1f GB", b/1073741824;
                 else if (b>=1048576) printf "%.1f MB", b/1048576;
                 else if (b>=1024)    printf "%.1f KB", b/1024;
                 else                 printf "%d B", b;
             }')

echo
echo "[clean] Pairs visited:     $n_pairs"
echo "[clean] Method dirs found: $n_dirs_deleted"
echo "[clean] Size freed:        $SIZE_HUMAN"
[[ "$DRY_RUN" -eq 1 ]] && echo "[clean] (dry-run; nothing actually deleted)"
