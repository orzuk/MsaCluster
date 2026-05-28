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
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY_RUN=1; shift ;;
        --pairs)   PAIRS="$2"; shift 2 ;;
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
        if [[ -d "$D" ]]; then
            SZ=$(du -sb "$D" 2>/dev/null | awk '{print $1}')
            SZ_HUMAN=$(du -sh "$D" 2>/dev/null | awk '{print $1}')
            if [[ "$DRY_RUN" -eq 1 ]]; then
                echo "[DRY-RUN] rm -rf $D  (${SZ_HUMAN})"
            else
                rm -rf "$D"
                echo "[clean] deleted $D  (${SZ_HUMAN})"
            fi
            n_dirs_deleted=$((n_dirs_deleted+1))
            total_size_freed=$((total_size_freed+SZ))
        fi
    done
done

# Convert bytes -> human
human() {
    local b=$1
    if   [[ $b -gt 1073741824 ]]; then awk "BEGIN{printf \"%.1f GB\", $b/1073741824}"
    elif [[ $b -gt 1048576    ]]; then awk "BEGIN{printf \"%.1f MB\", $b/1048576}"
    elif [[ $b -gt 1024       ]]; then awk "BEGIN{printf \"%.1f KB\", $b/1024}"
    else echo "${b} B"; fi
}
SIZE_HUMAN=$(human "$total_size_freed")

echo
echo "[clean] Pairs visited:     $n_pairs"
echo "[clean] Method dirs found: $n_dirs_deleted"
echo "[clean] Size freed:        $SIZE_HUMAN"
[[ "$DRY_RUN" -eq 1 ]] && echo "[clean] (dry-run; nothing actually deleted)"
