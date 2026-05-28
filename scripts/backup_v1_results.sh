#!/bin/bash
# ============================================================================
# backup_v1_results.sh -- back up v1 method results before v2 method runs.
#
# Uses rsync with --ignore-existing so re-running this script will pick up
# only what's missing (interrupted backups can be resumed safely).
#
# Backs up these per-pair subdirectories:
#   output_AF        output_DDG        output_s4pred
#   output_ccmpred   output_msa_transformer
#   output_ESMFold   output_esmfold    output_esm
#   output_Boltz2    output_AlphaFold  output_AlphaFold3
#
# Backup location: backups/results_v1_<YYYYMMDD>/<pair>__<subdir>/
#
# Usage:
#   bash scripts/backup_v1_results.sh                  # full or resume
#   bash scripts/backup_v1_results.sh --tag custom     # custom suffix
#   bash scripts/backup_v1_results.sh --dry-run
# ============================================================================
set -euo pipefail

TAG="$(date +%Y%m%d)"
DRY_RUN=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --tag)     TAG="$2"; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
        -h|--help) sed -n '2,22p' "$0"; exit 0 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

BACKUP_DIR="backups/results_v1_${TAG}"
mkdir -p "$BACKUP_DIR"

SUBDIRS=(
    output_AF
    output_DDG
    output_s4pred
    output_ccmpred
    output_msa_transformer
    output_ESMFold
    output_esmfold
    output_esm
    output_Boltz2
    output_AlphaFold
    output_AlphaFold3
)

# rsync flags:
#   -a   archive (preserve perms, timestamps, recursion)
#   --ignore-existing  skip files already at the destination (resume-safe)
#   --info=progress2   single-line live progress; comment out for quiet
RSYNC_FLAGS="-a --ignore-existing"
if [[ "$DRY_RUN" -eq 1 ]]; then
    RSYNC_FLAGS="$RSYNC_FLAGS --dry-run -v"
fi

echo "[backup] Target: $BACKUP_DIR"
echo "[backup] Dry-run: $DRY_RUN"
echo

n_pair=0
n_copied=0
for D in Pipeline/FoldPairs/*/; do
    P=$(basename "$D")
    [[ "$P" == "_s4pred_work" || "$P" == "jobs" ]] && continue
    n_pair=$((n_pair+1))
    for SUB in "${SUBDIRS[@]}"; do
        SRC="${D}${SUB}"
        if [[ -d "$SRC" ]]; then
            DST="${BACKUP_DIR}/${P}__${SUB}"
            # rsync needs trailing slashes on src to copy contents
            mkdir -p "$DST"
            rsync $RSYNC_FLAGS "${SRC}/" "${DST}/"
            n_copied=$((n_copied+1))
        fi
    done
done

echo
echo "[backup] Visited $n_pair pairs, $n_copied subdir copies"
if [[ "$DRY_RUN" -eq 0 ]]; then
    du -sh "$BACKUP_DIR"
    echo "[backup] Done. To resume, just rerun this script (idempotent)."
fi
