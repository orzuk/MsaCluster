#!/bin/bash
# ============================================================================
# compress_v1_backup.sh -- per-pair tar.gz the v1 results backup.
#
# The backup at backups/results_v1_<date>/ has this layout:
#     <pair>__output_AF/
#     <pair>__output_DDG/
#     <pair>__output_s4pred/
#     ... 11 method subdirs per pair, 93 pairs total ...
#
# This script groups all <pair>__* directories into a single archive:
#     <pair>.tar.gz
# stored inside the SAME backup root. Per-pair granularity = easy
# to list, extract, or copy individual pairs; not one giant archive.
#
# Steps per pair:
#   1. tar -czf <pair>.tar.gz <pair>__*
#   2. tar -tzf <pair>.tar.gz > /dev/null   (verify integrity)
#   3. (optional, with --delete) rm -rf <pair>__*
#
# Expected compression: 3-5x (PDB/JSON text compresses well, NPY/NPZ
# mildly). 61 GB -> ~15-20 GB.
#
# Usage:
#   bash scripts/compress_v1_backup.sh                       # compress, keep originals
#   bash scripts/compress_v1_backup.sh --delete              # compress then delete originals
#   bash scripts/compress_v1_backup.sh --backup-dir <path>   # custom backup root
#   bash scripts/compress_v1_backup.sh --dry-run             # show what would happen
# ============================================================================
set -euo pipefail

BACKUP_DIR=""
DELETE_ORIG=0
DRY_RUN=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --delete)     DELETE_ORIG=1; shift ;;
        --backup-dir) BACKUP_DIR="$2"; shift 2 ;;
        --dry-run)    DRY_RUN=1; shift ;;
        -h|--help)    sed -n '2,28p' "$0"; exit 0 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

# Auto-find most recent backup dir if not specified
if [[ -z "$BACKUP_DIR" ]]; then
    BACKUP_DIR=$(ls -dt backups/results_v1_* 2>/dev/null | head -1)
    if [[ -z "$BACKUP_DIR" ]]; then
        echo "No backups/results_v1_* directory found"
        exit 1
    fi
fi
[[ ! -d "$BACKUP_DIR" ]] && { echo "Backup dir not found: $BACKUP_DIR"; exit 1; }

cd "$BACKUP_DIR"
echo "[compress] BACKUP_DIR: $(pwd)"
echo "[compress] DELETE original subdirs after compress: $DELETE_ORIG"
echo "[compress] Dry-run: $DRY_RUN"
echo

# Gather unique pair names from <pair>__<subdir> directories
PAIRS=$(ls -1 -d */ 2>/dev/null \
        | sed 's#/$##' \
        | grep '__' \
        | sed 's/__[^_]*\(_[^_]*\)*$//' \
        | sort -u)

n=$(echo "$PAIRS" | wc -w)
echo "[compress] Found $n unique pairs"
echo

SIZE_BEFORE=$(du -sb . 2>/dev/null | awk '{print $1}')
SIZE_HUMAN_BEFORE=$(du -sh . 2>/dev/null | awk '{print $1}')

ok=0; skip=0; fail=0
for P in $PAIRS; do
    # Skip if archive already exists and is non-empty
    if [[ -f "${P}.tar.gz" ]] && [[ -s "${P}.tar.gz" ]]; then
        skip=$((skip+1))
        continue
    fi

    # Glob the subdirs for this pair
    shopt -s nullglob
    SUBS=( "${P}"__* )
    shopt -u nullglob
    if [[ ${#SUBS[@]} -eq 0 ]]; then
        echo "[compress] $P: no subdirs found; skip"
        skip=$((skip+1))
        continue
    fi

    if [[ "$DRY_RUN" -eq 1 ]]; then
        echo "[DRY-RUN] tar -czf ${P}.tar.gz ${SUBS[*]}"
        if [[ "$DELETE_ORIG" -eq 1 ]]; then
            echo "[DRY-RUN]   then rm -rf ${SUBS[*]}"
        fi
        ok=$((ok+1))
        continue
    fi

    # Create archive
    if tar -czf "${P}.tar.gz" "${SUBS[@]}" 2>"${P}.tar.gz.err"; then
        # Verify integrity (cheap: list, redirect to /dev/null)
        if tar -tzf "${P}.tar.gz" > /dev/null 2>>"${P}.tar.gz.err"; then
            rm -f "${P}.tar.gz.err"
            sz=$(du -sh "${P}.tar.gz" | awk '{print $1}')
            echo "[compress] $P  ok  (${#SUBS[@]} subdirs -> ${sz})"
            if [[ "$DELETE_ORIG" -eq 1 ]]; then
                rm -rf "${SUBS[@]}"
            fi
            ok=$((ok+1))
        else
            echo "[compress] $P  FAIL verify: see ${P}.tar.gz.err"
            fail=$((fail+1))
        fi
    else
        echo "[compress] $P  FAIL tar: see ${P}.tar.gz.err"
        fail=$((fail+1))
    fi
done

echo
echo "[compress] Summary: ok=$ok skip=$skip fail=$fail"

if [[ "$DRY_RUN" -eq 0 ]]; then
    SIZE_AFTER=$(du -sb . 2>/dev/null | awk '{print $1}')
    SIZE_HUMAN_AFTER=$(du -sh . 2>/dev/null | awk '{print $1}')
    echo "[compress] Size before: ${SIZE_HUMAN_BEFORE}"
    echo "[compress] Size after:  ${SIZE_HUMAN_AFTER}"
    if [[ "$SIZE_BEFORE" -gt 0 ]]; then
        ratio=$(awk "BEGIN{printf \"%.1f\", $SIZE_BEFORE/$SIZE_AFTER}")
        echo "[compress] Compression ratio: ${ratio}x"
    fi
fi

cat <<EOF

To extract a single pair's data back:
  tar -xzf <BACKUP_DIR>/<pair>.tar.gz -C /tmp/restore/
  ls /tmp/restore/<pair>__*/
EOF
