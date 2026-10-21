#!/bin/bash
# ===========================================================================
# One-off: MOVE (never delete) legacy MSA-Transformer / CCMpred contact outputs
# and their cached scores out of the way, so fine-cluster contacts can be
# recomputed with provenance manifests (utils/contact_utils.py).
#
# Legacy outputs were computed on an OLDER fine clustering: their ShallowMsa_k
# labels refer to different sequence sets than the current
# output_msa_cluster/ShallowMsa_k.a3m. The scorer already ignores them (no
# manifest), this just makes the pair directories unambiguous.
#
# Per pair, moves into <pair>/_legacy_contacts_<YYYYmmdd>/ (same relative paths):
#   output_cmaps/msa_transformer/   (only if it has NO _manifest.json)
#   output_cmaps/ccmpred/           (only if it has NO _manifest.json)
#   Analysis/df_cmap.csv, Analysis/df_cmap_ccmpred.csv
#   tmp_ccmpred/
#   output_msa_cluster_pairtrim/
#
# Usage (from the repo root; dry run by default, prints what it would move):
#   bash scripts/shell/archive_legacy_contacts.sh <pair_id> [pair_id ...]
#   bash scripts/shell/archive_legacy_contacts.sh ALL
#   bash scripts/shell/archive_legacy_contacts.sh --apply ALL
# Env: PIPELINE_ROOT (default Pipeline/FoldPairs), ARCHIVE_DATE (default today)
# ===========================================================================
set -euo pipefail

APPLY=0
PAIRS=()
for a in "$@"; do
    case "$a" in
        --apply) APPLY=1 ;;
        --dry-run) APPLY=0 ;;
        -h|--help) sed -n '2,25p' "$0"; exit 0 ;;
        *) PAIRS+=("$a") ;;
    esac
done
[[ ${#PAIRS[@]} -gt 0 ]] || { echo "usage: $0 [--apply] <pair_id ...|ALL>"; exit 2; }

ROOT="${PIPELINE_ROOT:-Pipeline/FoldPairs}"
DATE="${ARCHIVE_DATE:-$(date +%Y%m%d)}"
[[ -d "$ROOT" ]] || { echo "no such dir: $ROOT (run from the repo root or set PIPELINE_ROOT)"; exit 2; }

if [[ "${PAIRS[0]}" == "ALL" ]]; then
    mapfile -t PAIRS < <(ls -d "$ROOT"/*/ | xargs -n1 basename \
                         | grep -E '^[0-9a-zA-Z]{5}_[0-9a-zA-Z]{5}$')
fi

MODE=$([[ $APPLY -eq 1 ]] && echo "APPLY" || echo "DRY-RUN")
echo "[archive] mode=$MODE root=$ROOT pairs=${#PAIRS[@]} dest=<pair>/_legacy_contacts_$DATE"

n_items=0; n_skip=0
for pair in "${PAIRS[@]}"; do
    pdir="$ROOT/$pair"
    [[ -d "$pdir" ]] || { echo "[archive] WARN missing pair dir $pdir"; continue; }
    dest="$pdir/_legacy_contacts_$DATE"
    for rel in output_cmaps/msa_transformer output_cmaps/ccmpred \
               Analysis/df_cmap.csv Analysis/df_cmap_ccmpred.csv \
               tmp_ccmpred output_msa_cluster_pairtrim; do
        src="$pdir/$rel"
        [[ -e "$src" ]] || continue
        if [[ -d "$src" && -f "$src/_manifest.json" ]]; then
            echo "[archive] SKIP $src (has _manifest.json: already provenance-tracked)"
            n_skip=$((n_skip + 1)); continue
        fi
        if [[ -e "$dest/$rel" ]]; then
            echo "[archive] SKIP $src (destination $dest/$rel exists)"
            n_skip=$((n_skip + 1)); continue
        fi
        if [[ -d "$src" ]]; then
            what="dir, $(find "$src" -type f | wc -l) files, $(du -sh "$src" | cut -f1)"
        else
            what="file, $(du -h "$src" | cut -f1)"
        fi
        echo "[archive] $MODE mv $src -> $dest/$rel ($what)"
        n_items=$((n_items + 1))
        if [[ $APPLY -eq 1 ]]; then
            mkdir -p "$(dirname "$dest/$rel")"
            mv -n "$src" "$dest/$rel"
        fi
    done
done
echo "[archive] $MODE done: $n_items item(s) $([[ $APPLY -eq 1 ]] && echo moved || echo would be moved), $n_skip skipped"
[[ $APPLY -eq 1 ]] || echo "[archive] dry run only; re-run with --apply to move"
