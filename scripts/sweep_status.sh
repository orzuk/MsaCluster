#!/bin/bash
# ============================================================================
# sweep_status.sh -- per-pair completion status for the v2 8-method sweep.
#
# Walks Pipeline/FoldPairs/<pair>/ and reports, for each pair, which of the
# 8 v2 methods have produced their canonical output. Also lists any
# currently running msaclust jobs from squeue.
#
# A method is considered DONE if its output file/dir exists and is non-empty.
# Methods + their completion markers:
#   AF2     -> output_AF/AF2/ShallowMsa_*/         (non-empty subdirs)
#   AF3     -> output_AF/AF3/ShallowMsa_*/         (non-empty subdirs)
#   ESM2    -> output_esmfold2/*.{pdb,cif}         (any predicted structure)
#   Boltz   -> output_boltz2/ShallowMsa_*/         (any cluster dir)
#   DDG     -> Analysis/df_ddg.csv                 (aggregator CSV)
#   S4P     -> Analysis/df_s4pred.csv              (aggregator CSV)
#   CCM     -> output_cmaps/ccmpred/Coarse*        (coarse-resolution v2)
#   MSAT    -> output_cmaps/msa_transformer/Coarse* (coarse-resolution v2)
#
# Usage:
#   bash scripts/sweep_status.sh                       # all pairs
#   bash scripts/sweep_status.sh --pairs "PAIR1 PAIR2" # specific pairs
#   bash scripts/sweep_status.sh --summary             # counts only, no per-pair
#   bash scripts/sweep_status.sh --incomplete-only     # hide pairs where all 8 are done
# ============================================================================
set -eo pipefail   # NOT -u: empty associative-array dereference would trip it

PAIRS=""
SUMMARY_ONLY=0
INCOMPLETE_ONLY=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --pairs)            PAIRS="$2"; shift 2 ;;
        --summary)          SUMMARY_ONLY=1; shift ;;
        --incomplete-only)  INCOMPLETE_ONLY=1; shift ;;
        -h|--help)          sed -n '2,30p' "$0"; exit 0 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

if [[ -n "$PAIRS" ]]; then
    PAIR_LIST="$PAIRS"
else
    PAIR_LIST=$(ls -d Pipeline/FoldPairs/*/ 2>/dev/null \
                | sed -E 's#.*FoldPairs/##; s#/##' \
                | grep -vE '^(_s4pred_work|jobs)$')
fi

# ---- Currently running msaclust jobs (per-pair via job name) ----
declare -A RUNNING
while IFS=$'\t' read -r jname jid; do
    [[ -z "$jname" ]] && continue
    # job name is msaclust_<pair>
    pair="${jname#msaclust_}"
    RUNNING["$pair"]="$jid"
done < <(squeue -u "$USER" -h -o '%j	%i' 2>/dev/null | grep '^msaclust_')

# ---- Check each pair's 8 methods ----
n_total=0; n_done=0; n_partial=0; n_empty=0
declare -A METHOD_TOTAL
for M in AF2 AF3 ESM2 Boltz DDG S4P CCM MSAT; do METHOD_TOTAL["$M"]=0; done

check_method() {
    # args: pair_dir, method_name -> echoes ✓ or -
    local pdir="$1" method="$2"
    case "$method" in
        AF2)    compgen -G "${pdir}/output_AF/AF2/ShallowMsa_*"            > /dev/null 2>&1 ;;
        AF3)    compgen -G "${pdir}/output_AF/AF3/ShallowMsa_*"            > /dev/null 2>&1 ;;
        ESM2)   ls "${pdir}/output_esmfold2/"*.{pdb,cif} 2>/dev/null | grep -q .  ;;
        Boltz)  compgen -G "${pdir}/output_boltz2/ShallowMsa_*"            > /dev/null 2>&1 ;;
        DDG)    [[ -s "${pdir}/Analysis/df_ddg.csv" ]] ;;
        S4P)    [[ -s "${pdir}/Analysis/df_s4pred.csv" ]] ;;
        CCM)    compgen -G "${pdir}/output_cmaps/ccmpred/"*Coarse*         > /dev/null 2>&1 ;;
        MSAT)   compgen -G "${pdir}/output_cmaps/msa_transformer/"*Coarse* > /dev/null 2>&1 ;;
    esac
}

ROWS=()
for P in $PAIR_LIST; do
    PDIR="Pipeline/FoldPairs/$P"
    [[ ! -d "$PDIR" ]] && continue
    n_total=$((n_total+1))
    row="$P"
    done_count=0
    for M in AF2 AF3 ESM2 Boltz DDG S4P CCM MSAT; do
        if check_method "$PDIR" "$M"; then
            row+="\t✓"
            done_count=$((done_count+1))
            METHOD_TOTAL["$M"]=$((${METHOD_TOTAL["$M"]}+1))
        else
            row+="\t-"
        fi
    done
    jobinfo="-"
    if [[ -n "${RUNNING[$P]:-}" ]]; then
        jobinfo="R(${RUNNING[$P]})"
    fi
    row+="\t$jobinfo"
    if [[ "$done_count" -eq 8 ]]; then
        row+="\tDONE"
        n_done=$((n_done+1))
    elif [[ "$done_count" -eq 0 ]]; then
        row+="\tempty"
        n_empty=$((n_empty+1))
    else
        row+="\t${done_count}/8"
        n_partial=$((n_partial+1))
    fi
    # Hide complete rows if --incomplete-only
    if [[ "$INCOMPLETE_ONLY" -eq 1 ]] && [[ "$done_count" -eq 8 ]]; then
        :
    else
        ROWS+=("$row")
    fi
done

# ---- Print per-pair table ----
if [[ "$SUMMARY_ONLY" -eq 0 ]]; then
    {
        echo -e "Pair\tAF2\tAF3\tESM2\tBoltz\tDDG\tS4P\tCCM\tMSAT\tJob\tStatus"
        for r in "${ROWS[@]}"; do echo -e "$r"; done
    } | column -t -s $'\t'
    echo
fi

# ---- Summary ----
echo "==== Summary ===="
echo "  Pairs visited:  $n_total"
echo "  DONE (8/8):     $n_done"
echo "  Partial:        $n_partial"
echo "  Empty:          $n_empty"
echo "  Running jobs:   ${#RUNNING[@]}"
echo
echo "Per-method completion:"
for M in AF2 AF3 ESM2 Boltz DDG S4P CCM MSAT; do
    printf "  %-6s  %d / %d\n" "$M" "${METHOD_TOTAL[$M]}" "$n_total"
done
