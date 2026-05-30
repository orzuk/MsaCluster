#!/bin/bash
# ============================================================================
# run_v2_subset_analysis.sh -- run downstream analysis on DONE pairs only.
#
# When the sweep finishes some-but-not-all pairs, this gives you a quick
# preview of the v2 numbers without waiting for everything. Auto-detects
# DONE pairs from sweep_status, then runs:
#   1. fold_diversity_survey  (cross-method per-cluster matrix)
#   2. aggregate_s4pred       (add S4PRED rows)
#   3. permutation_null       (per-pair concordance p-values)
#   4. BH FDR full pool       (does any pair survive?)
#   5. BH FDR stratified by trigger class
#   6. cross-method correlation matrix (figure)
#   7. (optional) HTML regen for each DONE pair
#
# Output: everything under docs/. The BH FDR result is in
# docs/fold_diversity_concordance_perm_bh.csv and
# docs/stratified_concordance.csv.
#
# Usage:
#   bash scripts/run_v2_subset_analysis.sh                # auto-detect DONE
#   bash scripts/run_v2_subset_analysis.sh --pairs "P1 P2 P3"
#   bash scripts/run_v2_subset_analysis.sh --no-html      # skip HTML step
#   bash scripts/run_v2_subset_analysis.sh --tag v2_32   # add suffix to outputs
# ============================================================================
set -eo pipefail

PAIRS_OVERRIDE=""
DO_HTML=1
TAG=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --pairs)    PAIRS_OVERRIDE="$2"; shift 2 ;;
        --no-html)  DO_HTML=0; shift ;;
        --tag)      TAG="$2"; shift 2 ;;
        -h|--help)  sed -n '2,30p' "$0"; exit 0 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

REPO_DIR=$(pwd)
[[ -f run_foldswitch_pipeline.py ]] || { echo "Run from MsaCluster repo root"; exit 1; }
mkdir -p docs

# ---- 0. Identify DONE pairs ----
if [[ -n "$PAIRS_OVERRIDE" ]]; then
    PAIR_LIST="$PAIRS_OVERRIDE"
    echo "[v2-analysis] Using --pairs override: $(echo $PAIR_LIST | wc -w) pairs"
else
    # Auto-detect from sweep_status: any pair whose row ends in "DONE"
    PAIR_LIST=$(bash scripts/sweep_status.sh 2>/dev/null \
                | awk '$NF == "DONE" {print $1}' \
                | tr '\n' ' ')
    N_DONE=$(echo $PAIR_LIST | wc -w)
    if [[ "$N_DONE" -eq 0 ]]; then
        echo "[v2-analysis] No DONE pairs detected. Run a sweep first."
        exit 1
    fi
    echo "[v2-analysis] Auto-detected $N_DONE DONE pairs"
fi
PAIR_LIST_COMMA=$(echo $PAIR_LIST | tr ' ' ',')

# Save the input list for reproducibility
OUT_SUFFIX=""
[[ -n "$TAG" ]] && OUT_SUFFIX=".$TAG"
PAIRS_FILE="docs/v2_subset_pairs${OUT_SUFFIX}.txt"
echo "$PAIR_LIST" | tr ' ' '\n' > "$PAIRS_FILE"
echo "[v2-analysis] Saved pair list to $PAIRS_FILE"

echo
echo "==================================================="
echo "STEP 1/7: fold_diversity_survey (cross-method matrix)"
echo "==================================================="
python3 scripts/fold_diversity_survey.py --pairs "$PAIR_LIST_COMMA" || echo "(continuing)"

echo
echo "==================================================="
echo "STEP 2/7: aggregate S4PRED into survey"
echo "==================================================="
python3 scripts/aggregate_s4pred_into_survey.py --pairs "$PAIR_LIST_COMMA" \
    --write_mode replace_method 2>&1 | tail -10 || echo "(continuing)"

echo
echo "==================================================="
echo "STEP 3/7: permutation null for concordance p-values"
echo "==================================================="
python3 scripts/permutation_null_concordance.py --n-perm 1000 --seed 12345 \
    2>&1 | tail -10 || echo "(continuing)"

echo
echo "==================================================="
echo "STEP 4/7: BH FDR (full pool)"
echo "==================================================="
python3 scripts/bh_correct_concordance_p.py --alpha 0.05 2>&1 | tail -10 || \
    echo "(continuing)"

echo
echo "==================================================="
echo "STEP 5/7: BH FDR (stratified by trigger class)"
echo "==================================================="
python3 scripts/stratified_concordance_by_trigger.py --alpha 0.05 \
    2>&1 | tail -10 || echo "(continuing)"

echo
echo "==================================================="
echo "STEP 6/7: cross-method correlation matrix"
echo "==================================================="
python3 scripts/cross_method_correlation.py 2>&1 | tail -10 || echo "(continuing)"

if [[ "$DO_HTML" -eq 1 ]]; then
    echo
    echo "==================================================="
    echo "STEP 7/7: HTML regen for each DONE pair (slow)"
    echo "==================================================="
    python3 TableResults/gen_pair_html.py --pairs "$PAIR_LIST_COMMA" --mode inline \
        2>&1 | tail -10 || echo "(continuing)"
else
    echo "[v2-analysis] Skipping HTML step (--no-html)"
fi

echo
echo "==================================================="
echo "DONE. Headline files:"
echo "==================================================="
echo "  docs/fold_diversity_survey.csv             (per-cluster matrix, all methods)"
echo "  docs/fold_diversity_concordance_perm.csv   (per-pair p-values)"
echo "  docs/fold_diversity_concordance_perm_bh.csv (BH FDR full pool)"
echo "  docs/stratified_concordance.csv            (BH FDR stratified)"
echo "  docs/figs/cross_method_correlation.png      (method correlation)"
echo "  $PAIRS_FILE     (pair list used)"
[[ "$DO_HTML" -eq 1 ]] && echo "  docs/HTML/figs/*.html                       (per-pair pages)"

echo
echo "Quick view of BH results:"
echo "--- Full pool ---"
[[ -f docs/fold_diversity_concordance_perm_bh.csv ]] && \
    head -5 docs/fold_diversity_concordance_perm_bh.csv

echo "--- Stratified ---"
[[ -f docs/stratified_concordance.csv ]] && head docs/stratified_concordance.csv
