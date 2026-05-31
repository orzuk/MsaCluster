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
echo "STEP 0/7: score ESMFold2 predictions vs folds"
echo "==================================================="
# Writes docs/esmfold2_tmscores_<pair>.csv (one medoid prediction/cluster,
# scored vs both truth folds). The survey's 'ESM' slot now reads these
# (ESMFold2 = v2 single-sequence baseline; replaces ESMFold v1).
# --skip_existing: cache successfully-scored pairs so re-runs are fast
# (failed/all-NaN CSVs are still re-scored). Force a full re-score by
# deleting docs/esmfold2_tmscores_*.csv first.
python3 scripts/score_esmfold2.py --foldpair_ids "$PAIR_LIST_COMMA" --skip_existing \
    2>&1 | tail -15 || echo "(continuing)"

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
echo "STEP 6a/7: sequence-divergence regression correction"
echo "==================================================="
# Strips the per-cluster cluster_to_query seq-identity gradient from
# per-cluster TM-diffs (and ΔΔG). Required before residual-corrected
# phylo placement — v1's KaiB hit (p_NN1=0.005) was on CORRECTED labels.
python3 scripts/seq_divergence_correction.py 2>&1 | tail -10 \
    || echo "(continuing)"

echo
echo "==================================================="
echo "STEP 6b/7: phylogenetic placement (p_NN1 / Fitch) — raw + corrected"
echo "==================================================="
# Per-(pair,method) NN1 + Fitch parsimony permutation tests on F1c/F2c
# labels projected onto the cluster-level phylo tree. Output:
#   docs/phylo_placement.csv (raw) and docs/phylo_placement_corrected.csv
echo "--- raw labels -> docs/phylo_placement.csv ---"
python3 scripts/phylo_placement.py --pairs "$PAIR_LIST_COMMA" --no-d-statistic \
    2>&1 | tail -15 || echo "(continuing)"
echo "--- corrected labels (seq-divergence residuals) -> docs/phylo_placement_corrected.csv ---"
python3 scripts/phylo_placement.py --pairs "$PAIR_LIST_COMMA" --labels corrected --no-d-statistic \
    2>&1 | tail -15 || echo "(continuing)"

# BH within method, across pairs — standard v1 approach
echo
echo "--- BH FDR per method on phylo_placement results ---"
python3 - <<'PY' 2>&1 | tail -40
import pandas as pd
from pathlib import Path
import sys

def bh_within_method(df, p_col, alpha=0.05):
    """BH within each method's column, across pairs."""
    results = []
    for method in df['method'].unique():
        m = df[df['method'] == method].copy()
        m = m.dropna(subset=[p_col]).sort_values(p_col).reset_index(drop=True)
        n = len(m)
        if n == 0:
            continue
        m['bh_thresh'] = (m.index + 1) / n * alpha
        passing = m[m[p_col] <= m['bh_thresh']]
        results.append({
            'method': method,
            'n_pairs': n,
            'n_BH_sig': len(passing),
            'min_raw_p': m[p_col].min(),
            'passing_pairs': ', '.join(passing['pair_id'].tolist()) if len(passing) else '-',
        })
    return pd.DataFrame(results)

for src in ['docs/phylo_placement.csv', 'docs/phylo_placement_corrected.csv']:
    if not Path(src).is_file():
        continue
    df = pd.read_csv(src)
    pcols = [c for c in ('morans_p','nn_concordance_p','parsimony_p','nn3_concordance_p') if c in df.columns]
    print(f"\n=== {src} ===")
    print(f"columns available: {pcols}")
    for pcol in pcols:
        out = bh_within_method(df, pcol, alpha=0.05)
        print(f"\n  BH per method on {pcol}:")
        print(out.to_string(index=False))
PY


echo
echo "==================================================="
echo "STEP 6c/7: cross-method correlation matrix"
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
echo "  docs/phylo_placement.csv                   (per-(pair,method) p_NN1 + Fitch)"
echo "  $PAIRS_FILE     (pair list used)"
[[ "$DO_HTML" -eq 1 ]] && echo "  docs/HTML/figs/*.html                       (per-pair pages)"

echo
echo "Quick view of BH results:"
echo "--- Full pool ---"
[[ -f docs/fold_diversity_concordance_perm_bh.csv ]] && \
    head -5 docs/fold_diversity_concordance_perm_bh.csv

echo "--- Stratified ---"
[[ -f docs/stratified_concordance.csv ]] && head docs/stratified_concordance.csv

echo
echo "--- Phylo placement: top hits (best p_NN1 per pair) ---"
if [[ -f docs/phylo_placement.csv ]]; then
    python3 -c "
import pandas as pd
try:
    df = pd.read_csv('docs/phylo_placement.csv')
    # Find best (lowest) p_NN1 per pair across methods
    pcol = [c for c in df.columns if 'p_NN1' in c or 'p_nn1' in c]
    if pcol:
        p = pcol[0]
        best = df.loc[df.groupby('pair_id')[p].idxmin()][['pair_id', 'method', p]].sort_values(p).head(10)
        print(best.to_string(index=False))
    else:
        print('p_NN1 column not found; columns:', list(df.columns))
except Exception as e:
    print(f'failed: {e}')
"
fi
