#!/bin/bash
#SBATCH --job-name=final_chain
#SBATCH --account=orzuk
#SBATCH --partition=glacier
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --chdir=/sci/labs/orzuk/orzuk/github/MsaCluster
#SBATCH --output=Pipeline/FoldPairs/jobs/final_chain_%j.out
# ===========================================================================
# Final CPU-only re-run of the downstream analysis (all 9 methods, truth-length
# TM normalisation, DDG threshold 2.0 kcal/mol). Every flag is explicit so the
# run does not depend on script defaults. Stops on the first error.
#
# Prerequisite: Boltz-2 / BioEmu per-pair CSVs re-scored
# (run_foldswitch_pipeline.py --run_mode run_boltz2|run_bioemu --rescore_only TRUE).
# scripts/shell/run_final_rerun.sh submits those and then this job with an
# afterany dependency. Standalone:  sbatch scripts/run_final_chain.sbatch
#
# Steps
#   0  back up docs/*.csv -> docs/backup_<date>/
#   a  ESMFold2 re-score (existing PDBs, no prediction)
#   b  survey + S4PRED aggregation
#   c  sequence-divergence correction (PID cache rebuilt)
#   d  corrected concordance permutation (+ BH)
#   e  raw concordance permutation (+ BH)
#   f  stratified-by-trigger BH + Chakravarty Fisher test
#   g  Moran's I, raw and corrected
#   h  cross-method correlation + global plots
#   i  summary HTML table
# ===========================================================================
set -euo pipefail

REPO=/sci/labs/orzuk/orzuk/github/MsaCluster
VENV=/sci/labs/orzuk/orzuk/venvs/my-python-venv
cd "$REPO"
# shellcheck disable=SC1091
source "$VENV/bin/activate"
export PYTHONUNBUFFERED=1

ALL9="AF2,AF3,Boltz2,BioEmu,ESM,MSAT,CCMpred,DDG,S4PRED"
DELTA_TM=0.05
DELTA_DDG=2.0
SEED=12345
SUF=_9method

SURVEY=docs/fold_diversity_survey.csv
SURVEY_COR=docs/fold_diversity_survey_corrected${SUF}.csv
CONCORD_COR=docs/fold_diversity_concordance_corrected${SUF}.csv
PERM_COR=docs/fold_diversity_concordance_corrected_perm${SUF}.csv
PERM_COR_BH=docs/fold_diversity_concordance_corrected_perm${SUF}_bh.csv
PERM_RAW=docs/fold_diversity_concordance_perm${SUF}.csv
PERM_RAW_BH=docs/fold_diversity_concordance_perm${SUF}_bh.csv
STRAT_COR=docs/stratified_by_trigger_corrected${SUF}.csv
PHYLO_RAW=docs/phylo_placement.csv
PHYLO_COR=docs/phylo_placement_corrected.csv

step() { echo; echo "=== [$(date +%H:%M:%S)] $* ==="; }

# Print the number of distinct pairs in a CSV vs the expected count (a
# mismatch is reported, not fatal). Optional 3rd arg: also report how many
# pairs have a non-null value in that column.
count_pairs() {
    python3 - "$1" "$2" "${3:-}" <<'PY'
import sys
import pandas as pd
path, expected, col = sys.argv[1], int(sys.argv[2]), sys.argv[3]
df = pd.read_csv(path)
n = df["pair_id"].nunique()
tag = "OK" if n == expected else "WARN MISMATCH"
extra = ""
if col:
    k = df[pd.to_numeric(df[col], errors="coerce").notna()]["pair_id"].nunique()
    extra = f"; {k} with non-null {col}"
print(f"[pair-count] {path}: {n} pairs (expected {expected}) {tag}{extra}")
PY
}

echo "[final-chain] host=$(hostname) start=$(date) commit=$(git rev-parse --short HEAD)"
echo "[final-chain] methods=$ALL9 delta_tm=$DELTA_TM delta_ddg=$DELTA_DDG"

step "0: back up docs/*.csv"
BACKUP=docs/backup_$(date +%Y%m%d_%H%M%S)
mkdir -p "$BACKUP"
cp -p docs/*.csv "$BACKUP"/
echo "backed up $(ls "$BACKUP" | wc -l) CSVs to $BACKUP"

step "a: ESMFold2 re-score (truth-length TM; no --skip_existing => every pair re-scored)"
python3 scripts/score_esmfold2.py --foldpair_ids ALL \
    --mode_label ESMFold2_newmethod --pred_subdir output_esmfold2 --output_dir docs
echo "esmfold2 score CSVs: $(ls docs/esmfold2_tmscores_*.csv | wc -l)"

step "b1: fold-diversity survey (all pairs)"
python3 scripts/fold_diversity_survey.py --pairs ALL --delta "$DELTA_TM" \
    --delta-ddg "$DELTA_DDG" --msat-metric t_recall --output-dir docs
step "b2: aggregate S4PRED into survey (replace S4PRED rows)"
python3 scripts/aggregate_s4pred_into_survey.py --pairs ALL \
    --fpairs_root Pipeline/FoldPairs --survey_csv "$SURVEY" --write_mode replace_method
count_pairs "$SURVEY" 93
python3 - "$SURVEY" "$ALL9" <<'PY'
import sys
import pandas as pd
df = pd.read_csv(sys.argv[1])
want = sys.argv[2].split(",")
per = df.groupby("method")["pair_id"].nunique()
print("[survey] pairs per method:", per.reindex(want).fillna(0).astype(int).to_dict())
missing = [m for m in want if m not in per.index]
if missing:
    sys.exit(f"[survey] ERROR: methods missing from survey: {missing}")
PY

step "c: sequence-divergence correction (rebuilds docs/per_cluster_pid_symmetric.csv)"
python3 scripts/seq_divergence_correction.py --rebuild-pid-cache \
    --methods "$ALL9" --out-suffix "$SUF" --delta "$DELTA_TM" --delta-ddg "$DELTA_DDG"
count_pairs "$SURVEY_COR" 93
count_pairs "$CONCORD_COR" 92 mean_concordance_corrected
python3 - "$SURVEY_COR" <<'PY'
import sys
import pandas as pd
df = pd.read_csv(sys.argv[1])
n_f2 = int(pd.to_numeric(df["mean_pid_to_F2"], errors="coerce").notna().sum())
print(f"[correction] rows with non-null mean_pid_to_F2: {n_f2}/{len(df)}")
fits = df.groupby(["pair_id", "method"])["regression_n_predictors"].max()
print(f"[correction] (pair,method) fits: {len(fits)} total, "
      f"{int((fits == 2).sum())} with 2 predictors, "
      f"{int((fits == 1).sum())} with 1, {int((fits == 0).sum())} uncorrected")
if n_f2 == 0:
    sys.exit("[correction] ERROR: mean_pid_to_F2 is empty -- PID cache rebuild failed")
PY

step "d: corrected concordance permutation (100000 perms, + BH)"
python3 scripts/permutation_null_concordance.py --mode corrected --methods "$ALL9" \
    --input "$SURVEY_COR" --output "$PERM_COR" --n-perm 100000 --seed "$SEED" --alpha 0.05
count_pairs "$PERM_COR_BH" 92 perm_mean_concordance_p

step "e: raw concordance permutation (10000 perms, + BH)"
python3 scripts/permutation_null_concordance.py --mode raw --methods "$ALL9" \
    --input "$SURVEY" --output "$PERM_RAW" --n-perm 10000 --seed "$SEED" \
    --delta "$DELTA_TM" --delta-ddg "$DELTA_DDG"
python3 scripts/bh_correct_concordance_p.py --input "$PERM_RAW" --output "$PERM_RAW_BH" --alpha 0.05
count_pairs "$PERM_RAW_BH" 92 perm_mean_concordance_p

step "f: stratified-by-trigger BH on corrected perms + Chakravarty Fisher test"
python3 scripts/stratified_concordance_by_trigger.py --perm-csv "$PERM_COR" \
    --triggers-csv docs/triggers_from_pdb.csv --output "$STRAT_COR" \
    --alpha 0.05 --min-class-size 10

step "g1: Moran's I, raw labels"
python3 scripts/phylo_placement.py --pairs ALL --labels centered --survey-csv "$SURVEY" \
    --output "$PHYLO_RAW" --n-perm 1000 --seed "$SEED" \
    --delta-tm "$DELTA_TM" --delta-ddg "$DELTA_DDG" --no-d-statistic --make-figures 0
count_pairs "$PHYLO_RAW" 89 morans_p
step "g2: Moran's I, corrected labels"
python3 scripts/phylo_placement.py --pairs ALL --labels corrected --survey-csv "$SURVEY_COR" \
    --output "$PHYLO_COR" --n-perm 1000 --seed "$SEED" \
    --delta-tm "$DELTA_TM" --delta-ddg "$DELTA_DDG" --no-d-statistic --make-figures 0
count_pairs "$PHYLO_COR" 89 morans_p

step "h1: cross-method correlation (all 9 methods via utils.method_labels.METHOD_ORDER)"
python3 scripts/cross_method_correlation.py
step "h2: global plots"
python3 run_foldswitch_pipeline.py --run_mode plot --foldpair_ids ALL \
    --run_job_mode inline --plot_scope global

step "i: summary HTML table (reads $PERM_COR_BH)"
python3 TableResults/gen_html_table.py

echo
echo "[final-chain] done $(date). Backup of previous CSVs: $BACKUP"
