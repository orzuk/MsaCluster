#!/bin/bash
#SBATCH --job-name=xcl1_pilot_analysis
#SBATCH --account=orzuk
#SBATCH --partition=glacier
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --chdir=/sci/labs/orzuk/orzuk/github/MsaCluster
#SBATCH --output=Pipeline/Pilots/xcl1_dishman/jobs/analysis_%j.out
# ===========================================================================
# Final CPU step of the XCL1 / Dishman curated-MSA pilot: score, survey and
# phylo-test the pilot clusters. Submitted with an afterok dependency on the
# per-method jobs by scripts/shell/run_xcl1_dishman_pilot.sh.
#
#   sbatch scripts/xcl1_dishman_analysis.sbatch [pilot_root] [pair_id]
#
# Everything is written under <pilot_root>/ (per-pair Analysis CSVs) and
# <pilot_root>/results/ (survey + phylo CSVs). docs/ (the production tables)
# is never written: postprocess runs with --per_pair_only and every analysis
# script gets an explicit pilot output path.
# ===========================================================================
set -euo pipefail

PILOT_ROOT="${1:-/sci/labs/orzuk/orzuk/github/MsaCluster/Pipeline/Pilots/xcl1_dishman}"
PAIR="${2:-2n54B_2hdmA}"
REPO=/sci/labs/orzuk/orzuk/github/MsaCluster
VENV=/sci/labs/orzuk/orzuk/venvs/my-python-venv
MANIFEST="$REPO/data/pilots/xcl1_dishman/manifest.tsv"

cd "$REPO"
# shellcheck disable=SC1091
source "$VENV/bin/activate"
export MSACLUSTER_DATA_DIR="$PILOT_ROOT"
export PYTHONUNBUFFERED=1

RES="$PILOT_ROOT/results"
mkdir -p "$RES" "$PILOT_ROOT/jobs"
SURVEY="$RES/fold_diversity_survey.csv"
PHYLO="$RES/phylo_placement.csv"

echo "[xcl1-analysis] pair=$PAIR root=$PILOT_ROOT host=$(hostname) start=$(date)"
echo "[xcl1-analysis] clusters: $(ls "$PILOT_ROOT/$PAIR/output_msa_cluster"/ShallowMsa_*.a3m 2>/dev/null | wc -l)"

echo "=== a: per-pair Analysis CSVs (df_af / df_esm / df_cmap; pilot only) ==="
python3 Analysis/postprocess_unified.py --pairs "$PAIR" --mode full --force_rerun --per_pair_only

echo "=== b: ESMFold2 re-score vs the two truth folds ==="
python3 scripts/score_esmfold2.py --pair "$PAIR" \
    --mode_label ESMFold2_dishman --pred_subdir output_esmfold2 --output_dir "$RES"

echo "=== c: fold-diversity survey (pilot clusters) ==="
python3 scripts/fold_diversity_survey.py --pairs "$PAIR" --delta 0.05 --delta-ddg 2.0 \
    --msat-metric t_recall --esmfold2-dir "$RES" --output-dir "$RES"

echo "=== d: aggregate S4PRED rows into the pilot survey ==="
python3 scripts/aggregate_s4pred_into_survey.py --pairs "$PAIR" \
    --fpairs_root "$PILOT_ROOT" --survey_csv "$SURVEY" --write_mode replace_method

echo "=== e: phylogenetic tests on the pilot tree ==="
# The clade-shift scan (utils/phylo_tests.py + phylo_placement --tests) lives on
# branch feature/clade-shift. Use it when it is merged; otherwise fall back to
# Moran's I only.
# TODO(clade-shift): once feature/clade-shift is merged into master, drop the
# fallback branch below and always pass
#   --tests morans,clade_shift,clade_shift_combined --clade-n-perm 10000
if python3 -c "import utils.phylo_tests" 2>/dev/null \
   && python3 scripts/phylo_placement.py --help 2>/dev/null | grep -q -- "--tests"; then
    echo "[xcl1-analysis] clade-shift available"
    python3 scripts/phylo_placement.py --pairs "$PAIR" --labels centered \
        --survey-csv "$SURVEY" --output "$PHYLO" \
        --tests morans,clade_shift,clade_shift_combined --clade-n-perm 10000 \
        --n-perm 1000 --seed 12345 --delta-tm 0.05 --delta-ddg 2.0 \
        --no-d-statistic --make-figures 0 --figure-dir "$RES/figures"
else
    echo "[xcl1-analysis] TODO: feature/clade-shift is NOT merged -- running Moran's I only;"
    echo "[xcl1-analysis] TODO: re-run with --tests morans,clade_shift,clade_shift_combined after the merge."
    python3 scripts/phylo_placement.py --pairs "$PAIR" --labels centered \
        --survey-csv "$SURVEY" --output "$PHYLO" \
        --n-perm 1000 --seed 12345 --delta-tm 0.05 --delta-ddg 2.0 \
        --no-d-statistic --make-figures 0 --figure-dir "$RES/figures"
fi

echo "=== f: cluster family composition + clade verdict vs human XCL1 ==="
PHYLO_ARG=()
[[ -f "$PHYLO" ]] && PHYLO_ARG=(--phylo_csv "$PHYLO")
python3 scripts/build_curated_chemokine_msa.py annotate \
    --pair_dir "$PILOT_ROOT/$PAIR" --manifest "$MANIFEST" \
    --output "$RES/cluster_family_composition.csv" "${PHYLO_ARG[@]}"

echo "[xcl1-analysis] done=$(date). Results in $RES"
