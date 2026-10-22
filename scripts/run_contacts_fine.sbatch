#!/bin/bash
#SBATCH --job-name=contacts_fine
#SBATCH --account=orzuk
#SBATCH --partition=glacier
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --chdir=/sci/labs/orzuk/orzuk/github/MsaCluster
#SBATCH --output=Pipeline/FoldPairs/jobs/contacts_fine_%j.out
# ===========================================================================
# (Re)compute MSA-Transformer or CCMpred contacts on the CURRENT fine clusters
# (output_msa_cluster/ShallowMsa_*.a3m) for many pairs inside ONE job (the
# Moriah ~60 jobs/hour submission limit rules out per-pair fan-out), then
# re-score that method's per-pair CSV (Analysis/df_cmap.csv or
# df_cmap_ccmpred.csv) with the provenance guard (+ top-L for CCMpred).
#
# Only clusters without a matching provenance manifest entry are run, so the
# job is resumable: resubmit after a timeout and it continues.
#
#   CCMpred (CPU; defaults above; pairs in parallel, THREADS threads each):
#     sbatch scripts/run_contacts_fine.sbatch ccmpred ALL
#   MSA-Transformer (GPU; pairs sequential on one GPU):
#     sbatch -p salmon --gres=gpu:l40s:1 --cpus-per-task=8 --mem=48G \
#            --time=24:00:00 scripts/run_contacts_fine.sbatch msat ALL
#   Specific pairs instead of ALL:  ... ccmpred 1dzlA_5keqF 2qkeE_5jytA
#
# Env knobs: THREADS (CCMpred threads per pair, default 4; pairs in parallel =
# cpus / THREADS), RESCORE=0 to skip the re-score step.
# Per-pair logs: Pipeline/FoldPairs/<pair>/jobs/contacts_fine_<method>_<pair>.out
# ===========================================================================
set -euo pipefail
source /sci/labs/orzuk/orzuk/venvs/my-python-venv/bin/activate
export PYTHONUNBUFFERED=1

METHOD="${1:-}"
shift || true
case "$METHOD" in
    msat)    RUN_MODE=run_cmap_msa_transformer; CMAP_METHOD=msa_transformer ;;
    ccmpred) RUN_MODE=run_cmap_ccmpred;         CMAP_METHOD=ccmpred ;;
    *) echo "usage: run_contacts_fine.sbatch msat|ccmpred <pairs...|ALL>"; exit 2 ;;
esac

if [[ "${1:-}" == "ALL" || $# -eq 0 ]]; then
    mapfile -t PAIRS < <(ls -d Pipeline/FoldPairs/*/ | xargs -n1 basename \
                         | grep -E '^[0-9a-zA-Z]{5}_[0-9a-zA-Z]{5}$')
else
    PAIRS=("$@")
fi

NCPU="${SLURM_CPUS_PER_TASK:-8}"
if [[ "$METHOD" == "msat" ]]; then
    python3 -c "import torch, sys; sys.exit(0 if torch.cuda.is_available() else 1)" \
        || { echo "[contacts-fine] ERROR: no CUDA device; submit msat with -p salmon --gres=gpu:l40s:1"; exit 3; }
    NPAR=1                       # one GPU: pairs sequentially
    EXTRA=()
else
    THREADS="${THREADS:-4}"
    NPAR=$(( NCPU / THREADS )); (( NPAR >= 1 )) || NPAR=1
    EXTRA=(--ccmpred_threads "$THREADS")
fi
echo "[contacts-fine] method=$METHOD pairs=${#PAIRS[@]} parallel=$NPAR host=$(hostname) commit=$(git rev-parse --short HEAD) start=$(date)"

run_one() {
    local run_mode="$1" method="$2" pair="$3"; shift 3
    mkdir -p "Pipeline/FoldPairs/$pair/jobs"
    local log="Pipeline/FoldPairs/$pair/jobs/contacts_fine_${method}_${pair}.out"
    local t0=$SECONDS
    if python3 run_foldswitch_pipeline.py --run_mode "$run_mode" --foldpair_ids "$pair" \
            --run_job_mode inline --cluster_resolution fine "$@" > "$log" 2>&1; then
        echo "OK   $pair $((SECONDS - t0))s $(grep -hE '^\[(ccmpred|msat)\] ' "$log" | tr '\n' ' ' | cut -c1-200)"
    else
        echo "FAIL $pair $((SECONDS - t0))s (see $log)"
    fi
}
export -f run_one

printf '%s\n' "${PAIRS[@]}" | xargs -P "$NPAR" -I{} bash -c 'run_one "$@"' _ \
    "$RUN_MODE" "$METHOD" {} "${EXTRA[@]}"
echo "[contacts-fine] contacts done=$(date)"

if [[ "${RESCORE:-1}" == "1" ]]; then
    # Re-score only this method's CSV (provenance-checked; CCMpred top-L=1.0*L).
    # Equivalent for both methods at once:
    #   python3 run_foldswitch_pipeline.py --run_mode postprocess --foldpair_ids ALL \
    #       --run_job_mode inline --force_cmap TRUE
    echo "[contacts-fine] re-scoring $CMAP_METHOD (parallel=$NCPU)"
    printf '%s\n' "${PAIRS[@]}" | xargs -P "$NCPU" -n 4 \
        python3 Analysis/cmap_analysis.py --method "$CMAP_METHOD" --resolution fine \
        > "Pipeline/FoldPairs/jobs/contacts_fine_rescore_${METHOD}_${SLURM_JOB_ID:-local}.log" 2>&1 \
        || echo "[contacts-fine] WARN some re-score calls failed (see log)"
    python3 - "$CMAP_METHOD" "${PAIRS[@]}" <<'PY'
import os, sys
import pandas as pd
method, pairs = sys.argv[1], sys.argv[2:]
csv = "df_cmap_ccmpred.csv" if method == "ccmpred" else "df_cmap.csv"
n_pairs = n_rows = n_shallow = 0
for p in pairs:
    f = os.path.join("Pipeline/FoldPairs", p, "Analysis", csv)
    try:
        df = pd.read_csv(f)
    except Exception:
        continue
    if df.empty or "cluster_hash" not in df.columns:
        continue
    n_pairs += 1; n_rows += len(df)
    n_shallow += int(df["cluster"].astype(str).str.contains("ShallowMsa_").sum())
print(f"[contacts-fine] {csv}: {n_pairs}/{len(pairs)} pairs with provenance-tagged rows; "
      f"{n_rows} rows ({n_shallow} fine-cluster rows)")
PY
fi
echo "[contacts-fine] done=$(date)"
