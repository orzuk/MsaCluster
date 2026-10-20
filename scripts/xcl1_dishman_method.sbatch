#!/bin/bash
#SBATCH --job-name=xcl1_pilot
#SBATCH --account=orzuk
#SBATCH --chdir=/sci/labs/orzuk/orzuk/github/MsaCluster
#SBATCH --output=Pipeline/Pilots/xcl1_dishman/jobs/method_%x_%j.out
# ===========================================================================
# One prediction method for the XCL1 / Dishman curated-MSA pilot, run INLINE
# over all clusters inside a SINGLE job (the cluster allows ~60 submissions
# per hour, so do not fan out per cluster).
#
#   sbatch [-p salmon --gres=gpu:l40s:1 ...] scripts/xcl1_dishman_method.sbatch \
#          <af2|af3|af|esmfold2|boltz2|bioemu|ddg|s4pred> [pilot_root] [pair_id]
#
# Partition / GPU / time are chosen by the launcher
# (scripts/shell/run_xcl1_dishman_pilot.sh) on the sbatch command line.
#
# All pipeline paths are redirected to <pilot_root>/<pair_id>/ by
# MSACLUSTER_DATA_DIR (= run_foldswitch_pipeline.py --data_root), so the
# production Pipeline/FoldPairs/2n54B_2hdmA outputs are never touched.
# ===========================================================================
set -euo pipefail

METHOD="${1:-}"
PILOT_ROOT="${2:-/sci/labs/orzuk/orzuk/github/MsaCluster/Pipeline/Pilots/xcl1_dishman}"
PAIR="${3:-2n54B_2hdmA}"
[[ -n "$METHOD" ]] || { echo "usage: $0 <method> [pilot_root] [pair_id]"; exit 2; }

REPO=/sci/labs/orzuk/orzuk/github/MsaCluster
VENV=/sci/labs/orzuk/orzuk/venvs/my-python-venv
ESMFOLD2_VENV=/sci/labs/orzuk/orzuk/venvs/esmfold2-venv
S4PRED_DIR=/sci/labs/orzuk/orzuk/github/s4pred
THERMOMPNN_DIR=/sci/labs/orzuk/orzuk/github/ThermoMPNN
HF_CACHE=/sci/labs/orzuk/orzuk/tmp/xdg-cache

cd "$REPO"
export MSACLUSTER_DATA_DIR="$PILOT_ROOT"
export PYTHONUNBUFFERED=1
mkdir -p "$PILOT_ROOT/jobs"

echo "[xcl1-pilot] method=$METHOD pair=$PAIR root=$PILOT_ROOT host=$(hostname) start=$(date)"

PIPE=(python3 run_foldswitch_pipeline.py --foldpair_ids "$PAIR"
      --data_root "$PILOT_ROOT" --run_job_mode inline)

case "$METHOD" in
  af2|af3|af)
    # shellcheck disable=SC1091
    source "$VENV/bin/activate"
    case "$METHOD" in af2) VER=2 ;; af3) VER=3 ;; *) VER=both ;; esac
    "${PIPE[@]}" --run_mode run_AF --af_ver "$VER" \
        --query_type medoid --af_msa_top_n 10 \
        --af_output_suffix dishman --allow_inline_af
    ;;
  esmfold2)
    # Separate venv (own transformers/esm pins); run_ESMFold2.py reads
    # MSACLUSTER_DATA_DIR for the pair root, so it writes into the pilot dir.
    export XDG_CACHE_HOME="$HF_CACHE" HF_HOME="$HF_CACHE/huggingface" HF_HUB_OFFLINE=1
    # shellcheck disable=SC1091
    source "$ESMFOLD2_VENV/bin/activate"
    python3 run_ESMFold2.py --foldpair_ids "$PAIR" \
        --representative_method medoid --skip_existing
    ;;
  boltz2)
    # shellcheck disable=SC1091
    source "$VENV/bin/activate"
    "${PIPE[@]}" --run_mode run_boltz2 --boltz2_mode apo
    ;;
  bioemu)
    # shellcheck disable=SC1091
    source "$VENV/bin/activate"
    "${PIPE[@]}" --run_mode run_bioemu --query_type medoid --af_msa_top_n 10
    ;;
  ddg)
    # shellcheck disable=SC1091
    source "$VENV/bin/activate"
    "${PIPE[@]}" --run_mode run_ddg --thermompnn_dir "$THERMOMPNN_DIR" --ddg_device cuda
    ;;
  s4pred)
    # shellcheck disable=SC1091
    source "$VENV/bin/activate"
    "${PIPE[@]}" --run_mode run_s4pred --s4pred_dir "$S4PRED_DIR" --s4pred_device cpu \
        --s4pred_work_root "$PILOT_ROOT/_s4pred_work"
    ;;
  *)
    echo "unknown method: $METHOD" >&2
    exit 2
    ;;
esac

echo "[xcl1-pilot] method=$METHOD done=$(date)"
