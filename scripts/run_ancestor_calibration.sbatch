#!/bin/bash
#SBATCH --job-name=anc_calib
#SBATCH --account=orzuk
#SBATCH --partition=salmon
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --chdir=/sci/labs/orzuk/orzuk/github/MsaCluster
#SBATCH --output=ancestor_calibration.out
#
# Calibration experiment: do our per-sequence fold-preference predictors
# reproduce the NMR-measured fold populations of the reconstructed XCL1
# ancestors (Dishman et al., Science 2021, 371:86-90)?
#
# Outputs
#   docs/ancestor_calibration.csv            one row per (sequence, method, mode)
#   docs/ancestor_calibration_stats.csv      rho + exact permutation p + control null
#   docs/figures/global/ancestor_calibration.png (+ .pdf)
#
# Only 5 ancestral sequences + 20 random family controls, so this is a small
# job. Rough costs on one L40S:
#   ESMFold2  25 single-sequence folds of ~70 aa      ~5 min
#   DDG       2 ThermoMPNN backbones + table lookups  ~1 min
#   S4PRED    25 single sequences                     seconds (CPU)
#   Boltz-2   25 short single-sequence predictions    ~20-40 min
#   AF2       25 x 5 models of a 70-mer, per mode     ~25-40 min per mode
# The default method list is everything but AF2 and BioEmu, which are
# submitted separately (below) so the common case finishes in well under
# an hour.
#
# Knobs (environment, all optional):
#   ANC_METHODS   comma list from esmfold2,ddg,s4pred,boltz2,bioemu,af2,af2_single
#                 (default: esmfold2,ddg,s4pred,boltz2)
#   ANC_DEVICE    cuda | cpu           (default cuda)
#   ANC_WORK      staging root         (default Pipeline/AncestorCalibration)
#   ANC_PAIR      reference pair       (default 2n54B_2hdmA)
#   ANC_EXTRA     extra flags passed straight through
#
# SLURM splits --export on commas, so set these in the submitting shell and
# pass a bare --export=ALL rather than putting a comma list inside --export.
#
# CPU-only variant -- DDG and S4PRED both run fine on CPU for 25 short
# sequences, so the cheap half needs no GPU at all:
#
#   ANC_METHODS=ddg,s4pred ANC_DEVICE=cpu \
#   sbatch -p glacier --gres=NONE --cpus-per-task=8 --mem=16G --time=00:30:00 \
#          --export=ALL scripts/run_ancestor_calibration.sbatch
#
# AF2 (MSA mode and single-sequence mode) is 50 ColabFold runs, so give it
# its own job and its own output files, then merge the CSVs if you want one
# table. --n_controls 0 drops the negative control for the expensive method
# (the controls are what make it 50 runs instead of 10):
#
#   ANC_METHODS=af2,af2_single \
#   ANC_EXTRA='--n_controls 0 --out docs/ancestor_calibration_af2.csv --stats_out docs/ancestor_calibration_af2_stats.csv --fig docs/figures/global/ancestor_calibration_af2.png' \
#   sbatch --time=04:00:00 --export=ALL scripts/run_ancestor_calibration.sbatch
#
# BioEmu needs its own venv (see INSTALL_bioemu.md) and ~100 samples per
# sequence, so it is NOT in the default method list either:
#
#   ANC_METHODS=bioemu \
#   sbatch --time=08:00:00 --export=ALL scripts/run_ancestor_calibration.sbatch
#
# Run against a different working tree (e.g. a worktree) with
#   sbatch --chdir=<worktree> scripts/run_ancestor_calibration.sbatch
# The tree needs Pipeline/ (a symlink is fine) so that
# Pipeline/FoldPairs/2n54B_2hdmA/{chain_pdb_files,fasta_chain_files,output_get_msa}
# are readable. That pair dir is only ever READ -- the ancestors are staged
# into their own FoldPairs-shaped trees under $ANC_WORK.

set -e
source /sci/labs/orzuk/orzuk/venvs/my-python-venv/bin/activate
export MPLBACKEND=Agg QT_QPA_PLATFORM=offscreen PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4

# Tool locations (override in the environment if yours differ).
export THERMOMPNN_DIR="${THERMOMPNN_DIR:-/sci/labs/orzuk/orzuk/github/ThermoMPNN}"
export S4PRED_DIR="${S4PRED_DIR:-/sci/labs/orzuk/orzuk/github/s4pred}"
ESMFOLD2_VENV="${ESMFOLD2_VENV:-/sci/labs/orzuk/orzuk/venvs/esmfold2-venv}"
ESMFOLD2_HF_CACHE="${ESMFOLD2_HF_CACHE:-/sci/labs/orzuk/orzuk/tmp/xdg-cache}"

METHODS="${ANC_METHODS:-esmfold2,ddg,s4pred,boltz2}"
DEVICE="${ANC_DEVICE:-cuda}"
WORK="${ANC_WORK:-Pipeline/AncestorCalibration}"
PAIR="${ANC_PAIR:-2n54B_2hdmA}"

echo "host $(hostname)  cwd $(pwd)  start $(date)"
echo "methods=$METHODS  device=$DEVICE  pair=$PAIR  work=$WORK"
nvidia-smi || true

echo "=== unit tests (delta conventions, Spearman, exact 5! permutation p) ==="
python3 tests/test_ancestor_calibration.py

echo "=== planned commands (dry run) ==="
python3 scripts/ancestor_calibration.py --dry_run --methods "$METHODS" \
    --pair_id "$PAIR" --work_root "$WORK" --device "$DEVICE"

echo "=== calibration ==="
python3 scripts/ancestor_calibration.py \
    --methods "$METHODS" \
    --pair_id "$PAIR" \
    --work_root "$WORK" \
    --device "$DEVICE" \
    --n_controls "${ANC_N_CONTROLS:-20}" \
    --msa_top_n "${ANC_MSA_TOP_N:-10}" \
    --thermompnn_dir "$THERMOMPNN_DIR" \
    --s4pred_dir "$S4PRED_DIR" \
    --esmfold2_venv "$ESMFOLD2_VENV" \
    --esmfold2_hf_cache "$ESMFOLD2_HF_CACHE" \
    --out docs/ancestor_calibration.csv \
    --stats_out docs/ancestor_calibration_stats.csv \
    --fig docs/figures/global/ancestor_calibration.png \
    ${ANC_EXTRA:-}

echo "=== results ==="
column -s, -t < docs/ancestor_calibration_stats.csv || cat docs/ancestor_calibration_stats.csv

echo "=== DONE $(date) ==="
