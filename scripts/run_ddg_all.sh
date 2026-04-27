#!/bin/bash
# ============================================================================
# run_ddg_all.sh — submit ThermoMPNN DDG scoring jobs for fold-switching pairs.
#
# Submits one sbatch job per pair on the HUJI/Moriah cluster. Each job loads
# ThermoMPNN, computes (or reuses cached) per-backbone L x 20 DDG matrices,
# scores per-cluster homologs against both backbones, and writes:
#   Pipeline/FoldPairs/<pair>/output_ddg/ddg_matrix_<tag>.npz   (cached)
#   Pipeline/FoldPairs/<pair>/output_ddg/cluster_<NNN>.json     (per-seq detail)
#   Pipeline/FoldPairs/<pair>/Analysis/df_ddg.csv               (per-cluster summary)
#
# Usage
# -----
#   bash scripts/run_ddg_all.sh                 # all 93 pairs
#   bash scripts/run_ddg_all.sh 4dxrA_4dxtA     # one pair (smoke test)
#   bash scripts/run_ddg_all.sh "1dzlA_5keqF 4dxrA_4dxtA"   # multiple pairs
#
# After all jobs finish, fold DDG into the global tables:
#   python3 run_foldswitch_pipeline.py --run_mode postprocess \
#       --foldpair_ids ALL --run_job_mode inline
#
# ----------------------------------------------------------------------------
# Lessons learned during the debugging session that built this recipe
# (2026-04-26 / 2026-04-27):
#
# Environment:
#   - The `my-python-venv` venv has all DDG deps (tmtools, parasail, omegaconf,
#     pytorch-lightning, wandb, biopython, biotite, etc.). Newly added: those
#     last four were missing and pip-installed.
#   - The venvs were moved into `/sci/labs/orzuk/orzuk/venvs/`. The activate
#     scripts hardcoded the old `/sci/labs/orzuk/orzuk/<v>` path; they were
#     patched with sed to point at the new location. If activation prompts
#     appear but `which python3` still shows miniforge, the activate script
#     is pointing to a non-existent dir.
#   - Conda's auto-activated `(base)` prepends miniforge to PATH on shell
#     start. Most venvs win because their `bin/python3` is found first; if a
#     venv only has `python` (no `python3` symlink) you fall through to
#     conda. None of the project venvs have this issue today.
#
# SLURM partitions on Moriah (as of 2026-04-27):
#   - catfish  : 5 nodes, gpu:l4:8     (24 GB) — usually busy
#   - salmon   : 10 nodes, gpu:l40s:8  (48 GB) — best capacity, used here
#   - dogfish  : 2 nodes, gpu:a100:8   (80 GB) — was drained
#   - goldfish : 1 node,  gpu:h200:8   (141 GB) — newest, single-node
#   - puffin / a30 (referenced in old scripts) NO LONGER EXISTS.
#   - Each partition has a `<name>-k` preemptible variant for faster scheduling.
#
# Pipeline orchestrator quirks:
#   - `--run_job_mode sbatch` submits one outer sbatch per pair. The outer
#     sbatch flags come from `--sbatch_*` args at the top level. The
#     SBATCH_HINTS table in run_foldswitch_pipeline.py only flows into the
#     INNER inline rerun, NOT the outer sbatch — so always pass --sbatch_*.
#   - Bug fixed in commit 9ce6dd2: `_submit_pair_job` was clobbering --gres
#     for non-esmfold modes; now it only does so for non-GPU modes.
#   - Bug fixed in commit 65beca0: `seed_residue_indices` couldn't handle
#     deep MSAs with generic `>S1` / `>S2` headers; now accepts a fallback
#     positional index.
#
# Account / fair-share:
#   - The orzuk personal account has low FairShare (~0.04) due to high recent
#     usage; jobs queue for hours-to-days. The course-52017-25 account has
#     FairShare ~0.95 and gets jobs scheduled in seconds. Use it for bursty
#     submissions like this 93-pair batch.
#
# Time / resources:
#   - 00:30:00 wall is plenty (per-pair runtime: 10s with cached matrices, up
#     to a few minutes when computing both backbones from scratch). Lower
#     time requests pack 4x denser per node and reduce queue wait.
#   - 16 GB RAM, 2 CPUs, 1 GPU is sufficient for ThermoMPNN inference.
# ============================================================================

set -euo pipefail

# ---- arguments / defaults ---------------------------------------------------
PAIRS="${1:-ALL}"     # space-separated pair ids, or the literal "ALL"

# ---- paths ------------------------------------------------------------------
REPO_DIR="/sci/labs/orzuk/orzuk/github/MsaCluster"
VENV_DIR="/sci/labs/orzuk/orzuk/venvs/my-python-venv"
THERMOMPNN_DIR="/sci/labs/orzuk/orzuk/github/ThermoMPNN"

# ---- SLURM submission knobs (the proven combination) -----------------------
PARTITION="salmon"
GRES="gpu:l40s:1"
ACCOUNT="course-52017-25"
CPUS="2"
MEM="16G"
TIME="00:30:00"

# ---- run --------------------------------------------------------------------
cd "$REPO_DIR"
# shellcheck disable=SC1091
source "$VENV_DIR/bin/activate"

echo "[run_ddg_all] pairs:     $PAIRS"
echo "[run_ddg_all] partition: $PARTITION  gres=$GRES  account=$ACCOUNT"
echo "[run_ddg_all] resources: cpus=$CPUS  mem=$MEM  time=$TIME"
echo

# shellcheck disable=SC2086  # we want word-splitting on $PAIRS
python3 run_foldswitch_pipeline.py \
    --run_mode run_ddg \
    --foldpair_ids $PAIRS \
    --run_job_mode sbatch \
    --thermompnn_dir "$THERMOMPNN_DIR" \
    --sbatch_partition "$PARTITION" \
    --sbatch_gres "$GRES" \
    --sbatch_cpus "$CPUS" \
    --sbatch_mem "$MEM" \
    --sbatch_time "$TIME" \
    --sbatch_account "$ACCOUNT"

cat <<'EOF'

Monitor progress:
  squeue -u $USER -h | wc -l                                        # queued+running
  squeue -u $USER -t R -h | wc -l                                   # running
  ls Pipeline/FoldPairs/*/Analysis/df_ddg.csv 2>/dev/null | wc -l   # finished (target: 93)

Tail the most recent job log:
  ls -t Pipeline/FoldPairs/*/jobs/run_ddg_*.out | head -1 | xargs tail -30

After all 93 finish, fold DDG into global tables:
  python3 run_foldswitch_pipeline.py --run_mode postprocess \
      --foldpair_ids ALL --run_job_mode inline
EOF
