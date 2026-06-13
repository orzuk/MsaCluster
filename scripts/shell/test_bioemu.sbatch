#!/bin/bash
#SBATCH --job-name=bioemu_test
#SBATCH --gres=gpu:1
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --output=bioemu_test_%j.out
# If your GPU queue needs an explicit partition/account, add e.g.:
#   #SBATCH --partition=gpu

set -e
. /sci/labs/orzuk/orzuk/venvs/bioemu-venv/bin/activate
export BIOEMU_CACHE=/sci/labs/orzuk/orzuk/bioemu_cache
export HF_HOME=/sci/labs/orzuk/orzuk/hf_cache
export BIOEMU_SO3_CACHE="$BIOEMU_CACHE/so3"
mkdir -p "$BIOEMU_CACHE" "$HF_HOME" "$BIOEMU_SO3_CACHE"

# Home has ~2 GB quota; keep ALL caches in lab storage. bioemu/colabfold
# downloads AF2 params (~3.5 GB) to ~/.cache/colabfold, which is a symlink to
# lab -- ensure the symlink target exists so the download doesn't fail on a
# dangling link. (Writes go to lab via the symlink; nothing lands in home.)
CF_DIR="$(readlink -f "$HOME/.cache/colabfold" 2>/dev/null || true)"
[ -n "$CF_DIR" ] && mkdir -p "$CF_DIR/params"

# jaxlib + tensorflow both register xla_data.proto into protobuf's C++ pool;
# jaxlib aborts when it is the SECOND registrant. sitecustomize.py in jaxfirst/
# forces jax to init first (see that file). Put it at the FRONT of PYTHONPATH.
export PYTHONPATH="${SLURM_SUBMIT_DIR:-$PWD}/scripts/shell/jaxfirst${PYTHONPATH:+:$PYTHONPATH}"
# Quiet the TF/absl + httpx/httpcore DEBUG firehose.
export TF_CPP_MIN_LOG_LEVEL=3

echo "=== GPU check ==="
python -c "import torch; print('cuda', torch.cuda.is_available(), torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'NO GPU')"

echo "=== bioemu.sample signature (the real CLI flag names) ==="
python -c "import bioemu.sample as s, inspect; fn=getattr(s,'main',None) or getattr(s,'sample',None); print(fn.__name__, inspect.signature(fn))"

echo "=== tiny standalone sample (chignolin, 10 samples) ==="
OUT=/sci/labs/orzuk/orzuk/bioemu_test_out
rm -rf "$OUT"; mkdir -p "$OUT"
python -m bioemu.sample --sequence GYDPETGTWG --num_samples 10 --output_dir "$OUT" \
    --cache_embeds_dir "$BIOEMU_CACHE" --cache_so3_dir "$BIOEMU_SO3_CACHE"

echo "=== output files ==="
ls -la "$OUT"
echo "DONE"
