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
mkdir -p "$BIOEMU_CACHE" "$HF_HOME"

echo "=== GPU check ==="
python -c "import torch; print('cuda', torch.cuda.is_available(), torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'NO GPU')"

echo "=== bioemu.sample signature (the real CLI flag names) ==="
python -c "import bioemu.sample as s, inspect; fn=getattr(s,'main',None) or getattr(s,'sample',None); print(fn.__name__, inspect.signature(fn))"

echo "=== tiny standalone sample (chignolin, 10 samples) ==="
OUT=/sci/labs/orzuk/orzuk/bioemu_test_out
rm -rf "$OUT"; mkdir -p "$OUT"
python -m bioemu.sample --sequence GYDPETGTWG --num_samples 10 --output_dir "$OUT"

echo "=== output files ==="
ls -la "$OUT"
echo "DONE"
