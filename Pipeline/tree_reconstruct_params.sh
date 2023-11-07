#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=8
#SBATCH --mem=32G

#SBATCH --partition=puffin
#SBATCH --gres=gpu:a30


OUTPUT_NAME_DIR="$1"


. /sci/labs/orzuk/orzuk/my-python-venv/bin/activate

module load torch/1.3
mkdir -p ./Pipeline/$OUTPUT_NAME_DIR/output_phytree
python3  ./run_tree_reconstruct.py  --method distance  --input_msas ./Pipeline/$OUTPUT_NAME_DIR/output_get_msa/DeepMsa.a3m -o ./Pipeline/$OUTPUT_NAME_DIR/output_phytree
