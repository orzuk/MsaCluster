#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10G
#SBATCH --gres=gpu:a100-1-10

OUTPUT_NAME_DIR="$1"


. /sci/labs/orzuk/orzuk/my-python-venv/bin/activate

source /sci/labs/orzuk/steveabecassis/colabfold_new/bin/activate.csh
. /sci/labs/orzuk/steveabecassis/colabfold_new/bin/activate
module load torch/1.3
mkdir -p ./Pipeline/$OUTPUT_NAME_DIR/esm_cmap_output
python3  ./runESM.py  --input_msas ./Pipeline/$OUTPUT_NAME_DIR/output_msa_cluster -o ./Pipeline/$OUTPUT_NAME_DIR/output_cmap_esm