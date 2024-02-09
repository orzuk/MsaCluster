#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=16
#SBATCH --mem=256G
#SBATCH --partition=puffin
#SBATCH --gres=gpu:a100-3-40

OUTPUT_NAME_DIR="$1"


source /sci/labs/orzuk/steveabecassis/colabfold_new/bin/activate.csh

module load torch/1.3
mkdir -p ./Pipeline/$OUTPUT_NAME_DIR/output_esm_fold

python3  ./ESMFoldHF.py  -input $OUTPUT_NAME_DIR



