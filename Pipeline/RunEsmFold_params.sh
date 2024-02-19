#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --ntasks=16
#SBATCH --mem=32G
#SBATCH --partition=puffin
#SBATCH --gres=gpu:a100:8

OUTPUT_NAME_DIR="$1"


source /sci/labs/orzuk/steveabecassis/colabfold_new/bin/activate.csh

module load torch/1.3
mkdir -p ./$OUTPUT_NAME_DIR/output_esm_fold

python3  ./ESMFoldHF.py  -input $OUTPUT_NAME_DIR



