#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --ntasks=8
#SBATCH --mem=24G


OUTPUT_NAME_DIR="$1"

./sci/labs/orzuk/steveabecassis/colabfold_new/bin/activate.csh

#module load torch/1.3
mkdir -p ./Pipeline/$OUTPUT_NAME_DIR/output_esm_fold
export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:256
python3  ./ESMFoldHF.py  -input $OUTPUT_NAME_DIR



