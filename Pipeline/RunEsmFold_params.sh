#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --ntasks=8
#SBATCH --mem=24G




#. /sci/labs/orzuk/steveabecassis/colabfold_new/bin/activate
#source /sci/labs/orzuk/steveabecassis/colabfold_new/bin/activate.csh
#module load torch/1.3
#mkdir -p ./Pipeline/2pbkB_3njqA/output_esm_fold
#export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:256
python3  ./ESMFoldHF.py  -input 2pbkB_3njqA


