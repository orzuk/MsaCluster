#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10G


source /sci/labs/orzuk/steveabecassis/colabfold_new/bin/activate.csh
module load torch/1.3
mkdir output/esm_cmap_output
python3  ./runESM.py  -i ./cluster_msa_output_1jfk/ -o ./output/esm_cmap_output
