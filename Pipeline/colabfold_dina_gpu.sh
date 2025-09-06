#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10G
#SBATCH --gres=gpu:a100-1-10

module load cuda/11.1
module load cudnn/8.0.5

source /sci/labs/dina/dina/collabFold_phoenix/bin/activate.csh

mkdir ./output/af_preds
"$COLABFOLD_BIN"  ./output/output_msa_cluster ./output/af_preds/  ${CF_DATA_FLAGS[@]}
