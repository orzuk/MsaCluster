#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10G
#SBATCH --gres=gpu:a100-1-10

module load cuda/11.1
module load cudnn/8.0.5

source /sci/labs/dina/dina/collabFold_phoenix/bin/activate.csh
/sci/labs/dina/dina/localcolabfold/colabfold-conda/bin/colabfold_batch  /sci/labs/orzuk/steveabecassis/MsaCluster/Pipeline/output/output_msa_cluster /sci/labs/orzuk/steveabecassis/MsaCluster/Pipeline/output/AF_preds/  --data=/cs/labs/dina/seanco/colabfold/weights/















