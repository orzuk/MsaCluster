#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10G
#SBATCH --gres=gpu:a100-1-10

OUTPUT_NAME_DIR="$1"


module load torch/1.3
module load cuda/11.1
module load cudnn/8.0.5


source /sci/labs/dina/dina/collabFold_phoenix/bin/activate.csh
/sci/labs/dina/dina/localcolabfold/colabfold-conda/bin/colabfold_batch  ./Pipeline/$OUTPUT_NAME_DIR/output_msa_cluster ./Pipeline/$OUTPUT_NAME_DIR/AF_preds/  --data=/cs/labs/dina/seanco/colabfold/weights/











