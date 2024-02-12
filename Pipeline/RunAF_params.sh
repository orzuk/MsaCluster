#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10G
#SBATCH --gres=gpu:a100-1-10

OUTPUT_NAME_DIR="$1"


module load torch/1.3
module load cuda/11.1
module load cudnn/8.0.5


source /sci/labs/dina/dina/collabFold_phoenix/bin/activate.csh
# Run full alignment and also each cluster separately!!! 
/sci/labs/dina/dina/localcolabfold/colabfold-conda/bin/colabfold_batch  $OUTPUT_NAME_DIR/output_get_msa  $OUTPUT_NAME_DIR/AF_preds/  --data=/cs/labs/dina/seanco/colabfold/weights/
/sci/labs/dina/dina/localcolabfold/colabfold-conda/bin/colabfold_batch  $OUTPUT_NAME_DIR/output_msa_cluster $OUTPUT_NAME_DIR/AF_preds/  --data=/cs/labs/dina/seanco/colabfold/weights/

# Zip output json+pdb files to save space, also convert png to jpg 
gzip $OUTPUT_NAME_DIR/AF_preds/*.json
gzip $OUTPUT_NAME_DIR/AF_preds/*.pdb

ls -1 $OUTPUT_NAME_DIR/AF_preds/*.png | xargs -n 1 bash -c 'convert "$0" "${0%.*}.jpg"'
rm  $OUTPUT_NAME_DIR/AF_preds/*.png







