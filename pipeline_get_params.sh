#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10G

FASTA_FILE_INPUT="$1"
OUTPUT_NAME_DIR="$2"

mkdir -p Pipeline
mkdir -p Pipeline/$OUTPUT_NAME_DIR
mkdir -p Pipeline/$OUTPUT_NAME_DIR/output_get_msa
python3 ./get_msa.py $FASTA_FILE_INPUT ./Pipeline/$OUTPUT_NAME_DIR/output_get_msa  -name 'DeepMsa'
mkdir -p Pipeline/$OUTPUT_NAME_DIR/output_msa_cluster
python3  ./ClusterMSA_moriah.py --keyword ShallowMsa -i ./Pipeline/$OUTPUT_NAME_DIR/output_get_msa/DeepMsa.a3m  -o ./Pipeline/$OUTPUT_NAME_DIR/output_msa_cluster
mkdir -p ./Pipeline/$OUTPUT_NAME_DIR/AF_preds
sbatch -o ./Pipeline/$OUTPUT_NAME_DIR/RunAF.out  ./Pipeline/RunAF_params.sh $OUTPUT_NAME_DIR
#mkdir -p ./Pipeline/output/esm_fold_output
#sbatch -o ./Pipeline/$OUTPUT_NAME_DIR/ESMFold.out ./Pipeline/RunEsmFold.sh
mkdir -p ./Pipeline/$OUTPUT_NAME_DIR/output_cmap_esm
sbatch -o ./Pipeline/$OUTPUT_NAME_DIR/CmapESM.out  ./CmapESM_params.sh $OUTPUT_NAME_DIR
