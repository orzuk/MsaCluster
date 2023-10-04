#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10G

OUTPUT_NAME_DIR="$1"

module load torch/1.3
mkdir -p ./Pipeline/$OUTPUT_NAME_DIR/esm_cmap_output
python3  ./runESM.py  --input_msas ./Pipeline/$OUTPUT_NAME_DIR/output_msa_cluster -o ./Pipeline/$OUTPUT_NAME_DIR/output_cmap_esm
