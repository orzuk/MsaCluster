#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10G


FASTA_FILE_INPUT="$1"
OUTPUT_NAME_DIR="$2"

. /sci/labs/orzuk/orzuk/my-python-venv/bin/activate

mkdir -p Pipeline/$OUTPUT_NAME_DIR/output_get_msa
python3 ./get_msa.py $FASTA_FILE_INPUT ./Pipeline/$OUTPUT_NAME_DIR/output_get_msa  -name 'DeepMsa'

