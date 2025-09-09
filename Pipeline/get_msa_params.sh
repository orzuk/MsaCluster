#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10G


FASTA_FILE_INPUT="$1"
OUTPUT_NAME_DIR="$2"

. /sci/labs/orzuk/orzuk/my-python-venv/bin/activate

mkdir -p $OUTPUT_NAME_DIR/output_get_msa
rm -r $OUTPUT_NAME_DIR/output_get_msa/DeepMsa_env
python3 ./get_msa.py "Pipeline/<pair>/_seed_both.a3m" "Pipeline/<pair>/output_get_msa" --name DeepMsa

