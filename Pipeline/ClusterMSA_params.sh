#!/bin/bash

#SBATCH --time=20:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10G

OUTPUT_NAME_DIR="$1"

. /sci/labs/orzuk/orzuk/my-python-venv/bin/activate
# source /sci/labs/orzuk/steveabecassis/colabfold_new/bin/activate.csh
mkdir -p $OUTPUT_NAME_DIR/output_msa_cluster
python3  ./ClusterMSA_moriah.py --keyword ShallowMsa -i $OUTPUT_NAME_DIR/output_get_msa/DeepMsa.a3m  -o  $OUTPUT_NAME_DIR/output_msa_cluster
