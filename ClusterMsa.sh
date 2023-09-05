#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10G


source /sci/labs/orzuk/steveabecassis/colabfold_new/bin/activate.csh
python3  ./ClusterMSA_moriah.py --keyword MSA -i ./output/output_get_msa/test.a3m  -o ./output/output_msa_cluster
