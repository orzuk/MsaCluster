#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=24
#SBATCH --mem=50G


mkdir output/esm_fold_output
python3 ./ESMFoldHF.py -input ./msas4esm/cluster_008.a3m -output ./output/esm_fold_output -name 'Cluster008_'