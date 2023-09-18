#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10G
#SBATCH --gres=gpu:a100-3-40


mkdir output/esm_fold_output
python3 ./ESMFoldHF.py -input ./msas4esm/cluster_008.a3m -output ./output/esm_fold_output -name 'Cluster008_'
