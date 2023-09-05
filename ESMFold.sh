#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10G

mkdir output/esm_fold_output
python3 ./get_sample_msa_algn.py ./output/output_get_msa  ./output/esm_fold_output/ -name 'test'
