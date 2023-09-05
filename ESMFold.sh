#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10G

mkdir output/esm_fold_output
python3 ./ESMFoldHF.py -input ./2QKEE_002.a3m  -output ./output/esm_fold_output/ -name 'test'
