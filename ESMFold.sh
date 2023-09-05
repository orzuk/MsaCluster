#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10G

mkdir output/esm_fold_output
python3 ./ESMFoldHF.py ./2QKEE_002.a3m  ./output/esm_fold_output/ -name 'test'
