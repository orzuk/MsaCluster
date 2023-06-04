#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=8
#SBATCH --mem=4G


python3  ./RunAF2_moriah.py --input_msa_folder ./msas
python3  ./RunAF2_moriah.py --input_msa_folder ./msas1

