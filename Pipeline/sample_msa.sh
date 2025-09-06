#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10G


python3 ./get_sample_msa_algn.py ./output/output_get_msa  ./output/output_sample_msa/
