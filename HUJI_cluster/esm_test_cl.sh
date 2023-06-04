#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=2
#SBATCH --mem=1G

module load torch/1.3
python3  ./runESM.py --input_msa ./msas/2QKEE_002.a3m -o ./esm_output

