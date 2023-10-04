#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10G

. /sci/labs/orzuk/steveabecassis/colabfold_new/bin/activate

python3 ./get_msa.py ./2QKEE_002.a3m    ./output/output_get_msa  -name 'test'
