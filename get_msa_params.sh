#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10G

arg1="$1"
arg2="$2"


mkdir $arg1
python3 ./get_msa.py  $arg2   ./$arg1  -name 'test'
