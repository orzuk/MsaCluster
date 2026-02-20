#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=1
#SBATCH --mem=2G

./sci/labs/orzuk/steveabecassis/colabfold_new/bin/activate.csh

python3 ./test.py
