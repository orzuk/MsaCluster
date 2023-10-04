#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10G
. /sci/labs/orzuk/orzuk/my-python-venv/bin/activate
 pip install -r requirements.txt
