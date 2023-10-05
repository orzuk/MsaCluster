#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10G

module load torch/1.3

python3 run_foldswitch_pipeline.py
