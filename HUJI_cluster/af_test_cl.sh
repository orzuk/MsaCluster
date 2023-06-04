#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4G


python3  ./RunAF2_moriah.py --input_msa ./ClusterMSA/output/test_pred_002.a3m


