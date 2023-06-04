#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=2
#SBATCH --mem=1G



python3  ./ClusterMSA_moriah.py --keyword test_pred -i ./AF_cluster_73815.a3m -o ./output 

