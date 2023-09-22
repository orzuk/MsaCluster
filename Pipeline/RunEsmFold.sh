#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=16
#SBATCH --mem=256G

source /sci/labs/orzuk/steveabecassis/colabfold_new/bin/activate.csh
mkdir output/esm_fold_output
python3 /sci/labs/orzuk/steveabecassis/MsaCluster/ESMFoldHF.py -input /sci/labs/orzuk/steveabecassis/MsaCluster/Pipeline/output/output_msa_cluster -output /sci/labs/orzuk/steveabecassis/MsaCluster/Pipeline/output/esm_fold_output -name 'Cluster'

