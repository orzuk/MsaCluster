#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10G


source /sci/labs/orzuk/steveabecassis/colabfold_new/bin/activate.csh
module load torch/1.3
mkdir /sci/labs/orzuk/steveabecassis/MsaCluster/Pipeline/output/esm_cmap_output
python3  ./runESM.py  -input_msas /sci/labs/orzuk/steveabecassis/MsaCluster/Pipeline/output/output_msa_cluster -o /sci/labs/orzuk/steveabecassis/MsaCluster/Pipeline/output/output_cmap_esm
