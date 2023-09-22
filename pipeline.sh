#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10G

mkdir Pipeline
mkdir Pipeline/output
mkdir Pipeline/output/output_get_msa
python3 /sci/labs/orzuk/steveabecassis/MsaCluster/get_msa.py ./input/FASTA_FILE   /sci/labs/orzuk/steveabecassis/MsaCluster/Pipeline/output/output_get_msa  -name 'DeepMsa'
mkdir Pipeline/output/output_msa_cluster
python3  /sci/labs/orzuk/steveabecassis/MsaCluster/ClusterMSA_moriah.py --keyword ShallowMsa -i /sci/labs/orzuk/steveabecassis/MsaCluster/Pipeline/output/output_get_msa/DeepMsa.a3m  -o /sci/labs/orzuk/steveabecassis/MsaCluster/Pipeline/output/output_msa_cluster
mkdir /sci/labs/orzuk/steveabecassis/MsaCluster/Pipeline/output/AF_preds
sbatch /sci/labs/orzuk/steveabecassis/MsaCluster/Pipeline/RunAF.sh
mkdir /sci/labs/orzuk/steveabecassis/MsaCluster/Pipeline/output/esm_fold_output
sbatch /sci/labs/orzuk/steveabecassis/MsaCluster/Pipeline/RunEsmFold.sh
mkdir  /sci/labs/orzuk/steveabecassis/MsaCluster/Pipeline/output/output_cmap_esm
sbatch /sci/labs/orzuk/steveabecassis/MsaCluster/CmapESM.sh



