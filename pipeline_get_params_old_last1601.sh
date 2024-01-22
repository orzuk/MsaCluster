#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10G

#FASTA_FILE_INPUT="$1"
#OUTPUT_NAME_DIR="$2"

FOLD1 = "$1"
FOLD2 = "$2"

mkdir Pipeline
mkdir Pipeline/$FOLD1_$FOLD2
mkdir Pipeline/$FOLD1_$FOLD2/pdb_file_path
mkdir Pipeline/$FOLD1_$FOLD2/chain_pdb_file_path
mdir  Pipeline/$FOLD1_$FOLD2/fasta_files
mkdir Pipeline/$FOLD1_$FOLD2/org_cmaps
mkdir Pipeline/$FOLD1_$FOLD2/output_get_msa
python3 ./get_msa.py ./Pipeline/fasta_files/$FOLD1.fasta ./Pipeline/$FOLD1_$FOLD2/output_get_msa  -name 'DeepMsa'
mkdir Pipeline/$FOLD1_$FOLD2/output_msa_cluster
python3  ./ClusterMSA_moriah.py --keyword ShallowMsa -i ./Pipeline/$FOLD1_$FOLD2/output_get_msa/DeepMsa.a3m  -o ./Pipeline/$FOLD1_$FOLD2/output_msa_cluster
mkdir ./Pipeline/$FOLD1_$FOLD2/AF_preds
sbatch ./Pipeline/RunAF_params.sh $FOLD1_$FOLD2
mkdir Pipeline/$FOLD1_$FOLD2/esm_fold_output
sbatch /sci/labs/orzuk/steveabecassis/MsaCluster/Pipeline/RunEsmFold_params.sh
mkdir  ./Pipeline/$FOLD1_$FOLD2/output_cmap_esm
sbatch ./CmapESM_params.sh $FOLD1_$FOLD2





