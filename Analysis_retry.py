from argparse import  ArgumentParser
from Bio.PDB import PDBParser, PDBIO, Select
from contact_map import ContactFrequency, ContactDifference
import numpy as np
import mdtraj as md
import os
import pandas as pd
from scipy.stats import spearmanr
import gzip
import subprocess
import re
from Bio import pairwise2
from Bio import PDB
from tqdm import tqdm

def save_org_cmaps(path,fold):
    traj = md.load(filename_or_filenames=f'{path}/chain_pdb_files/{fold}.pdb')
    frame_contacts = ContactFrequency(traj[0])
    cmap = frame_contacts.residue_contacts.df.to_numpy()
    np.save(f'{path}/cmaps_pairs/{fold}.npy',cmap)
    return frame_contacts

def get_tmscore_align(path_fold1,path_fold2):
    command = f"./TMalign {path_fold1} {path_fold2}"
    output = subprocess.check_output(command, shell=True)
    match = re.search(r"TM-score=\s+(\d+\.\d+)", str(output))
    if match:
        result = match.group(1)
        return float(result)
    else:
        return None

def get_jaccard_index_score(cmap1,cmap2,th=0.4):
    # Calculate intersection and union
    intersection = np.logical_and((cmap1> th).astype(int), (cmap2 > th).astype(int))
    union = np.logical_or((cmap1> th).astype(int), (cmap2 > th).astype(int))
    # Calculate the Jaccard Index
    jaccard_index = round(np.sum(intersection) / np.sum(union),3)
    print(f"Jaccard Index: {jaccard_index}")
    return jaccard_index
def spearman_score(cmap1,cmap2):
    # Calculate Spearman Correlation
    correlation,_ = spearmanr(cmap1.flatten(), cmap2.flatten())
    print(f"Spearman Correlation Coefficient: {correlation}")
    return correlation

def get_TP_TN_FP_FN_metrics(cmap1,cmap2,th=0.4):
    cmap1 = (cmap1 > th).astype(int)
    cmap2 = (cmap2 > th).astype(int)

    TP = np.sum((cmap1 == 1) & (cmap2 == 1))
    TN = np.sum((cmap1 == 0) & (cmap2 == 0))
    FP = np.sum((cmap1 == 0) & (cmap2 == 1))
    FN = np.sum((cmap1 == 1) & (cmap2 == 0))
    print(f"True Positives (TP): {TP}")
    print(f"True Negatives (TN): {TN}")
    print(f"False Positives (FP): {FP}")
    print(f"False Negatives (FN): {FN}")

    size = max(cmap1.size,cmap2.size)
    TP_norm = np.sum((cmap1 == 1) & (cmap2 == 1))/size
    TN_norm = np.sum((cmap1 == 0) & (cmap2 == 0))/size
    FP_norm = np.sum((cmap1 == 0) & (cmap2 == 1))/size
    FN_norm = np.sum((cmap1 == 1) & (cmap2 == 0))/size
    print(f"True Positives Normalized  (TP): {TP_norm}")
    print(f"True Negatives Normalized  (TN): {TN_norm}")
    print(f"False Positives Normalized (FP): {FP_norm}")
    print(f"False Negatives Normalized (FN): {FN_norm}")
    return TP,TN,FP,FN,TP_norm,TN_norm,FP_norm,FN_norm

def unzip_files(folder_path):
    files_in_folder = os.listdir(folder_path)
    # Iterate through the files and unzip .gz files
    for file_name in files_in_folder:
        file_path = os.path.join(folder_path, file_name)
        # Check if the file is a .gz file
        if file_name.endswith('.gz'):
            # Create the output file name by removing the '.gz' extension
            output_file_name = os.path.splitext(file_name)[0]
            output_file_path = os.path.join(folder_path, output_file_name)

            # Open and unzip the .gz file
            with gzip.open(file_path, 'rb') as gz_file, open(output_file_path, 'wb') as output_file:
                output_file.write(gz_file.read())

            print(f"Unzipped {file_name} to {output_file_name}")




def get_align_indexes(seqA,seqB):
    alignments = pairwise2.align.globalxx(seqA,seqB,one_alignment_only=True)
    best_align = alignments[0]
    seqA = best_align.seqA
    seqB = best_align.seqB
    cursA = 0
    cursB = 0
    seqA_idxs = []
    seqB_idxs = []
    for aa in range(len(seqA)):
        if (seqA[aa] != '-')&(seqB[aa] != '-'):
            seqA_idxs.append(cursA)
            seqB_idxs.append(cursB)
            cursA +=1
            cursB +=1
        if (seqA[aa] == '-')&(seqB[aa] != '-'):
            cursB +=1
        if (seqA[aa] != '-')&(seqB[aa] == '-'):
            cursA +=1
    return seqA_idxs,seqB_idxs

def extract_protein_sequence(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure("protein_structure", pdb_file)

    residue_sequence = ""

    # Iterate through the structure and extract the residue sequence
    for model in structure:
        for chain in model:
            for residue in chain:
                if PDB.is_aa(residue):
                    residue_sequence += PDB.Polypeptide.three_to_one(residue.get_resname())

    return residue_sequence

if __name__ == '__main__':

    '''
    Cmap analysis
    '''
    for fold_pair in tqdm(os.listdir('/Users/steveabecassis/Desktop/Pipeline')):
        try:
            if('df_af.csv' in os.listdir(f'/Users/steveabecassis/Desktop/Pipeline/{fold_pair}/Analysis')):
                continue
            if fold_pair == '4cmqB_4zt0C':
                continue
            print(f'Fold_pair: {fold_pair}')
            path = f'/Users/steveabecassis/Desktop/Pipeline/{fold_pair}'
            print(f'Path: {path}')
            if not os.path.exists(f'{path}/cmaps_pairs'):
                print("Mkdir: " + f'{path}/cmaps_pairs')
                os.mkdir(f'{path}/cmaps_pairs')
            if not os.path.exists(f'{path}/Analysis'):
                print("Mkdir: " + f'{path}/Analysis')
                os.mkdir(f'{path}/Analysis')

            folds = fold_pair.split("_")
            fold1 = folds[0]
            fold2 = folds[1]
            save_org_cmaps(f'{path}', fold1)
            save_org_cmaps(f'{path}', fold2)

            unzip_files(folder_path=f'{path}/AF_preds')

            pdb_files = os.listdir(f'{path}/AF_preds')
            res = []
            pdb_files = [i for i in pdb_files if 'pdb' in str(i)]
            pdb_fold1_path = f'{path}/chain_pdb_files/{fold1}.pdb'
            pdb_fold2_path = f'{path}/chain_pdb_files/{fold2}.pdb'
            for pdb_file in pdb_files:
                try:
                    if ('pdb' in str(pdb_file)) & ('gz' not in str(pdb_file)):
                        # score_pdb1 = tmscoring.get_tm(f'{path}/AF_preds/{pdb_file}', pdb_fold1_path)
                        # score_pdb2 = tmscoring.get_tm(f'{path}/AF_preds/{pdb_file}', pdb_fold2_path)

                        score_pdb1 = get_tmscore_align(f'{path}/AF_preds/{pdb_file}', pdb_fold1_path)
                        score_pdb2 = get_tmscore_align(f'{path}/AF_preds/{pdb_file}', pdb_fold2_path)

                        temp = {'pdb_file':pdb_file,'score_pdb1':score_pdb1,'score_pdb2':score_pdb2}
                        res.append(temp)
                except:
                    continue
            df_af = pd.DataFrame(res)
            df_af.loc[df_af.score_pdb1 > df_af.score_pdb2, 'Fold'] = fold1[:-1]
            df_af.loc[df_af.score_pdb1 < df_af.score_pdb2, 'Fold'] = fold2[:-1]
            df_af.sort_values(by='pdb_file', inplace=True)
            df_af['cluster_num'] = df_af['pdb_file'].apply(lambda x: x[11:14]).replace('ela','Query')
            df_af['is_fold1'] = (df_af['Fold'] == fold1[:-1]).astype(int)
            df_af['unique_fold_score'] = df_af.groupby('cluster_num')['is_fold1'].transform('mean')
            df_af.drop(columns=['is_fold1'],inplace=True)
            df_af.to_csv(f'{path}/Analysis/df_af.csv')
            print('Finish !')
        except:
            continue












