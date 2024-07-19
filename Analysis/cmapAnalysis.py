from tqdm import tqdm
import os
import mdtraj as md
from scipy.ndimage import binary_dilation
from Bio.PDB import PDBParser
from Bio import PDB
from Bio import pairwise2
import pandas as pd
import numpy as np

def extract_protein_sequence(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure("protein_structure", pdb_file)

    residue_sequence = ""
    flag = 0
    # Iterate through the structure and extract the residue sequence
    for model in structure:
        for chain in model:
            if flag == 1:
                break
            flag = 1
            for residue in chain:
                if PDB.is_aa(residue):
                    residue_sequence += PDB.Polypeptide.three_to_one(residue.get_resname())

    return residue_sequence


def get_align_indexes(seqA, seqB):
    alignments = pairwise2.align.globalxx(seqA, seqB, one_alignment_only=True)
    best_alignment = alignments[0]
    aligned_seq1, aligned_seq2, score, begin, end = best_alignment
    index_pairs = []
    seq1_index = seq2_index = 0

    for a, b in zip(aligned_seq1, aligned_seq2):
        if a != '-' and b != '-':
            index_pairs.append((seq1_index, seq2_index))
            seq1_index += 1
            seq2_index += 1
        elif a != '-':
            seq1_index += 1
        elif b != '-':
            seq2_index += 1

    idx_seqA = [i[0] for i in index_pairs]
    idx_seqB = [i[-1] for i in index_pairs]
    return idx_seqA, idx_seqB


def get_contact_map_from_pdb(pdb_file, size=None, start_pos=0, end_pos=-1):
    pdb_obj = md.load_pdb(pdb_file)
    distances, pairs = md.compute_contacts(pdb_obj)
    contacts = md.geometry.squareform(distances, pairs)[0]
    if size is None:
        size = contacts.shape[0]
    arr = np.zeros([size, size])
    arr[np.where(contacts[start_pos:end_pos, start_pos:end_pos] < 0.5)] = 1
    return arr


def dilate_with_tolerance(array, tolerance):
    structure = np.ones((2 * tolerance + 1, 2 * tolerance + 1), dtype=int)
    return binary_dilation(array, structure=structure).astype(int)



if __name__ == '__main__':
    folder    = '/Users/steveabecassis/Desktop/Pipeline'
    tolerance = 2
    for fold_pair in tqdm(os.listdir(folder)):
        print(fold_pair)
        fold_pair = '1qlnA_1h38D'
        try:
            fold1     = fold_pair.split('_')[0]
            fold2     = fold_pair.split('_')[-1]
            seq_fold1 = extract_protein_sequence(f'{folder}/{fold_pair}/chain_pdb_files/{fold1}.pdb')
            seq_fold2 = extract_protein_sequence(f'{folder}/{fold_pair}/chain_pdb_files/{fold2}.pdb')

            fold1_idxs ,fold2_idxs = get_align_indexes(seq_fold1, seq_fold2)
            cmap_pdb1             = get_contact_map_from_pdb(f'{folder}/{fold_pair}/chain_pdb_files/{fold1}.pdb')[fold1_idxs][:, fold1_idxs]
            cmap_pdb2             = get_contact_map_from_pdb(f'{folder}/{fold_pair}/chain_pdb_files/{fold2}.pdb')[fold2_idxs][:, fold2_idxs]

            dilated_cmap1 = dilate_with_tolerance(cmap_pdb1, tolerance)
            dilated_cmap2 = dilate_with_tolerance(cmap_pdb2, tolerance)

            overlap = np.multiply(dilated_cmap1, cmap_pdb2).astype(int)
            overlap_dilated = dilate_with_tolerance(overlap, tolerance)
            only_cmap1      = np.clip(cmap_pdb1 - overlap_dilated,0,1).astype(int)
            only_cmap2 = np.clip(cmap_pdb2 - overlap_dilated, 0, 1).astype(int)

            cmap_pred_files = os.listdir(f'{folder}/{fold_pair}/output_cmap_esm')
            res = []
            for cmap_pred_file in cmap_pred_files:
                if 'Shall' in cmap_pred_file:
                    cmap_pred                   = np.load(f'{folder}/{fold_pair}/output_cmap_esm/{cmap_pred_file}')[fold1_idxs][:, fold1_idxs]
                    cmap_pred[cmap_pred > 0.4]  = 1
                    cmap_pred[cmap_pred <= 0.4] = 0
                    cmap_pred_without_overlap = np.clip(cmap_pred - overlap_dilated, 0, 1).astype(int)

                    cmap_pred_dilated  = dilate_with_tolerance(cmap_pred, tolerance)
                    only_cmap1_dilated = dilate_with_tolerance(only_cmap1, tolerance)
                    only_cmap2_dilated = dilate_with_tolerance(only_cmap2, tolerance)

                    cmap_pred_without_overlap_dilated = dilate_with_tolerance(cmap_pred_without_overlap, tolerance)
                    cmap_pred_1_true = np.multiply(cmap_pred_without_overlap_dilated,only_cmap1).astype(int)
                    precision_1      = np.sum(np.multiply(cmap_pred_without_overlap,cmap_pred_1_true))/np.sum(cmap_pred_without_overlap)
                    recall_1         = np.sum(cmap_pred_1_true)/np.sum(only_cmap1)

                    cmap_pred_2_true = np.multiply(cmap_pred_without_overlap_dilated,only_cmap2).astype(int)
                    precision_2      = np.sum(np.multiply(cmap_pred_without_overlap,cmap_pred_2_true))/np.sum(cmap_pred_without_overlap)
                    recall_2         = np.sum(cmap_pred_2_true)/np.sum(only_cmap2)

                    cluster_num      = cmap_pred_file[-7:-4]
                    res.append({'fileName' :cmap_pred_file[:-4],'cluster_num' :cluster_num,'precision_1' :precision_1,'recall_1' :recall_1 ,
                                'precision_2' :precision_2 ,'recall_2' :recall_2})

            temp_df = pd.DataFrame(res)
            temp_df.to_csv(f'{folder}/{fold_pair}/Analysis/new_cmap_res.csv' ,index=False)
        except Exception as e:
            print(e)
            continue








