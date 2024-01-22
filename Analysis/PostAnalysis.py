import os

import numpy as np
import pandas as pd
from contact_map import ContactFrequency, ContactDifference
from scipy.stats import spearmanr


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

if __name__ == '__main__':
    path = '/Users/steveabecassis/Desktop/pdb_file'
    fold_1 = '1jfk.npy'
    fold_2 = '1jfk.npy'
    cmap_fold_1 = np.loadtxt(f'{path}/{fold_1}')
    cmap_fold_2 = np.loadtxt(f'{path}/{fold_2}')

    cmap_res_analysis = []
    for cmap_cluster in os.listdir(f'{path}/esm_cmap_output'):
        if ('npy' in str(cmap_cluster)):
            file_name = str(cmap_cluster[:-4])
            cmap_cluster_loaded = np.loadtxt(f'{path}/esm_cmap_output/{cmap_cluster}')
            jaccard_index_score_fold1 = get_jaccard_index_score(cmap_fold_1,cmap_cluster_loaded)
            jaccard_index_score_fold2 = get_jaccard_index_score(cmap_fold_2,cmap_cluster_loaded)
            spearmanr_score_fold1 = spearman_score(cmap_fold_1, cmap_cluster_loaded)
            spearmanr_score_fold2 = spearman_score(cmap_fold_2, cmap_cluster_loaded)
            TP_f1,TN_f1,FP_f1,FN_f1,TP_norm_f1,TN_norm_f1,FP_norm_f1,FN_norm_f1 = get_TP_TN_FP_FN_metrics(cmap_fold_1,cmap_cluster_loaded,th=0.4)
            TP_f2,TN_f2,FP_f2,FN_f2,TP_norm_f2,TN_norm_f2,FP_norm_f2,FN_norm_f2 = get_TP_TN_FP_FN_metrics(cmap_fold_2, cmap_cluster_loaded, th=0.4)

            scores_dict = {'file_name':file_name,'jaccard_index_score_fold1':jaccard_index_score_fold1,'jaccard_index_score_fold2':jaccard_index_score_fold2,
                           'spearmanr_score_fold1':spearmanr_score_fold1,'spearmanr_score_fold2':spearmanr_score_fold2,
                           'TP_f1':TP_f1,'TN_f1':TN_f1,'FP_f1':FP_f1,'FN_f1':FN_f1,'TP_norm_f1':TP_norm_f1,'TN_norm_f1':TN_norm_f1,'FP_norm_f1':FP_norm_f1,'FN_norm_f1':FN_norm_f1,
            'TP_fold2': TP_f2, 'TN_f2': TN_f2, 'FP_f2':FP_f2, 'FN_f2': FN_f2, 'TP_norm_f2': TP_norm_f2, 'TN_norm_f2': TN_norm_f2, 'FP_norm_f2': FP_norm_f2, 'FN_norm_f2': FN_norm_f2
            }
            cmap_res_analysis.append(scores_dict)

    final_res_df = pd.DataFrame(cmap_res_analysis)








