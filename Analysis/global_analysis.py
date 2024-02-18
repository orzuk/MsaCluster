import pandas as pd
import os

final_res = []
pipeline_folder = '/Users/steveabecassis/Desktop/Pipeline_res/Pipeline'
fold_pairs = os.listdir(pipeline_folder)
# fold_pair = '1fzpD_2frhA'

for fold_pair in fold_pairs:
    try:
        df_cmap = pd.read_csv(f'{pipeline_folder}/{fold_pair}/Analysis/cmap_scores_df.csv')
        df_cmap = df_cmap[df_cmap.file_name.str.contains('Sha')]
        df_cmap['cluster_num'] = df_cmap['file_name'].apply(lambda x: int(x[-3:]))

        df_cmap['JI_mean_cluster_fold1'] = df_cmap.groupby('cluster_num')['jaccard_index_score_fold1'].transform('mean')
        df_cmap['JI_mean_cluster_fold2'] = df_cmap.groupby('cluster_num')['jaccard_index_score_fold2'].transform('mean')
        df_cmap['Spearman_mean_cluster_fold1'] = df_cmap.groupby('cluster_num')['spearmanr_score_fold1'].transform('mean')
        df_cmap['Spearman_mean_cluster_fold2'] = df_cmap.groupby('cluster_num')['spearmanr_score_fold2'].transform('mean')

        JI_fold1_ratio = round((df_cmap['jaccard_index_score_fold1'] >= df_cmap['jaccard_index_score_fold2'] ).mean(),2)
        JI_max_mean_cluster_fold_1 = df_cmap['JI_mean_cluster_fold1'].max()

        JI_fold2_ratio = round((df_cmap['jaccard_index_score_fold1'] < df_cmap['jaccard_index_score_fold2'] ).mean(),2)
        JI_max_mean_cluster_fold_2 = df_cmap['JI_mean_cluster_fold2'].max()

        Spearman_fold1_ratio = round((df_cmap['spearmanr_score_fold1'] >= df_cmap['spearmanr_score_fold2'] ).mean(),2)
        Spearman_max_mean_cluster_fold_1 = df_cmap['Spearman_mean_cluster_fold1'].max()

        Spearman_fold2_ratio = round((df_cmap['spearmanr_score_fold1'] < df_cmap['spearmanr_score_fold2'] ).mean(),2)
        Spearman_max_mean_cluster_fold_2 = df_cmap['Spearman_mean_cluster_fold2'].max()

        df_af = pd.read_csv(f'{pipeline_folder}/{fold_pair}/Analysis/df_af.csv')
        df_af = df_af[df_af.pdb_file.str.contains('Sha')]
        df_af['TM_mean_cluster_pdb1'] = df_af.groupby('cluster_num')['score_pdb1'].transform('mean')
        df_af['TM_mean_cluster_pdb2'] = df_af.groupby('cluster_num')['score_pdb2'].transform('mean')

        AF_max_ratio_fold1 = df_af['unique_fold_score'].max()
        AF_max_mean_cluster_fold_1 = df_af['TM_mean_cluster_pdb1'].max()
        AF_max_ratio_fold2 = 1 - df_af['unique_fold_score'].min()
        AF_max_mean_cluster_fold_2 = df_af['TM_mean_cluster_pdb2'].max()
        final_res.append({'fold_pair':fold_pair,'JI_fold1_ratio':JI_fold1_ratio,'JI_max_mean_cluster_fold_1':JI_max_mean_cluster_fold_1,
                          'JI_fold2_ratio':JI_fold2_ratio,'JI_max_mean_cluster_fold_2':JI_max_mean_cluster_fold_2,
                          'Spearman_fold1_ratio':Spearman_fold1_ratio,'Spearman_max_mean_cluster_fold_1':Spearman_max_mean_cluster_fold_1,
                          'Spearman_fold2_ratio':Spearman_fold2_ratio,'Spearman_max_mean_cluster_fold_2':Spearman_max_mean_cluster_fold_2,
                          'AF_max_ratio_fold1':AF_max_ratio_fold1,'AF_max_mean_cluster_fold_1':AF_max_mean_cluster_fold_1,
                          'AF_max_ratio_fold2':AF_max_ratio_fold2,'AF_max_mean_cluster_fold_2':AF_max_mean_cluster_fold_2})
    except:
        continue


final_res_df = pd.DataFrame(final_res)

