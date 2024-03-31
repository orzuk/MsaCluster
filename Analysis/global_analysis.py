import pandas as pd
import os

final_res = []
pipeline_folder = '/Users/steveabecassis/Desktop/Pipeline_res/Pipeline'
fold_pairs = os.listdir(pipeline_folder)
# fold_pair = '1fzpD_2frhA'


if __name__=='__main__':

    for fold_pair in fold_pairs:
        if 'sh' in fold_pair:
            continue
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
        except:
            JI_fold1_ratio=None
            JI_max_mean_cluster_fold_1=None
            JI_fold2_ratio=None
            JI_max_mean_cluster_fold_2=None
            Spearman_fold1_ratio=None
            Spearman_max_mean_cluster_fold_1=None
            Spearman_fold2_ratio=None
            Spearman_max_mean_cluster_fold_2=None

        try:
            df_af = pd.read_csv(f'{pipeline_folder}/{fold_pair}/Analysis/df_af.csv')
            df_af = df_af[df_af.pdb_file.str.contains('Sha')]
            df_af['TM_mean_cluster_pdb1'] = df_af.groupby('cluster_num')['score_pdb1'].transform('mean')
            df_af['TM_mean_cluster_pdb2'] = df_af.groupby('cluster_num')['score_pdb2'].transform('mean')

            # AF_max_ratio_fold1 = df_af['unique_fold_score'].max()
            AF_max_mean_cluster_fold_1 = df_af['TM_mean_cluster_pdb1'].max()
            # AF_max_ratio_fold2 = 1 - df_af['unique_fold_score'].min()
            AF_max_mean_cluster_fold_2 = df_af['TM_mean_cluster_pdb2'].max()
            df_af['is_fold_1']      = ((df_af.score_pdb1) >= (df_af.score_pdb2)).astype(int)
            df_af['is_fold_2']      = ((df_af.score_pdb1) < (df_af.score_pdb2)).astype(int)
            df_af['cluster_fold_1'] = df_af.groupby('cluster_num')['is_fold_1'].transform('mean')
            df_af['cluster_fold_2'] = df_af.groupby('cluster_num')['is_fold_2'].transform('mean')
            AF_max_ratio_fold1 = df_af['cluster_fold_1'].max()
            AF_max_ratio_fold2 = df_af['cluster_fold_2'].max()

        except:
            AF_max_ratio_fold1 = None
            AF_max_mean_cluster_fold_1 = None
            AF_max_ratio_fold2 = None
            AF_max_mean_cluster_fold_2 = None

        final_res.append({'fold_pair':fold_pair,'JI_fold1_ratio':JI_fold1_ratio,'JI_max_mean_cluster_fold_1':JI_max_mean_cluster_fold_1,
                          'JI_fold2_ratio':JI_fold2_ratio,'JI_max_mean_cluster_fold_2':JI_max_mean_cluster_fold_2,
                          'Spearman_fold1_ratio':Spearman_fold1_ratio,'Spearman_max_mean_cluster_fold_1':Spearman_max_mean_cluster_fold_1,
                          'Spearman_fold2_ratio':Spearman_fold2_ratio,'Spearman_max_mean_cluster_fold_2':Spearman_max_mean_cluster_fold_2,
                          'AF_max_ratio_fold1':AF_max_ratio_fold1,'AF_max_mean_cluster_fold_1':AF_max_mean_cluster_fold_1,
                          'AF_max_ratio_fold2':AF_max_ratio_fold2,'AF_max_mean_cluster_fold_2':AF_max_mean_cluster_fold_2})
        # except:
        #     continue


    final_res_df = pd.DataFrame(final_res)
    final_res_df = final_res_df[(~final_res_df.JI_fold1_ratio.isna())&(~final_res_df.AF_max_ratio_fold2.isna())]
    df_interesting = final_res_df[final_res_df.AF_max_ratio_fold1==1][final_res_df.AF_max_ratio_fold2==1][~final_res_df.AF_max_ratio_fold1.isna()]
    df_interesting_cmap =  final_res_df[~final_res_df.fold_pair.isin(df_interesting.fold_pair.tolist())][(final_res_df.JI_fold1_ratio>0)&(final_res_df.JI_fold1_ratio<1)][(final_res_df.Spearman_fold1_ratio>0)&(final_res_df.Spearman_fold1_ratio<1)]

    fold_pair  = '3g0hA_3ewsB'
    # df_cmap = pd.read_csv(f'{pipeline_folder}/{fold_pair}/Analysis/cmap_scores_df.csv')
    df_af = pd.read_csv(f'{pipeline_folder}/{fold_pair}/Analysis/df_af.csv').iloc[:,1:]
    df_af = df_af[df_af.pdb_file.str.contains('Sha')]
    df_af['TM_mean_cluster_pdb1'] = df_af.groupby('cluster_num')['score_pdb1'].transform('mean')
    df_af['TM_mean_cluster_pdb2'] = df_af.groupby('cluster_num')['score_pdb2'].transform('mean')

    df_af['is_fold_1'] = ((df_af.score_pdb1) >= (df_af.score_pdb2)).astype(int)
    df_af['is_fold_2'] = ((df_af.score_pdb1) < (df_af.score_pdb2)).astype(int)
    df_af['cluster_fold_1'] = df_af.groupby('cluster_num')['is_fold_1'].transform('mean')
    df_af['cluster_fold_2'] = df_af.groupby('cluster_num')['is_fold_2'].transform('mean')
    df_af_interesting = df_af[(df_af['cluster_fold_1']==1)|(df_af['cluster_fold_2']==1)]

    df_af_interesting_ = pd.concat([df_af_interesting.sort_values(by=['cluster_fold_1','TM_mean_cluster_pdb1']).iloc[:5,:],df_af_interesting.sort_values(by=['cluster_fold_2','TM_mean_cluster_pdb2']).iloc[:5,:]])