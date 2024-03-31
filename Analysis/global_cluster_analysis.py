import pandas as pd
import os

final_res = []
pipeline_folder = '/Users/steveabecassis/Desktop/Pipeline_res/Pipeline'
fold_pairs = os.listdir(pipeline_folder)
# fold_pair = '1fzpD_2frhA'


if __name__=='__main__':

    final_df = pd.DataFrame()

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

            df_cmap['cmap_fold1_JI'] = (df_cmap['jaccard_index_score_fold1'] > df_cmap['jaccard_index_score_fold2'])
            df_cmap['cmap_fold1_spearman'] = (df_cmap['spearmanr_score_fold1'] > df_cmap['spearmanr_score_fold2'])

            if ((len(df_cmap[df_cmap['cmap_fold1_JI']==1])>0)&(len(df_cmap[df_cmap['cmap_fold1_spearman']==1])>0)):
                df_cmap['CMAP2FOLDS'] = 1
            else:
                df_cmap['CMAP2FOLDS'] = 0


            df_af = pd.read_csv(f'{pipeline_folder}/{fold_pair}/Analysis/df_af.csv')
            df_af = df_af[df_af.pdb_file.str.contains('Sha')]
            df_af['TM_mean_cluster_pdb1'] = df_af.groupby('cluster_num')['score_pdb1'].transform('mean')
            df_af['TM_mean_cluster_pdb2'] = df_af.groupby('cluster_num')['score_pdb2'].transform('mean')

            df_af['is_fold_1']        = ((df_af.score_pdb1) >= (df_af.score_pdb2)).astype(int)
            df_af['is_fold_2']        = ((df_af.score_pdb1) < (df_af.score_pdb2)).astype(int)
            df_af['cluster_fold_1']   = df_af.groupby('cluster_num')['is_fold_1'].transform('mean')
            df_af['cluster_fold_2']   = df_af.groupby('cluster_num')['is_fold_2'].transform('mean')
            if ((len(df_af[df_af['cluster_fold_1']==1])>0)&(len(df_af[df_af['cluster_fold_2']==1])>0)):
                df_af['AF2FOLDS'] = 1
            else:
                df_af['AF2FOLDS'] = 0
            df_af = df_af[(df_af['cluster_fold_1']==1)|(df_af['cluster_fold_2']==1)]
            df_cmap['cluster_num'] = df_cmap['cluster_num'].astype(int)
            df_af['cluster_num']   = df_af['cluster_num'].astype(int)
            df = pd.merge(df_cmap,df_af,on='cluster_num')

            df['fold_pair'] = fold_pair
            cols2keep = ['fold_pair','cluster_num','jaccard_index_score_fold1','jaccard_index_score_fold2', 'spearmanr_score_fold1',
                          'spearmanr_score_fold2','cmap_fold1_JI','cmap_fold1_spearman','TM_mean_cluster_pdb1','TM_mean_cluster_pdb2',
                          'cluster_fold_1','cluster_fold_2','AF2FOLDS','CMAP2FOLDS']
            df = df[cols2keep]
            final_df = final_df.append(df)
        except:
            continue

    final_df.drop_duplicates(inplace=True)
    final_df.fold_pair.nunique()
    df_2foldsAF = final_df[final_df.AF2FOLDS==1]
    df_2foldsAF.fold_pair.nunique()
    df_2folds_only_CMAP = final_df[(final_df.AF2FOLDS==0)&(final_df.CMAP2FOLDS==1)]
    df_2folds_only_CMAP.fold_pair.nunique()


    print('Finish !')
    # df_af = df_af[(df_af['cluster_fold_1']==1)|(df_af['cluster_fold_2']==1)]
    # df_interesting =