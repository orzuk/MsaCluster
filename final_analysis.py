import numpy as np
import pandas as pd
import os
import warnings
from tqdm import tqdm

warnings.filterwarnings('ignore')

final_res = []
# pipeline_folder = '/Users/steveabecassis/Desktop/Pipeline_res/Pipeline'
pipeline_folder = '/Users/steveabecassis/Desktop/Pipeline'
fold_pairs = os.listdir(pipeline_folder)

final_df = pd.DataFrame()

df_esmfold_analysis = pd.read_csv(f'{pipeline_folder}_res/df_esmfold_analysis.csv')


for fold_pair in tqdm(fold_pairs):
    if 'sh' in fold_pair:
        continue
    try:
        df_cmap = pd.read_csv(f'{pipeline_folder}/{fold_pair}/Analysis/cmap_scores_df.csv')
        df_cmap = df_cmap[df_cmap.file_name.str.contains('Sha')]
        df_cmap['cluster_num'] = df_cmap['file_name'].apply(lambda x: int(x[-3:]))
        #         df_cmap['cmap_fold1_JI'] = (df_cmap['jaccard_index_score_fold1'] > df_cmap['jaccard_index_score_fold2'])
        df_cmap['CMAP_fold1'] = (df_cmap['spearmanr_score_fold1'] > df_cmap['spearmanr_score_fold2']).astype(int)
        df_cmap['CMAP_fold2'] = (df_cmap['spearmanr_score_fold1'] < df_cmap['spearmanr_score_fold2']).astype(int)
        df_cmap['Spearman4th'] = df_cmap.apply(lambda x: max(x.spearmanr_score_fold1, x.spearmanr_score_fold2), axis=1)
        df_cmap['CMAPdiff'] = abs(df_cmap['spearmanr_score_fold1'] - df_cmap['spearmanr_score_fold2'])


        df_af = pd.read_csv(f'{pipeline_folder}/{fold_pair}/Analysis/df_af.csv')
        df_af = df_af[df_af.pdb_file.str.contains('Sha')]
        df_af['TM_mean_cluster_pdb1'] = df_af.groupby('cluster_num')['score_pdb1'].transform('mean')
        df_af['TM_mean_cluster_pdb2'] = df_af.groupby('cluster_num')['score_pdb2'].transform('mean')
        df_af['AFTMdiff'] = abs(df_af['TM_mean_cluster_pdb1'] - df_af['TM_mean_cluster_pdb2'])
        df_af['AFTM4th'] = df_af.apply(lambda x: max(x.TM_mean_cluster_pdb1, x.TM_mean_cluster_pdb2), axis=1)

        df_af['is_fold_1']      = ((df_af.score_pdb1) >= (df_af.score_pdb2)).astype(int)
        df_af['is_fold_2']      = ((df_af.score_pdb1) < (df_af.score_pdb2)).astype(int)
        df_af['cluster_fold_1'] = df_af.groupby('cluster_num')['is_fold_1'].transform('mean')
        df_af['cluster_fold_2'] = df_af.groupby('cluster_num')['is_fold_2'].transform('mean')

        df_af['AF_fold1'] = df_af.groupby('cluster_num')['is_fold_1'].transform('mean')
        df_af['AF_fold2'] = df_af.groupby('cluster_num')['is_fold_2'].transform('mean')

        df_cmap['cluster_num'] = df_cmap['cluster_num'].astype(int)
        df_af['cluster_num'] = df_af['cluster_num'].astype(int)
        df = pd.merge(df_cmap, df_af, on='cluster_num')

        df['fold_pair'] = fold_pair

        df_esmfold = df_esmfold_analysis[df_esmfold_analysis.fold_pair == fold_pair]
        df_esmfold['EsmTM4th'] = df_esmfold.apply(lambda x: max(x.TM_mean_cluster_pdb1, x.TM_mean_cluster_pdb2), axis=1)
        df_esmfold['EsmTMdiff'] = abs(df_esmfold.TM_mean_cluster_pdb1 - df_esmfold.TM_mean_cluster_pdb2)
        df_esmfold = df_esmfold[df_esmfold['cluster_num'].apply(lambda x: len(x) <= 3)]
        df_esmfold['cluster_num'] = df_esmfold['cluster_num'].astype(int)

        df_esmfold = df_esmfold[['cluster_num','is_fold_1','is_fold_2','EsmTMdiff','EsmTM4th']]
        df_esmfold.columns = ['cluster_num','ESM_fold1','ESM_fold2','EsmTMdiff','EsmTM4th']

        cols2keep = ['fold_pair','cluster_num','CMAP_fold1','CMAP_fold2','Spearman4th','CMAPdiff','AF_fold1','AF_fold2','AFTMdiff','AFTM4th','ESM_fold1','ESM_fold2','EsmTMdiff','EsmTM4th']

        df = pd.merge(df,df_esmfold,on='cluster_num')
        df = df[cols2keep].drop_duplicates()

        final_df = final_df.append(df)
    except:
        continue


final_df.sort_values(by=['fold_pair','cluster_num'],inplace=True)
final_df.fold_pair.nunique()

af_cmap_esm_2FOLDS = []
af_cmap_2FOLDS     = []
cmap_2FOLDS        = []
af_esm_2FOLDS      = []
AF_2FOLDS          = []
cmap_esm_2FOLDS    = []

TM_th = 0.4
AFTMdiff_th= 0.005
spearman_th = 0.1
CMAPdiff_th = 0.001

for fold_pair in final_df.fold_pair.unique():
    df_temp = final_df[final_df.fold_pair == fold_pair]
    df_temp = df_temp[df_temp.Spearman4th > spearman_th][df_temp.AFTM4th > TM_th][df_temp.AFTMdiff > AFTMdiff_th][df_temp.CMAPdiff > CMAPdiff_th]
    f1 = len(df_temp[(df_temp.CMAP_fold1 == 1)&(df_temp.AF_fold1 == 1)&(df_temp.ESM_fold1 == 1)]) > 0
    f2 = len(df_temp[(df_temp.CMAP_fold2 == 1)&(df_temp.AF_fold2 == 1)&(df_temp.ESM_fold2 == 1)]) > 0
    if (f1 & f2):
        af_cmap_esm_2FOLDS.append(fold_pair)
    df_temp = final_df[final_df.fold_pair == fold_pair]
    df_temp = df_temp[df_temp.AFTM4th > TM_th][df_temp.AFTMdiff > AFTMdiff_th]
    f1 = len(df_temp[(df_temp.AF_fold1 == 1)&(df_temp.ESM_fold1 == 1)]) > 0
    f2 = len(df_temp[(df_temp.AF_fold2 == 1)&(df_temp.ESM_fold2 == 1)]) > 0
    if (f1 & f2):
        af_esm_2FOLDS.append(fold_pair)
    df_temp = final_df[final_df.fold_pair == fold_pair]
    df_temp = df_temp[df_temp.Spearman4th > spearman_th][df_temp.CMAPdiff > CMAPdiff_th]
    f1 = len(df_temp[(df_temp.CMAP_fold1 == 1)&(df_temp.ESM_fold1 == 1)]) > 0
    f2 = len(df_temp[(df_temp.CMAP_fold2 == 1)&(df_temp.ESM_fold2 == 1)]) > 0
    if (f1 & f2):
        cmap_esm_2FOLDS.append(fold_pair)
    df_temp = final_df[final_df.fold_pair == fold_pair]
    df_temp = df_temp[df_temp.Spearman4th > spearman_th][df_temp.CMAPdiff > CMAPdiff_th][df_temp.AFTM4th > TM_th][df_temp.AFTMdiff > AFTMdiff_th]
    f1 = len(df_temp[(df_temp.CMAP_fold1 == 1)&(df_temp.AF_fold1 == 1)]) > 0
    f2 = len(df_temp[(df_temp.CMAP_fold2 == 1)&(df_temp.AF_fold2 == 1)]) > 0
    if (f1 & f2):
        af_cmap_2FOLDS.append(fold_pair)
    df_temp = final_df[final_df.fold_pair == fold_pair]
    df_temp = df_temp[df_temp.Spearman4th > spearman_th][df_temp.CMAPdiff > CMAPdiff_th]
    f1 = len(df_temp[(df_temp.CMAP_fold1 == 1)]) > 0
    f2 = len(df_temp[(df_temp.CMAP_fold2 == 1)]) > 0
    if (f1 & f2):
        cmap_2FOLDS.append(fold_pair)
    df_temp = final_df[final_df.fold_pair == fold_pair]
    df_temp = df_temp[df_temp.AFTM4th > TM_th][df_temp.AFTMdiff > AFTMdiff_th]
    f1 = len(df_temp[(df_temp.AF_fold1 == 1)]) > 0
    f2 = len(df_temp[(df_temp.AF_fold2 == 1)]) > 0
    if (f1 & f2):
        AF_2FOLDS.append(fold_pair)




len(af_cmap_esm_2FOLDS)
len(af_cmap_2FOLDS)
len(cmap_2FOLDS)
len(af_esm_2FOLDS)
len(AF_2FOLDS)
len(cmap_esm_2FOLDS)

df_interesting = final_df[final_df.fold_pair.isin(af_cmap_esm_2FOLDS)]


all_int_pair = af_cmap_esm_2FOLDS+af_cmap_2FOLDS+cmap_2FOLDS+af_esm_2FOLDS+AF_2FOLDS+cmap_esm_2FOLDS
import numpy as np
len(np.unique(all_int_pair))

np.quantile(final_df['CMAPdiff'],0.05)


print('Finish !')