import os
import pandas as pd
from tqdm import tqdm
files = os.listdir('/Users/steveabecassis/Desktop/Pipeline_res/HTMLs_new_1507')
folder = '/Users/steveabecassis/Desktop/Pipeline'
files_notebook = os.listdir('/Users/steveabecassis/Desktop/Pipeline_res/HTMLs_new_2110')
fold_pairs = [i[:-5] for i in files_notebook]
res = []
for fold_pair in tqdm(fold_pairs):
    df_cmap_analysis = pd.read_parquet('/Users/steveabecassis/Desktop/Pipeline/cmap_exact_analysis_tol0_2510.parq')
    df_esmfold_analysis = pd.read_csv(f'/Users/steveabecassis/Desktop/Pipeline_res/df_esmfold_analysis.csv')
    try:
        if (len(os.listdir(f'/Users/steveabecassis/Desktop/Pipeline/{fold_pair}/output_msa_cluster')) < 2):
            continue
        df_af = pd.read_csv(f'{folder}/{fold_pair}/Analysis/df_af.csv')
        df_af = df_af[df_af.cluster_num != 'Query'].iloc[:, 1:-1]
        if(len(df_af[df_af.score_pdb1 > df_af.score_pdb2])>0):
            BEST_AF1         = round(df_af[df_af.score_pdb1 > df_af.score_pdb2].sort_values(by='score_pdb1',ascending=False).score_pdb1.iloc[0],3)
            BEST_AF1_CLUSTER = df_af[df_af.score_pdb1 > df_af.score_pdb2].sort_values(by='score_pdb1',ascending=False).iloc[0].cluster_num
        else:
            BEST_AF1 = '-'
            BEST_AF1_CLUSTER = ''
        if(len(df_af[df_af.score_pdb1 < df_af.score_pdb2])>0):
            BEST_AF2         = round(df_af[df_af.score_pdb1 < df_af.score_pdb2].sort_values(by='score_pdb2', ascending=False).score_pdb2.iloc[0],3)
            BEST_AF2_CLUSTER = df_af[df_af.score_pdb1 < df_af.score_pdb2].sort_values(by='score_pdb2',ascending=False).iloc[0].cluster_num

        else:
            BEST_AF2 = '-'
            BEST_AF2_CLUSTER = ''

        df_esmfold = df_esmfold_analysis[df_esmfold_analysis.fold_pair == fold_pair]
        df_esmfold = df_esmfold[df_esmfold.fold.str.contains('ShallowMsa')]

        if len(df_esmfold[df_esmfold.TM_mean_cluster_pdb1 > df_esmfold.TM_mean_cluster_pdb2]) > 0:
            BEST_ESM1 = round(df_esmfold[df_esmfold.TM_mean_cluster_pdb1 > df_esmfold.TM_mean_cluster_pdb2].sort_values(by='TMscore_fold1',ascending=False).TMscore_fold1.iloc[0],3)
            try:
                BEST_ESM1_AF = round(df_esmfold[df_esmfold.TM_mean_cluster_pdb1 > df_esmfold.TM_mean_cluster_pdb2][df_esmfold.cluster_num == BEST_AF1_CLUSTER].sort_values(by='TMscore_fold1',ascending=False).TMscore_fold1.iloc[0],3)
            except:
                BEST_ESM1_AF = '-'
        else:
            BEST_ESM1 = '-'
            BEST_ESM1_AF = '-'
        if len(df_esmfold[df_esmfold.TM_mean_cluster_pdb1 < df_esmfold.TM_mean_cluster_pdb2]) > 0:
            BEST_ESM2 = round(df_esmfold[df_esmfold.TM_mean_cluster_pdb1 < df_esmfold.TM_mean_cluster_pdb2].sort_values(by='TMscore_fold2',ascending=False).TMscore_fold2.iloc[0],3)
            try:
                BEST_ESM2_AF = round(df_esmfold[df_esmfold.TM_mean_cluster_pdb1 < df_esmfold.TM_mean_cluster_pdb2][df_esmfold.cluster_num == BEST_AF2_CLUSTER].sort_values(by='TMscore_fold2',ascending=False).TMscore_fold2.iloc[0], 3)
            except:
                BEST_ESM2_AF = '-'
        else:
            BEST_ESM2 = '-'
            BEST_ESM2_AF = '-'

        df_cmap_analysis = df_cmap_analysis[df_cmap_analysis.FoldPair == fold_pair]
        if len(df_cmap_analysis[df_cmap_analysis.recall_only_fold1 > df_cmap_analysis.recall_only_fold2])>0:
            BEST_RECALL_FOLD1 = df_cmap_analysis[df_cmap_analysis.recall_only_fold1 > df_cmap_analysis.recall_only_fold2].sort_values(by='recall_only_fold1',ascending=False).recall_only_fold1.iloc[0]
        else:
            BEST_RECALL_FOLD1 = '-'
        if len(df_cmap_analysis[df_cmap_analysis.recall_only_fold1 < df_cmap_analysis.recall_only_fold2])>0:
            BEST_RECALL_FOLD2 = df_cmap_analysis[df_cmap_analysis.recall_only_fold1 < df_cmap_analysis.recall_only_fold2].sort_values(by='recall_only_fold2',ascending=False).recall_only_fold2.iloc[0]
        else:
            BEST_RECALL_FOLD2 = '-'

        try:
            BEST_RECALL_FOLD1_AF = df_cmap_analysis[df_cmap_analysis.cluster_num ==BEST_AF1_CLUSTER].recall_only_fold1.iloc[0]
        except:
            BEST_RECALL_FOLD1_AF = '-'

        try:
            BEST_RECALL_FOLD2_AF = df_cmap_analysis[df_cmap_analysis.cluster_num ==BEST_AF2_CLUSTER].recall_only_fold2.iloc[0]
        except:
            BEST_RECALL_FOLD2_AF ='-'




        res.append({'fold_pair':fold_pair,'BEST_AF1':BEST_AF1,'BEST_AF2':BEST_AF2,'BEST_ESM1':BEST_ESM1,'BEST_ESM1_AF':BEST_ESM1_AF,
                    'BEST_ESM2':BEST_ESM2,'BEST_ESM2_AF':BEST_ESM2_AF,'BEST_RECALL_FOLD1':BEST_RECALL_FOLD1,'BEST_RECALL_FOLD1_AF':BEST_RECALL_FOLD1_AF,
                    'BEST_RECALL_FOLD2':BEST_RECALL_FOLD2,'BEST_RECALL_FOLD2_AF':BEST_RECALL_FOLD2_AF})

    except Exception as e:
        print(f'Fold Pair:{fold_pair}')
        print(e)
        continue


final_res_df = pd.DataFrame(res)
# final_res_df.astype(str).to_parquet('/Users/steveabecassis/Desktop/Pipeline_res/final_res_df_2510.parq')
# '4rmbB_4rmbA' , 1fzpD_2frhA 2lepA_4hddA
df = pd.read_parquet('/Users/steveabecassis/Desktop/Pipeline_res/final_res_df_2510.parq')
