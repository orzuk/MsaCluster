import os
import pandas as pd
from tqdm import tqdm
import re
# Folder containing all the output results of the fold pairs
RES_PATH = '/Users/steveabecassis/Desktop/Pipeline'
# File create by the script cmap_analysis.py (in the analysis folder)
CMAP_ANALYSIS_FILE = '/Users/steveabecassis/Desktop/Pipeline/cmap_exact_analysis_tol0_2510.parq'
# File create by the script esmfold_analysis.py (in the analysis folder)
ESMFOLD_ANALYSIS_FILE = '/Users/steveabecassis/Desktop/Pipeline_res/df_esmfold_analysis.csv'

files = os.listdir(RES_PATH)
pattern = r'^[0-9a-zA-Z]{5}_[0-9a-zA-Z]{5}$'
fold_pairs = [i for i in files if re.match(pattern, i) ]
res = []

for fold_pair in tqdm(fold_pairs):
    df_cmap_analysis = pd.read_parquet(CMAP_ANALYSIS_FILE)
    df_esmfold_analysis = pd.read_csv(ESMFOLD_ANALYSIS_FILE)
    try:
        if (len(os.listdir(f'{RES_PATH}/{fold_pair}/output_msa_cluster')) < 2):
            continue
        df_af = pd.read_csv(f'{RES_PATH}/{fold_pair}/Analysis/df_af.csv')
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
df = pd.read_parquet('/Users/steveabecassis/Desktop/Pipeline_res/final_res_df_2510.parq')
