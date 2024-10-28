import os
import pandas as pd
from tqdm import tqdm
files = os.listdir('/Users/steveabecassis/Desktop/Pipeline_res/HTMLs_new_1507')
folder = '/Users/steveabecassis/Desktop/Pipeline'
for file in files:
    try:
        fold_pair = file[:-5]
        if (len(os.listdir(f'/Users/steveabecassis/Desktop/Pipeline/{fold_pair}/output_msa_cluster')) < 2):
            continue
        df_af = pd.read_csv(f'{folder}/{fold_pair}/Analysis/df_af.csv')
        df_af = df_af[df_af.cluster_num != 'Query'].iloc[:, 1:-1]
        if(len(df_af[df_af.score_pdb1 > df_af.score_pdb2])>0):
            BEST_AF1  = round(df_af[df_af.score_pdb1 > df_af.score_pdb2].sort_values(by='score_pdb1',ascending=False).score_pdb1.iloc[0],3)
        else:
            BEST_AF1 = '-'
        if(len(df_af[df_af.score_pdb1 < df_af.score_pdb2])>0):
            BEST_AF2 = round(df_af[df_af.score_pdb1 < df_af.score_pdb2].sort_values(by='score_pdb2', ascending=False).score_pdb2.iloc[0],3)
        else:
            BEST_AF2 = '-'

        df_esmfold_analysis = pd.read_csv(f'/Users/steveabecassis/Desktop/Pipeline_res/df_esmfold_analysis.csv')
        df_esmfold = df_esmfold_analysis[df_esmfold_analysis.fold_pair == fold_pair]
        df_esmfold = df_esmfold[df_esmfold.fold.str.contains('ShallowMsa')]

        if len(df_esmfold[df_esmfold.TM_mean_cluster_pdb1 > df_esmfold.TM_mean_cluster_pdb2]) > 0:
            BEST_ESM1 = round(df_esmfold[df_esmfold.TM_mean_cluster_pdb1 > df_esmfold.TM_mean_cluster_pdb2].sort_values(by='TMscore_fold1',ascending=False).TMscore_fold1.iloc[0],3)
        else:
            BEST_ESM1 = '-'
        if len(df_esmfold[df_esmfold.TM_mean_cluster_pdb1 < df_esmfold.TM_mean_cluster_pdb2]) > 0:
            BEST_ESM2 = round(df_esmfold[df_esmfold.TM_mean_cluster_pdb1 < df_esmfold.TM_mean_cluster_pdb2].sort_values(by='TMscore_fold2',ascending=False).TMscore_fold2.iloc[0],3)
        else:
            BEST_ESM2 = '-'

        try:
            a = ''
        except:
            MAX_RECALL_1 = '-'
            MAX_RECALL_2 = '-'


        print(f'<tr>')
        print(f'    <td><a href="https://steveabecassis.github.io/MsaCluster/HTML/{file}" target="_blank">{fold_pair}</a></td>')
        print(f'    <td>{BEST_AF1}</td>')
        print(f'    <td>{BEST_AF2}</td>')
        print(f'    <td>{BEST_ESM1}</td>')
        print(f'    <td>{BEST_ESM2}</td>')
        print(f'    <td>{MAX_RECALL_1}</td>')
        print(f'    <td>{MAX_RECALL_2}</td>')
        print(f'</tr>')
    except:
        continue
