import os
import pandas as pd
from tqdm import tqdm
files = os.listdir('./HTML')
folder = '/Users/steveabecassis/Desktop/Pipeline'
for file in files:
    try:
        fold_pair = file[:-5]
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

        if len(df_esmfold[df_esmfold.TM_mean_cluster_pdb1 > df_esmfold.TM_mean_cluster_pdb2]) > 0:
            BEST_ESM1 = round(df_esmfold[df_esmfold.TM_mean_cluster_pdb1 > df_esmfold.TM_mean_cluster_pdb2].sort_values(by='TMscore_fold1',ascending=False).TMscore_fold1.iloc[0],3)
        else:
            BEST_ESM1 = '-'
        if len(df_esmfold[df_esmfold.TM_mean_cluster_pdb1 < df_esmfold.TM_mean_cluster_pdb2]) > 0:
            BEST_ESM2 = round(df_esmfold[df_esmfold.TM_mean_cluster_pdb1 < df_esmfold.TM_mean_cluster_pdb2].sort_values(by='TMscore_fold2',ascending=False).TMscore_fold2.iloc[0],3)
        else:
            BEST_ESM2 = '-'

        try:
            df_cmap = pd.read_csv(f'{folder}/{fold_pair}/Analysis/new_cmap_res.csv')
            df_cmap['f1Score_1'] = (df_cmap['precision_1'] * df_cmap['recall_1']) / (df_cmap['precision_1'] + df_cmap['recall_1'])
            df_cmap['f1Score_2'] = (df_cmap['precision_2'] * df_cmap['recall_2']) / (df_cmap['precision_2'] + df_cmap['recall_2'])
            if len(df_cmap[df_cmap.f1Score_1 > df_cmap.f1Score_2])>0:
                MAX_F1_1 = round(df_cmap[df_cmap.f1Score_1 > df_cmap.f1Score_2].sort_values(by='f1Score_1',ascending=False).f1Score_1.iloc[0],3)
            else:
                MAX_F1_1 = '-'
            if len(df_cmap[df_cmap.f1Score_1 < df_cmap.f1Score_2])>0:
                MAX_F1_2 = round(df_cmap[df_cmap.f1Score_1 < df_cmap.f1Score_2].sort_values(by='f1Score_2',ascending=False).f1Score_2.iloc[0],3)
            else:
                MAX_F1_2 = '-'
            if len(df_cmap[df_cmap.recall_1 > df_cmap.recall_2])>0:
                MAX_RECALL_1 = round(df_cmap[df_cmap.recall_1 > df_cmap.recall_2].sort_values(by='recall_1',ascending=False).recall_1.iloc[0],3)
            else:
                MAX_RECALL_1 = '-'
            if len(df_cmap[df_cmap.recall_1 < df_cmap.recall_2])>0:
                MAX_RECALL_2 = round(df_cmap[df_cmap.recall_1 < df_cmap.recall_2].sort_values(by='recall_2',ascending=False).recall_2.iloc[0],3)
            else:
                MAX_RECALL_2 = '-'
        except:
            MAX_F1_1 = '-'
            MAX_F1_2 = '-'
            MAX_RECALL_1 = '-'
            MAX_RECALL_2 = '-'


        print(f'<tr>')
        print(f'    <td><a href="https://steveabecassis.github.io/MsaCluster/HTML/{file}" target="_blank">{fold_pair}</a></td>')
        print(f'    <td>{BEST_AF1}</td>')
        print(f'    <td>{BEST_AF2}</td>')
        print(f'    <td>{BEST_ESM1}</td>')
        print(f'    <td>{BEST_ESM2}</td>')
        print(f'    <td>{MAX_F1_1}</td>')
        print(f'    <td>{MAX_F1_2}</td>')
        print(f'    <td>{MAX_RECALL_1}</td>')
        print(f'    <td>{MAX_RECALL_2}</td>')
        print(f'</tr>')
    except:
        continue
