import os
import pandas as pd
from tqdm import tqdm
import re

folder = '/Users/steveabecassis/Desktop/Pipeline'
df_cmap_analysis = pd.read_parquet('/Users/steveabecassis/Desktop/Pipeline/cmap_exact_analysis_tol0_2510.parq')
df_esmfold_analysis = pd.read_csv(f'/Users/steveabecassis/Desktop/Pipeline_res/df_esmfold_analysis.csv')
pattern = r'^[0-9a-zA-Z]{5}_[0-9a-zA-Z]{5}$'
res = []
df_af_all = pd.DataFrame()
fold_pairs = os.listdir('/Users/steveabecassis/Desktop/Pipeline')
for fold_pair in tqdm(fold_pairs):
   try:
        df_af = pd.read_csv(f'{folder}/{fold_pair}/Analysis/df_af.csv')
        df_af['fold_pair'] = fold_pair
        df_af_all = pd.concat([df_af_all,df_af])
   except:
        print(fold_pair)
        continue

df_af_all = df_af_all.iloc[:,1:]
# df_af_all.to_csv('/Users/steveabecassis/Desktop/df_af_all.csv',index=False)
# df_cmap_analysis.to_csv('/Users/steveabecassis/Desktop/df_cmap_all.csv',index=False)
# df_esmfold_analysis.to_csv('/Users/steveabecassis/Desktop/df_esmfold_all.csv',index=False)


'''
final_res_df = pd.DataFrame(res)
# final_res_df.astype(str).to_parquet('/Users/steveabecassis/Desktop/Pipeline_res/final_res_df_2510.parq')
# '4rmbB_4rmbA' , 1fzpD_2frhA 2lepA_4hddA
df = pd.read_parquet('/Users/steveabecassis/Desktop/Pipeline_res/final_res_df_2510.parq')
df.columns

af_success = df[(df['BEST_AF1'] != '-') & (df['BEST_AF2'] != '-')]
perc_AF_success = 100 - round(100*len(af_success)/len(df),2)
perc_AF_success


af1 = af_success[af_success['BEST_RECALL_FOLD1_AF'] >  af_success['BEST_RECALL_FOLD2_AF']][(af_success['BEST_AF1'] >  af_success['BEST_AF2']) | ( af_success['BEST_AF2'] == '-')]
af2 = af_success[af_success['BEST_RECALL_FOLD1_AF'] <  af_success['BEST_RECALL_FOLD2_AF']][(af_success['BEST_AF1'] <  af_success['BEST_AF2']) | ( af_success['BEST_AF1'] == '-')]
perc_cmap_af_success = 100 - round(100*(len(af1) + len(af2))/len(af_success),2)
perc_cmap_af_success

cmap_success = df[(df['BEST_RECALL_FOLD1'] == '-') | (df['BEST_RECALL_FOLD2'] == '-')]
perc_cmap_success_success = 100 - round(100*len(cmap_success)/len(df),2)
perc_cmap_success_success
'''









'''
AF SUCCESS --> 62/85
AF FAIL    --> 23/85
'''
