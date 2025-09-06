import pandas as pd
import os
import re
from config import *
'''
check get msa cluster
get af prediction
get msat pred
get esmfold pred
get df_af file analysis
get esmfold analysis
get msat analysis
'''
def filter_pdb_pairs(file_list):
    # Pattern: 4 characters + 1 letter + underscore + 4 characters + 1 letter
    pattern = r'^[A-Za-z0-9]{4}[A-Za-z]_[A-Za-z0-9]{4}[A-Za-z]$'

    # Filter the list
    filtered_list = [item for item in file_list if re.match(pattern, item)]
    return filtered_list

fold_pairs =  filter_pdb_pairs(os.listdir(DATA_DIR))

monitoring = []
for fold_pair in fold_pairs:
    try:
        cluster_count = len(os.listdir(f'{DATA_DIR}/{fold_pair}/output_msa_cluster'))
        af_pdb_files = [i for i in os.listdir(f'{DATA_DIR}/{fold_pair}/AF_preds') if '.pdb.gz' in i]
        af_pred_cnt = (len(af_pdb_files)//5) - 1
        if af_pred_cnt < cluster_count:
            af_pdb_files = [i for i in os.listdir(f'{DATA_DIR}/{fold_pair}/AF_preds') if '.pdb' in i]
            af_pred_cnt = (len(af_pdb_files) // 5) - 1
        af_analysis = 'df_af.csv' in os.listdir(f'{DATA_DIR}/{fold_pair}/Analysis')
        esm_fold_pdb_files = [i for i in os.listdir(f'{DATA_DIR}/{fold_pair}/output_esm_fold') if '.pdb' in i]
        esm_fold_cnt =  (len(esm_fold_pdb_files) - 2)//10
        cmap_msa_transformer_pdb_files = [i for i in os.listdir(f'{DATA_DIR}/{fold_pair}/output_cmaps/msa_transformer')]
        cmap_cnt = len(cmap_msa_transformer_pdb_files) - 2
        monitoring.append({'fold_pair':fold_pair,'cluster_count':cluster_count,'af_pred_cnt':af_pred_cnt,'af_analysis':af_analysis,'cmap_cnt':cmap_cnt,'esm_fold_cnt':esm_fold_cnt})
    except Exception as e:
        print(e)
        continue

monitoring_df = pd.DataFrame(monitoring)
monitoring_df.shape


cluster_fail = monitoring_df[monitoring_df.cluster_count == 0].fold_pair.tolist()
len(cluster_fail)
monitoring_df = monitoring_df[monitoring_df.cluster_count != 0]

error_df = monitoring_df[(monitoring_df.cluster_count != monitoring_df.af_pred_cnt)|(monitoring_df.cluster_count > monitoring_df.cmap_cnt)|
                         (monitoring_df.cluster_count != monitoring_df.esm_fold_cnt) | (monitoring_df.af_analysis == False)]


