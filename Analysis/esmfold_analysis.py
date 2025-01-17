import pandas as pd
import os
import subprocess
import re
from tqdm import tqdm


'''
To run this script on your local computer you have to specify where the TM align software is installed on your computer and adapt the function
get_tmscore_align changing your path.
And to set in pipeline_folder your path where the results are on your machine
'''

final_res       = []
# pipeline_folder = '/Users/steveabecassis/Desktop/Pipeline_res/Pipeline'
pipeline_folder = '/Users/steveabecassis/Desktop/Pipeline'
fold_pairs      = os.listdir(pipeline_folder)

# fold_pair = '1eboE_5fhcJ'
# fold1     = fold_pair.split('_')[0]
# fold2     = fold_pair.split('_')[1]

def get_tmscore_align(path_fold1,path_fold2):
    command = f"/Users/steveabecassis/Desktop/TMalign {path_fold1} {path_fold2}"
    output = subprocess.check_output(command, shell=True)
    match = re.search(r"TM-score=\s+(\d+\.\d+)", str(output))
    if match:
        result = match.group(1)
        return float(result)
    else:
        return None


if __name__=='__main__':

    for fold_pair in tqdm(fold_pairs):
        print(fold_pair)
        try:
            if 'sh' in fold_pair:
                continue
            for fold_pred in tqdm(os.listdir(f'{pipeline_folder}/{fold_pair}/output_esm_fold')):
                if len(os.listdir(f'{pipeline_folder}/{fold_pair}/output_esm_fold')) == 0:
                    continue
                fold1     = fold_pair.split('_')[0]
                fold2     = fold_pair.split('_')[1]
                path_fold_1   = f'{pipeline_folder}/{fold_pair}/chain_pdb_files/{fold1}.pdb'
                path_fold_2   = f'{pipeline_folder}/{fold_pair}/chain_pdb_files/{fold2}.pdb'
                TMscore_fold1 = get_tmscore_align(path_fold_1, f'{pipeline_folder}/{fold_pair}/output_esm_fold/{fold_pred}')
                TMscore_fold2 = get_tmscore_align(path_fold_2, f'{pipeline_folder}/{fold_pair}/output_esm_fold/{fold_pred}')
                final_res.append({'fold_pair':fold_pair,'fold':fold_pred,'TMscore_fold1':TMscore_fold1,'TMscore_fold2':TMscore_fold2})
        except Exception as e:
            print(e)
            continue

    df                         = pd.DataFrame(final_res)
    df['cluster_num']          = df.fold.apply(lambda x: x.split('_')[1] if 'Shallow' in x else x.split('_')[0])
    df['TM_mean_cluster_pdb1'] = df.groupby(['fold_pair','cluster_num'])['TMscore_fold1'].transform('mean')
    df['TM_mean_cluster_pdb2'] = df.groupby(['fold_pair','cluster_num'])['TMscore_fold2'].transform('mean')
    df['sample_count']         = df.groupby(['fold_pair','cluster_num'])['TMscore_fold1'].transform('count')
    df['is_fold_1']            = ((df.TM_mean_cluster_pdb1) >= (df.TM_mean_cluster_pdb2)).astype(int)
    df['is_fold_2']            = ((df.TM_mean_cluster_pdb1) < (df.TM_mean_cluster_pdb2)).astype(int)
    df['cluster_fold_1']       = df.groupby(['fold_pair','cluster_num'])['is_fold_1'].transform('mean')
    df['cluster_fold_2']       = df.groupby(['fold_pair','cluster_num'])['is_fold_2'].transform('mean')
    df.sort_values(by=['fold_pair','cluster_num'],inplace=True)
    df_summary                 = df[['fold_pair','TM_mean_cluster_pdb1', 'TM_mean_cluster_pdb2', 'cluster_num','is_fold_1', 'is_fold_2', 'cluster_fold_1', 'cluster_fold_2','sample_count']].drop_duplicates(subset=['fold_pair','cluster_num','TM_mean_cluster_pdb1', 'TM_mean_cluster_pdb2'])

    # df.to_csv('/Users/steveabecassis/Desktop/Pipeline_res/df_esmfold_analysis.csv',index=False)









    # for fold_pair in fold_pairs:
    #     if 'sh' in fold_pair:
    #         continue
    #     try:
    #         a=1
    #     except Exception as e:
    #         print(e)
