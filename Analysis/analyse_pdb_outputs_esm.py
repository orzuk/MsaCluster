import tmscoring
import os
import pandas as pd
from tqdm import tqdm
import re

#output_pipeline_1iyt
#1iyt
#2nao


if __name__ == '__main__':
    input_path = '/Users/steveabecassis/Desktop/PipelineTest/output_pipeline_1jfk'
    # input_path = '/Users/steveabecassis/Desktop/PipelineTest/output_pipeline_1iyt'
    pdb1 = f'{input_path}/1jfk.pdb'
    # pdb1 = f'{input_path}/1iyt.pdb'
    pdb1_name = pdb1[-8:-4]
    pdb2 =f'{input_path}/2nxq.pdb'
    # pdb2 = f'{input_path}/2nao.pdb'
    pdb2_name = pdb2[-8:-4]


    # pdb_folder = '/Users/steveabecassis/Desktop/pdb_output/2LEP'
    pdb_folder = '/Users/steveabecassis/Desktop/PipelineTest/output_pipeline_1jfk/esm_output/esm_fold_output'
    pdb_files = os.listdir(pdb_folder)
    res = []
    pdb_files = [i for i in pdb_files if 'pdb' in str(i)]
    for pdb_file in tqdm(pdb_files):
        if 'pdb' in str(pdb_file):
            pdb_file_path = f'{pdb_folder}/{pdb_file}'
            score_pdb1 = tmscoring.get_tm(pdb_file_path, pdb1)
            score_pdb2 = tmscoring.get_tm(pdb_file_path, pdb2)
            temp = {'pdb_file':pdb_file,'score_pdb1':score_pdb1,'score_pdb2':score_pdb2}
            res.append(temp)
    print('Finish !')


    df = pd.DataFrame(res)
    df.loc[df.score_pdb1 > df.score_pdb2 ,'Fold'] = pdb1_name
    df.loc[df.score_pdb1 < df.score_pdb2, 'Fold'] = pdb2_name
    df.sort_values(by='pdb_file',inplace=True)
    df['cluster_num'] = df.pdb_file.apply(lambda x : x[7:10])
    df['seq_num']     = df.pdb_file.apply(lambda x : re.search(r'_(\d+)\.pdb', x).group(1))
    df['seq_num'] = df['seq_num'].astype(int)
    df.sort_values(by=['cluster_num','seq_num'], inplace=True)
    df = df[df['pdb_file'].str.contains('Cl')]
    fold1             = df.Fold.iloc[0]
    df['is_fold1']    = (df['Fold'] == fold1).astype(int)
    df['unique_fold_score']       = df.groupby('cluster_num')['is_fold1'].transform('mean')

    # df.to_csv('/Users/steveabecassis/Desktop/PipelineTest/output_pipeline_1jfk/esm_output_analysis.csv',index=False)

    print(f"High score for fold {pdb1_name} is {round(df[df['Fold'] == pdb1_name].score_pdb1.max(),2)}")
    print(f"High score for fold {pdb2_name} is {round(df[df['Fold'] == pdb2_name].score_pdb2.max(), 2)}")
    print(f"Count for fold {pdb1_name} is {len(df[df['Fold'] == pdb1_name])}")
    print(f"Count for fold {pdb2_name} is {len(df[df['Fold'] == pdb2_name])}")

    print('Finish to run !')

    # df.to_parquet(f'/Users/steveabecassis/Desktop/pdb_output/tm_scores_{pdb1_name}_{pdb2_name}.parq')
    # df = pd.read_parquet('/Users/steveabecassis/Desktop/Kaib/AF_cluster/HDBSCAN15/tm_scores.parq')

    # df = pd.read_parquet('/Users/steveabecassis/Desktop/pdb_output/tm_scores_2JP1A_4QHF.parq')




