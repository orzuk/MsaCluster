import tmscoring
import os
import pandas as pd
from tqdm import tqdm


if __name__ == '__main__':
    pdb1 ='/Users/steveabecassis/PycharmProjects/AF_Cluster/data_sep2022/03_OtherFoldswitchers/02_selecase/2JP1A.pdb'
    # pdb1 = '/Users/steveabecassis/Desktop/pdb_output/4hdd.pdb'
    pdb1_name = pdb1[-8:-4]
    pdb2 ='/Users/steveabecassis/PycharmProjects/AF_Cluster/data_sep2022/03_OtherFoldswitchers/02_selecase/4QHF.pdb'
    # pdb2 = '/Users/steveabecassis/Desktop/pdb_output/2lep.pdb'
    pdb2_name = pdb2[-8:-4]

    # pdb_folder = '/Users/steveabecassis/Desktop/pdb_output/2LEP'
    pdb_folder = '/Users/steveabecassis/Desktop/4QHF'
    pdb_files = os.listdir(pdb_folder)
    res = []
    for pdb_file in tqdm(pdb_files):
        if 'pdb' in str(pdb_file):
            pdb_file_path = f'{pdb_folder}/{pdb_file}'
            score_pdb1 = tmscoring.get_tm(pdb_file_path, pdb1)
            score_pdb2 = tmscoring.get_tm(pdb_file_path, pdb2)
            pdb_file_path = '/Users/steveabecassis/Desktop/test.pdb'
            temp = {'pdb_file':pdb_file,'score_pdb1':score_pdb1,'score_pdb2':score_pdb2}
            res.append(temp)
    print('Finish !')


    df = pd.DataFrame(res)
    df.loc[df.score_pdb1 > df.score_pdb2 ,'Fold'] = pdb1_name
    df.loc[df.score_pdb1 < df.score_pdb2, 'Fold'] = pdb2_name
    print(f"High score for fold {pdb1_name} is {round(df[df['Fold'] == pdb1_name].score_pdb1.max(),2)}")
    print(f"High score for fold {pdb2_name} is {round(df[df['Fold'] == pdb2_name].score_pdb2.max(), 2)}")
    print(f"Count for fold {pdb1_name} is {len(df[df['Fold'] == pdb1_name])}")
    print(f"Count for fold {pdb2_name} is {len(df[df['Fold'] == pdb2_name])}")


    # df.to_parquet(f'/Users/steveabecassis/Desktop/pdb_output/tm_scores_{pdb1_name}_{pdb2_name}.parq')
    # df = pd.read_parquet('/Users/steveabecassis/Desktop/Kaib/AF_cluster/HDBSCAN15/tm_scores.parq')

    # df = pd.read_parquet('/Users/steveabecassis/Desktop/pdb_output/tm_scores_2JP1A_4QHF.parq')




