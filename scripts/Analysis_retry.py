from argparse import  ArgumentParser
from Bio.PDB import PDBParser, PDBIO, Select
from contact_map import ContactFrequency, ContactDifference
import numpy as np
import mdtraj as md
import os
import pandas as pd
from scipy.stats import spearmanr
import gzip
import subprocess
import re
from Bio import pairwise2
from Bio import PDB
from tqdm import tqdm
from Analysis import *



if __name__ == '__main__':

    '''
    Cmap analysis
    '''
    for fold_pair in tqdm(os.listdir('/Users/steveabecassis/Desktop/Pipeline')):
        try:
            if('df_af.csv' in os.listdir(f'/Users/steveabecassis/Desktop/Pipeline/{fold_pair}/Analysis')):
                continue
            if fold_pair == '4cmqB_4zt0C':
                continue
            print(f'Fold_pair: {fold_pair}')
            path = f'/Users/steveabecassis/Desktop/Pipeline/{fold_pair}'
            print(f'Path: {path}')
            if not os.path.exists(f'{path}/cmaps_pairs'):
                print("Mkdir: " + f'{path}/cmaps_pairs')
                os.mkdir(f'{path}/cmaps_pairs')
            if not os.path.exists(f'{path}/Analysis'):
                print("Mkdir: " + f'{path}/Analysis')
                os.mkdir(f'{path}/Analysis')

            folds = fold_pair.split("_")
            fold1 = folds[0]
            fold2 = folds[1]
            save_org_cmaps(f'{path}', fold1)
            save_org_cmaps(f'{path}', fold2)

            unzip_files(folder_path=f'{path}/AF_preds')

            pdb_files = os.listdir(f'{path}/AF_preds')
            res = []
            pdb_files = [i for i in pdb_files if 'pdb' in str(i)]
            pdb_fold1_path = f'{path}/chain_pdb_files/{fold1}.pdb'
            pdb_fold2_path = f'{path}/chain_pdb_files/{fold2}.pdb'
            for pdb_file in pdb_files:
                try:
                    if ('pdb' in str(pdb_file)) & ('gz' not in str(pdb_file)):
                        # score_pdb1 = tmscoring.get_tm(f'{path}/AF_preds/{pdb_file}', pdb_fold1_path)
                        # score_pdb2 = tmscoring.get_tm(f'{path}/AF_preds/{pdb_file}', pdb_fold2_path)

                        score_pdb1 = get_tmscore_align(f'{path}/AF_preds/{pdb_file}', pdb_fold1_path)
                        score_pdb2 = get_tmscore_align(f'{path}/AF_preds/{pdb_file}', pdb_fold2_path)

                        temp = {'pdb_file':pdb_file,'score_pdb1':score_pdb1,'score_pdb2':score_pdb2}
                        res.append(temp)
                except:
                    continue
            df_af = pd.DataFrame(res)
            df_af.loc[df_af.score_pdb1 > df_af.score_pdb2, 'Fold'] = fold1[:-1]
            df_af.loc[df_af.score_pdb1 < df_af.score_pdb2, 'Fold'] = fold2[:-1]
            df_af.sort_values(by='pdb_file', inplace=True)
            df_af['cluster_num'] = df_af['pdb_file'].apply(lambda x: x[11:14]).replace('ela','Query')
            df_af['is_fold1'] = (df_af['Fold'] == fold1[:-1]).astype(int)
            df_af['unique_fold_score'] = df_af.groupby('cluster_num')['is_fold1'].transform('mean')
            df_af.drop(columns=['is_fold1'],inplace=True)
            df_af.to_csv(f'{path}/Analysis/df_af.csv')
            print('Finish !')
        except:
            continue












