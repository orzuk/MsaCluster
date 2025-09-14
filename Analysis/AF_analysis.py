import os, sys
import pandas as pd
from tqdm import tqdm
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, ROOT)

import glob
from utils.align_utils import *
from utils.utils import list_protein_pairs
from config import *


# Compute TMscores of AF predictions for all clusters in a fold pair
def compute_AF_pred_tmscores(fold_pair):
#    fold_pair = pair_output_path.replace("Pipeline/", "", 1)
    print(f'Fold_pair: {fold_pair}')
    path = f'{DATA_DIR}/{fold_pair}'
    print(f'Path: {path}')
    if not os.path.exists(f'{path}/Analysis'):
        print("Mkdir: " + f'{path}/Analysis')
        os.mkdir(f'{path}/Analysis')

    folds = fold_pair.split("_")
    chains = [folds[0][-1],folds[1][-1]]
    AF_OUTPUT_DIR = path + '/' + PAIR_DIRS['AF']
#    af_pdb_files = os.listdir(AF_OUTPUT_DIR)
    af_pdb_files = [os.path.basename(p) for p in glob.glob(f"{AF_OUTPUT_DIR}/**/*.pdb", recursive=True)]

    res = []
    af_pdb_files = [i for i in af_pdb_files if 'pdb' in str(i)]
    for af_pdb_file in tqdm(af_pdb_files):
        if ('pdb' in str(af_pdb_file)) & ('gz' not in str(af_pdb_file)):
            score_pdb1 = compute_tmscore_align(f'{AF_OUTPUT_DIR}/{af_pdb_file}',
                                               f'{path}/{folds[0][:-1]}.pdb', chain2=chains[0])
            score_pdb2 = compute_tmscore_align(f'{AF_OUTPUT_DIR}/{af_pdb_file}',
                                               f'{path}/{folds[1][:-1]}.pdb', chain2=chains[1])
            res.append({'pdb_file':af_pdb_file,'score_pdb1':score_pdb1,'score_pdb2':score_pdb2})
    df_af = pd.DataFrame(res)
    df_af.loc[df_af.score_pdb1 > df_af.score_pdb2, 'Fold'] = folds[0][:-1]
    df_af.loc[df_af.score_pdb1 < df_af.score_pdb2, 'Fold'] = folds[1][:-1]
    df_af.sort_values(by='pdb_file', inplace=True)
    df_af['cluster_num'] = df_af['pdb_file'].apply(lambda x: x[11:14]).replace('ela','Query')
    df_af['is_fold1'] = (df_af['Fold'] == folds[0][:-1]).astype(int)
    df_af['unique_fold_score'] = df_af.groupby('cluster_num')['is_fold1'].transform('mean')
    df_af.drop(columns=['is_fold1'],inplace=True)
    df_af.to_csv(f'{path}/Analysis/df_af.csv')
    print('Finish AF TMscores computations!')


if __name__=='__main__':
    fold_pairs = list_protein_pairs()
    print("All pairs to compute AF TMscores: ", fold_pairs)

    for fold_pair in fold_pairs:
        compute_AF_pred_tmscores(fold_pair[0] + "_" + fold_pair[1])

    print("Finish all AF TMscores computations!")