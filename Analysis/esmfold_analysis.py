import pandas as pd
import os, sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, ROOT)

from utils.align_utils import *
from utils.utils import list_protein_pairs
from config import *
from tqdm import tqdm

'''
To run this script on your local computer you have to specify where the TM align software is installed on your computer and adapt the function
compute_tmscore_align changing your path.
And to set in pipeline_folder your path where the results are on your machine
'''

# Compute TMscores of AF predictions for all clusters in a fold pair
def compute_esmfold_pred_tmscores(fold_pair):
    final_res       = []
    print(f'Fold_pair: {fold_pair}')
    path = f'{DATA_DIR}/{fold_pair}'
    print(f'Path: {path}')
    if not os.path.exists(f'{path}/Analysis'):
        print("Mkdir: " + f'{path}/Analysis')
        os.mkdir(f'{path}/Analysis')

    folds = fold_pair.split("_")
    chains = [folds[0][-1],folds[1][-1]]

    esm_pdb_files =  tqdm(os.listdir(f'{DATA_DIR}/{fold_pair}/output_esm_fold'))


    for fold_pred in esm_pdb_files:
        TMscore_fold1 = compute_tmscore_align(f'{DATA_DIR}/{fold_pair}/chain_pdb_files/{folds[0]}.pdb',
                                              f'{DATA_DIR}/{fold_pair}/output_esm_fold/{fold_pred}')
        TMscore_fold2 = compute_tmscore_align(f'{DATA_DIR}/{fold_pair}/chain_pdb_files/{folds[1]}.pdb',
                                              f'{DATA_DIR}/{fold_pair}/output_esm_fold/{fold_pred}')
        final_res.append({'fold_pair': fold_pair, 'fold': fold_pred,
                          'TMscore_fold1': TMscore_fold1, 'TMscore_fold2': TMscore_fold2})

    print("final_res: ", final_res)
    df                         = pd.DataFrame(final_res)
    print("df columns: " , df.columns)
    df['cluster_num']          = df.fold.apply(lambda x: x.split('_')[1] if 'Shallow' in x else x.split('_')[0])
    df['TM_mean_cluster_pdb1'] = df.groupby(['fold_pair','cluster_num'])['TMscore_fold1'].transform('mean')
    df['TM_mean_cluster_pdb2'] = df.groupby(['fold_pair','cluster_num'])['TMscore_fold2'].transform('mean')
    df['sample_count']         = df.groupby(['fold_pair','cluster_num'])['TMscore_fold1'].transform('count')
    df['is_fold_1']            = ((df.TM_mean_cluster_pdb1) >= (df.TM_mean_cluster_pdb2)).astype(int)
    df['is_fold_2']            = ((df.TM_mean_cluster_pdb1) < (df.TM_mean_cluster_pdb2)).astype(int)
    df['cluster_fold_1']       = df.groupby(['fold_pair','cluster_num'])['is_fold_1'].transform('mean')
    df['cluster_fold_2']       = df.groupby(['fold_pair','cluster_num'])['is_fold_2'].transform('mean')
    df.sort_values(by=['fold_pair','cluster_num'],inplace=True)
#    df.to_csv('./data/df_esmfold_analysis.csv',index=False)
    df.to_csv(f'{path}/Analysis/df_af.csv')


if __name__=='__main__':
    fold_pairs = list_protein_pairs()
    print("All pairs to compute ESMFold TMscores: ", fold_pairs)

    for fold_pair in fold_pairs:
        fold_pair_subdir = fold_pair[0] + "_" + fold_pair[1]
        if 'sh' in fold_pair_subdir:
            continue
        out_dir = f"{DATA_DIR}/{fold_pair_subdir}/output_esm_fold"
        if not os.path.isdir(out_dir) or len(os.listdir(out_dir)) == 0:
            print(f"[skip] No ESMFold outputs for {fold_pair_subdir} at {out_dir}")
            continue

        compute_esmfold_pred_tmscores(fold_pair_subdir)

    print("Finish all ESMFold TMscores computations!")


