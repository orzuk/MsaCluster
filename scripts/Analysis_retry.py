from argparse import ArgumentParser
from Bio.PDB import PDBParser, PDBIO, Select
from contact_map import ContactFrequency, ContactDifference
import numpy as np
import mdtraj as md
import os
import sys
import pandas as pd
from scipy.stats import spearmanr
import gzip
import subprocess
import re
from Bio import pairwise2
from Bio import PDB
from tqdm import tqdm

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, ROOT)

from config import DATA_DIR
from utils.align_utils import compute_tmscore_align
from utils.utils import save_org_cmaps


def unzip_files(folder_path):
    """Unzip all .gz files in folder_path."""
    for f in os.listdir(folder_path):
        if f.endswith('.gz'):
            gz_path = os.path.join(folder_path, f)
            out_path = gz_path[:-3]
            if not os.path.exists(out_path):
                with gzip.open(gz_path, 'rb') as f_in:
                    with open(out_path, 'wb') as f_out:
                        f_out.write(f_in.read())


if __name__ == '__main__':

    parser = ArgumentParser(description="Retry AF analysis for fold pairs")
    parser.add_argument("--pipeline_dir", default=None,
                        help="Path to Pipeline directory (default: DATA_DIR from config)")
    args = parser.parse_args()

    pipeline_dir = args.pipeline_dir or DATA_DIR

    '''
    Cmap analysis
    '''
    for fold_pair in tqdm(os.listdir(pipeline_dir)):
        try:
            pair_path = os.path.join(pipeline_dir, fold_pair)
            analysis_dir = os.path.join(pair_path, 'Analysis')
            if os.path.isdir(analysis_dir) and 'df_af.csv' in os.listdir(analysis_dir):
                continue
            if fold_pair == '4cmqB_4zt0C':
                continue
            print(f'Fold_pair: {fold_pair}')
            print(f'Path: {pair_path}')
            if not os.path.exists(os.path.join(pair_path, 'cmaps_pairs')):
                os.mkdir(os.path.join(pair_path, 'cmaps_pairs'))
            if not os.path.exists(analysis_dir):
                os.mkdir(analysis_dir)

            folds = fold_pair.split("_")
            fold1 = folds[0]
            fold2 = folds[1]
            save_org_cmaps(pair_path, fold1)
            save_org_cmaps(pair_path, fold2)

            unzip_files(folder_path=os.path.join(pair_path, 'AF_preds'))

            pdb_files = os.listdir(os.path.join(pair_path, 'AF_preds'))
            res = []
            pdb_files = [i for i in pdb_files if 'pdb' in str(i)]
            pdb_fold1_path = os.path.join(pair_path, 'chain_pdb_files', f'{fold1}.pdb')
            pdb_fold2_path = os.path.join(pair_path, 'chain_pdb_files', f'{fold2}.pdb')
            for pdb_file in pdb_files:
                try:
                    if ('pdb' in str(pdb_file)) and ('gz' not in str(pdb_file)):
                        score_pdb1 = compute_tmscore_align(
                            os.path.join(pair_path, 'AF_preds', pdb_file), pdb_fold1_path)
                        score_pdb2 = compute_tmscore_align(
                            os.path.join(pair_path, 'AF_preds', pdb_file), pdb_fold2_path)

                        temp = {'pdb_file':pdb_file,'score_pdb1':score_pdb1,'score_pdb2':score_pdb2}
                        res.append(temp)
                except Exception as e:
                    print(f"[retry] WARN: skipping {pdb_file}: {e}")
                    continue
            df_af = pd.DataFrame(res)
            df_af.loc[df_af.score_pdb1 > df_af.score_pdb2, 'Fold'] = fold1[:-1]
            df_af.loc[df_af.score_pdb1 < df_af.score_pdb2, 'Fold'] = fold2[:-1]
            df_af.sort_values(by='pdb_file', inplace=True)
            df_af['cluster_num'] = df_af['pdb_file'].apply(lambda x: x[11:14]).replace('ela','Query')
            df_af['is_fold1'] = (df_af['Fold'] == fold1[:-1]).astype(int)
            df_af['unique_fold_score'] = df_af.groupby('cluster_num')['is_fold1'].transform('mean')
            df_af.drop(columns=['is_fold1'],inplace=True)
            df_af.to_csv(os.path.join(analysis_dir, 'df_af.csv'))
            print('Finish !')
        except Exception as e:
            print(f"[retry] ERROR: {fold_pair}: {e}")
            continue









