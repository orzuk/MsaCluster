#from Analysis import *
from argparse import ArgumentParser
from tqdm import tqdm
from config import *
from utils.protein_utils import *
from utils.align_utils import *
import os
import pandas as pd
import gzip

if __name__ == '__main__':

    '''
    Cmap analysis
    '''

    parser = ArgumentParser()
    parser.add_argument("-input",help="Should be the pdb pair including the chain like 1eboE_5fhcJ")
    args = parser.parse_args()
    pair_output_path = args.input
    fold_pair = pair_output_path.replace("Pipeline/", "", 1)
    print(f'Fold_pair: {fold_pair}')
    path = f'{DATA_DIR}/{fold_pair}'
    print(f'Path: {path}')
    if not os.path.exists(f'{path}/Analysis'):
        print("Mkdir: " + f'{path}/Analysis')
        os.mkdir(f'{path}/Analysis')


    folds = pair_output_path.split("_")
    fold1 = folds[0]
    fold2 = folds[1]
    chain1 = fold1[-1]
    chain2 = fold2[-1]

#    unzip_files(folder_path=f'{path}/AF_preds')

    pdb_files = os.listdir(f'{path}/AF_preds')
    res = []
    pdb_files = [i for i in pdb_files if 'pdb' in str(i)]
#    pdb_fold1_path = f'{path}/chain_pdb_files/{fold1}.pdb'
#    pdb_fold2_path = f'{path}/chain_pdb_files/{fold2}.pdb'
    pdb_fold1_path = f'{path}/{fold1[:-1]}.pdb'  # exclude chain name
    pdb_fold2_path = f'{path}/{fold2[:-1]}.pdb'  # exclude chain name
    for pdb_file in tqdm(pdb_files):
        if ('pdb' in str(pdb_file)) & ('gz' not in str(pdb_file)):
            # score_pdb1 = tmscoring.get_tm(f'{path}/AF_preds/{pdb_file}', pdb_fold1_path)
            # score_pdb2 = tmscoring.get_tm(f'{path}/AF_preds/{pdb_file}', pdb_fold2_path)

            score_pdb1 = compute_tmscore_align(f'{path}/AF_preds/{pdb_file}', pdb_fold1_path, chain2=chain1)
            score_pdb2 = compute_tmscore_align(f'{path}/AF_preds/{pdb_file}', pdb_fold2_path, chain2=chain2)

            temp = {'pdb_file':pdb_file,'score_pdb1':score_pdb1,'score_pdb2':score_pdb2}
            res.append(temp)
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












