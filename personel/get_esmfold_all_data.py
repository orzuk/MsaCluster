import os
import pandas as pd
from tqdm import tqdm

from Bio.PDB import PDBParser, PPBuilder


def extract_sequence_from_pdb(file_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", file_path)
    # Extract sequence using Polypeptide Builder
    pp_builder = PPBuilder()
    sequence = ""
    for pp in pp_builder.build_peptides(structure):
        sequence += pp.get_sequence()
    return str(sequence)


if __name__ == "__main__":
    pdb_file = '/Users/steveabecassis/Desktop/Pipeline/3r9jC_2kxoA/output_esm_fold/ShallowMsa_000_3.pdb'
    df_esmfold_analysis = pd.read_csv(f'/Users/steveabecassis/Desktop/Pipeline_res/df_esmfold_analysis.csv')
    cols2keep =['fold_pair', 'fold', 'TMscore_fold1', 'TMscore_fold2', 'cluster_num','TM_mean_cluster_pdb1', 'TM_mean_cluster_pdb2']
    df_esmfold_analysis = df_esmfold_analysis[cols2keep]
    s = extract_sequence_from_pdb(pdb_file)
    from tqdm import tqdm
    tqdm.pandas()
    df_esmfold_analysis['sequence'] = df_esmfold_analysis.progress_apply(lambda x: extract_sequence_from_pdb(f'/Users/steveabecassis/Desktop/Pipeline/{x.fold_pair}/output_esm_fold/{x.fold}'),axis=1)
    fold_pairs = os.listdir('/Users/steveabecassis/Desktop/Pipeline')

    df_esmfold_analysis.to_csv('/Users/steveabecassis/Desktop/df_esmfold_all.csv',index=False)
