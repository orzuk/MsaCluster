import os
import subprocess
import re
import pandas as pd
from tqdm import tqdm
from utils.protein_utils import *
from utils.align_utils import *


if __name__ == '__main__':
    res = []
    for foldpair in tqdm(os.listdir(DATA_DIR)):
        try:
            pdb_files = os.listdir(f'{DATA_DIR}/{foldpair}/chain_pdb_files/')
            path_fold1 = f'{DATA_DIR}/{foldpair}/chain_pdb_files/{pdb_files[0]}'
            path_fold2 = f'{DATA_DIR}/{foldpair}/chain_pdb_files/{pdb_files[1]}'

            tmscore = compute_tmscore_align(path_fold1,path_fold2)
            res.append({'fold_pair':foldpair,'tm_score':tmscore})
        except:
            continue

    import pandas as pd
df = pd.DataFrame(res)
