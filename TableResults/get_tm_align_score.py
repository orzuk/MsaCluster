import os
import subprocess
import re
import pandas as pd
from tqdm import tqdm
from utils.protein_utils import *
from utils.align_utils import *


# path_fold1 = '/Users/steveabecassis/Desktop/Pipeline/3t5oA_4a5wB/chain_pdb_files/3t5oA.pdb'
# path_fold2 = '/Users/steveabecassis/Desktop/Pipeline/3t5oA_4a5wB/chain_pdb_files/4a5wB.pdb'

if __name__ == '__main__':
    res = []
    for foldpair in tqdm(os.listdir(f'/Users/steveabecassis/Desktop/Pipeline/')):
        try:
            pdb_files = os.listdir(f'/Users/steveabecassis/Desktop/Pipeline/{foldpair}/chain_pdb_files/')
            path_fold1 = f'/Users/steveabecassis/Desktop/Pipeline/{foldpair}/chain_pdb_files/{pdb_files[0]}'
            path_fold2 = f'/Users/steveabecassis/Desktop/Pipeline/{foldpair}/chain_pdb_files/{pdb_files[1]}'

            tmscore = compute_tmscore_align(path_fold1,path_fold2)
            res.append({'fold_pair':foldpair,'tm_score':tmscore})
        except:
            continue

    import pandas as pd
df = pd.DataFrame(res)
