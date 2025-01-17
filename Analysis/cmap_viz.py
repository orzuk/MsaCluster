import os

import prody as pr
import py3Dmol
import tempfile
import pandas as pd
from PlotUtils import *
from scipy.ndimage import binary_dilation
import matplotlib.patches as mpatches
from IPython.display import display, HTML
from Bio import PDB
from Bio import pairwise2
from Bio.PDB import Superimposer
import py3Dmol
import ipywidgets as widgets
from Bio.PDB import PDBParser
from Bio import pairwise2
import py3Dmol
from Bio import PDB, Align
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Polypeptide import three_to_one
from io import StringIO
import numpy as np
from tqdm import tqdm


def align_and_resize_contact_maps(cmap1, cmap2, window_size=10, step_size=1):
    """
    Align two contact maps and resize them to match the smaller map's dimensions.

    :param cmap1: First contact map (2D symmetric numpy array)
    :param cmap2: Second contact map (2D symmetric numpy array)
    :param window_size: Size of the sliding window for comparison
    :param step_size: Step size for sliding the window
    :return: Tuple of (aligned_cmap1, aligned_cmap2), both with dimensions of the smaller input map
    """

    # Determine which map is larger
    if cmap1.shape[0] * cmap1.shape[1] >= cmap2.shape[0] * cmap2.shape[1]:
        larger_cmap, smaller_cmap = cmap1, cmap2
    else:
        larger_cmap, smaller_cmap = cmap2, cmap1

    best_score = float('-inf')
    best_offset = (0, 0)

    # Find the best alignment
    for i in range(0, larger_cmap.shape[0] - smaller_cmap.shape[0] + 1, step_size):
        for j in range(0, larger_cmap.shape[1] - smaller_cmap.shape[1] + 1, step_size):
            window = larger_cmap[i:i + smaller_cmap.shape[0], j:j + smaller_cmap.shape[1]]
            score = np.sum(window * smaller_cmap)  # Simple dot product for similarity

            if score > best_score:
                best_score = score
                best_offset = (i, j)

    # Extract the aligned portion from the larger map
    aligned_larger = larger_cmap[best_offset[0]:best_offset[0] + smaller_cmap.shape[0],
                     best_offset[1]:best_offset[1] + smaller_cmap.shape[1]]

    # Ensure we return cmap1 and cmap2 in the correct order
    if cmap1.shape[0] * cmap1.shape[1] >= cmap2.shape[0] * cmap2.shape[1]:
        return aligned_larger, smaller_cmap
    else:
        return smaller_cmap, aligned_larger

def load_pred_cmap(fileName):
    cmap = np.load(f'{plot_tool.folder}/{plot_tool.fold_pair}/output_cmap_esm/{fileName}.npy')
    cmap[cmap > 0.3] = 1
    cmap[cmap <= 0.3] = 0
    return cmap

def get_only_cmaps(cmap1,cmap2):
    diff_folds = cmap1 - cmap2
    only_fold1 = diff_folds.copy()
    only_fold1[only_fold1 == -1] = 0
    only_fold2 = diff_folds.copy()
    only_fold2[only_fold2 == 1] = 0
    only_fold2[only_fold2 == -1] = 1
    return only_fold1,only_fold2


def process_array_tolerance(arr, tolerance=1):
    result = arr.copy()
    rows, cols = arr.shape
    # Find locations of 1.25
    locations_125 = np.where(arr == 1.25)
    for x, y in zip(*locations_125):
        # Define the neighborhood based on tolerance
        x_start = max(0, x - tolerance)
        x_end = min(rows, x + tolerance + 1)
        y_start = max(0, y - tolerance)
        y_end = min(cols, y + tolerance + 1)

        # If 1.5 is in the neighborhood, set the current cell to 1.75
        for i in range(x_start, x_end):
            for j in range(y_start, y_end):
                if arr[i, j] == 1.5:
                    result[i, j] = 1.75

    return result





if __name__ == '__main__':
    '''
    Define the models outputs path on your local computer
    '''
    folder    = '/Users/steveabecassis/Desktop/Pipeline'
    tolerance = 0
    final_df = pd.read_parquet(f'{folder}/cmap_exact_analysis.parq')

    res = []
    for fold_pair in tqdm(os.listdir(folder)):
        fold_pair = '4cmqB_4zt0C'
        if '.sh' in fold_pair:
            continue
        # fold_pair = '1svfC_4wsgC'
        try:
            fold1     = fold_pair.split('_')[0]
            fold2     = fold_pair.split('_')[-1]
            plot_tool = PlotTool(folder=folder, fold_pair=fold_pair)

            cmap_pdb1 = plot_tool.get_contact_map_from_pdb(f'{plot_tool.folder}/{plot_tool.fold_pair}/chain_pdb_files/{plot_tool.fold1}.pdb')
            cmap_pdb2 = plot_tool.get_contact_map_from_pdb(f'{plot_tool.folder}/{plot_tool.fold_pair}/chain_pdb_files/{plot_tool.fold2}.pdb')
            cmap_aligned_pdb1, cmap_aligned_pdb2 = align_and_resize_contact_maps(cmap_pdb1, cmap_pdb2, window_size=1, step_size=1)
            cmap_only_pdb1,cmap_only_pdb2 = get_only_cmaps(cmap_aligned_pdb1,cmap_aligned_pdb2)
            cmaps_path = f'{plot_tool.folder}/{plot_tool.fold_pair}/output_cmap_esm'
            cmaps = os.listdir(cmaps_path)
            os.makedirs(f'{cmaps_path}/VizCmaps', exist_ok=True)
            for cmap in tqdm(cmaps):
                try:
                    if 'Shallow' in cmap:
                        cmap_pred = load_pred_cmap(f'{cmap[:-4]}')
                    else:
                        continue
                    #     current_shape = cmap_pred.shape
                    #     target_shape = min(cmap_pdb1.shape[0], cmap_pdb2.shape[0])
                    #     if current_shape[0] != target_shape:
                    #         diff = target_shape - current_shape[0]  # 323 - 320 = 3
                    #         pad_before = diff // 2  # 1
                    #         pad_after = diff - pad_before  # 2
                    #         # Pad the array symmetrically
                    #         cmap_pred = np.pad(cmap_pred,
                    #                                   pad_width=((pad_before, pad_after), (pad_before, pad_after)),
                    #                                   mode='constant',
                    #                                   constant_values=0)


                        cmap_predicted_1_aligned,_  = align_and_resize_contact_maps(cmap_pred,cmap_only_pdb1,window_size=1, step_size=1)
                        cmap_predicted_2_aligned,_  = align_and_resize_contact_maps(cmap_pred,cmap_only_pdb2,window_size=1, step_size=1)

                        visualization_map_1 = cmap_aligned_pdb1 + 0.5*cmap_predicted_1_aligned + 0.25*cmap_only_pdb1
                        visualization_map_1_tol = process_array_tolerance(visualization_map_1, tolerance=0)

                        visualization_map_2 = cmap_aligned_pdb2 + 0.5*cmap_predicted_2_aligned + 0.25*cmap_only_pdb2
                        visualization_map_2_tol = process_array_tolerance(visualization_map_2, tolerance=0)

                        recall_only_fold1 = round(np.count_nonzero(visualization_map_1_tol == 1.75) / np.count_nonzero(cmap_only_pdb1 == 1),2)
                        recall_only_fold2 = round(np.count_nonzero(visualization_map_2_tol == 1.75) / np.count_nonzero(cmap_only_pdb2 == 1), 2)

                        np.save(f'{cmaps_path}/VizCmaps/{cmap[:-4]}_visualization_map_1_tol_0.npy', visualization_map_1_tol)
                        np.save(f'{cmaps_path}/VizCmaps/{cmap[:-4]}_visualization_map_2_tol_0.npy', visualization_map_2_tol)

                        res.append({'FoldPair':fold_pair,'File':cmap,'recall_only_fold1':recall_only_fold1,'recall_only_fold2':recall_only_fold2})
                except Exception as e:
                    print(e)
                    continue

        except Exception as e:
            print(f'Error for {fold_pair}: {e}')
            continue

    final_df = pd.DataFrame(res)

    # final_df.to_parquet(f'{folder}/cmap_exact_analysis_tol0_2310.parq')
