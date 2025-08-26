import os, sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, ROOT)

from tqdm import tqdm
from Analysis.PlotUtils import *
from config import *
import pandas as pd
from utils.utils import *
from utils.align_utils import *
import argparse


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


def load_cmap(path):
    cmap = np.load(path)
    return cmap

def get_only_cmaps(cmap1,cmap2):
    diff_folds = cmap1 - cmap2
    only_fold1 = diff_folds.copy()
    only_fold1[only_fold1 == -1] = 0
    only_fold2 = diff_folds.copy()
    only_fold2[only_fold2 == 1] = 0
    only_fold2[only_fold2 == -1] = 1
    return only_fold1,only_fold2



# --- replace your existing main with this ---
if __name__ == "__main__":
    # Parse optional CLI pair
    parser = argparse.ArgumentParser(description="CMAP analysis for one pair or all pairs")
    parser.add_argument("pair", nargs="?", help="Optional pair like 1dzlA_5keqF")
    args = parser.parse_args()

    DEFAULT_PYCHARM_PAIR = "1dzlA_5keqF"

    if args.pair:  # explicit CLI pair
        fold_pairs = [pair_str_to_tuple(args.pair)]
        mode = "single-arg"
    elif is_pycharm():  # running from PyCharm, no CLI arg
        fold_pairs = [pair_str_to_tuple(DEFAULT_PYCHARM_PAIR)]
        mode = "pycharm-default"
    else:  # no arg, non-PyCharm → run all pairs
        fold_pairs = list_protein_pairs()
        mode = "all"

    print(f"[cmap_analysis] Mode={mode} | total pairs={len(fold_pairs)}")

    fold_pair_errors = []
    for fold_pair in fold_pairs:
        chains = (fold_pair[0][-1], fold_pair[1][-1])
        # ensure tuple form like ('1dzlA','5keqF')
        if isinstance(fold_pair, str):
            fold_pair = pair_str_to_tuple(fold_pair)

        fold_pair_subdir = f"{fold_pair[0]}_{fold_pair[1]}"
        res = []
        try:
            print("Set plot_tool for: ", fold_pair_subdir)
            plot_tool = PlotTool(folder=DATA_DIR, fold_pair=fold_pair_subdir)

            print("Extract sequences for: ", fold_pair_subdir)
            seq1 = extract_protein_sequence(f"{DATA_DIR}/{fold_pair_subdir}/{fold_pair[0][:-1]}.pdb",
                                            chain=chains[0], ca_only=True)
            seq2 = extract_protein_sequence(f"{DATA_DIR}/{fold_pair_subdir}/{fold_pair[1][:-1]}.pdb",
                                            chain=chains[1], ca_only=True)  # Add chain!

            print("Compute contact maps: ", fold_pair_subdir)
            cmap_pdb1, _ = pdb_to_contact_map(f"{DATA_DIR}/{fold_pair_subdir}/{fold_pair[0][:-1]}.pdb",
                                              chain=chains[0])
            cmap_pdb2, _ = pdb_to_contact_map(f"{DATA_DIR}/{fold_pair_subdir}/{fold_pair[1][:-1]}.pdb",
                                              chain=chains[1])

            print(f"[debug] cmap sizes: n1={cmap_pdb1.shape[0]}, n2={cmap_pdb2.shape[0]}")
            print(f"[debug] seq lengths: L1={len(seq1)}, L2={len(seq2)}")

            print("NEW!!! Align cmaps for: ", fold_pair_subdir)
#            cmap_aligned_pdb1, cmap_aligned_pdb2 = align_and_resize_contact_maps(
#                cmap_pdb1, cmap_pdb2, window_size=1, step_size=1)  # OLD!!!
            cmap_aligned_pdb1, cmap_aligned_pdb2, (idx1, idx2) = align_cmaps_by_sequence(
                cmap_pdb1, seq1, cmap_pdb2, seq2, mode="standard") #"standard") # 'blosum')

            print("Get only cmaps for: ", fold_pair_subdir)
            cmap_only_pdb1, cmap_only_pdb2 = get_only_cmaps(cmap_aligned_pdb1, cmap_aligned_pdb2)

            cmaps = os.listdir(f"{DATA_DIR}/{fold_pair_subdir}/output_cmap_esm")
            path_viz_maps = f"{DATA_DIR}/{fold_pair_subdir}/output_cmap_esm/VizCmaps"

            print("Loop on cmaps for: ", fold_pair_subdir)
            for cmap in tqdm(cmaps):
                try:
                    if (cmap == "VizCmaps") or "deep" in cmap:
                        continue
                    if "Shallow" in cmap:
                        cmap_pred = load_pred_cmap(f"{cmap[:-4]}")
                        cluster = cmap[-7:-4]
                    else:
                        continue

                    visualization_map_1_tol = load_cmap(
                        f"{path_viz_maps}/msa_t__ShallowMsa_{cluster}_visualization_map_1_tol_0.npy")
                    visualization_map_2_tol = load_cmap(
                        f"{path_viz_maps}/msa_t__ShallowMsa_{cluster}_visualization_map_2_tol_0.npy")
                    recall_only_fold1 = round(
                        np.count_nonzero(visualization_map_1_tol == 1.75)
                        / np.count_nonzero(cmap_only_pdb1 == 1), 2)
                    recall_only_fold2 = round(
                        np.count_nonzero(visualization_map_2_tol == 1.75)
                        / np.count_nonzero(cmap_only_pdb2 == 1), 2)

                    res.append(
                        {
                            "fold_pair": fold_pair,
                            "File": cmap,
                            "recall_only_fold1": recall_only_fold1,
                            "recall_only_fold2": recall_only_fold2,
                        }
                    )
                except Exception as e:
                    print(e)
                    continue

        except Exception as e:
            fold_pair_errors.append(fold_pair)
            print(f"Error for {fold_pair}: {e}")
            continue

        # Save data-frame for each protein pair
        df_cmap = pd.DataFrame(res)
        df_cmap.to_csv(f"{DATA_DIR}/{fold_pair_subdir}/Analysis/df_cmap.csv")

    print("Finish all MSA-Transformer CMAP Similarity computations!")

