from tqdm import tqdm
from olds.Analysis import *
from PlotUtils import *


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



if __name__ == '__main__':
    '''
    Define the models outputs path on your local computer
    '''
    folder    = '/Users/steveabecassis/Desktop/Pipeline'
    files = os.listdir(folder)
    pattern = r'^[0-9a-zA-Z]{5}_[0-9a-zA-Z]{5}$'
    fold_pairs = [i for i in files if re.match(pattern, i)]
    fold_pair_errors  = []
    res = []
    for fold_pair in fold_pairs:
        if '.sh' in fold_pair:
            continue
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
            path_viz_maps = f'{folder}/{fold_pair}/output_cmap_esm/VizCmaps'

            for cmap in tqdm(cmaps):
                try:
                    if 'Shallow' in cmap:
                        cmap_pred = load_pred_cmap(f'{cmap[:-4]}')
                        cluster = cmap[-7:-4]
                    else:
                        continue
                    visualization_map_1_tol = load_cmap(f'{path_viz_maps}/msa_t__ShallowMsa_{cluster}_visualization_map_1_tol_0.npy')
                    visualization_map_2_tol = load_cmap(f'{path_viz_maps}/msa_t__ShallowMsa_{cluster}_visualization_map_2_tol_0.npy')
                    recall_only_fold1 = round(np.count_nonzero(visualization_map_1_tol == 1.75) / np.count_nonzero(cmap_only_pdb1 == 1),2)
                    recall_only_fold2 = round(np.count_nonzero(visualization_map_2_tol == 1.75) / np.count_nonzero(cmap_only_pdb2 == 1), 2)

                    # np.save(f'{cmaps_path}/VizCmaps/{cmap[:-4]}_visualization_map_1_tol_0.npy', visualization_map_1_tol)
                    # np.save(f'{cmaps_path}/VizCmaps/{cmap[:-4]}_visualization_map_2_tol_0.npy', visualization_map_2_tol)

                    res.append({'fold_pair':fold_pair,'File':cmap,'recall_only_fold1':recall_only_fold1,'recall_only_fold2':recall_only_fold2})
                except Exception as e:
                    print(e)
                    continue

        except Exception as e:
            fold_pair_errors.append(fold_pair)
            print(f'Error for {fold_pair}: {e}')
            continue

    final_df = pd.DataFrame(res)
    success_cmap_score = final_df.FoldPair.unique()
    # final_df.to_parquet(f'{folder}/cmap_exact_analysis_tol0_2510.parq')
    final_df['cluster_num'] = final_df.File.apply(lambda x:x[-7:-4])
    print('steve')
