from config import *
# import nglview as nv
import re

# sys.path.append(os.path.join(os.path.dirname(__file__), '.'))

#if not platform.system() == "Linux":  # Plotting doesn't work on unix
import pymol
from pymol import cmd  # , stored

from Bio import Align
from utils.phytree_utils import *
from scripts.MSA_Clust import *
from utils.utils import *
from utils.energy_utils import *
from matplotlib.colors import ListedColormap
import math
import pandas as pd


# Make all plots
# Input:
# pdbids - IDs of pdb proteins
# fasta_dir - directory of fasta file
# foldpair_id - ID of protein pair
# pdbchains - IDs of pdb proteins chains
# plot_tree_clusters
#
# Output:
# Three figures for each fold-switch pair:
# 1. Phylogenetic tree with matching scores to each of the fold switches
# 2. Cmap of each cluster and its match to the two folds
# 3. Two folds aligned
def make_foldswitch_all_plots(pdbids, fasta_dir, foldpair_id, pdbchains,
                              plot_tree_clusters= False, plot_contacts = True, global_plots = False):

    print("Plot for foldpair_id: " + foldpair_id)
    fasta_file_names = {pdbids[fold] + pdbchains[fold]: fasta_dir + "/" + foldpair_id + "/" + \
                        pdbids[fold] + pdbchains[fold] + '.fasta' for fold in range(2)}  # Added chain to file ID
    msa_pred_files = glob(fasta_dir + "/" + foldpair_id + "/output_cmap_esm/*.npy")
    n_cmaps = len(msa_pred_files)
    msa_files = glob(fasta_dir + "/" + foldpair_id + "/output_msa_cluster/*.a3m")
    msa_clusters = {file.split("\\")[-1][:-4]: read_msa(file) for file in msa_files}

    # Filter 'bad' families: too shallow alignments (no clusters), same PDB ID, other reasons??
    if len(msa_files) == 0:
        print("Shallow alignment! No MSA Clusters! Skipping family")
        return [None]*3
    if pdbids[0] == pdbids[1]:
        print("Same PDB-ID for both folds! Might be buggy!")  # Not supported yet!  Skipping family")

    # First load files
    seqs = {}
    for fold in range(2):
        print("Read fold-seq=", fasta_file_names[pdbids[fold] + pdbchains[fold]])
        with open(fasta_file_names[pdbids[fold] + pdbchains[fold]], "r") as text_file:
            seqs[pdbids[fold] + pdbchains[fold]] = text_file.read().split("\n")[1]
    pairwise_alignment = Align.PairwiseAligner().align(seqs[pdbids[0] + pdbchains[0]],
                                                       seqs[pdbids[1] + pdbchains[1]])

    if max([len(seqs[fold]) for fold in seqs]) > 1024:
        print( "Lengths: ", [len(seqs[fold]) for fold in seqs])
        print("Long sequence! Length = " + " > maximum supported length of 1024")
        return [None]*4

    try:  # read in text format or python format the pairwise predicted contact maps of msa transformer for each cluster
        msa_transformer_pred = {file.split("\\")[-1][9:-4]: np.genfromtxt(file) for file in msa_pred_files}
    except:
        msa_transformer_pred = {file.split("\\")[-1][9:-4]: np.load(file, allow_pickle=True) for file in msa_pred_files}  # remove 'shallow' string
#    print("msa_pred_files: ", msa_pred_files)
    try:
        print("Current dir: ", os.getcwd())
#        print("Try to extract true cmap from: ",
#              fasta_dir + "/" + foldpair_id + "/" + pdbids[fold] + pdbchains[fold] + "_pdb_contacts.npy")
        true_cmap = {pdbids[fold] + pdbchains[fold]: np.load(fasta_dir +  # problem with first !! # genfromtxt
                    "/" + foldpair_id + "/" + pdbids[fold] + pdbchains[fold] + "_pdb_contacts.npy").astype(int)
                     for fold in range(2)}
        print("Got true cmap!!!")
    except:
        print("Couldn't extract true cmap !!! ")
        true_cmap = {pdbids[fold] + pdbchains[fold]: np.load(fasta_dir +  # problem with first !!
                    "/" + foldpair_id + "/" + pdbids[fold] + pdbchains[fold] + "_pdb_contacts.npy",
                        allow_pickle=True).astype(int) for fold in range(2)}

    print("Get matching indices: pdbids", pdbids, "pdbchains", pdbchains)
    match_true_cmap, match_predicted_cmaps = \
        get_matching_indices_two_cmaps(pairwise_alignment, true_cmap, msa_transformer_pred)

    if plot_contacts:
        print("Plot Array Contact Map")
        plot_array_contacts_and_predictions(match_predicted_cmaps, match_true_cmap,
                                    fasta_dir + "/Results/Figures/Cmap_MSA/" + foldpair_id + '_all_clusters_cmap')

    shared_unique_contacts, shared_unique_contacts_metrics, contacts_united = match_predicted_and_true_contact_maps(
        match_predicted_cmaps, match_true_cmap)  # here number of cmaps is #clusters + 1

#                           match_predicted_cmaps}  # Why only shared?
    cluster_node_values = {ctype: (shared_unique_contacts_metrics["shared"][ctype]['long_P@L5'],
                                   shared_unique_contacts_metrics[pdbids[0] + pdbchains[0]][ctype]['long_P@L5'],
                                   shared_unique_contacts_metrics[pdbids[1] + pdbchains[1]][ctype]['long_P@L5']) for ctype in
                           match_predicted_cmaps}  # Why only shared?
    # load tree
    #        phytree_msa_str = "sbatch -o './Pipeline/" + foldpair_id + "/tree_reconstruct_for_" + foldpair_id + ".out' ./Pipeline/tree_reconstruct_params.sh " + foldpair_id  # Take one of the two !!! # ""./input/2qke.fasta 2qke
    #        print(phytree_msa_str)
    phytree_file = './Pipeline/' + foldpair_id + '/output_phytree/DeepMsa_tree.nwk'
    print("Load Biopython treefile: " + phytree_file)
    bio_tree = Phylo.read(phytree_file, "newick")  # This is different from write_newick_with_quotes !!!!
    print("Convert to ete3 tree:")
    ete_tree = convert_biopython_to_ete3(bio_tree)
#    print("Load treefile: " + phytree_file)
#    ete_tree = read_tree_ete(phytree_file)
#    for node in ete_tree.traverse():
#        print(node.name)
#        node.name = node.name.strip("'")
#    print("Set node cluster ids values and draw:")
#    print(fasta_dir + "/" + foldpair_id + "/output_msa_cluster/*.a3m")
    ete_leaves_cluster_ids = seqs_ids_to_cluster_ids(fasta_dir + "/" + foldpair_id + "/output_msa_cluster/*.a3m",
                                                     [n.name for n in ete_tree])
    print("Converted seq ids to cluster ids:")

    ete_leaves_node_values = {n.name: cluster_node_values[foldpair_id + '/output_cmap_esm/msa_t__Shallow' + ete_leaves_cluster_ids[n.name]]
                              for n in ete_tree if ete_leaves_cluster_ids[n.name] != 'p'}  # update to include matching two folds !!
    ete_leaves_node_values = pd.DataFrame(ete_leaves_node_values).T
    ete_leaves_node_values.columns = ["shared", pdbids[0] + pdbchains[0], pdbids[1] + pdbchains[1]]
#    print("Dump pickle:")
#    with open('tree_draw.pkl', 'wb') as f:  # Python 3: open(..., 'wb')
#        pickle.dump([phytree_file, fasta_dir + "/Results/Figures/PhyTree/" + foldpair_id + "_phytree", ete_leaves_node_values], f)
#    print("Node Values: ")
#    print(ete_leaves_node_values)
#    print(type(ete_leaves_node_values))
#    with open('bad_tree_and_msa.pkl', 'wb') as f:  # Python 3: open(..., 'rb')
#        pickle.dump([phytree_file, ete_leaves_node_values, fasta_dir + "/Results/Figures/PhyTree/" + foldpair_id + "_phytree"], f)
    print("plot tree cluster flag = ", plot_tree_clusters)
    if plot_tree_clusters:  # plot only clusters
        print("Plot Tree Clusters:")
##        cluster_node_values.pop('p')  # remove nodes without cluster
        cluster_node_values = pd.DataFrame(cluster_node_values).T  # convert to 3*[#clusters] pandas data-frame
        representative_cluster_leaves = unique_values_dict({n.name: ete_leaves_cluster_ids[n.name] for n in ete_tree if ete_leaves_cluster_ids[n.name] != 'p'} )

        # Compute tm scores of the predicted models of AF, ESM fold and the two structures
#        print("cluster_node_values", cluster_node_values, " Index: ", cluster_node_values.index)

        new_indices = [match.group() for s in cluster_node_values.index for match in re.finditer(r'M[sS][aA][a-zA-Z0-9_].*', s)]
        print("Parsed shortened indices: ", new_indices)

        tmscores_df = pd.DataFrame(index=new_indices,
                                   columns=['AF_TMscore_fold1', 'AF_TMscore_fold2', 'ESMF_TMscore_fold1', 'ESMF_TMscore_fold2'])  # modify index to exclude directory
        print("total cmaps: " + str(n_cmaps))
#        print("total node values: ", cluster_node_values.shape)
        AF_model_files = glob('Pipeline/' + foldpair_id + "/AF_preds/*Msa*model_1_*pdb")
        af_df = pd.read_csv(AF_MODEL_FILE, dtype=str)
        msa_trans_df = pd.read_csv(MSA_TRANS_MODEL_FILE, dtype=str)
        esmf_df = pd.read_csv(ESMF_MODEL_FILE, dtype=str)
#        print("read results files!")

        print("ALL AF MODEL FILES: ", AF_model_files)
        print("All of TMScors indices: ", tmscores_df.index)
        print("Their intersection with AF_model_files: ", set(tmscores_df.index).intersection(set(AF_model_files)))
        print("Missing from TMScors indices: ", set(AF_model_files) - set(tmscores_df.index))
        print("Missing from AF_model_files indices: ", set(tmscores_df.index) - set(AF_model_files))

#        print("AF_model_files= ", AF_model_files)
        for fold in range(2):
 #           ctr = 0
            for c in range(len(tmscores_df.index)): # cluster_node_values.index:  # loop over cluster names
                cur_AF_file = next((element for element in AF_model_files if tmscores_df.index[c] in element
                                    or (tmscores_df.index[c][-4:] == 'deep' and 'Deep' in element)), None)

                if cur_AF_file is None:
                    tmscores_df.iloc[c, fold] = 0
                    continue
#                    cluster_node_values.index)[c] # AF_model_files[AF_model_files == cluster_node_values.index[c]]
                print("Cluster: ", tmscores_df.index[c], " ; Cur AF file: ",  cur_AF_file, " fold-pair id=", foldpair_id)
                print("pdb_file1: ", 'Pipeline/' + foldpair_id + "/" + pdbids[fold] + '.pdb', 
                      " ; pdb_file2: ", cur_AF_file)
#                tmscores_df.iloc[c, fold] = compute_tmscore_align('Pipeline/' + foldpair_id + "/" + pdbids[fold] + '.pdb',
#                    cur_AF_file, pdbchains[fold], pdbchains[0])  # AF PREDICITON ALWAYS THE FIRST!!! # what chain to give the prediction? of first or second??

                if not 'eep' in cur_AF_file:
                    # Update AF Fold values
                    filtered_af_df = af_df[(af_df['fold_pair'] == foldpair_id) &
                                        (af_df['cluster_num'] == tmscores_df.index[c][4:])] #  &
#                                        ('unrelaxed_rank' in af_df['pdb_file'])]
#                    good_inds = 'unrelaxed_rank_001' in af_df['pdb_file']
#                    print("good inds:", good_inds)
#                    print("Head of filtered_af_df: ", filtered_af_df.head())
#                    tmscores_df.iloc[c, fold] = float(filtered_af_df['score_pdb' + str(fold+1)].iloc[0])
                    tmscores_df.iloc[c, fold]  = float(filtered_af_df.loc[filtered_af_df['pdb_file'].str.contains(
                        'unrelaxed_rank_00' + str(AF2_MODEL)), 'score_pdb' + str(fold+1)].values[0])
#                    print("tm-score=", tmscores_df.iloc[c, fold])
                    # Update ESM Fold values
                    filtered_esm_df = esmf_df[(esmf_df['fold_pair'] == foldpair_id) &
                                          (esmf_df['cluster_num'] == tmscores_df.index[c][4:])]
                    tmscores_df.iloc[c, fold+2] = float(filtered_esm_df['TM_mean_cluster_pdb' + str(fold+1)].iloc[0])


        # Get induced subtree
        clusters_subtree = extract_induced_subtree(phytree_file, representative_cluster_leaves)
        cluster_node_values.index = new_indices
        for n in clusters_subtree.iter_leaves(): # change name
            n.name = ete_leaves_cluster_ids[n.name]
        concat_scores = pd.concat([tmscores_df, cluster_node_values], ignore_index= True, axis=1)
        concat_scores.columns = ['TM-AF1', 'TM-AF2', 'TM-ESM1', 'TM-ESM2', 'RE-MSAT-COM', 'RE-MSAT1', 'RE-MSAT2']
 #        print("Concat scores:", concat_scores)
#        print("Clusters subtree", clusters_subtree)
        # save to pickle:
        with open('tree_clusters.pkl', 'wb') as f:  # Python 3: open(..., 'rb')
                pickle.dump([clusters_subtree, concat_scores,
                             fasta_dir + "/Results/Figures/PhyTreeCluster/" + foldpair_id + "_phytree_cluster",
                             tmscores_df, phytree_file, representative_cluster_leaves,
                             ete_leaves_node_values, ete_leaves_cluster_ids], f) #  true_cmap, msa_transformer_pred], f)
        print("Call tree visialization !!!! foldpair_id=",
              foldpair_id, " output dir: ", fasta_dir + "/Results/Figures/PhyTreeCluster/" + foldpair_id + "_phytree_cluster")
        visualize_tree_with_heatmap(clusters_subtree, concat_scores, fasta_dir + "/Results/Figures/PhyTreeCluster/" + foldpair_id + "_phytree_cluster")
    else:  # plot entire tree
        visualize_tree_with_heatmap(phytree_file, ete_leaves_node_values, fasta_dir + "/Results/Figures/PhyTree/" + foldpair_id + "_phytree")
        concat_scores = [] # temp

    # Collect :
##    print("Match predicted before compute_cmap_distances:")
##    print(match_predicted_cmaps)
    cmap_dists_vec = compute_cmap_distances(match_predicted_cmaps)  # msa_transformer_pred)
    seqs_dists_vec = np.mean(compute_seq_distances(msa_clusters))  # can take the entire sequence !
    num_seqs_msa_vec = len(seqs)

    if not platform.system() == "Linux":  # Plot two structures aligned (doesn't work in unix)
        align_and_visualize_proteins('Pipeline/' + foldpair_id + "/" + pdbids[0] + '.pdb',
                                     'Pipeline/' + foldpair_id + "/" + pdbids[1] + '.pdb',
                                     fasta_dir + "/Results/Figures/3d_struct/" + foldpair_id + "_3d_aligned.png", False)

    if global_plots:
        print("Make global plots!")
        global_pairs_statistics_plots(output_file=FIGURE_RES_DIR + "/fold_pair_scatter_plot.png")

    return cmap_dists_vec, seqs_dists_vec, num_seqs_msa_vec, concat_scores

#        print("Cmap dist: " + str(cmap_dists_vec[i]) + ", seq dist:" + str(seqs_dists_vec[i]))
#        break
# next plotP


# Align two pdb structures from two pdb files
def align_and_visualize_proteins(pdb_file1, pdb_file2, output_file, open_environment=True):
    """
    Align two protein structures and save the visualization as an image.

    Args:
    pdb_file1 (str): Path to the first PDB file.
    pdb_file2 (str): Path to the second PDB file.
    output_file (str): Path to save the output image.
    """

    if open_environment:  # Initialize PyMOL
#        pymol.cmd.set("quiet", 1)
        pymol.finish_launching(['pymol', '-cq'])  # '-c' for command line (no GUI)

    # Delete existing objects
    cmd.delete('all')

    # Load the PDB files
    cmd.load(pdb_file1, 'protein1')
    cmd.load(pdb_file2, 'protein2')

    # Align the structures
    cmd.align('protein2', 'protein1')

    # Set different colors for visualization
    cmd.color('red', 'protein1')
    cmd.color('blue', 'protein2')

    # Set the view (optional, you can customize this)
    cmd.zoom('all', buffer=10)

    # Save the image
    cmd.png(output_file)

    # Quit PyMOL
    if open_environment:  # Initialize PyMOL
        cmd.quit()

# Example usage
# align_and_visualize_proteins('path/to/pdb1.pdb', 'path/to/pdb2.pdb', 'output.png')


# Plot multiple contacts and predictions together
def plot_array_contacts_and_predictions(predictions, contacts, save_file=[]):
    """
    Plot multiple contacts and predictions together

    Parameters:
    predictions: Contact map predictions
    contacts: True Contact MAps (pair)
    save_file (str): Path to save the output image.
    """
    n_pred = len(predictions)
    n_row = math.ceil(math.sqrt(n_pred))  # *2
    if n_row * (n_row - 1) >= n_pred:  # *2
        n_col = n_row - 1
    else:
        n_col = n_row

    n_AA_aligned = len(contacts[next(iter(contacts))])  # number of aligned amino-acids in contacts
    fig, axes = plt.subplots(figsize=(18, 18), nrows=n_row, ncols=n_col, layout="compressed")
    #    print("Num cmaps: " + str(n_pred))
    #    print(axes.shape)
    #    fig, axes = plt.subplots(figsize=(18, 6), ncols=n_pred)
    ctr = 0
    #    for ax, name in zip(axes, PDB_IDS):
#    print("Contact lens:" + str(len(contacts)))
    recall = {}
    for name in predictions.keys():  # loop over predictions
        if n_col == 1:
            ax = axes[ctr]
        else:
            ax = axes[ctr // n_col, ctr % n_col]
        ctr = ctr + 1
#        print("Plotting prediction: " + name)  # + " -> true: " + true_name)
        recall[name] = plot_foldswitch_contacts_and_predictions(
            predictions[name], contacts, ax=ax, title=name, show_legend= ctr == 1)


    if len(save_file) > 0:  # save and close plot (enable automatic saving of multiple plots)
        plt.savefig(save_file + '.png')
        print("Save cmap fig: " + save_file + '.png')
    else:
        plt.show()

    # Find maximum and display plot for the best clusters:
#    print("Recall all clusters:", recall)
#    print("Recall no deep: ", {x:recall[x] for x in recall.keys() if "deep" not in x})
    best_recall_clusters = find_max_keys({x:recall[x] for x in recall.keys() if "deep" not in x})  # Need to exclude the deep alignment here !!!
    print("Best recall clusters: ", best_recall_clusters)
    best_cluster_ids = {x : int(best_recall_clusters[x][0][-3:]) for x in best_recall_clusters}
    print("Best recall clusters with IDS: ", best_recall_clusters, best_cluster_ids)
    fold_ids = list(contacts.keys())
    # Create a new figure for the second plot
    plt.figure(figsize=(10, 8))  # Adjust size as needed

    # Load energies:
    energy_dir = "Pipeline/output_deltaG"
    print("Load: ", os.path.join(energy_dir, f"deltaG_{fold_ids[0][:4]}.txt"))
    residue_energies_0 = read_energy_tuples(os.path.join(energy_dir, f"deltaG_{fold_ids[0][:4]}.txt"))
    residue_energies_1 = read_energy_tuples(os.path.join(energy_dir, f"deltaG_{fold_ids[1][:4]}.txt"))

    delta_energies, delta_energies_filtered = align_and_compare_residues(residue_energies_0, residue_energies_1, fold_ids[0][:4], fold_ids[1][:4])
#    print("delta_energies", delta_energies, " len=", len(delta_energies))
#    print("delta_energies_filtered", delta_energies_filtered, " len=", len(delta_energies_filtered))
    print("n_AA_aligned=", n_AA_aligned)

    delta_energies_filtered = np.array(delta_energies_filtered[:n_AA_aligned])  # Temp: need to fix alignment here!!!

    best_recall = plot_foldswitch_contacts_and_predictions( predictions=(predictions[best_recall_clusters[fold_ids[0]][0]],
                                                            predictions[best_recall_clusters[fold_ids[1]][0]]),
                                                            contacts=contacts, title="Best clusters", show_legend=True,
                                                            cluster_names= (str(best_cluster_ids[fold_ids[0]]), str(best_cluster_ids[fold_ids[1]])),
                                                            x_vector = delta_energies_filtered,
                                                            y_vector = delta_energies_filtered)
    print("best recall: ", best_recall)
    plt.savefig(save_file.replace('all', 'best'))
    plt.close()  # Close the figure to avoid reuse issues

"""Adapted from: https://github.com/rmrao/evo/blob/main/evo/visualize.py"""
def plot_contacts_and_predictions(
        predictions: Union[torch.Tensor, np.ndarray],
        contacts: Union[torch.Tensor, np.ndarray],
        ax: Optional[mpl.axes.Axes] = None,
        # artists: Optional[ContactAndPredictionArtists] = None,
        cmap: str = "Blues",
        ms: float = 0.5,  # size of circles
        title: Union[bool, str, Callable[[float], str]] = True,
        animated: bool = False,
) -> None:
    if isinstance(predictions, torch.Tensor):
        predictions = predictions.detach().cpu().numpy()
    if isinstance(contacts, torch.Tensor):
        contacts = contacts.detach().cpu().numpy()
    if ax is None:
        ax = plt.gca()

    seqlen = contacts.shape[0]
    relative_distance = np.add.outer(-np.arange(seqlen), np.arange(seqlen))
    bottom_mask = relative_distance < 0
    masked_image = np.ma.masked_where(bottom_mask, predictions)
    invalid_mask = np.abs(np.add.outer(np.arange(seqlen), -np.arange(seqlen))) < 6
    predictions = predictions.copy()
    predictions[invalid_mask] = float("-inf")

    topl_val = np.sort(predictions.reshape(-1))[-seqlen]
    pred_contacts = predictions >= topl_val
    true_positives = contacts & pred_contacts & ~bottom_mask
    false_positives = ~contacts & pred_contacts & ~bottom_mask
    other_contacts = contacts & ~pred_contacts & ~bottom_mask

    if isinstance(title, str):
        title_text: Optional[str] = title
    elif title:
        long_range_pl = compute_precisions(predictions, contacts, minsep=24)["P@L"].item()
        if callable(title):
            title_text = title(long_range_pl)
        else:
            title_text = f"Long Range P@L: {100 * long_range_pl:0.1f}"
    else:
        title_text = None

    img = ax.imshow(masked_image, cmap=cmap, animated=animated)  # Show main image
    oc = ax.plot(*np.where(other_contacts), "o", c="grey", ms=ms, label="other")[0]
    fp = ax.plot(*np.where(false_positives), "o", c="r", ms=ms, label="FP")[0]
    tp = ax.plot(*np.where(true_positives), "o", c="b", ms=ms, label="TP")[0]
    ti = ax.set_title(title_text) if title_text is not None else None
    # artists = ContactAndPredictionArtists(img, oc, fp, tp, ti)

    # Show second structure here!
    ax.legend(loc="upper left")
    ax.axis("square")
    ax.set_xlim([0, seqlen])
    ax.set_ylim([0, seqlen])
    save_flag = False  # add as input
#    print("Now plotting cmap: save-flag=", save_flag)
    if save_flag:
        plt.savefig('%s.pdf' % title, bbox_inches='tight')


# Plot contacts and predictions for TWO folds !!!
def plot_foldswitch_contacts_and_predictions(
        predictions: Union[torch.Tensor, np.ndarray],
        contacts: Union[torch.Tensor, np.ndarray],
        x_vector=None,
        y_vector=None,
        cluster_names=None,
        vector_cmap="viridis",
        ax: Optional[mpl.axes.Axes] = None,
        # artists: Optional[ContactAndPredictionArtists] = None,
        cmap: str = "gray_r",  # "Blues",
        ms: float = 3,
        hit_sign: str = 'x',
        title: Union[bool, str, Callable[[float], str]] = True,
        animated: bool = False,
        show_legend: bool = False,
        save_flag: bool = False):

    fold_ids = list(contacts.keys())

    # Handle single or double predictions
    if not isinstance(predictions, tuple):
        predictions = [predictions, predictions]
        two_predictions = False
    else:
        two_predictions = True
    for p in range(2):
        if isinstance(predictions[p], torch.Tensor):
            predictions[p] = predictions[p].detach().cpu().numpy()
    predictions = tuple(predictions)

    for fold in fold_ids: # always two contact maps
        if isinstance(contacts[fold], torch.Tensor):
            contacts[fold] = contacts[fold].detach().cpu().numpy()
    if ax is None:
        ax = plt.gca()

    if len(fold_ids) == 1:  # same PDB ID, duplicate
        fold_ids[1] = fold_ids[0]
    seqlen = contacts[fold_ids[0]].shape[0]  # take length of the first
    relative_distance = np.add.outer(-np.arange(seqlen), np.arange(seqlen))
    top_bottom_mask = {fold_ids[0]: relative_distance < 0, fold_ids[1]: relative_distance > 0}  # first is lower left, second is upper right
    #    masked_image = np.ma.masked_where(bottom_mask, predictions)
#    masked_image = np.ma.masked_where(top_bottom_mask[list(fold_ids)[0]], predictions[0])
    invalid_mask = np.abs(np.add.outer(np.arange(seqlen), -np.arange(seqlen))) < 6  # stripe around diagonal
    topl_val = [[], []]
    pred_contacts = [[], []]
    _, _, contacts_united = match_predicted_and_true_contact_maps({title: predictions[0]}, contacts)  # only contacts matter, not predictions
    predictions_copy = copy.deepcopy(predictions)
    for p in range(2):
        predictions_copy[p][invalid_mask] = float("-inf") # don't show predictions near diagonal
        topl_val[p] = np.sort(predictions_copy[p].reshape(-1))[-seqlen]
        pred_contacts[p] = predictions_copy[p] >= topl_val[p]
#        print("fold ", p , fold_ids[p], " toplval=", topl_val[p], " num prediction=", sum(pred_contacts[p]))

    true_positives, true_positives_unique, false_positives, other_contacts = {}, {}, {}, {}
    recall = {}

    p = 0
    for fold in fold_ids:
        true_positives[fold] = contacts[fold] & pred_contacts[p] & top_bottom_mask[fold]
        true_positives_unique[fold] = (np.transpose(contacts_united) == 2) & pred_contacts[p] & top_bottom_mask[fold]
        false_positives[fold] = ~contacts[fold] & pred_contacts[p] & top_bottom_mask[fold]
        other_contacts[fold] = contacts[fold] & ~pred_contacts[p] & top_bottom_mask[fold]
        recall[fold] = sum(true_positives[fold]) / ( sum(true_positives[fold]) + sum(other_contacts[fold] ) )
        p += 1

    if isinstance(title, str):
        title_text: Optional[str] = title
    elif title:
        long_range_pl0 = compute_precisions(predictions_copy[0], contacts, minsep=24)["P@L"].item()
        long_range_pl1 = compute_precisions(predictions_copy[1], contacts, minsep=24)["P@L"].item()
        if callable(title):
            title_text = title(long_range_pl0 + long_range_pl1)
        else:
            title_text = f"Long Range P@L: {100 * long_range_pl0:0.1f} , {100 * long_range_pl1:0.1f}"
    else:
        title_text = None

    # Start drawing
    # Check if vectors are provided
    include_vectors = x_vector is not None or y_vector is not None

    # Create the figure and axes
    if include_vectors:
        fig = plt.figure(figsize=(10, 10))
        gs = GridSpec(2, 2, width_ratios=[1, 0.03], height_ratios=[0.03, 1], figure=fig, wspace=0.05, hspace=0.05)
        ax = fig.add_subplot(gs[1, 0])

    colors = ['white', 'lightgray', 'darkgray']  # background, unique, shared
    custom_cmap = ListedColormap(colors)
    img = ax.imshow(contacts_united, cmap=custom_cmap, animated=animated)

    ms = ms * 50 / seqlen
    # Create legend entries with square markers for contacts
    shared_contacts = ax.scatter([], [], marker='s', c='lightgray', s=ms * 50, label='Shared Contacts')
    unique_contacts = ax.scatter([], [], marker='s', c='darkgray', s=ms * 50, label='Unique Contacts')

    categories = ["false_positives", "true_positives", "true_positives_unique"]
    colors = ["r", "b", "g"]
    relative_size = [1,1,1.5]
    labels = ["False Positives", "True Shared Positives", "True Unique Positives"]
    plots = []
    offset = -0.15  # Amount to shift circles to the left
    for i, category in enumerate(categories):
        x_coords, y_coords = np.where(locals()[category][fold_ids[0]])
        plots.append(ax.plot(x_coords + offset, y_coords, hit_sign, c=colors[i], ms=ms*relative_size[i], label=labels[i])[0])
        x_coords, y_coords = np.where(locals()[category][fold_ids[1]])
        plots.append(ax.plot(x_coords + offset, y_coords, hit_sign, c=colors[i], ms=ms*relative_size[i])[0])

    # Add a single colorbar for x_vector and y_vector
    if include_vectors:
        ax_xvec = fig.add_subplot(gs[0, 0], sharex=ax)
        ax_yvec = fig.add_subplot(gs[1, 1], sharey=ax)
        # Normalize vector values
        combined_vector = np.concatenate([x_vector.flatten(), y_vector.flatten()])
        norm = plt.Normalize(vmin=combined_vector.min(), vmax=combined_vector.max())

        # Plot x_vector heatmap

        ax_xvec.imshow(x_vector[np.newaxis, :], aspect="auto", cmap=vector_cmap, norm=norm)
#        original_aspect_x = ax_xvec.get_data_ratio()  # Get the current aspect ratio
#        ax_xvec.imshow(x_vector[np.newaxis, :], aspect=original_aspect_x *0.7, cmap=vector_cmap, norm=norm)  # Reduce width by 30%
        ax_xvec.axis("off")

        # Plot y_vector heatmap
        ax_yvec.imshow(y_vector[:, np.newaxis], aspect="auto", cmap=vector_cmap, norm=norm)
#        original_aspect_y = ax_yvec.get_data_ratio()
#        ax_yvec.imshow(y_vector[:, np.newaxis], aspect=original_aspect_y *0.7, cmap=vector_cmap, norm=norm)  # Reduce height by 30%
        ax_yvec.axis("off")

        # Add colorbar for vectors
        cbar_ax = fig.add_axes([0.93, 0.3, 0.02, 0.4])  # [left, bottom, width, height]
        cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=vector_cmap), cax=cbar_ax)
        cbar.set_label("ΔΔG (kcal/mol)", rotation=90, labelpad=-45, va='bottom')

        plt.tight_layout(rect=[0, 0, 0.9, 1])  # Leave space for the colorbar

    if show_legend:
        legend_elements = [
            Line2D([0], [0], marker='s', color='lightgray', label='Shared Contacts', markersize=10, linestyle='None'),
            Line2D([0], [0], marker='s', color='darkgray', label='Unique Contacts', markersize=10, linestyle='None'),
            Line2D([0], [0], marker='o', color='red', label='False Positives', markersize=6, linestyle='None'),
            Line2D([0], [0], marker='o', color='blue', label='True Shared Positives', markersize=6, linestyle='None'),
            Line2D([0], [0], marker='o', color='green', label='True Unique Positives', markersize=6, linestyle='None')]
        # Add the legend to the plot
        ax.legend(
            handles=legend_elements,
            loc="upper center",  # Place legend at the top center
            bbox_to_anchor=(0.5, -0.075),  # Adjust position just above the plot
            ncol=5,  # Reduce the number of columns
            frameon=False,  # Remove legend box
            fontsize=7,  # Smaller font size
            columnspacing=0.9,  # Reduce spacing between columns
            handletextpad=0.5,  # Reduce spacing between handles and text
            borderaxespad=0.2,  # Adjust padding between the legend and axes
        )

    ax.axis("square")
#    print("Recall is: ", {k:round(recall[k], 4) for k in recall}, " fold ids is: ", fold_ids)
    if cluster_names is None:
        ax.set_xlabel(fold_ids[0] + ", recall=" + str(round(recall[fold_ids[0]], 4)) , fontsize=14)
        ax.set_ylabel(fold_ids[1] + ", recall=" + str(round(recall[fold_ids[1]], 4)), fontsize=14)
    else:
        ax.set_xlabel(fold_ids[0] + ", recall=" + str(round(recall[fold_ids[0]], 4)) + ", cluster " + cluster_names[0], fontsize=14)
        ax.set_ylabel(fold_ids[1] + ", recall=" + str(round(recall[fold_ids[1]], 4)) + ", cluster " + cluster_names[0], fontsize=14)

    ax.set_xlim([0, seqlen])
    ax.set_ylim([0, seqlen])

#    predictions = save_predictions
    if save_flag:
        plt.savefig('%s.pdf' % title, bbox_inches='tight')
    return recall


# MAke global plots for all families
def global_pairs_statistics_plots(file_path=None, output_file="fold_pair_scatter_plot.png"):
    """
    Reads a CSV file and creates a scatter plot with specific values.

    Args:
        file_path (str): Path to the input CSV file.
        output_file (str): Path to save the scatter plot image.
    """
    # Load the data
    # load three output files
    af_df = pd.read_csv(AF_MODEL_FILE)
    msa_trans_df = pd.read_csv(MSA_TRANS_MODEL_FILE)
    esmf_df = pd.read_csv(ESMF_MODEL_FILE)

    print(af_df.columns.to_list())
    print(msa_trans_df.columns.to_list())
    print(esmf_df.columns.to_list())

#    esmf_df = pd.read_csv(file_path)

    # Group data by the first column
    grouped_esmf = esmf_df.groupby(esmf_df.columns[0])
    grouped_af = af_df.groupby(af_df.columns[0])
    grouped_msa_trans = msa_trans_df.groupby(msa_trans_df.columns[0])

    # Initialize the plot
    for model in ['AF', 'MSA_TRANS', 'ESMFold']:
        grouped = grouped_esmf if model == 'ESMFold' else grouped_af if model == 'AF' else grouped_msa_trans
        plt.figure(figsize=(10, 8))
        for fold_pair, group in grouped:


            if model == 'ESMFold':
                mean_x = group['TM_mean_cluster_pdb1'].mean()
                mean_y = group['TM_mean_cluster_pdb2'].mean()

                # Compute maximums
                max_x = group['TMscore_fold1'].max()
                max_y = group['TMscore_fold2'].max()
            elif model == 'AF':
                mean_x = group['score_pdb1'].mean()
                mean_y = group['score_pdb2'].mean()

                # Compute maximums
                max_x = group['score_pdb1'].max()
                max_y = group['score_pdb2'].max()
            else:    
                mean_x = group['recall_only_fold1'].mean()
                mean_y = group['recall_only_fold2'].mean()

                # Compute maximums
                max_x = group['recall_only_fold1'].max()
                max_y = group['recall_only_fold2'].max()

            # Plot the mean values with '+'
            # Convert dict_keys to a list for comparison
            first_key = list(grouped.groups.keys())[0]

            # Plot the mean values with '+'
            plt.scatter(mean_x, mean_y, marker='+', color='blue', label='Mean' if fold_pair == first_key else "")

            # Plot the maximum values with '*'
            plt.scatter(max_x, max_y, marker='*', color='red', label='Max' if fold_pair == first_key else "")

            plt.plot([mean_x, max_x], [mean_y, max_y], linestyle='--', color='gray', alpha=0.5)

            # Add fold-pair text for max values
            plt.text(max_x, max_y, fold_pair, fontsize=8, ha='right', color='black')

        # Customize the plot
        if model == 'cmap':
            plt.xlabel('Recall Fold 1')
            plt.ylabel('Recall Fold 2')
        else:
            plt.xlabel('TMscore fold1')
            plt.ylabel('TMscore fold2')
        plt.title(f'Scatter Plot of {model} vs. Fold Pairs')
        plt.legend(loc='upper left')
        plt.grid(True)

        # Save the plot to a file
        plt.tight_layout()
        plt.savefig(output_file[:-4] + f"_{model}.png")
        plt.close()
        print(f"{model} Scatter plot saved to {output_file}")

    # Now compare AF to ESMFold
    plt.figure(figsize=(10, 8))
    for fold_pair, group in grouped_esmf:
        mean_x = group['TM_mean_cluster_pdb1'].mean()
        mean_y = grouped_af.get_group(fold_pair)['score_pdb1'].mean()


        plt.scatter(mean_x, mean_y, marker='+', color='blue', label='Mean' if fold_pair == first_key else "")

    plt.xlabel("TMScore ESMF")
    plt.ylabel("TMScore AF")
    plt.title("TMScore Comparison ESMF vs. AF")
    plt.legend(loc='upper left')
    plt.tight_layout()


    # Plot cmap recalls
    grouped = msa_trans_df.groupby('fold_pair')
    plt.figure(figsize=(10, 8))

    for fold_pair, group in grouped:
        # Calculate mean values for recall_only_fold1 and recall_only_fold2
        mean_x = group['recall_only_fold1'].mean()
        mean_y = group['recall_only_fold2'].mean()

        # Find the max values and corresponding row
        max_row = group.loc[group[['recall_only_fold1', 'recall_only_fold2']].sum(axis=1).idxmax()]
        max_x = max_row['recall_only_fold1']
        max_y = max_row['recall_only_fold2']

        first_key = list(grouped.groups.keys())[0]

        # Plot the mean with '+'
        plt.scatter(mean_x, mean_y, marker='+', color='blue', label='Mean' if fold_pair == first_key else "")

        # Plot the max with '*' and add text annotation
        plt.scatter(max_x, max_y, marker='*', color='red', label='Max' if fold_pair == first_key else "")
        plt.text(max_x, max_y, fold_pair, fontsize=9, ha='right')

        # Add a faint line connecting mean and max points
        plt.plot([mean_x, max_x], [mean_y, max_y], linestyle='--', color='gray', alpha=0.5)

    plt.xlabel("Recall Fold 1")
    plt.ylabel("Recall Fold 2")
    plt.title("Recall Comparison Unique Contacts")
    plt.axhline(0, color='black', linewidth=0.8, linestyle='--')
    plt.axvline(0, color='black', linewidth=0.8, linestyle='--')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(loc='upper left')
    plt.tight_layout()

    # Save the plot to the specified file
    plt.savefig(output_file[:-4] + "_msa_trans.png")
    plt.close()



