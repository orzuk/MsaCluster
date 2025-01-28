from config import *
# import nglview as nv

import platform
import re

# sys.path.append(os.path.join(os.path.dirname(__file__), '.'))

#if not platform.system() == "Linux":  # Plotting doesn't work on unix
import pymol
from pymol import cmd  # , stored

from Bio import Align
# from scripts import phytree_utils
from utils.phytree_utils import *
from scripts.MSA_Clust import *
from utils.utils import *
from utils.energy_utils import *
from matplotlib.colors import ListedColormap
import math


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
                              plot_tree_clusters= True, plot_contacts = False, global_plots = False):
    #    i = foldpair_ids.index(foldpair_id)
    #    cur_family_dir = fasta_dir + "/" + foldpair_id
#    plot("Start plotting inside make_foldswith_all_plots2!!!")
 #   pymol.finish_launching(['pymol', '-cq'])
#    plot("Start plotting inside make_foldswith_all_plots!!!")

     # Temp decide what to plot here

    print("Plot for foldpair_id: " + foldpair_id)
    fasta_file_names = {pdbids[fold] + pdbchains[fold]: fasta_dir + "/" + foldpair_id + "/" + \
                        pdbids[fold] + pdbchains[fold] + '.fasta' for fold in range(2)}  # Added chain to file ID
    #    msa_file = fasta_dir + "/" + foldpair_id + "/output_get_msa/DeepMsa.a3m"
    #    MSA = read_msa(msa_file)  # AlignIO.read(open(msa_file), "fasta")
    msa_pred_files = glob(fasta_dir + "/" + foldpair_id + "/output_cmap_esm/*.npy")
    n_cmaps = len(msa_pred_files)
#    n_cmaps = min(3, n_cmaps)  # temp for debug !!
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

#    print("Full msa_trans keys: ")
#    print({file.split("\\")[-1] for file in msa_pred_files})
    try:  # read in text format or python format the pairwise predicted contact maps of msa transformer for each cluster
        msa_transformer_pred = {file.split("\\")[-1][9:-4]: np.genfromtxt(file) for file in msa_pred_files}
    except:
        msa_transformer_pred = {file.split("\\")[-1][9:-4]: np.load(file, allow_pickle=True) for file in msa_pred_files}  # remove 'shallow' string
#    print("msa_pred_files: ", msa_pred_files)
    try:
        print("Current dir: ", os.getcwd())
        print("Try to extract true cmap from: ",
              fasta_dir + "/" + foldpair_id + "/" + pdbids[fold] + pdbchains[fold] + "_pdb_contacts.npy")
        true_cmap = {pdbids[fold] + pdbchains[fold]: np.load(fasta_dir +  # problem with first !! # genfromtxt
                    "/" + foldpair_id + "/" + pdbids[fold] + pdbchains[fold] + "_pdb_contacts.npy").astype(int)
                     for fold in range(2)}
#        print(true_cmap[pdbids[0] + pdbchains[0]])
#        print(true_cmap[pdbids[1] + pdbchains[1]])
    except:
        print("Couldn't extract true cmap !!! ")
        true_cmap = {pdbids[fold] + pdbchains[fold]: np.load(fasta_dir +  # problem with first !!
                    "/" + foldpair_id + "/" + pdbids[fold] + pdbchains[fold] + "_pdb_contacts.npy",
                        allow_pickle=True).astype(int) for fold in range(2)}
#    print("All predicted MSA transformer files:")
#    print(msa_transformer_pred.keys())
#    print(msa_transformer_pred.values())
#    print("Exit function! shapes:")
#    print(msa_transformer_pred['p'].shape)
#    print("Exit function! shapes2:")
#    print(msa_transformer_pred['Msa_000'].shape)  # Show shape. Should be square form: n*n residues

#    print("Dump pickle for debug:")
#    with open('debug_match_index.pkl', 'wb') as f:  # Python 3: open(..., 'wb')
#        pickle.dump([pairwise_alignment], f) #  true_cmap, msa_transformer_pred], f)

    print("Get matching indices:")
    match_true_cmap, match_predicted_cmaps = \
        get_matching_indices_two_cmaps(pairwise_alignment, true_cmap, msa_transformer_pred)

#    print("Pairwise alignment: ", pairwise_alignment)
#    print("True cmap: ", true_cmap)
#    print("msa transfomers keys: ", msa_transformer_pred.keys())
#    print("Exit function! match true:")
#    print(match_true_cmap)
#    print("Match predicted: ", match_predicted_cmaps)
    print("pdbids", pdbids, "pdbchains", pdbchains)
#    for f in range(2):
#        plt.imshow(1-match_true_cmap[pdbids[f]+pdbchains[f]].astype(int),  cmap='gray', vmin=0, vmax=1)
#        plt.axis('off')  # Turn off axes for a cleaner image
#        plt.show()
#
#        # Save the array to an image file
#        output_file = "temp_binary_" + pdbids[f]+pdbchains[f] + "_cmap.png"
#        plt.imsave(output_file, 1-match_true_cmap[pdbids[f]+pdbchains[f]].astype(int), cmap='gray')
#        print(f"Image saved to {output_file}")

#    differences = np.argwhere(match_true_cmap[pdbids[0]+pdbchains[0]] != match_true_cmap[pdbids[1]+pdbchains[1]])
#    print("Indices where the arrays are different:")
#    print(differences)
#    print("Means of fold differences: ", np.mean(match_true_cmap[pdbids[0]+pdbchains[0]] ), np.mean(match_true_cmap[pdbids[1]+pdbchains[1]] ))


    if plot_contacts:
        print("Plot Array")
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
    if plot_tree_clusters:  # plot only clusters
##        cluster_node_values.pop('p')  # remove nodes without cluster
        cluster_node_values = pd.DataFrame(cluster_node_values).T  # convert to 3*[#clusters] pandas data-frame
        representative_cluster_leaves = unique_values_dict({n.name: ete_leaves_cluster_ids[n.name] for n in ete_tree if ete_leaves_cluster_ids[n.name] != 'p'} )

        # Compute tm scores of the predicted models of AF, ESM fold and the two structures
#        print("cluster_node_values", cluster_node_values, " Index: ", cluster_node_values.index)

        new_indices = [match.group() for s in cluster_node_values.index for match in re.finditer(r'M[sS][aA][a-zA-Z0-9_].*', s)]
        print("Parsed shortenede indices: ", new_indices)

        tmscores_df = pd.DataFrame(index=new_indices,
                                   columns=['AF_TMscore_fold1', 'AF_TMscore_fold2', 'ESMF_TMscore_fold1', 'ESMF_TMscore_fold2'])  # modify index to exclude directory
        print("total cmaps: " + str(n_cmaps))
        print("total node values: ", cluster_node_values.shape)
        AF_model_files = glob('Pipeline/' + foldpair_id + "/AF_preds/*Msa*model_1_*pdb")
        ESM_model_files = ''  # don't have them currently
        af_df = pd.read_csv(AF_MODEL_FILE, dtype=str)
        msa_trans_df = pd.read_csv(MSA_TRANS_MODEL_FILE, dtype=str)
        esmf_df = pd.read_csv(ESMF_MODEL_FILE, dtype=str)
#        print("read results files!")

#        print("AF_model_files= ", AF_model_files)
        for fold in range(2):
 #           ctr = 0
            for c in range(len(tmscores_df.index)): # cluster_node_values.index:  # loop over cluster names
                cur_AF_file = next((element for element in AF_model_files if tmscores_df.index[c] in element
                                    or (tmscores_df.index[c][-4:] == 'deep' and 'Deep' in element)), None)

#                    cluster_node_values.index)[c] # AF_model_files[AF_model_files == cluster_node_values.index[c]]
                print("Cluster: ", tmscores_df.index[c], " ; Cur AF file: ",  cur_AF_file, " fold-pair id=", foldpair_id)

                # Save to debug just the function "compute_tmscore"
#                true_pdb_file = 'Pipeline/' + foldpair_id + "/" + pdbids[fold] + '.pdb'
#                with open('compute_tmscore.pkl', 'wb') as f:  # Python 3: open(..., 'rb')
#                    pickle.dump([fold, true_pdb_file, cur_AF_file, pdbchains], f)
#                with open('compute_tmscore.pkl', 'rb') as f:  # Python 3: open(..., 'rb')
#                    fold, true_pdb_file, cur_AF_file, pdbchains = pickle.load(f)
#                tmscores_mat[fold, ctr] = compute_tmscore(true_pdb_file, cur_AF_file, pdbchains[fold], pdbchains[0])
#                tmscores_mat[fold, ctr] = compute_tmscore('Pipeline/' + foldpair_id + "/" + pdbids[fold] + '.pdb',
#                print("Compute tm score for: ", 'Pipeline/' + foldpair_id + "/" + pdbids[fold] + '.pdb',
#                      cur_AF_file, pdbchains[fold], pdbchains[0])
#                print("first chain: ", pdbchains[fold])
#                print("second chain: ", pdbchains[0], " ; finished")
                # Why recalculate tm score each time? could load precomputed scores
                tmscores_df.iloc[c, fold] = compute_tmscore('Pipeline/' + foldpair_id + "/" + pdbids[fold] + '.pdb',
                    cur_AF_file, pdbchains[fold], pdbchains[0])  # AF PREDICITON ALWAYS THE FIRST!!! # what chain to give the prediction? of first or second??

                if not 'eep' in cur_AF_file:
                    # Update ESM Fold values
#                    print("esm clus num", esmf_df['cluster_num'], "tmscoresindex", tmscores_df.index[c][4:])
                    filtered_df = esmf_df[(esmf_df['fold_pair'] == foldpair_id) &
                                          (esmf_df['cluster_num'] == tmscores_df.index[c][4:])]
#                    print("ESMf Filtered df: ", filtered_df, " key: ", 'TM_mean_cluster_pdb' + str(fold+1))
                    tmscores_df.iloc[c, fold+2] = float(filtered_df['TM_mean_cluster_pdb' + str(fold+1)].iloc[0])

                    # Update MSA-Transformer recall values

#                    tmp_clust_num = str(int(tmscores_df.index[c][4:]))
#                    print("Cluster num: ", tmscores_df.index[c][4:], " and ", tmp_clust_num)
#                    filtered_df = msa_trans_df[(msa_trans_df['FoldPair'] == foldpair_id) &
#                                               (msa_trans_df['cluster_num'] == tmscores_df.index[c][4:])]
#                    print("MSA_Trans Filtered df: ", filtered_df, " key: ", 'recall_only_fold' + str(fold+1))

#                    print("Filtered values: ", filtered_df['recall_only_fold' + str(fold+1)])
#                    print("Filtered values again: ", filtered_df['recall_only_fold1'])
#                    print("Filtered values again zero: ", filtered_df['recall_only_fold1'].iloc[0])
#                    print("Filtered values again zero zero: ", filtered_df['recall_only_fold1'][0])

#                    tmscores_df.iloc[c, fold+5] = filtered_df['recall_only_fold' + str(fold+1)].iloc[0]

        # Get induced subtree
#        print("phytree_file=", phytree_file)
#        print("representative_cluster_leaves", representative_cluster_leaves)
        clusters_subtree = extract_induced_subtree(phytree_file, representative_cluster_leaves)
#        print("New cluster_node_values index: ", cluster_node_values.index)

        cluster_node_values.index = new_indices
#        print("New cluster_node_values:", cluster_node_values)
#        print("cluster_subtree: ", clusters_subtree)
        for n in clusters_subtree.iter_leaves(): # change name
            n.name = ete_leaves_cluster_ids[n.name]
#        print("Now renamed cluster_subtree:")
#        print(clusters_subtree)

#        print("Cluster node values:")
#        print(cluster_node_values)
#        print("TMScores mat: ")
#        print(tmscores_df)
#        print("Types node_values, tmscores_mat:")
#        print(type(cluster_node_values))
#        print(type(tmscores_df))
#        print("Now concatenate shapes: cluster_node, tmscores: " + str(cluster_node_values.shape) + ", " + str(tmscores_df.shape))
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
        visualize_tree_with_heatmap(clusters_subtree, concat_scores, fasta_dir + "/Results/Figures/PhyTreeCluster/" + foldpair_id + "_phytree_cluster")
    else:  # plot entire tree
        visualize_tree_with_heatmap(phytree_file, ete_leaves_node_values, fasta_dir + "/Results/Figures/PhyTree/" + foldpair_id + "_phytree")
        concat_scores = [] # temp

#    with open('bad_tree_and_msa.pkl', 'rb') as f:  # Python 3: open(..., 'rb')
#        phytree_file, ete_leaves_node_values, outfile = pickle.load(f)

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
        global_pairs_statistics_plots(output_file="Pipeline/Results/Figures/fold_pair_scatter_plot.png")

    return cmap_dists_vec, seqs_dists_vec, num_seqs_msa_vec, concat_scores

#        print("Cmap dist: " + str(cmap_dists_vec[i]) + ", seq dist:" + str(seqs_dists_vec[i]))
#        break
# next plotP

import pymolPy3


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
    contacts_ids = contacts.keys()
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



    ##        for true_name in contacts_ids: # loop over two folds
    ##            print("Plotting prediction: " + name + " -> true: " + true_name)
    ##           plot_contacts_and_predictions(
    ##                predictions[name], contacts[true_name], ax=ax, title = name)
    #            prediction, target, ax=ax, title = lambda prec: f"{name}: Long Range P@L: {100 * prec:0.1f}")
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
#    fig = plt.gcf()  # Get current figure
#    ax = plt.gca()  # Get current Axes
#    print("New labels: ", ax.get_xlabel() + ' cluster ' + str(best_cluster_ids[fold_ids[0]]),
#          ax.get_ylabel() + ' cluster ' + str(best_cluster_ids[fold_ids[1]]))
#    ax.set_xlabel(ax.get_xlabel() + ' cluster ' + str(best_cluster_ids[fold_ids[0]]))
#    ax.set_ylabel(ax.get_ylabel() + ' cluster ' + str(best_cluster_ids[fold_ids[1]]))
#    print("New NEW labels: ", ax.get_xlabel() ,  ax.get_ylabel() )
#    plt.draw()
#    ax.figure.canvas.draw()
    print("best recall: ", best_recall)
    plt.savefig(save_file.replace('all', 'best'))
    plt.close()  # Close the figure to avoid reuse issues

#    print("Save best cmap fig: " + save_file + '_best_clusters_cmap.png')

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
        hit_sign: str = 'o',
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

    # Plot the contact maps
#    maps = [("True Contacts " + fold_ids[0], contacts[fold_ids[0]]),
#            ("True Contacts " + fold_ids[1], contacts[fold_ids[1]]),
#            ("Predicted Contacts 1", pred_contacts[0]),
#            ("Predicted Contacts 2", pred_contacts[1])]
#
#    plot_2by2 = False
#    if plot_2by2:
#        # Create a 2x2 figure
#        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 10), layout="compressed")
#        for ax, (tmp_title, contact_map) in zip(axes.ravel(), maps):
#            ax.imshow(contact_map, cmap="gray_r", animated=animated)  # Black (1) and White (0)
#            #        ax.plot([0, contact_map.shape[0]], [0, contact_map.shape[1]], c="red", linestyle="-", linewidth=0.5)
#            #        ax.plot(*np.where(false_positives[fold_ids[0]]), hit_sign, c="r", ms=ms, label="FP")[0]
#            ax.set_title(tmp_title)
#            ax.axis("off")
##            ax.invert_xaxis()  # Flip the x-axis direction
#            ax.invert_yaxis()  # Flip the x-axis direction
#
#        ax.axis("square")
#
#        # Adjust layout and save to file
#        plt.tight_layout()
#        plt.savefig('cmap_vs_predictions.png')
#        plt.close()

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
#        gs = GridSpec(2, 2, width_ratios=[1, 0.03], height_ratios=[0.03, 1], figure=fig)
        gs = GridSpec(2, 2, width_ratios=[1, 0.03], height_ratios=[0.03, 1], figure=fig, wspace=0.05, hspace=0.05)

        ax = fig.add_subplot(gs[1, 0])
        ax_xvec = fig.add_subplot(gs[0, 0], sharex=ax)
        ax_yvec = fig.add_subplot(gs[1, 1], sharey=ax)

#        fig = plt.figure(figsize=(8, 8))
#        gs = GridSpec(2, 2, width_ratios=[1, 30], height_ratios=[30, 1], figure=fig)
#        ax = fig.add_subplot(gs[0, 1])

    colors = ['white', 'lightgray', 'darkgray']  # background, unique, shared
    custom_cmap = ListedColormap(colors)
    img = ax.imshow(contacts_united, cmap=custom_cmap, animated=animated)

    ms = ms * 50 / seqlen
    # Create legend entries with square markers for contacts
    shared_contacts = ax.scatter([], [], marker='s', c='lightgray', s=ms * 50, label='Shared Contacts')
    unique_contacts = ax.scatter([], [], marker='s', c='darkgray', s=ms * 50, label='Unique Contacts')

    categories = ["false_positives", "true_positives", "true_positives_unique"]
    colors = ["r", "b", "g"]
    labels = ["False Positives", "True Shared Positives", "True Unique Positives"]
    plots = []
    offset = -0.15  # Amount to shift circles to the left
    for i, category in enumerate(categories):
        x_coords, y_coords = np.where(locals()[category][fold_ids[0]])
        plots.append(ax.plot(x_coords + offset, y_coords, hit_sign, c=colors[i], ms=ms, label=labels[i])[0])
        x_coords, y_coords = np.where(locals()[category][fold_ids[1]])
        plots.append(ax.plot(x_coords + offset, y_coords, hit_sign, c=colors[i], ms=ms)[0])

#    if x_vector is not None:
#        ax_xvec = fig.add_subplot(gs[1, 1])
#        x_vector = np.array(x_vector).reshape(1, -1)  # Ensure it's a 1-row heatmap
#        img_xvec = ax_xvec.imshow(x_vector, cmap=vector_cmap, aspect='auto')
#        ax_xvec.axis("off")
#
#        # Add y_vector heatmap
#    if y_vector is not None:
#        ax_yvec = fig.add_subplot(gs[0, 0])
#        y_vector = np.array(y_vector).reshape(-1, 1)  # Ensure it's a 1-column heatmap
#        img_yvec = ax_yvec.imshow(y_vector, cmap=vector_cmap, aspect='auto')
#        ax_yvec.axis("off")

    # Add a single colorbar for x_vector and y_vector
    if include_vectors:
        # Normalize vector values
        combined_vector = np.concatenate([x_vector.flatten(), y_vector.flatten()])
        norm = plt.Normalize(vmin=combined_vector.min(), vmax=combined_vector.max())

        # Plot x_vector heatmap
        # x_vector = np.array(x_vector).reshape(1, -1)  # Ensure it's a 1-row heatmap
        ax_xvec.imshow(x_vector[np.newaxis, :], aspect="auto", cmap=vector_cmap, norm=norm)
        ax_xvec.axis("off")

        # Plot y_vector heatmap
        # y_vector = np.array(y_vector).reshape(-1, 1)  # Ensure it's a 1-column heatmap
        ax_yvec.imshow(y_vector[:, np.newaxis], aspect="auto", cmap=vector_cmap, norm=norm)
        ax_yvec.axis("off")

        # Add colorbar for vectors
        cbar_ax = fig.add_axes([0.93, 0.3, 0.02, 0.4])  # [left, bottom, width, height]
        cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=vector_cmap), cax=cbar_ax)
#        cbar.set_label("ΔΔG (kcal/mol)", fontsize=10)
        cbar.set_label("ΔΔG (kcal/mol)", rotation=90, labelpad=-40, va='bottom')

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
    grouped = esmf_df.groupby(esmf_df.columns[0])

    # Initialize the plot
    plt.figure(figsize=(10, 8))
    for fold_pair, group in grouped:
        # Compute means
        mean_x = group['TM_mean_cluster_pdb1'].mean()
        mean_y = group['TM_mean_cluster_pdb2'].mean()

        # Compute maximums
        max_x = group['TMscore_fold1'].max()
        max_y = group['TMscore_fold2'].max()

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
    plt.xlabel('TMscore fold1')
    plt.ylabel('TMscore fold2')
    plt.title('Scatter Plot of ESMFold vs. Fold Pairs')
    plt.legend(loc='upper left')
    plt.grid(True)

    # Save the plot to a file
    plt.tight_layout()
    plt.savefig(output_file[:-4] + "_ESMFold.png")
    plt.close()
    print(f"ESMFold Scatter plot saved to {output_file}")


    # Plot cmap recalls
    grouped = msa_trans_df.groupby('FoldPair')
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



