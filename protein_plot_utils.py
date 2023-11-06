from protein_utils import *
import nglview as nv
import py3Dmol
from pyvirtualdisplay import Display
from IPython.display import display, Image
import requests
from PIL import Image
from io import BytesIO
from MSA_Clust import match_predicted_and_true_contact_maps


def make_foldswitch_all_plots(pdbids, fasta_dir, foldpair_id, pdbchains):
    fasta_file_names = {pdbids[i][fold] + pdbchains[i][fold]: fasta_dir + "/" + foldpair_id + "/" + \
                                                              pdbids[i][fold] + pdbchains[i][fold] + '.fasta' for fold
                        in range(2)}  # Added chain to file ID
    msa_file = fasta_dir + "/" + foldpair_id + "/output_get_msa/DeepMsa.a3m"
    MSA = read_msa(msa_file)  # AlignIO.read(open(msa_file), "fasta")
    msa_pred_files = glob(fasta_dir + "/" + foldpair_id + "/output_cmap_esm/*.npy")
    n_cmaps = len(msa_pred_files)
    n_cmaps = min(3, n_cmaps)  # temp for debug !!
    msa_files = glob(fasta_dir + "/" + foldpair_id + "/output_msa_cluster/*.a3m")
    msa_clusters = {file.split("\\")[-1][:-4]: read_msa(file) for file in msa_files}

    # Filter 'bad' families: too shallow alignments (no clusters), same PDB ID, other reasons??
    if len(msa_files) == 0:
        print("Shallow alignment! No MSA Clusters! Skipping family")
        continue
    if pdbids[i][0] == pdbids[i][1]:
        print("Same PDB-ID for both folds! Might be buggy!")  # Not supported yet!  Skipping family")
    #            continue
    #        if foldpair_id in ["1nqjB_1nqdA", "1qlnA_1h38D", "3l5nB_2a73B", "2bzyB_2lqwA", "4cmqB_4zt0C", "5tpnA_1g2cF"]:
    #            print("Special problematic family! " + foldpair_id)
    #            continue
    #       if i < 15:
    #            print("Ignore first!")
    #            continue

    # First load files
    seqs = {}
    for fold in range(2):
        with open(fasta_file_names[pdbids[i][fold] + pdbchains[i][fold]], "r") as text_file:
            seqs[pdbids[i][fold] + pdbchains[i][fold]] = text_file.read().split("\n")[1]
    pairwise_alignment = Align.PairwiseAligner().align(seqs[pdbids[i][0] + pdbchains[i][0]],
                                                       seqs[pdbids[i][1] + pdbchains[i][1]])

    if max([len(seqs[fold]) for fold in seqs]) > 1024:
        print("Long sequence! Length = " + " > maximum supported length of 1024")
        continue

    #        msa_transformer_pred = [None]*n_cmaps
    #        for file in msa_pred_files:
    #            print(file)
    #            msa_transformer_pred[i] = np.genfromtxt(file)
    try:  # read in text format or python format
        msa_transformer_pred = {file.split("\\")[-1][14:-4]: np.genfromtxt(file) for file in msa_pred_files}
    except:
        msa_transformer_pred = {file.split("\\")[-1][14:-4]: np.load(file, allow_pickle=True) for file in
                                msa_pred_files}

    #        print("Predicted cmap sizes for:" + foldpair_id)
    #        print([c.shape[0] for c in msa_transformer_pred.values()])

    #        true_cmap = [None]*2
    #        for fold in range(2):
    #            true_cmap[fold] = np.genfromtxt(fasta_dir +
    #                "/" + foldpair_id + "/" + pdbids[i][fold] + pdbchains[i][fold] + "_pdb_contacts.npy")  # load for pickle, genfromtxt for tab-delimited
    try:
        true_cmap = {pdbids[i][fold] + pdbchains[i][fold]: np.genfromtxt(fasta_dir +  # problem with first !!
                                                                         "/" + foldpair_id + "/" + pdbids[i][fold] +
                                                                         pdbchains[i][
                                                                             fold] + "_pdb_contacts.npy").astype(int)
                     for fold in range(2)}
    #            print("Loaded text true!")
    except:
        #            print("load binary true!")
        true_cmap = {pdbids[i][fold] + pdbchains[i][fold]: np.load(fasta_dir +  # problem with first !!
                                                                   "/" + foldpair_id + "/" + pdbids[i][fold] +
                                                                   pdbchains[i][fold] + "_pdb_contacts.npy",
                                                                   allow_pickle=True).astype(int) for fold in range(2)}

    #        print("seq lens:")
    #        print( [len(seqs[fold]) for fold in seqs])
    #        print("aligned lens:")
    #        print(len((pairwise_alignment[0])))

    #        print("True cmap sizes:")
    #        print([c.shape[0] for c in true_cmap.values()])
    #        print("plot")
    #        print("Predicted cmap sizes:")
    #        print([c.shape[0] for c in msa_transformer_pred.values()])

    #        ttt = copy.deepcopy(true_cmap)
    #        ttt['4rmbB'] = true_cmap['4rmbA']
    #        ttt['4rmbA'] = true_cmap['4rmbB']

    # Match indices of the two folds:
    #        print("Matching: ")
    #        print(len(msa_transformer_pred))
    #        print(true_cmap)
    #        print(len(true_cmap))
    #        match_true_cmap, match_predicted_cmaps = get_matching_indices_two_maps(pairwise_alignment, ttt, msa_transformer_pred)
    match_true_cmap, match_predicted_cmaps = get_matching_indices_two_maps(pairwise_alignment, true_cmap,
                                                                           msa_transformer_pred)

    #        for fold in match_true_cmap.keys():
    #            print(fold)
    #            print(true_cmap[fold].shape)
    #            print(match_true_cmap[fold].shape)
    #        for fold in match_predicted_cmaps.keys():
    #            print(fold)
    #            print(match_predicted_cmaps[fold].shape)
    #        plot_array_contacts_and_predictions(msa_transformer_pred, true_cmap)

    plot_array_contacts_and_predictions(match_predicted_cmaps, match_true_cmap,
                                        fasta_dir + "/Results/Figures/Cmap_MSA/" + foldpair_id)

    shared_unique_contacts, shared_unique_contacts_metrics, contacts_united = match_predicted_and_true_contact_maps(
        match_predicted_cmaps, match_true_cmap)
    print("Cmap metrics shared:")
    print(shared_unique_contacts_metrics["shared"])
    for fold in range(2):
        print("Cmap metrics for " + pdbids[i][fold] + pdbchains[i][fold] + ":")
        print(shared_unique_contacts_metrics[pdbids[i][fold] + pdbchains[i][fold]])
    cluster_node_values = {ctype: shared_unique_contacts_metrics["shared"][ctype]['long_P@L5'] for ctype in
                           match_predicted_cmaps}

    # load tree
    #        phytree_msa_str = "sbatch -o './Pipeline/" + foldpair_id + "/tree_reconstruct_for_" + foldpair_id + ".out' ./Pipeline/tree_reconstruct_params.sh " + foldpair_id  # Take one of the two !!! # ""./input/2qke.fasta 2qke
    #        print(phytree_msa_str)
    phytree_file = './Pipeline/' + foldpair_id + '/output_phytree/DeepMsa_tree.nwk'
    ete_tree = Tree(phytree_file, format=1)
    print("Set node cluster ids values and draw:")
    print(fasta_dir + "/" + foldpair_id + "/output_msa_cluster/*.a3m")
    ete_leaves_cluster_ids = seqs_ids_to_cluster_ids(fasta_dir + "/" + foldpair_id + "/output_msa_cluster/*.a3m",
                                                     [n.name for n in ete_tree])

    print("Set node values and draw:")
    #        nv = list(range(len(ete_tree)))
    ete_leaves_node_values = {n.name: cluster_node_values[ete_leaves_cluster_ids[n.name]] for n in
                              ete_tree}  # update these values to include matching to two folds !!
    draw_tree_with_values(phytree_file, fasta_dir + "/Results/Figures/PhyTree/" + foldpair_id + "_phytree",
                          ete_leaves_node_values)

    # Collect :
    cmap_dists_vec[i] = compute_cmap_distances(match_predicted_cmaps)  # msa_transformer_pred)
    seqs_dists_vec[i] = np.mean(compute_seq_distances(msa_clusters))  # can take the entire sequence !
    num_seqs_msa_vec[i] = len(seqs)

    return cmap_dists_vec, seqs_dists_vec, num_seqs_msa_vec


#        print("Cmap dist: " + str(cmap_dists_vec[i]) + ", seq dist:" + str(seqs_dists_vec[i]))
#        break
# next plotP


# Plot structure using nglview
def plot_pdb_struct(pdb_file, output_3d_image_file = []):
#    view = nv.show_file(pdb_file)
#    view

    with open(pdb_file, "r") as f:
        pdb_data = f.read()

    # Create a 3Dmol view
    viewer = py3Dmol.view(width=400, height=400)

    # Add the protein structure data to the viewer
    viewer.addModel(pdb_data, "pdb")

    # Style the visualization (optional)
    viewer.setStyle({"cartoon": {"color": "spectrum"}})

    # Zoom to fit the structure
    viewer.zoomTo()

    output_file = f"protein_structure.png"
    with open(output_file, "wb") as f:
        f.write(viewer.toImage("png"))
    print(f"Saved {output_file}")

#    display(viewer)
#    output_file = "protein_structure.png"  # Replace with your desired output file name

#    with Display():
#        viewer.savefig(output_file)

    # Show the 3D visualization
#    viewer.show()

    if len(output_3d_image_file) == 0:
        output_3d_image_file = "Pipeline/Results/Figures/3d_struct/" + os.path.basename(pdb_file) + ".png"


#    screenshot = viewer.png()
#    with open(output_3d_image_file, 'wb') as f:
#        f.write(screenshot)

#    viewer.png(output_3d_image_file) # , width=400, height=400)

    return 0


# Plot multiple contacts and predictions together
def plot_array_contacts_and_predictions(predictions, contacts, save_file=[]):
    n_pred = len(predictions)
    n_row = math.ceil(math.sqrt(n_pred))  # *2
    if n_row * (n_row - 1) >= n_pred:  # *2
        n_col = n_row - 1
    else:
        n_col = n_row
    PDB_IDS = predictions.keys()  # [p[name] for p in predictions]
    contacts_ids = contacts.keys()
    fig, axes = plt.subplots(figsize=(18, 18), nrows=n_row, ncols=n_col, layout="compressed")
    #    print("Num cmaps: " + str(n_pred))
    #    print(axes.shape)
    #    fig, axes = plt.subplots(figsize=(18, 6), ncols=n_pred)
    ctr = 0
    #    for ax, name in zip(axes, PDB_IDS):
#    print("n_col=" + str(n_col))
#    print("n_row=" + str(n_row))
#    print(PDB_IDS)
#    print("Contact lens:" + str(len(contacts)))
    for name in PDB_IDS:  # loop over predictions
        if n_col == 1:
            ax = axes[ctr]
        else:
            ax = axes[ctr // n_col, ctr % n_col]
        ctr = ctr + 1
#        print("Plotting prediction: " + name)  # + " -> true: " + true_name)
        plot_foldswitch_contacts_and_predictions(
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


"""Adapted from: https://github.com/rmrao/evo/blob/main/evo/visualize.py"""


def plot_contacts_and_predictions(
        predictions: Union[torch.Tensor, np.ndarray],
        contacts: Union[torch.Tensor, np.ndarray],
        ax: Optional[mpl.axes.Axes] = None,
        # artists: Optional[ContactAndPredictionArtists] = None,
        cmap: str = "Blues",
        ms: float = 1,
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
        long_range_pl = compute_precisions(predictions, contacts, minsep=24)[
            "P@L"
        ].item()
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
    if save_flag:
        plt.savefig('%s.pdf' % title, bbox_inches='tight')


# Plot contacts and predictions for TWO folds !!!
def plot_foldswitch_contacts_and_predictions(
        predictions: Union[torch.Tensor, np.ndarray],
        contacts: Union[torch.Tensor, np.ndarray],
        ax: Optional[mpl.axes.Axes] = None,
        # artists: Optional[ContactAndPredictionArtists] = None,
        cmap: str = "gray_r",  # "Blues",
        ms: float = 5,
        title: Union[bool, str, Callable[[float], str]] = True,
        animated: bool = False,
        show_legend: bool = False,
) -> None:
    fold_ids = list(contacts.keys())
    if isinstance(predictions, torch.Tensor):
        predictions = predictions.detach().cpu().numpy()
    for fold in fold_ids:
        if isinstance(contacts[fold], torch.Tensor):
            contacts[fold] = contacts[fold].detach().cpu().numpy()
    if ax is None:
        ax = plt.gca()

    if len(fold_ids) == 1: # same PDB ID, duplicate
        fold_ids[1] = fold_ids[0]
    seqlen = contacts[fold].shape[0]
#    print(seqlen)
#    for fold in fold_ids:
#        print(contacts[fold].shape)
#    print(fold_ids)
    relative_distance = np.add.outer(-np.arange(seqlen), np.arange(seqlen))
    top_bottom_mask = {fold_ids[0]: relative_distance < 0, fold_ids[1]: relative_distance > 0}
    #    masked_image = np.ma.masked_where(bottom_mask, predictions)
    masked_image = np.ma.masked_where(top_bottom_mask[list(fold_ids)[0]], predictions)
    invalid_mask = np.abs(np.add.outer(np.arange(seqlen), -np.arange(seqlen))) < 6
    predictions = predictions.copy()
    predictions[invalid_mask] = float("-inf")

##    contacts_united = (contacts[fold_ids[0]] + contacts[fold_ids[1]])  # 0: no contact, 1: contact in one, 2: contact in both
##    for fold in fold_ids:
##        contacts_united[np.where(contacts[fold] & (contacts_united == 1) & top_bottom_mask[fold])] = 0
    # Flip contact in one and both:
 ##   cc = copy.deepcopy(contacts_united)
##    contacts_united[cc == 1] = 2
##    contacts_united[cc == 2] = 1

    _, _, contacts_united = match_predicted_and_true_contact_maps({title:predictions}, contacts)

    topl_val = np.sort(predictions.reshape(-1))[-seqlen]
    pred_contacts = predictions >= topl_val
    true_positives, false_positives, other_contacts = {}, {}, {}  # [None]*2, [None]*2, [None]*2

    for fold in fold_ids:
        #        print(fold)
        #        print(true_positives[fold])
        #        print(contacts[fold])
        #        print("Top-Bottom")
        #       print(top_bottom_mask[fold])
        true_positives[fold] = contacts[fold] & pred_contacts & top_bottom_mask[fold]
        false_positives[fold] = ~contacts[fold] & pred_contacts & top_bottom_mask[fold]
        other_contacts[fold] = contacts[fold] & ~pred_contacts & top_bottom_mask[fold]

    if isinstance(title, str):
        title_text: Optional[str] = title
    elif title:
        long_range_pl = compute_precisions(predictions, contacts, minsep=24)[
            "P@L"
        ].item()
        if callable(title):
            title_text = title(long_range_pl)
        else:
            title_text = f"Long Range P@L: {100 * long_range_pl:0.1f}"
    else:
        title_text = None

    #    img = ax.imshow(masked_image, cmap=cmap, animated=animated)  # Show main image
    img = ax.imshow(contacts_united, cmap=cmap, animated=animated)  # Show main image
#    for fold in fold_ids:
#        oc = ax.plot(*np.where(other_contacts[fold]), "o", c="grey", ms=ms, label="other")[0]
    ms = ms * 50 / seqlen
#    print("ms: " + str(ms))
    fp = ax.plot(*np.where(false_positives[fold_ids[0]]), "o", c="r", ms=ms, label="FP")[0]
    tp = ax.plot(*np.where(true_positives[fold_ids[0]]), "o", c="b", ms=ms, label="TP")[0]
    fp = ax.plot(*np.where(false_positives[fold_ids[1]]), "o", c="r", ms=ms)[0]
    tp = ax.plot(*np.where(true_positives[fold_ids[1]]), "o", c="b", ms=ms)[0]
    ti = ax.set_title(title_text) if title_text is not None else None
    # artists = ContactAndPredictionArtists(img, oc, fp, tp, ti)

    # Show second structure here!
    if show_legend:
        ax.legend(loc="upper left")
    ax.axis("square")
    ax.set_xlabel(fold_ids[0], fontsize=14)
    ax.set_ylabel(fold_ids[1], fontsize=14)
    ax.set_xlim([0, seqlen])
    ax.set_ylim([0, seqlen])
    save_flag = False  # add as input
    if save_flag:
        plt.savefig('%s.pdf' % title, bbox_inches='tight')