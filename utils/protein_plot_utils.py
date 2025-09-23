from config import *
import re
import os

# sys.path.append(os.path.join(os.path.dirname(__file__), '.'))

#if not platform.system() == "Linux":  # Plotting doesn't work on unix

PYMOL_AVAILABLE = False
try:
    import pymol
    from pymol import cmd
    PYMOL_AVAILABLE = True
except Exception:
    # Keep plotting utils importable on machines without PyMOL
    PYMOL_AVAILABLE = False

from Bio import Align
from utils.phytree_utils import *
from scripts.MSA_Clust import *
from utils.utils import *
from utils.energy_utils import *
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
import matplotlib.pyplot as plt

import math
import pandas as pd


# Map ete leaves to cluster metrics; accept several cluster-key aliases
def _resolve_cluster_key(raw_key: str, cluster_node_values: dict) -> str | None:
    # Examples of incoming raw_key from seqs_ids_to_cluster_ids:
    #   "ShallowMsa_000", "MSA_deep", "DeepMsa", "000", etc.
    cands = []
    k = str(raw_key)

    # as-is
    cands.append(k)

    # if it's just a number, prefer ShallowMsa_### (zero-padded)
    m = re.fullmatch(r"\d+", k)
    if m:
        cands.append(f"ShallowMsa_{int(k):03d}")

    # normalize deep naming
    if "deep" in k.lower():
        cands.extend(["DeepMsa", "MSA_deep"])

    # normalize "Shallow..." variants
    if k.lower().startswith("shallow") and not k.lower().startswith("shallowmsa_"):
        # turn e.g. "Shallow_007" / "ShallowMsa007" into "ShallowMsa_007"
        nn = re.sub(r"(?i)^shallow(msa)?_?", "ShallowMsa_", k)
        # ensure zero-padding
        nn = re.sub(r"(\D)(\d{1,3})$", lambda m: f"{m.group(1)}{int(m.group(2)):03d}", nn)
        cands.append(nn)

    # last: strip any file/prefix like "msa_t__"
    cands.append(k.replace("msa_t__", ""))

    # return the first candidate that exists
    for c in cands:
        if c in cluster_node_values:
            return c
    return None


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
def make_foldswitch_all_plots(
    pdbids, fasta_dir, foldpair_id, pdbchains,
    plot_tree_clusters: bool = False, plot_contacts: bool = True, global_plots: bool = False
):
    """
    Create per-pair plots (contact maps, tree heatmaps, etc.)
    Returns (cmap_dists_vec, seqs_dists_vec, num_seqs_msa_vec, concat_scores)
    """
    print("Plot for foldpair_id: " + foldpair_id)

    # ---------- Paths & ensure a single per-pair figure dir ----------
    fig_dir_root = os.path.join("Pipeline", foldpair_id, "output_figs")
    os.makedirs(fig_dir_root, exist_ok=True)

    # ---------- Inputs ----------
    fasta_file_names = {
        pdbids[i] + pdbchains[i]:
            os.path.join(fasta_dir, foldpair_id, f"{pdbids[i]}{pdbchains[i]}.fasta")
        for i in range(2)
    }
    msa_pred_files = glob(os.path.join(fasta_dir, foldpair_id, "output_cmaps", "msa_transformer", "*.npy"))
    msa_files = glob(os.path.join(fasta_dir, foldpair_id, "output_msa_cluster", "*.a3m"))
    msa_clusters = {os.path.basename(fp)[:-4]: read_msa(fp) for fp in msa_files}

    # Filter 'bad' families: too shallow alignments (no clusters), same PDB ID, etc.
    if len(msa_files) == 0:
        print("Shallow alignment! No MSA Clusters! Skipping family")
        return [None] * 4
    if pdbids[0] == pdbids[1]:
        print("Same PDB-ID for both folds! Might be buggy!")  # not supported fully yet

    # ---------- Load FASTA sequences ----------
    seqs = {}
    for i in range(2):
        key = pdbids[i] + pdbchains[i]
        with open(fasta_file_names[key], "r") as fh:
            seqs[key] = fh.read().splitlines()[1].strip()

    # Quick length guard (MSA-Transformer plotting assumptions)
    if max(len(s) for s in seqs.values()) > 1024:
        print("Long sequence! Length > 1024 (not supported by this plotting path)")
        return [None] * 4

    # ---------- Load predicted CMAPs (shallow + deep) ----------
    # Keys: "ShallowMsa_XXX" or "MSA_deep"
    msa_transformer_pred = {}
    for path in msa_pred_files:
        base = os.path.splitext(os.path.basename(path))[0]  # e.g., "msa_t__ShallowMsa_008"
        key = base[len("msa_t__"):] if base.startswith("msa_t__") else base
        try:
            msa_transformer_pred[key] = np.load(path, allow_pickle=True)
        except Exception:
            msa_transformer_pred[key] = np.genfromtxt(path)

    # ---------- Truth contacts ----------
    try:
        true_cmap = {
            pdbids[i] + pdbchains[i]:
                np.load(os.path.join(fasta_dir, foldpair_id, f"{pdbids[i]}{pdbchains[i]}_pdb_contacts.npy")).astype(int)
            for i in range(2)
        }
        print("Got true cmap!!!")
    except Exception:
        true_cmap = {
            pdbids[i] + pdbchains[i]:
                np.load(os.path.join(fasta_dir, foldpair_id, f"{pdbids[i]}{pdbchains[i]}_pdb_contacts.npy"),
                        allow_pickle=True).astype(int)
            for i in range(2)
        }

    # ---------- Align truth/pred indices ----------
    pairwise_alignment = Align.PairwiseAligner().align(
        seqs[pdbids[0] + pdbchains[0]],
        seqs[pdbids[1] + pdbchains[1]],
    )
    print("Get matching indices: pdbids", pdbids, "pdbchains", pdbchains)
    match_true_cmap, match_predicted_cmaps = get_matching_indices_two_cmaps(
        pairwise_alignment, true_cmap, msa_transformer_pred
    )

    # ---------- Plot CMAPs ----------
    if plot_contacts:
        print("Plot Array Contact Map")
        save_root = os.path.join(fig_dir_root, f"{foldpair_id}_all_clusters_cmap")
        # Ensure directory exists (done above), pass root (plot util appends .png)
        plot_array_contacts_and_predictions(match_predicted_cmaps, match_true_cmap, save_root, foldpair_id=foldpair_id)


    # ---------- Metrics on shared/unique contacts ----------
    shared_unique_contacts, shared_unique_contacts_metrics, contacts_united = \
        match_predicted_and_true_contact_maps(match_predicted_cmaps, match_true_cmap)

    cluster_node_values = {
        ctype: (
            shared_unique_contacts_metrics["shared"][ctype]['long_P@L5'],
            shared_unique_contacts_metrics[pdbids[0] + pdbchains[0]][ctype]['long_P@L5'],
            shared_unique_contacts_metrics[pdbids[1] + pdbchains[1]][ctype]['long_P@L5'],
        )
        for ctype in match_predicted_cmaps
    }

    # ---------- Tree & overlays ----------
    phytree_file = os.path.join('Pipeline', foldpair_id, 'output_phytree', 'DeepMsa_tree.nwk')
    print("Load Biopython treefile: " + phytree_file)
    bio_tree = Phylo.read(phytree_file, "newick")
    print("Convert to ete3 tree:")
    ete_tree = convert_biopython_to_ete3(bio_tree)

    ete_leaves_cluster_ids = seqs_ids_to_cluster_ids(
        os.path.join(fasta_dir, foldpair_id, "output_msa_cluster", "*.a3m"),
        [n.name for n in ete_tree])
    print("Converted seq ids to cluster ids:")

    # Build a dict of leaf-name -> (shared, foldA, foldB) metrics
    entries = {}
    for n in ete_tree.iter_leaves():  # leaves only
        cid = ete_leaves_cluster_ids.get(n.name)
        if not cid or cid == 'p':
            continue
        key = _resolve_cluster_key(cid, cluster_node_values)
        if key is None:
            print(f"[tree] WARN: no metrics for leaf '{n.name}' with cid='{cid}' → skipping")
            continue
        entries[n.name] = cluster_node_values[key]

    if not entries:
        # nothing to plot; keep a minimal frame to avoid downstream crashes
        ete_leaves_node_values = pd.DataFrame(columns=["shared", pdbids[0] + pdbchains[0], pdbids[1] + pdbchains[1]])
    else:
        ete_leaves_node_values = pd.DataFrame(entries).T
        ete_leaves_node_values.columns = ["shared", pdbids[0] + pdbchains[0], pdbids[1] + pdbchains[1]]


    if plot_tree_clusters:
        print("Plot Tree Clusters:")
        cluster_node_values_df = pd.DataFrame(cluster_node_values).T

        representative_cluster_leaves = unique_values_dict({
            n.name: ete_leaves_cluster_ids[n.name] for n in ete_tree if ete_leaves_cluster_ids[n.name] != 'p'
        })

        # Parse compact cluster keys for indexing
        new_indices = [m.group() for s in cluster_node_values_df.index for m in re.finditer(r'M[sS][aA][a-zA-Z0-9_].*', s)]
        print("Parsed shortened indices: ", new_indices)

        tmscores_df = pd.DataFrame(index=new_indices,
                                   columns=['AF_TMscore_fold1', 'AF_TMscore_fold2', 'ESMF_TMscore_fold1', 'ESMF_TMscore_fold2'])


                # --- REPLACE the legacy CSV/glob-based filling of tmscores_df with this unified-table version ---

        # Fill TM scores per cluster from the unified detailed table
        from config import DETAILED_RESULTS_TABLE
        if os.path.isfile(DETAILED_RESULTS_TABLE) and os.path.getsize(DETAILED_RESULTS_TABLE) > 0:
            det = pd.read_csv(DETAILED_RESULTS_TABLE)
        else:
            det = pd.DataFrame(columns=["fold_pair","model","cluster_num","TMscore_fold1","TMscore_fold2"])

        # Accept legacy column names if present
        if "TMscore_fold1" not in det.columns and "score_pdb1" in det.columns:
            det = det.rename(columns={"score_pdb1": "TMscore_fold1", "score_pdb2": "TMscore_fold2"})
        if "fold_pair" not in det.columns and "pair_id" in det.columns:
            det = det.rename(columns={"pair_id": "fold_pair"})

        # Normalize cluster id to something like "ShallowMsa_007" or "DeepMsa"
        det["_cluster_tag"] = det.get("cluster_num", det.get("cluster", "")).astype(str)

        # Keep only rows for this pair
        det = det[det["fold_pair"] == foldpair_id].copy()
        det["model"] = det.get("model", "").astype(str).str.upper()

        # helper to map a displayed index like 'MSA_deep' or 'ShallowMsa_007' back to det._cluster_tag
        def _match_cluster(row_index: str) -> str | None:
            idx = row_index
            # our 'new_indices' already extracted compact names (e.g., 'ShallowMsa_007', maybe lower/upper)
            # make a case-insensitive match
            hits = det[det["_cluster_tag"].astype(str).str.lower() == idx.lower()]
            if len(hits):
                return str(hits["_cluster_tag"].iloc[0])
            # also tolerate bare numbers (e.g., '007')
            m = re.search(r"(\d+)$", idx)
            if m is not None:
                hits = det[det["_cluster_tag"].astype(str).str.contains(m.group(1))]
                if len(hits):
                    return str(hits["_cluster_tag"].iloc[0])
            # finally, try DeepMsa aliasing
            if "deep" in idx.lower():
                hits = det[det["_cluster_tag"].astype(str).str.contains("Deep", case=False, regex=False)]
                if len(hits):
                    return str(hits["_cluster_tag"].iloc[0])
            return None

        # Fill per cluster: AF_* from AF2/AF3 best (choose max over AF2/AF3), ESM_* from ESM2/ESM3 best (max)
        for c in range(len(tmscores_df.index)):
            clus_key = _match_cluster(tmscores_df.index[c])
            if clus_key is None:
                tmscores_df.iloc[c, :] = 0.0
                continue

            # AF best for this cluster (max over AF2/AF3)
            af_rows = det[(det["_cluster_tag"] == clus_key) & (det["model"].isin(["AF2","AF3"]))]
            if len(af_rows):
                tmscores_df.iloc[c, 0] = pd.to_numeric(af_rows["TMscore_fold1"], errors="coerce").max()
                tmscores_df.iloc[c, 1] = pd.to_numeric(af_rows["TMscore_fold2"], errors="coerce").max()
            else:
                tmscores_df.iloc[c, 0] = 0.0
                tmscores_df.iloc[c, 1] = 0.0

            # ESM best for this cluster (max over ESM2/ESM3)
            es_rows = det[(det["_cluster_tag"] == clus_key) & (det["model"].isin(["ESM2","ESM3"]))]
            if len(es_rows):
                tmscores_df.iloc[c, 2] = pd.to_numeric(es_rows["TMscore_fold1"], errors="coerce").max()
                tmscores_df.iloc[c, 3] = pd.to_numeric(es_rows["TMscore_fold2"], errors="coerce").max()
            else:
                tmscores_df.iloc[c, 2] = 0.0
                tmscores_df.iloc[c, 3] = 0.0


        concat_scores = pd.concat([tmscores_df, cluster_node_values_df], ignore_index=True, axis=1)
        concat_scores.columns = ['TM-AF1', 'TM-AF2', 'TM-ESM1', 'TM-ESM2', 'RE-MSAT-COM', 'RE-MSAT1', 'RE-MSAT2']

        # Save figure
        out_root = os.path.join(fig_dir_root, f"{foldpair_id}_phytree_cluster")
        visualize_tree_with_heatmap(ete_tree, concat_scores, out_root)
    else:
        out_root = os.path.join(fig_dir_root, f"{foldpair_id}_phytree")
        visualize_tree_with_heatmap(phytree_file, ete_leaves_node_values, out_root)
        concat_scores = []

    # ---------- Global summaries ----------
    cmap_dists_vec = compute_cmap_distances(match_predicted_cmaps)
    seqs_dists_vec = np.mean(compute_seq_distances(msa_clusters))
    num_seqs_msa_vec = len(seqs)

    # Optional 3D-align plot (only on non-Linux in your original code)
    if not platform.system() == "Linux":
        out_3d = os.path.join(fig_dir_root, f"{foldpair_id}_3d_aligned.png")
        align_and_visualize_proteins(
            os.path.join('Pipeline', foldpair_id, f"{pdbids[0]}.pdb"),
            os.path.join('Pipeline', foldpair_id, f"{pdbids[1]}.pdb"),
            out_3d, False)

    if global_plots:
        print("Make global plots!")
        global_pairs_statistics_plots(output_dir=FIGURE_RES_DIR)

    return cmap_dists_vec, seqs_dists_vec, num_seqs_msa_vec, concat_scores



# Align two pdb structures from two pdb files
def align_and_visualize_proteins(pdb_file1, pdb_file2, output_file, open_environment=True):
    """
    Align two protein structures and save the visualization as an image.

    Args:
    pdb_file1 (str): Path to the first PDB file.
    pdb_file2 (str): Path to the second PDB file.
    output_file (str): Path to save the output image.
    """

    if not PYMOL_AVAILABLE:
        raise RuntimeError("PyMOL is not installed/available. Install it or run this on a machine with PyMOL.")

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
def plot_array_contacts_and_predictions(
    predictions,
    contacts,
    save_file="",
    *,
    foldpair_id: str | None = None
):
    """
    Plot multiple contact-map predictions against truth, and then a separate
    panel for the best-recall clusters. If per-residue ΔG files exist for the
    pair, overlay them as vectors; otherwise skip that overlay.

    Parameters
    ----------
    predictions : dict[str, np.ndarray]
        Map: cluster-name -> predicted contact matrix (NxN).
    contacts : dict[str, np.ndarray]
        Map: fold_id -> true contacts (NxN) for the two folds.
    save_file : str
        Root path (without extension) to save the 'all clusters' figure.
    foldpair_id : str | None
        Pair id like '1dzlA_5keqF'; used to locate DeltaG files under
        Pipeline/<pair>/output_deltaG. If None, ΔG overlay is skipped.
    """
    n_pred = len(predictions)
    # grid shape for the 'all clusters' panel
    n_row = math.ceil(math.sqrt(n_pred))
    n_col = n_row - 1 if n_row * (n_row - 1) >= n_pred else n_row

    # ---- draw all clusters and collect per-cluster recall ----
    n_AA_aligned = len(contacts[next(iter(contacts))])
    fig, axes = plt.subplots(figsize=(18, 18), nrows=n_row, ncols=n_col, layout="compressed")

    ctr = 0
    recall: dict[str, dict[str, float]] = {}  # cluster -> {foldA: rA, foldB: rB}
    for name in predictions.keys():
        ax = axes if n_col == 1 else axes[ctr // n_col, ctr % n_col]
        if n_col == 1:
            ax = axes[ctr]
        ctr += 1

        # plot and receive per-fold recall for this cluster
        recall[name] = plot_foldswitch_contacts_and_predictions(
            predictions[name], contacts, ax=ax, title=name, show_legend=(ctr == 1)
        )

    if save_file:
        plt.savefig(save_file + ".png")
        print("Save cmap fig:", save_file + ".png")
    else:
        plt.show()

    # ---- choose the best-recall cluster for each fold and render a focused panel ----
    if not recall:
        # nothing rendered above; just return early
        return

    # exclude "deep" clusters from best-cluster selection (your original logic)
    filtered = {k: v for k, v in recall.items() if "deep" not in k.lower()}
    # if filtering makes it empty, fall back to all clusters so we don't crash
    pick_from = filtered if filtered else recall

    best_recall_clusters = find_max_keys(pick_from)  # {'foldA': (bestCluster, maxVal), 'foldB': (...)}
    print("Best recall clusters:", best_recall_clusters)

    fold_ids = list(contacts.keys())
    best_cluster_ids = {f: best_recall_clusters[f][0] for f in fold_ids}
    print("Best recall cluster names:", best_cluster_ids)

    # ---- Optional ΔG overlay: read per-pair energy files if available ----
    xvec = yvec = None
    if foldpair_id:
        energy_dir = os.path.join("Pipeline", foldpair_id, "output_deltaG")
        try:
            e0 = read_energy_tuples(os.path.join(energy_dir, f"deltaG_{fold_ids[0][:4]}.txt"))
            e1 = read_energy_tuples(os.path.join(energy_dir, f"deltaG_{fold_ids[1][:4]}.txt"))
            delta_energies, delta_energies_filtered = align_and_compare_residues(
                e0, e1, fold_ids[0][:4], fold_ids[1][:4]
            )
            print("n_AA_aligned=", n_AA_aligned)
            delta_energies_filtered = np.array(delta_energies_filtered[:n_AA_aligned])
            xvec = yvec = delta_energies_filtered
        except Exception as e:
            # Make this non-fatal; just skip ΔG overlay if files not found or alignment fails
            print(f"[plot] NOTE: ΔG overlay skipped for {foldpair_id}: {e}")

    # ---- render 'best clusters' focused panel ----
    plt.figure(figsize=(10, 8))
    best = plot_foldswitch_contacts_and_predictions(
        predictions=(predictions[best_cluster_ids[fold_ids[0]]], predictions[best_cluster_ids[fold_ids[1]]]),
        contacts=contacts,
        title="Best clusters",
        show_legend=True,
        cluster_names=(best_cluster_ids[fold_ids[0]], best_cluster_ids[fold_ids[1]]),
        x_vector=xvec,
        y_vector=yvec,
    )
    print("best recall:", best)

    if save_file:
        plt.savefig(save_file.replace("all", "best"))
    plt.close()



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
            borderaxespad=0.2)  # Adjust padding between the legend and axes

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



# Make global plots for all families
def global_pairs_statistics_plots(output_dir: str | None = None) -> None:
    """
    Build ALL global scatter plots from the unified CSV(s) and save them under output_dir:

      • AF2 / AF3 / ESM2 / ESM3:
           per-pair mean(+) and max(*) of TM1 vs TM2, with a dashed line connecting them.
      • MSA-Transformer:
           per-pair Max Recall (t1_recall vs t2_recall) from each pair's df_cmap.csv.

    Files written (if data present):
      output_dir/fold_pair_scatter_plot_AF2.png
      output_dir/fold_pair_scatter_plot_AF3.png
      output_dir/fold_pair_scatter_plot_ESM2.png
      output_dir/fold_pair_scatter_plot_ESM3.png
      output_dir/fold_pair_scatter_plot_MSA_TRANS.png
    """
    # Default directory
    if output_dir is None:
        output_dir = FIGURE_RES_DIR
    os.makedirs(output_dir, exist_ok=True)

    # -------- AF/ESM from unified detailed table --------
    if (not os.path.isfile(DETAILED_RESULTS_TABLE)) or os.path.getsize(DETAILED_RESULTS_TABLE) == 0:
        print(f"[warn] No detailed CSV at {DETAILED_RESULTS_TABLE}. Run postprocess first.")
        df_det = pd.DataFrame(columns=["fold_pair","model","TMscore_fold1","TMscore_fold2"])
    else:
        df_det = pd.read_csv(DETAILED_RESULTS_TABLE)

    # Legacy column compatibility
    if "TMscore_fold1" not in df_det.columns and "score_pdb1" in df_det.columns:
        df_det = df_det.rename(columns={"score_pdb1": "TMscore_fold1", "score_pdb2": "TMscore_fold2"})
    if "fold_pair" not in df_det.columns:
        if "pair_id" in df_det.columns:
            df_det = df_det.rename(columns={"pair_id": "fold_pair"})
        else:
            df_det["fold_pair"] = "UNK"

    # Normalize model tags
    df_det["model"] = df_det.get("model", "").astype(str).str.upper()

    def _save_tm_scatter(sub_df: pd.DataFrame, model_tag: str, out_path: str) -> None:
        if sub_df.empty:
            print(f"[info] No rows for {model_tag} in {DETAILED_RESULTS_TABLE}")
            return

        g = sub_df.groupby("fold_pair", as_index=False).agg(
            mean_TM1=("TMscore_fold1","mean"),
            mean_TM2=("TMscore_fold2","mean"),
            max_TM1 =("TMscore_fold1","max"),
            max_TM2 =("TMscore_fold2","max"),
        )

        plt.figure(figsize=(8,7))
        first = True
        for _, r in g.iterrows():
            # '+' = mean; '*' = max; dashed line between them
            plt.scatter(r["mean_TM1"], r["mean_TM2"], marker="+", label="Mean" if first else "")
            plt.scatter(r["max_TM1"],  r["max_TM2"],  marker="*", label="Max"  if first else "")
            plt.plot([r["mean_TM1"], r["max_TM1"]], [r["mean_TM2"], r["max_TM2"]],
                     linestyle="--", color="gray", alpha=0.5)
            # annotate the max point with the pair id
            plt.text(r["max_TM1"], r["max_TM2"], r["fold_pair"], fontsize=7, ha="right")
            first = False

        plt.xlabel("TMscore fold1")
        plt.ylabel("TMscore fold2")
        plt.title(f"Per-pair TM (mean vs max): {model_tag}")
        plt.legend(loc="upper left")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(out_path, dpi=160)
        plt.close()
        print(f"[plot] wrote {out_path}")

    # Emit AF/ESM plots
    for tag in ["AF2", "AF3", "ESM2", "ESM3"]:
        _save_tm_scatter(
            df_det[df_det["model"] == tag],
            model_tag=tag,
            out_path=os.path.join(output_dir, f"fold_pair_scatter_plot_{tag}.png"),
        )

    # -------- MSA-Transformer (recall) from per-pair df_cmap.csv --------
    pairs = list_protein_pairs(parsed=False)
    rows = []
    for pid in pairs:
        p_csv = os.path.join(DATA_DIR, pid, "Analysis", "df_cmap.csv")
        if not (os.path.isfile(p_csv) and os.path.getsize(p_csv) > 0):
            continue
        try:
            d = pd.read_csv(p_csv)
        except Exception:
            continue

        c1 = "t1_recall" if "t1_recall" in d.columns else None
        c2 = "t2_recall" if "t2_recall" in d.columns else None
        if not c1 or not c2:
            continue

        rows.append({
            "fold_pair": pid,
            "max_R1": pd.to_numeric(d[c1], errors="coerce").max(),
            "max_R2": pd.to_numeric(d[c2], errors="coerce").max(),
        })

    if rows:
        df_msat = pd.DataFrame(rows)
        out = os.path.join(output_dir, "fold_pair_scatter_plot_MSA_TRANS.png")
        plt.figure(figsize=(8,7))
        first = True
        for _, r in df_msat.iterrows():
            plt.scatter(r["max_R1"], r["max_R2"], marker="*", label="Max" if first else "")
            plt.text(r["max_R1"], r["max_R2"], r["fold_pair"], fontsize=7, ha="right")
            first = False
        plt.xlabel("MSA-Transformer Recall (fold1)")
        plt.ylabel("MSA-Transformer Recall (fold2)")
        plt.title("Per-pair max recall: MSA-Transformer")
        plt.legend(loc="upper left")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(out, dpi=160)
        plt.close()
        print(f"[plot] wrote {out}")
    else:
        print("[info] No df_cmap.csv files found; skipping MSA-Transformer plot.")


# Used in making notebooks
def plot_viz_cmap(file, legend_plot):
    if not file or not os.path.isfile(file):
        print(f"[skip] viz cmap (missing): {file}")
        return
    try:
        visualization_map = np.load(file)
    except Exception as e:
        print(f"[skip] np.load failed: {e}")
        return

    size = visualization_map.shape[0]

#    data = np.random.choice([0, 1, 1.25, 1.5, 1.75], size=(size, size))

    # Create a custom colormap
    colors = ['grey', 'lightblue', 'purple', 'blue', 'magenta']
#    values = [0, 1, 1.25, 1.5, 1.75]
    bounds = [0, 0.99, 1.24, 1.49, 1.74, 2]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(bounds, cmap.N)

    # Create the plot
    plt.figure(figsize=(10, 8))
    im = plt.imshow(visualization_map, cmap=cmap, norm=norm, interpolation='nearest',origin='lower')

    # Create legend elements
    legend_elements = [
        Patch(facecolor='lightblue', edgecolor='black', label='Experiment contact'),
        Patch(facecolor='purple', edgecolor='black', label='Unique State Experiment contact'),
        Patch(facecolor='blue', edgecolor='black', label='Experiment contact predicted'),
        Patch(facecolor='magenta', edgecolor='black', label='Unique State Experiment contact predicted')
    ]

    # Add the legend
    plt.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5))

    # Set title and labels
    plt.title(legend_plot)
    plt.xlabel('Residue Index')
    plt.ylabel('Residue Index')

    # Adjust layout to make room for the legend
    plt.tight_layout()

    # Show the plot
    plt.show()



# Temp: why two versions??
def plot_viz_cmap2(file1, file2, legend_plot):

    print("Input to plot_viz_cmap2: ", file1, file2, legend_plot, flush=True)
    # Load the two 2D numpy arrays
    visualization_map1 = np.load(file1)
    visualization_map2 = np.load(file2)

    # Ensure both arrays have the same shape
    assert visualization_map1.shape == visualization_map2.shape, "Both input arrays must have the same shape"

    # Create a combined array
    combined_map = np.tril(visualization_map1) + np.triu(visualization_map2, k=1)

    # Create a custom colormap
    colors = ['grey', 'blue', 'lightblue', 'purple', 'blue', 'magenta']
    bounds = [0, 0.49, 0.99, 1.24, 1.49, 1.74, 2]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(bounds, cmap.N)

    # Create the plot
    plt.figure(figsize=(10, 8))
    im = plt.imshow(combined_map, cmap=cmap, norm=norm, interpolation='nearest', origin='lower')

    # Create legend elements
    legend_elements = [
        Patch(facecolor='lightblue', edgecolor='black', label='Experiment contact'),
        Patch(facecolor='purple', edgecolor='black', label='Unique State Experiment contact'),
        Patch(facecolor='blue', edgecolor='black', label='Predicted contacts'),
        Patch(facecolor='magenta', edgecolor='black', label='Unique State Experiment contact predicted')
    ]

    # Add the legend
    plt.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5))

    # Set title and labels
    plt.title(legend_plot)
    plt.xlabel('Residue Index')
    plt.ylabel('Residue Index')

    # Adjust layout to make room for the legend
    plt.tight_layout()

    # Show the plot
    plt.show()
