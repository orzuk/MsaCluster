"""Contact map plotting and evaluation helpers for fold-switching analysis."""

import copy
import math
import os
import random
import re

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

import torch

from utils.protein_utils import compute_precisions, evaluate_prediction
from utils.msa_utils import load_fasta
from utils.energy_utils import read_energy_tuples, align_and_compare_residues
from utils.utils import find_max_keys

from glob import glob


# ---------------------------------------------------------------------------
# Restored from scripts/MSA_Clust.py (deleted in dead-code removal)
# ---------------------------------------------------------------------------
def compute_cmap_distances(cmaps, cmap_main=None):
    """Avg squared distance of cluster cmaps from centroid."""
    if cmap_main is None:
        cmap_main = []
    D = 0
    n_maps = len(cmaps)
    keys_list = list(cmaps.keys())
    if len(cmap_main) == 0:
        cmap_main = cmaps[keys_list[0]]
        for i in range(1, n_maps):
            cmap_main += cmaps[keys_list[i]]
        cmap_main = cmap_main / n_maps
    for cur_cmap in cmaps.keys():
        D += sum(cmaps[cur_cmap] - cmap_main) ** 2
    return D / n_maps


def compute_seq_distances(MSA_clust):
    """Pairwise Levenshtein distances between MSA clusters."""
    from polyleven import levenshtein

    n_clusters = len(MSA_clust)
    D = np.zeros((n_clusters, n_clusters))
    keys_list = list(MSA_clust.keys())
    max_seq_per_cluster = 10
    for i in range(n_clusters):
        for j in range(i, n_clusters):
            n_i = len(MSA_clust[keys_list[i]])
            n_j = len(MSA_clust[keys_list[j]])
            II = random.sample(range(n_i), min(n_i, max_seq_per_cluster))
            JJ = random.sample(range(n_j), min(n_j, max_seq_per_cluster))
            for seq_i in II:
                for seq_j in JJ:
                    D[i, j] += levenshtein(
                        str(MSA_clust[keys_list[i]][seq_i][1]),
                        str(MSA_clust[keys_list[j]][seq_j][1]),
                    )
            D[i, j] = D[i, j] / (len(II) * len(JJ))
            D[j, i] = D[i, j]
    return D


def match_predicted_and_true_contact_maps(cmap_clusters, cmap_true):
    """Split contacts into shared/unique and evaluate predictions against each."""
    fold_ids = list(cmap_true.keys())
    seqlen = cmap_true[fold_ids[0]].shape[0]

    relative_distance = np.add.outer(-np.arange(seqlen), np.arange(seqlen))
    top_bottom_mask = {fold_ids[0]: relative_distance < 0, fold_ids[1]: relative_distance > 0}

    contacts_united = (cmap_true[fold_ids[0]] + cmap_true[fold_ids[1]])
    for fold in fold_ids:
        contacts_united[np.where(cmap_true[fold] & (contacts_united == 1) & top_bottom_mask[fold])] = 0
    cc = copy.deepcopy(contacts_united)
    contacts_united[cc == 1] = 2
    contacts_united[cc == 2] = 1

    shared_unique_contacts = {"shared": (contacts_united == 1).astype(int), fold_ids[0]: None, fold_ids[1]: None}
    for fold in fold_ids:
        shared_unique_contacts[fold] = ((contacts_united == 2) & top_bottom_mask[fold]).astype(int)
        shared_unique_contacts[fold] = shared_unique_contacts[fold] + shared_unique_contacts[fold].transpose()

    shared_unique_contacts_metrics = {"shared": None, fold_ids[0]: None, fold_ids[1]: None}
    for ctype in shared_unique_contacts_metrics:
        shared_unique_contacts_metrics[ctype] = {}
        for clust in cmap_clusters:
            shared_unique_contacts_metrics[ctype][clust] = evaluate_prediction(cmap_clusters[clust], shared_unique_contacts[ctype])

    return shared_unique_contacts, shared_unique_contacts_metrics, contacts_united


def seqs_ids_to_cluster_ids(msa_cluster_dir, seqs_ids=None):
    """Map individual MSA sequence IDs to their cluster names."""
    if seqs_ids is None:
        seqs_ids = []
    if msa_cluster_dir[-1] == "/":
        msa_files = os.listdir(msa_cluster_dir)
    else:
        msa_files = [os.path.basename(f) for f in glob(msa_cluster_dir)]

    seqs_IDs = {
        msa_file_name.replace('.a3m', '')[7:]:
            load_fasta(os.path.dirname(msa_cluster_dir) + "/" + msa_file_name)[0]
        for msa_file_name in msa_files
    }
    cluster_ids = {}
    for cluster in seqs_IDs:
        cluster_ids.update({s: cluster for s in seqs_IDs[cluster]})

    if len(seqs_ids) == 0:
        return cluster_ids
    else:
        return {s: cluster_ids[s] if s in cluster_ids else 'p' for s in seqs_ids}


# ---------------------------------------------------------------------------
# Helper shared with plot_3d
# ---------------------------------------------------------------------------
def _normalize_cluster_tag(s: str) -> str:
    """Canonical cluster tag: 'ShallowMsa_019', 'DeepMsa', etc."""
    s = s.strip()
    if s.lower() in ("deep", "deepmsa", "msa_deep", "deep_msa"):
        return "DeepMsa"
    m = re.fullmatch(r"(\d+)", s)
    if m:
        return f"ShallowMsa_{int(m.group(1)):03d}"
    m2 = re.fullmatch(r"ShallowMsa_(\d+)", s)
    if m2:
        return f"ShallowMsa_{int(m2.group(1)):03d}"
    return s


def _cluster_short_disp(tag: str) -> str:
    """ShallowMsa_019 -> C#19 ; DeepMsa/MSA_deep -> D; else pass-through."""
    t = _normalize_cluster_tag(tag or "")
    if "deep" in t.lower():
        return "D"
    m = re.search(r"(\d+)$", t)
    return f"C#{int(m.group(1))}" if m else t


# ---------------------------------------------------------------------------
# Contact map plotting
# ---------------------------------------------------------------------------
def _plot_contacts_panel(match_predicted_cmaps, match_true_cmap, fig_dir_root, foldpair_id, pred_truth_slices=None):
    """Write the all-clusters + best-clusters contact-map panels."""
    all_save = os.path.join(fig_dir_root, f"{foldpair_id}_all_clusters_cmap")
    plot_array_contacts_and_predictions(
        predictions=match_predicted_cmaps,
        contacts=match_true_cmap,
        save_file=all_save,
        foldpair_id=foldpair_id,
        pred_truth_slices=pred_truth_slices,
    )

    # Deep MSA panel (single subplot)
    deep_key = None
    for k in match_predicted_cmaps:
        if "deep" in k.lower():
            deep_key = k
            break
    if deep_key is None:
        return

    plt.figure(figsize=(10, 8))
    plot_foldswitch_contacts_and_predictions(
        predictions=match_predicted_cmaps[deep_key],
        contacts=match_true_cmap,
        title=f"Deep-MSA ({deep_key})",
        show_legend=True,
        cluster_names=("D", "D"),
    )
    out = os.path.join(fig_dir_root, f"{foldpair_id}_deep_cmap.png")
    plt.savefig(out, dpi=200)
    plt.close()
    print("[plot] wrote cmaps", out)


def plot_array_contacts_and_predictions(
    predictions,
    contacts,
    save_file="",
    *,
    foldpair_id: str | None = None,
    pred_truth_slices: dict | None = None
):
    """
    Draw grid of ALL clusters and one-panel BEST clusters.
    """
    from typing import Union, Optional, Callable, Tuple, Dict, List, Any

    preds = dict(predictions or {})
    n_pred = len(preds)
    if n_pred == 0:
        print("[plot] no predictions to plot.")
        return

    fold_ids = list(contacts.keys())

    n_row = int(math.ceil(math.sqrt(n_pred)))
    n_col = n_row - 1 if n_row * (n_row - 1) >= n_pred else n_row

    fig, axes = plt.subplots(n_row, n_col, figsize=(12, 12), layout="compressed")
    axes = np.array(axes).reshape(-1)

    recall: dict[str, dict[str, float]] = {}
    names = list(preds.keys())
    for i, name in enumerate(names):
        ax = axes[i]
        pred_ab = preds[name]
        targets_override = pred_truth_slices.get(name) if pred_truth_slices else None

        recall[name] = plot_foldswitch_contacts_and_predictions(
            predictions=pred_ab,
            contacts=contacts,
            ax=ax,
            title=_cluster_short_disp(name),
            show_legend=False,
            targets_override=targets_override
        )

    for j in range(n_pred, len(axes)):
        axes[j].axis("off")

    legend_handles = [
        Patch(facecolor='lightgray', edgecolor='#444', label='Shared Contacts'),
        Patch(facecolor='darkgray',  edgecolor='#444', label='Unique Contacts'),
        Line2D([0], [0], marker='o', linestyle='None', color='red',  label='False Positives',        markersize=6),
        Line2D([0], [0], marker='o', linestyle='None', color='blue', label='True Shared Positives',  markersize=6),
        Line2D([0], [0], marker='o', linestyle='None', color='green',label='True Unique Positives',  markersize=6),
    ]
    fig.legend(handles=legend_handles, loc='lower center', bbox_to_anchor=(0.5, -0.02),
               ncol=5, frameon=False, fontsize=8)
    fig.subplots_adjust(bottom=0.18)

    if save_file:
        plt.savefig(save_file + ".png", dpi=200)
        print("Save cmap fig:", save_file + ".png")
    plt.close(fig)

    if not recall:
        return

    filtered = {k: v for k, v in recall.items() if "deep" not in k.lower()}
    pick_from = filtered if filtered else recall

    best = find_max_keys(pick_from)
    print("Best recall clusters:", best)
    best_names = {fold_ids[0]: best[fold_ids[0]][0], fold_ids[1]: best[fold_ids[1]][0]}
    print("Best recall cluster names:", best_names)

    pred_best = (preds[best_names[fold_ids[0]]][0],
                 preds[best_names[fold_ids[1]]][1])

    targets_override_best = None
    if pred_truth_slices:
        tA = pred_truth_slices.get(best_names[fold_ids[0]], {})
        tB = pred_truth_slices.get(best_names[fold_ids[1]], {})
        targets_override_best = {
            fold_ids[0]: tA.get(fold_ids[0], contacts[fold_ids[0]]),
            fold_ids[1]: tB.get(fold_ids[1], contacts[fold_ids[1]]),
        }

    xvec = yvec = None
    if foldpair_id:
        energy_dir = os.path.join("Pipeline", "FoldPairs", foldpair_id, "output_deltaG")
        try:
            e0 = read_energy_tuples(os.path.join(energy_dir, f"deltaG_{fold_ids[0][:4]}.txt"))
            e1 = read_energy_tuples(os.path.join(energy_dir, f"deltaG_{fold_ids[1][:4]}.txt"))
            tmp_delta = align_and_compare_residues(e0, e1, fold_ids[0][:4], fold_ids[1][:4])
            delta_energies_filtered = tmp_delta[1]
            xvec = yvec = np.array(delta_energies_filtered[:len(contacts[next(iter(contacts))])])
        except Exception as e:
            print(f"[plot] NOTE: DG overlay skipped for {foldpair_id}: {e}")

    plt.figure(figsize=(10, 8))
    plot_foldswitch_contacts_and_predictions(
        predictions=pred_best,
        contacts=contacts,
        title="Best clusters",
        show_legend=True,
        cluster_names=(
            _cluster_short_disp(best_names[fold_ids[0]]),
            _cluster_short_disp(best_names[fold_ids[1]]),
        ),
        x_vector=xvec, y_vector=yvec,
        targets_override=targets_override_best
    )
    if save_file:
        plt.savefig(save_file.replace("all", "best") + ".png", dpi=200)
    plt.close()


"""Adapted from: https://github.com/rmrao/evo/blob/main/evo/visualize.py"""
def plot_contacts_and_predictions(
        predictions,
        contacts,
        ax=None,
        cmap: str = "Blues",
        ms: float = 0.5,
        title=True,
        animated: bool = False,
) -> None:
    from typing import Union, Optional, Callable

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
        title_text = title
    elif title:
        long_range_pl = compute_precisions(predictions, contacts, minsep=24)["P@L"].item()
        if callable(title):
            title_text = title(long_range_pl)
        else:
            title_text = f"Long Range P@L: {100 * long_range_pl:0.1f}"
    else:
        title_text = None

    img = ax.imshow(masked_image, cmap=cmap, animated=animated)
    oc = ax.plot(*np.where(other_contacts), "o", c="grey", ms=ms, label="other")[0]
    fp = ax.plot(*np.where(false_positives), "o", c="r", ms=ms, label="FP")[0]
    tp = ax.plot(*np.where(true_positives), "o", c="b", ms=ms, label="TP")[0]
    ti = ax.set_title(title_text) if title_text is not None else None

    ax.legend(loc="upper left")
    ax.axis("square")
    ax.set_xlim([0, seqlen])
    ax.set_ylim([0, seqlen])
    save_flag = False
    if save_flag:
        plt.savefig('%s.pdf' % title, bbox_inches='tight')


def plot_foldswitch_contacts_and_predictions(
        predictions,
        contacts,
        x_vector=None,
        y_vector=None,
        cluster_names=None,
        vector_cmap="viridis",
        ax=None,
        cmap: str = "gray_r",
        ms: float = 3,
        hit_sign: str = 'x',
        title=True,
        animated: bool = False,
        show_legend: bool = False,
        save_flag: bool = False,
        targets_override=None):

    def _ensure_square_2d(arr, fallback_L: int, tag: str):
        if isinstance(arr, torch.Tensor):
            arr = arr.detach().cpu().numpy()
        arr = np.asarray(arr)
        if arr.ndim == 2 and arr.shape[0] == arr.shape[1]:
            return arr.astype(float, copy=False)
        if arr.ndim == 1:
            n = arr.size
            L = int(np.sqrt(n))
            if L * L == n and L > 0:
                return arr.reshape(L, L).astype(float, copy=False)
        if fallback_L > 0:
            print(f"[plot] WARN: prediction '{tag}' is not 2D square (shape={arr.shape}); "
                  f"using zeros of size {fallback_L}x{fallback_L}")
            return np.zeros((fallback_L, fallback_L), dtype=float)
        return None

    if targets_override is not None:
        contacts = {k: targets_override[k] for k in list(contacts.keys()) if k in targets_override}

    fold_ids = list(contacts.keys())

    for k in fold_ids:
        if isinstance(contacts[k], torch.Tensor):
            contacts[k] = contacts[k].detach().cpu().numpy()
        contacts[k] = np.asarray(contacts[k])

    if len(fold_ids) == 1:
        fold_ids = [fold_ids[0], fold_ids[0]]

    seqlen = int(contacts[fold_ids[0]].shape[0])

    if not isinstance(predictions, tuple):
        predictions = [predictions, predictions]

    preds_np = []
    for i in range(2):
        fixed = _ensure_square_2d(predictions[i], seqlen, f"fold{i}")
        if fixed is None:
            fixed = np.zeros((seqlen, seqlen), dtype=float)
        preds_np.append(fixed)

    if ax is None:
        ax = plt.gca()

    rel = np.add.outer(-np.arange(seqlen), np.arange(seqlen))
    top_bottom_mask = {
        fold_ids[0]: (rel < 0),
        fold_ids[1]: (rel > 0)
    }

    display_sep_min = 6
    predictions_copy = [p.copy() for p in preds_np]
    topl_val = [0.0, 0.0]
    pred_contacts = [None, None]

    for p in range(2):
        Lp = predictions_copy[p].shape[0]
        if Lp <= 1:
            predictions_copy[p] = np.zeros((seqlen, seqlen), dtype=float)
            pred_contacts[p] = np.zeros((seqlen, seqlen), dtype=bool)
            continue
        idx = np.arange(Lp)
        invalid_mask_p = (np.abs(idx[:, None] - idx[None, :]) < display_sep_min)
        predictions_copy[p][invalid_mask_p] = float("-inf")
        flat = predictions_copy[p].ravel()
        if np.all(~np.isfinite(flat)):
            topl_val[p] = np.inf
        else:
            topl_val[p] = np.sort(flat)[-Lp]
        pred_contacts[p] = (predictions_copy[p] >= topl_val[p])

    _, _, contacts_united = match_predicted_and_true_contact_maps({str(title): predictions_copy[0]}, contacts)

    true_positives, true_positives_unique, false_positives, other_contacts = {}, {}, {}, {}
    recall = {}
    p = 0
    for fold in fold_ids:
        C = np.asarray(contacts[fold], dtype=bool)
        P = np.asarray(pred_contacts[p], dtype=bool)
        top_mask = np.asarray(top_bottom_mask[fold], dtype=bool)
        if C.shape != P.shape:
            Lc = min(C.shape[0], P.shape[0])
            C = C[:Lc, :Lc]; P = P[:Lc, :Lc]; top_mask = top_mask[:Lc, :Lc]

        true_positives[fold]        = C  & P & top_mask
        tu = (np.transpose(contacts_united) == 2)
        if tu.shape != P.shape:
            Lc = min(tu.shape[0], P.shape[0])
            tu = tu[:Lc, :Lc]
        true_positives_unique[fold] = tu & P & top_mask
        false_positives[fold]       = (~C) & P & top_mask
        other_contacts[fold]        = C & (~P) & top_mask
        tp = np.count_nonzero(true_positives[fold])
        miss = np.count_nonzero(other_contacts[fold])
        recall[fold] = (tp / (tp + miss)) if (tp + miss) > 0 else 0.0
        p += 1

    if isinstance(title, str):
        title_text = title
    elif title:
        long_range_pl0 = compute_precisions(predictions_copy[0], contacts[fold_ids[0]], minsep=24)["P@L"].item()
        long_range_pl1 = compute_precisions(predictions_copy[1], contacts[fold_ids[1]], minsep=24)["P@L"].item()
        if callable(title):
            title_text = title(long_range_pl0 + long_range_pl1)
        else:
            title_text = f"Long Range P@L: {100 * long_range_pl0:0.1f} , {100 * long_range_pl1:0.1f}"
    else:
        title_text = None

    include_vectors = (x_vector is not None) or (y_vector is not None)
    if include_vectors:
        fig = plt.figure(figsize=(10, 10))
        gs = GridSpec(2, 2, width_ratios=[1, 0.03], height_ratios=[0.03, 1],
                      figure=fig, wspace=0.05, hspace=0.05)
        ax = fig.add_subplot(gs[1, 0])

    colors_bg = ['white', 'lightgray', 'darkgray']
    custom_cmap = ListedColormap(colors_bg)
    img = ax.imshow(contacts_united, cmap=custom_cmap, animated=animated)

    ms_scaled = ms * 50 / seqlen
    ax.scatter([], [], marker='s', c='lightgray', s=ms_scaled * 50, label='Shared Contacts')
    ax.scatter([], [], marker='s', c='darkgray',  s=ms_scaled * 50, label='Unique Contacts')

    categories = ["false_positives", "true_positives", "true_positives_unique"]
    colors_pt  = ["r", "b", "g"]
    relative_size = [1, 1, 1.5]
    labels = ["False Positives", "True Shared Positives", "True Unique Positives"]
    offset = -0.15

    cat_maps = {
        "false_positives": false_positives,
        "true_positives": true_positives,
        "true_positives_unique": true_positives_unique,
    }

    for i, category in enumerate(categories):
        x_coords, y_coords = np.where(cat_maps[category][fold_ids[0]])
        ax.plot(x_coords + offset, y_coords, hit_sign, c=colors_pt[i], ms=ms_scaled*relative_size[i], label=labels[i])
        x_coords, y_coords = np.where(cat_maps[category][fold_ids[1]])
        ax.plot(x_coords + offset, y_coords, hit_sign, c=colors_pt[i], ms=ms_scaled*relative_size[i])

    if include_vectors:
        ax_xvec = fig.add_subplot(gs[0, 0], sharex=ax)
        ax_yvec = fig.add_subplot(gs[1, 1], sharey=ax)
        combined_vector = np.concatenate([x_vector.flatten(), y_vector.flatten()])
        norm = plt.Normalize(vmin=combined_vector.min(), vmax=combined_vector.max())
        ax_xvec.imshow(x_vector[np.newaxis, :], aspect="auto", cmap=vector_cmap, norm=norm)
        ax_xvec.axis("off")
        ax_yvec.imshow(y_vector[:, np.newaxis], aspect="auto", cmap=vector_cmap, norm=norm)
        ax_yvec.axis("off")
        cbar_ax = fig.add_axes([0.93, 0.3, 0.02, 0.4])
        cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=vector_cmap), cax=cbar_ax)
        cbar.set_label("DDG (kcal/mol)", rotation=90, labelpad=-45, va='bottom')
        plt.tight_layout(rect=[0, 0, 0.9, 1])

    if show_legend:
        legend_elements = [
            Line2D([0], [0], marker='s', color='lightgray', label='Shared Contacts', markersize=10, linestyle='None'),
            Line2D([0], [0], marker='s', color='darkgray', label='Unique Contacts', markersize=10, linestyle='None'),
            Line2D([0], [0], marker='o', color='red',   label='False Positives', markersize=6, linestyle='None'),
            Line2D([0], [0], marker='o', color='blue',  label='True Shared Positives', markersize=6, linestyle='None'),
            Line2D([0], [0], marker='o', color='green', label='True Unique Positives', markersize=6, linestyle='None'),
        ]
        ax.legend(handles=legend_elements, loc="upper center", bbox_to_anchor=(0.5, -0.075),
                  ncol=5, frameon=False, fontsize=7, columnspacing=0.9, handletextpad=0.5, borderaxespad=0.2)

    ax.axis("square")

    if cluster_names is None:
        ax.set_xlabel(f"{fold_ids[0]}, recall={round(recall[fold_ids[0]], 4)}", fontsize=12)
        ax.set_ylabel(f"{fold_ids[1]}, recall={round(recall[fold_ids[1]], 4)}", fontsize=12)
    else:
        cx = _cluster_short_disp(cluster_names[0])
        cy = _cluster_short_disp(cluster_names[1])
        ax.set_xlabel(f"{fold_ids[0]}, recall={round(recall[fold_ids[0]], 4)} — {cx}", fontsize=12)
        ax.set_ylabel(f"{fold_ids[1]}, recall={round(recall[fold_ids[1]], 4)} — {cy}", fontsize=12)

    ax.set_xlim([0, seqlen])
    ax.set_ylim([0, seqlen])

    if isinstance(title, str):
        ax.set_title(title)
    elif title:
        ax.set_title(title_text)

    if save_flag:
        plt.savefig(f"{title}.pdf", bbox_inches='tight')

    return recall


def plot_viz_cmap(file, legend_plot):
    """Used in making notebooks."""
    if not file or not os.path.isfile(file):
        print(f"[skip] viz cmap (missing): {file}")
        return
    try:
        visualization_map = np.load(file)
    except Exception as e:
        print(f"[skip] np.load failed: {e}")
        return

    size = visualization_map.shape[0]

    colors = ['grey', 'lightblue', 'purple', 'blue', 'magenta']
    bounds = [0, 0.99, 1.24, 1.49, 1.74, 2]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(bounds, cmap.N)

    plt.figure(figsize=(10, 8))
    im = plt.imshow(visualization_map, cmap=cmap, norm=norm, interpolation='nearest', origin='lower')

    legend_elements = [
        Patch(facecolor='lightblue', edgecolor='black', label='Experiment contact'),
        Patch(facecolor='purple', edgecolor='black', label='Unique State Experiment contact'),
        Patch(facecolor='blue', edgecolor='black', label='Experiment contact predicted'),
        Patch(facecolor='magenta', edgecolor='black', label='Unique State Experiment contact predicted')
    ]

    plt.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title(legend_plot)
    plt.xlabel('Residue Index')
    plt.ylabel('Residue Index')
    plt.tight_layout()
    plt.show()
