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
    """Avg squared distance of cluster cmaps from their centroid.

    Aligned cluster cmaps can have slightly different sizes (a few columns), so
    restrict to the MODAL shape - the centroid is only defined over equal-sized
    maps - rather than padding to a common frame."""
    items = {k: np.asarray(v, dtype=float) for k, v in (cmaps or {}).items()}
    items = {k: v for k, v in items.items() if v.ndim == 2}
    if not items:
        return 0.0
    if cmap_main is None or np.size(cmap_main) == 0:
        from collections import Counter
        modal = Counter(v.shape for v in items.values()).most_common(1)[0][0]
        items = {k: v for k, v in items.items() if v.shape == modal}
        stack = np.stack(list(items.values()))
        cmap_main = stack.mean(axis=0)
    else:
        cmap_main = np.asarray(cmap_main, dtype=float)
        items = {k: v for k, v in items.items() if v.shape == cmap_main.shape}
    if not items:
        return 0.0
    return float(np.mean([np.sum((v - cmap_main) ** 2) for v in items.values()]))


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


def _decompose_shared_unique(tA, tB, fold_a, fold_b):
    """Shared / fold-unique contact decomposition for two EQUAL-SIZED truth cmaps."""
    n = tA.shape[0]
    rel = np.add.outer(-np.arange(n), np.arange(n))
    tbm = {fold_a: rel < 0, fold_b: rel > 0}
    tmap = {fold_a: tA, fold_b: tB}
    united = (tA + tB)
    for f in (fold_a, fold_b):
        united[np.where(tmap[f] & (united == 1) & tbm[f])] = 0
    cc = copy.deepcopy(united)
    united[cc == 1] = 2
    united[cc == 2] = 1
    out = {"shared": (united == 1).astype(int)}
    for f in (fold_a, fold_b):
        u = ((united == 2) & tbm[f]).astype(int)
        out[f] = u + u.transpose()
    return out, united


def match_predicted_and_true_contact_maps(cmap_clusters, cmap_true, pred_truth_slices=None):
    """Split contacts into shared/unique and evaluate predictions against each.

    Uses the canonical per-pair alignment computed ONCE upstream
    (align_truth_and_preds_via_msa_first2 / the single pairwise alignment): when
    `pred_truth_slices` is supplied, each cluster prediction is evaluated against
    the truth sliced to the SAME aligned columns as that prediction, so sizes
    always agree - no padding/cropping needed."""
    fold_ids = list(cmap_true.keys())
    fa, fb = fold_ids[0], fold_ids[1]

    # Whole-frame decomposition (returned for callers that want the global view).
    shared_unique_contacts, contacts_united = _decompose_shared_unique(
        cmap_true[fa], cmap_true[fb], fa, fb)

    shared_unique_contacts_metrics = {"shared": {}, fa: {}, fb: {}}
    for clust in cmap_clusters:
        pred = cmap_clusters[clust]
        if pred_truth_slices and clust in pred_truth_slices:
            sl = pred_truth_slices[clust]
            su, _ = _decompose_shared_unique(sl[fa], sl[fb], fa, fb)   # matched columns
        else:
            su = shared_unique_contacts
        for ctype in ("shared", fa, fb):
            shared_unique_contacts_metrics[ctype][clust] = evaluate_prediction(pred, su[ctype])

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
        title="Deep MSA",
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

    # The two folds' best clusters live on DIFFERENT aligned column sets. Slice
    # both (and their matched truths) to the COMMON columns - from the single
    # canonical alignment's per-prediction `cols` - so the F1-upper/F2-lower
    # single-matrix overlay is one consistent frame and never size-mismatches.
    nameA, nameB = best_names[fold_ids[0]], best_names[fold_ids[1]]
    predA, predB = preds[nameA], preds[nameB]   # whole 2D cmaps (not a row!)
    slA = pred_truth_slices.get(nameA, {}) if pred_truth_slices else {}
    slB = pred_truth_slices.get(nameB, {}) if pred_truth_slices else {}
    targets_override_best = None
    if "cols" in slA and "cols" in slB:
        colsA, colsB = np.asarray(slA["cols"]), np.asarray(slB["cols"])
        common = np.intersect1d(colsA, colsB)
        if common.size:
            posA = np.searchsorted(colsA, common)
            posB = np.searchsorted(colsB, common)
            predA = predA[np.ix_(posA, posA)]
            predB = predB[np.ix_(posB, posB)]
            targets_override_best = {
                fold_ids[0]: slA[fold_ids[0]][np.ix_(posA, posA)],
                fold_ids[1]: slB[fold_ids[1]][np.ix_(posB, posB)],
            }
    if targets_override_best is None and (slA or slB):
        targets_override_best = {
            fold_ids[0]: slA.get(fold_ids[0], contacts[fold_ids[0]]),
            fold_ids[1]: slB.get(fold_ids[1], contacts[fold_ids[1]]),
        }
    pred_best = (predA, predB)

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

    # Sparse top-L predicted contacts at |i-j| >= 6 (see select_topl_contacts:
    # the old ``predictions >= sort()[-L]`` rule produced a dense lattice on
    # constant/saturated score maps).
    pred_contacts = select_topl_contacts(predictions, num_contacts=seqlen, min_sep=6)
    contacts = np.asarray(contacts, dtype=bool)
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


def select_topl_contacts(score_map, num_contacts=None, min_sep=6, topk_factor=1.0):
    """Select a SPARSE set of predicted contacts from a continuous score matrix.

    Standard contact-map evaluation convention: keep only the highest-scoring
    off-diagonal residue pairs at sequence separation |i-j| >= ``min_sep``, taking
    the top ``L * topk_factor`` of them (L = sequence length), and return a
    symmetric boolean contact mask.

    This replaces the previous ``predictions >= sort(...)[-L]`` rule, which on a
    constant / all-zero / saturated matrix would select EVERY off-diagonal cell
    (the dense red "False Positive" lattice bug): with ``>=`` and many tied
    values the threshold is met everywhere. Here we instead rank the upper
    triangle and pick a fixed-size top-K set, so the result is always sparse.

    Parameters
    ----------
    score_map : (L, L) array of continuous prediction scores (higher = more likely).
    num_contacts : explicit number of contacts to keep; if None, uses round(L * topk_factor).
    min_sep : minimum sequence separation |i-j| for a pair to be eligible (default 6).
    topk_factor : multiplier on L when ``num_contacts`` is None (1.0 -> top-L, 0.5 -> top-L/2).

    Returns
    -------
    (L, L) boolean array, symmetric, True at selected predicted contacts.
    """
    score_map = np.asarray(score_map, dtype=float)
    L = score_map.shape[0]
    pred = np.zeros((L, L), dtype=bool)
    if L <= min_sep:
        return pred

    # Upper-triangle indices at the required minimum separation.
    iu, ju = np.triu_indices(L, k=min_sep)
    scores = score_map[iu, ju]

    # Only finite, non-degenerate scores are eligible (drops -inf / NaN cells).
    finite = np.isfinite(scores)
    if not np.any(finite):
        return pred
    iu, ju, scores = iu[finite], ju[finite], scores[finite]

    if num_contacts is None:
        num_contacts = int(round(L * topk_factor))
    num_contacts = int(max(0, min(num_contacts, scores.size)))
    if num_contacts == 0:
        return pred

    # Top-`num_contacts` highest scores (partition is O(n); ties resolved arbitrarily
    # but the SET stays size-num_contacts, so the map remains sparse).
    top_idx = np.argpartition(scores, -num_contacts)[-num_contacts:]
    pred[iu[top_idx], ju[top_idx]] = True
    pred |= pred.T  # symmetrize: contact maps are undirected
    return pred


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
    # ``predictions_copy`` keeps the continuous scores (near-diagonal masked to
    # -inf) for the precision metrics; ``pred_contacts`` is the SPARSE discrete
    # top-L contact set that actually gets plotted as colored dots.
    predictions_copy = [p.copy() for p in preds_np]
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
        # Top-L predicted contacts at |i-j| >= display_sep_min (sparse, symmetric).
        pred_contacts[p] = select_topl_contacts(
            predictions_copy[p], num_contacts=Lp, min_sep=display_sep_min)

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
        # Per-residue ThermoMPNN ΔΔG carries extreme outliers (raw sums reach
        # thousands of kcal/mol); scaling the strip to raw min/max blows the
        # colorbar out to ±20000 and washes everything to one colour. Use a
        # SYMMETRIC, robust (98th-percentile-of-|ΔΔG|) limit so the diverging
        # map is centred at 0: red = residue favours F1, blue = favours F2.
        finite = combined_vector[np.isfinite(combined_vector)]
        lim = float(np.nanpercentile(np.abs(finite), 98)) if finite.size else 1.0
        if not np.isfinite(lim) or lim <= 0:
            lim = 1.0
        norm = plt.Normalize(vmin=-lim, vmax=lim)
        strip_cmap = "coolwarm"
        ax_xvec.imshow(x_vector[np.newaxis, :], aspect="auto", cmap=strip_cmap, norm=norm)
        ax_xvec.axis("off")
        ax_yvec.imshow(y_vector[:, np.newaxis], aspect="auto", cmap=strip_cmap, norm=norm)
        ax_yvec.axis("off")
        cbar_ax = fig.add_axes([0.93, 0.3, 0.02, 0.4])
        cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=strip_cmap), cax=cbar_ax)
        cbar.set_label("ΔΔG (kcal/mol)", rotation=90, labelpad=-45, va='bottom')
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
        ax.set_xlabel(f"{fold_ids[0]} ({cx}), recall={round(recall[fold_ids[0]], 4)}", fontsize=12)
        ax.set_ylabel(f"{fold_ids[1]} ({cy}), recall={round(recall[fold_ids[1]], 4)}", fontsize=12)

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
