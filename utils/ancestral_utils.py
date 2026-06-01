"""Ancestral fold-preference reconstruction at the cluster level.

For each fold-switching pair:
  1. Build a coarse phylogenetic tree (one leaf per MSA cluster)
  2. Assign fold preference (F1/F2/Amb) per cluster from AF2 TM-scores
  3. Run ancestral state reconstruction (pastml MPPA or Fitch parsimony)
  4. Classify the evolutionary pattern and count transitions
  5. Permutation test for significance
"""

import csv
import os
import re
import copy
import random
import tempfile
from collections import Counter
from typing import Optional

import numpy as np
import pandas as pd
from Bio import Phylo

from ete3 import Tree, TreeStyle, NodeStyle, TextFace

from config import DATA_DIR
from utils.phytree_utils import convert_biopython_to_ete3, prune_tree_to_leaves
from utils.plot_contacts import seqs_ids_to_cluster_ids


# ---------------------------------------------------------------------------
# Cluster tag normalization
# ---------------------------------------------------------------------------

_MSA_TAG_RE = re.compile(r"(?:Shallow)?Msa_(\d+)", re.IGNORECASE)


def _normalize_tag(s: str) -> str:
    """Normalize any cluster tag to canonical form.

    Handles: ShallowMsa_019, Msa_019, 19, deepmsa, deep -> ShallowMsa_019 or DeepMsa
    """
    s = str(s).strip()
    if s.lower() in ("deep", "deepmsa", "msa_deep", "deep_msa"):
        return "DeepMsa"
    m = _MSA_TAG_RE.fullmatch(s)
    if m:
        return f"ShallowMsa_{int(m.group(1)):03d}"
    if s.isdigit():
        return f"ShallowMsa_{int(s):03d}"
    return s


# ---------------------------------------------------------------------------
# 1. Build coarse cluster-level tree
# ---------------------------------------------------------------------------

def _group_leaves_by_cluster(ete_tree, leaf_cluster_map):
    """Return {cluster_tag: [leaf_name, ...]} for leaves present in tree."""
    groups = {}
    for leaf in ete_tree.iter_leaves():
        raw = leaf_cluster_map.get(leaf.name)
        if not raw or raw == "p":
            continue
        tag = _normalize_tag(raw)
        groups.setdefault(tag, []).append(leaf.name)
    return groups


def _cluster_short_label(tag):
    """ShallowMsa_019 -> S19, DeepMsa -> D."""
    if not tag:
        return ""
    if tag.lower().startswith("deep"):
        return "D"
    m = re.search(r"(\d+)$", tag)
    return f"S{int(m.group(1)):02d}" if m else tag


def build_coarse_tree(pair_id: str) -> tuple:
    """Load raw tree, prune to one leaf per cluster, rename leaves.

    Returns
    -------
    tree : ete3.Tree
        Pruned tree with leaves named "D", "S00", "S01", ...
    tag_to_label : dict
        {"ShallowMsa_000": "S00", "DeepMsa": "D", ...}
    label_to_tag : dict
        Inverse of tag_to_label.
    """
    pair_dir = os.path.join(DATA_DIR, pair_id)
    tree_file = os.path.join(pair_dir, "output_phytree", "DeepMsa_tree.nwk")
    cluster_dir = os.path.join(pair_dir, "output_msa_cluster")

    bio_tree = Phylo.read(tree_file, "newick")
    ete_tree = convert_biopython_to_ete3(bio_tree)

    # Map sequence IDs to cluster tags
    leaf_names = [n.name for n in ete_tree.iter_leaves()]
    leaf_cluster_map = seqs_ids_to_cluster_ids(
        os.path.join(cluster_dir, "*.a3m"), leaf_names
    )

    # Group and pick one representative per cluster
    cluster_to_leaves = _group_leaves_by_cluster(ete_tree, leaf_cluster_map)
    reps = {tag: leaves[0] for tag, leaves in cluster_to_leaves.items() if leaves}

    # Prune
    prune_tree_to_leaves(ete_tree, list(reps.values()))

    # Rename leaves: seq_id -> short label
    tag_to_label = {tag: _cluster_short_label(tag) for tag in reps}
    rep_to_tag = {v: k for k, v in reps.items()}

    for leaf in ete_tree.iter_leaves():
        tag = rep_to_tag.get(leaf.name)
        if tag:
            leaf.name = tag_to_label[tag]

    label_to_tag = {v: k for k, v in tag_to_label.items()}
    return ete_tree, tag_to_label, label_to_tag


# ---------------------------------------------------------------------------
# 2. Load AF2 TM-scores per cluster
# ---------------------------------------------------------------------------

def _parse_cluster_tag(s):
    """Parse a cluster_num value into a normalized tag, or None if unparseable."""
    tag = _normalize_tag(str(s))
    if tag.startswith("ShallowMsa_") or tag == "DeepMsa":
        return tag
    return None


def load_cluster_af2_tm(pair_id: str) -> dict:
    """Load df_af.csv and compute per-cluster max AF2 TM1/TM2.

    Returns
    -------
    dict : {cluster_tag: (TM1, TM2)}
        e.g. {"ShallowMsa_000": (0.92, 0.85), "DeepMsa": (0.78, 0.91)}
    """
    csv_path = os.path.join(DATA_DIR, pair_id, "Analysis", "df_af.csv")
    if not os.path.isfile(csv_path):
        raise FileNotFoundError(f"No df_af.csv for {pair_id}: {csv_path}")

    df = pd.read_csv(csv_path)

    # Normalize column names (old format: score_pdb1/score_pdb2)
    if "TMscore_fold1" not in df.columns and "score_pdb1" in df.columns:
        df = df.rename(columns={"score_pdb1": "TMscore_fold1",
                                "score_pdb2": "TMscore_fold2"})

    # Filter AF2 only (new format has 'model' column)
    if "model" in df.columns:
        df = df[df["model"].astype(str).str.lower() == "af2"]

    # Parse cluster tags
    col = "cluster_num" if "cluster_num" in df.columns else "cluster"
    df["_tag"] = df[col].apply(_parse_cluster_tag)
    df = df.dropna(subset=["_tag"])

    # Numeric TM-scores
    df["TM1"] = pd.to_numeric(df["TMscore_fold1"], errors="coerce")
    df["TM2"] = pd.to_numeric(df["TMscore_fold2"], errors="coerce")

    # Per-cluster max
    result = {}
    for tag, grp in df.groupby("_tag"):
        result[tag] = (grp["TM1"].max(), grp["TM2"].max())
    return result


# ---------------------------------------------------------------------------
# 3. Assign fold preference
# ---------------------------------------------------------------------------

def assign_fold_preference(cluster_tms: dict, delta: float = 0.05,
                           centered: bool = True) -> dict:
    """Assign F1/F2/Amb to each cluster based on TM-score difference.

    If ``centered`` is True (default since 2026-05-15), the per-cluster
    TMdiff is centered by subtracting the per-pair median of TMdiff
    across ShallowMsa clusters before applying the threshold. This is
    the same per-pair-median-centering used by fold_diversity_survey.py
    and is necessary because AF2 (and the other methods) often have
    systematic bias toward one fold per family --- under absolute
    classification every cluster gets the same label and the ancestral
    reconstruction sees no within-pair variation, regardless of whether
    real clade-level structure exists.

    Parameters
    ----------
    cluster_tms : {cluster_tag: (TM1, TM2)}
        The DeepMsa entry, if present (tag == 'DeepMsa'), is excluded
        from the median calculation and from the per-leaf classification.
    delta : threshold on the (centered) TMdiff
    centered : whether to subtract per-pair median TMdiff before
        thresholding (default True)

    Returns
    -------
    dict : {cluster_tag: "F1"|"F2"|"Amb"}
    """
    diffs = {}
    for tag, (tm1, tm2) in cluster_tms.items():
        if np.isnan(tm1) or np.isnan(tm2):
            diffs[tag] = float("nan")
        else:
            diffs[tag] = float(tm1) - float(tm2)
    if centered:
        shallow_diffs = [d for tag, d in diffs.items()
                         if tag != "DeepMsa" and not np.isnan(d)]
        median_diff = float(np.median(shallow_diffs)) if len(shallow_diffs) >= 2 else 0.0
    else:
        median_diff = 0.0

    prefs = {}
    for tag, d in diffs.items():
        if np.isnan(d):
            prefs[tag] = "Amb"
        else:
            centered_d = d - (median_diff if tag != "DeepMsa" else 0.0)
            if centered_d > delta:
                prefs[tag] = "F1"
            elif centered_d < -delta:
                prefs[tag] = "F2"
            else:
                prefs[tag] = "Amb"
    return prefs


# ---------------------------------------------------------------------------
# 4. Ancestral state reconstruction
# ---------------------------------------------------------------------------

def _fitch_bottom_up(tree):
    """Fitch parsimony: bottom-up pass. Sets node._fitch_set for each node."""
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            node._fitch_set = {node._state}
        else:
            child_sets = [c._fitch_set for c in node.children]
            inter = set.intersection(*child_sets)
            if inter:
                node._fitch_set = inter
            else:
                node._fitch_set = set.union(*child_sets)


def _fitch_top_down(tree):
    """Fitch parsimony: top-down pass. Resolves ambiguities, sets node._state."""
    for node in tree.traverse("preorder"):
        if node.is_leaf():
            continue  # already set
        if node.is_root():
            # Pick one state from the set (prefer non-Amb)
            s = node._fitch_set
            for pref in ["F1", "F2", "Amb"]:
                if pref in s:
                    node._state = pref
                    break
        else:
            parent_state = node.up._state
            if parent_state in node._fitch_set:
                node._state = parent_state
            else:
                # Pick one from fitch set
                for pref in ["F1", "F2", "Amb"]:
                    if pref in node._fitch_set:
                        node._state = pref
                        break


def run_asr(tree, leaf_states: dict, work_dir: str = None) -> dict:
    """Run ancestral state reconstruction on the coarse tree.

    Tries pastml (MPPA method) first; falls back to Fitch parsimony.

    Parameters
    ----------
    tree : ete3.Tree
        Coarse tree with leaves named by short labels (D, S00, S01, ...)
    leaf_states : {leaf_name: "F1"|"F2"|"Amb"}
    work_dir : directory for temporary files (pastml output)

    Returns
    -------
    node_states : {node_name_or_id: state} for all nodes
    method : "pastml" or "fitch"
    posteriors : {node_name_or_id: {"F1": p, "F2": p, "Amb": p}} or None
    """
    # Assign leaf states to tree nodes
    tree_copy = tree.copy("newick")
    leaf_set = set()
    for leaf in tree_copy.iter_leaves():
        state = leaf_states.get(leaf.name)
        if state is None:
            state = "Amb"
        leaf._state = state
        leaf_set.add(leaf.name)

    # Try pastml
    try:
        return _run_pastml(tree_copy, leaf_states, work_dir)
    except Exception as e:
        print(f"  [info] pastml unavailable ({e}), using Fitch parsimony")

    # Fitch fallback
    _fitch_bottom_up(tree_copy)
    _fitch_top_down(tree_copy)

    node_states = {}
    for i, node in enumerate(tree_copy.traverse("preorder")):
        key = node.name if node.name else f"N{i}"
        node_states[key] = node._state

    return node_states, "fitch", None


def _run_pastml(tree, leaf_states, work_dir):
    """Run pastml MPPA ancestral reconstruction."""
    from pastml.acr import acr

    if work_dir is None:
        work_dir = tempfile.mkdtemp(prefix="pastml_")

    # Write tree to newick
    tree_path = os.path.join(work_dir, "tree.nwk")
    tree.write(outfile=tree_path, format=1)

    # Write annotation CSV
    annot_path = os.path.join(work_dir, "annotations.csv")
    with open(annot_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "fold_pref"])
        for name, state in leaf_states.items():
            w.writerow([name, state])

    # Run ACR
    result = acr(tree_path, annot_path,
                 prediction_method="MPPA",
                 columns=["fold_pref"])

    # Parse result: pastml returns an annotated tree
    result_tree = Tree(os.path.join(work_dir, "named.tree_tree.nwk"), format=1)

    node_states = {}
    posteriors = {}
    for node in result_tree.traverse("preorder"):
        key = node.name if node.name else ""
        # pastml stores predictions as node features
        state = getattr(node, "fold_pref", None)
        if state is None and hasattr(node, "fold_pref_marginal"):
            mp = node.fold_pref_marginal
            state = max(mp, key=mp.get)
            posteriors[key] = dict(mp)
        if state is not None:
            node_states[key] = state

    return node_states, "pastml", posteriors


# ---------------------------------------------------------------------------
# 5. Count transitions and classify pattern
# ---------------------------------------------------------------------------

def count_transitions(tree, node_states: dict) -> int:
    """Count all parent->child state changes."""
    n = 0
    for node in tree.traverse("preorder"):
        if node.is_root():
            continue
        parent_key = node.up.name if node.up.name else ""
        child_key = node.name if node.name else ""
        ps = node_states.get(parent_key)
        cs = node_states.get(child_key)
        if ps and cs and ps != cs:
            n += 1
    return n


def count_fold_transitions(tree, node_states: dict) -> int:
    """Count only F1<->F2 transitions (ignoring Amb edges)."""
    n = 0
    for node in tree.traverse("preorder"):
        if node.is_root():
            continue
        parent_key = node.up.name if node.up.name else ""
        child_key = node.name if node.name else ""
        ps = node_states.get(parent_key)
        cs = node_states.get(child_key)
        if ps and cs and {ps, cs} == {"F1", "F2"}:
            n += 1
    return n


def classify_pattern(tree, node_states: dict) -> str:
    """Classify evolutionary pattern.

    Returns
    -------
    str : one of
        'ancient_divergence' - root children include both F1 and F2 lineages
        'recent_emergence'   - root is one fold, the other appears only in derived clade
        'uniform'            - all leaves same fold preference
        'complex'            - multiple transitions, no clear clade structure
    """
    leaf_states = set()
    for leaf in tree.iter_leaves():
        key = leaf.name if leaf.name else ""
        s = node_states.get(key)
        if s and s != "Amb":
            leaf_states.add(s)

    if len(leaf_states) <= 1:
        return "uniform"

    # Check root children
    root = tree.get_tree_root()
    root_key = root.name if root.name else ""
    root_state = node_states.get(root_key, "Amb")

    child_states = set()
    for child in root.children:
        ck = child.name if child.name else ""
        cs = node_states.get(ck, "Amb")
        if cs != "Amb":
            child_states.add(cs)
        # Also look at leaf descendants
        for leaf in child.iter_leaves():
            lk = leaf.name if leaf.name else ""
            ls = node_states.get(lk, "Amb")
            if ls != "Amb":
                child_states.add(ls)

    # If root's direct children lead to different folds
    if len(root.children) >= 2:
        subtree_folds = []
        for child in root.children:
            folds_in_subtree = set()
            for leaf in child.iter_leaves():
                lk = leaf.name if leaf.name else ""
                ls = node_states.get(lk, "Amb")
                if ls != "Amb":
                    folds_in_subtree.add(ls)
            subtree_folds.append(folds_in_subtree)

        # Ancient divergence: each major subtree has a different dominant fold
        if (len(subtree_folds) >= 2 and
                any("F1" in s and "F2" not in s for s in subtree_folds) and
                any("F2" in s and "F1" not in s for s in subtree_folds)):
            return "ancient_divergence"

    # Count transitions
    n_trans = count_fold_transitions(tree, node_states)
    if n_trans == 1:
        return "recent_emergence"
    elif n_trans >= 2:
        return "complex"

    return "uniform"


# ---------------------------------------------------------------------------
# 6. Permutation test
# ---------------------------------------------------------------------------

def run_permutation_test(tree, leaf_states: dict, n_perm: int = 1000,
                         seed: int = 42) -> float:
    """Shuffle leaf labels, recount transitions via Fitch, compute p-value.

    Returns p-value: fraction of permutations with <= observed transitions.
    (Low p-value = significant clustering of fold preferences on tree.)
    """
    # Observed transitions via Fitch
    obs_tree = tree.copy("newick")
    for leaf in obs_tree.iter_leaves():
        leaf._state = leaf_states.get(leaf.name, "Amb")
    _fitch_bottom_up(obs_tree)
    _fitch_top_down(obs_tree)
    obs_trans = 0
    for node in obs_tree.traverse("preorder"):
        if not node.is_root() and hasattr(node, "_state") and hasattr(node.up, "_state"):
            if node._state != node.up._state:
                obs_trans += 1

    # Permutations
    rng = random.Random(seed)
    labels = list(leaf_states.keys())
    states = list(leaf_states.values())
    n_leq = 0

    for _ in range(n_perm):
        perm_states = states.copy()
        rng.shuffle(perm_states)
        perm_map = dict(zip(labels, perm_states))

        perm_tree = tree.copy("newick")
        for leaf in perm_tree.iter_leaves():
            leaf._state = perm_map.get(leaf.name, "Amb")
        _fitch_bottom_up(perm_tree)
        _fitch_top_down(perm_tree)

        perm_trans = 0
        for node in perm_tree.traverse("preorder"):
            if not node.is_root() and hasattr(node, "_state") and hasattr(node.up, "_state"):
                if node._state != node.up._state:
                    perm_trans += 1

        if perm_trans <= obs_trans:
            n_leq += 1

    return n_leq / n_perm


# ---------------------------------------------------------------------------
# 7. Render ancestral tree figure
# ---------------------------------------------------------------------------

_STATE_COLORS = {"F1": "#d62728", "F2": "#1f77b4", "Amb": "#999999"}


def render_ancestral_tree(tree, node_states: dict, pair_id: str,
                          output_png: str, posteriors: dict = None) -> None:
    """Render the coarse tree with nodes colored by fold preference.

    Saves a PNG to output_png. Uses matplotlib for broad compatibility.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch
    import matplotlib.patches as mpatches

    # Collect layout info
    leaves = list(tree.iter_leaves())
    n_leaves = len(leaves)
    if n_leaves == 0:
        return

    # Compute x-positions (distance from root)
    root = tree.get_tree_root()
    node_x = {}
    node_y = {}

    def _assign_x(node, x):
        key = node.name if node.name else id(node)
        dist = node.dist if node.dist and node.dist > 0 else 0.1
        node_x[key] = x + dist
        for child in node.children:
            _assign_x(child, x + dist)

    _assign_x(root, 0)
    # Root at x=0
    root_key = root.name if root.name else id(root)
    node_x[root_key] = 0

    # Assign y-positions (leaves evenly spaced)
    leaf_idx = 0
    def _assign_y(node):
        nonlocal leaf_idx
        key = node.name if node.name else id(node)
        if node.is_leaf():
            node_y[key] = leaf_idx
            leaf_idx += 1
        else:
            child_ys = []
            for child in node.children:
                _assign_y(child)
                ck = child.name if child.name else id(child)
                child_ys.append(node_y[ck])
            node_y[key] = np.mean(child_ys)

    _assign_y(root)

    # Normalize x
    max_x = max(node_x.values()) if node_x else 1
    if max_x == 0:
        max_x = 1

    fig, ax = plt.subplots(figsize=(6, max(3, 0.5 * n_leaves + 1)))

    # Draw edges
    for node in tree.traverse("preorder"):
        if node.is_root():
            continue
        key = node.name if node.name else id(node)
        pk = node.up.name if node.up.name else id(node.up)
        x0, y0 = node_x[pk], node_y[pk]
        x1, y1 = node_x[key], node_y[key]
        # L-shaped lines
        ax.plot([x0, x0], [y0, y1], color="#444444", linewidth=1.2, zorder=1)
        ax.plot([x0, x1], [y1, y1], color="#444444", linewidth=1.2, zorder=1)

    # Draw nodes
    for node in tree.traverse("preorder"):
        key = node.name if node.name else id(node)
        state = node_states.get(node.name if node.name else "", "Amb")
        color = _STATE_COLORS.get(state, "#999999")
        x, y = node_x[key], node_y[key]

        if node.is_leaf():
            ax.scatter(x, y, c=color, s=120, zorder=3, edgecolors="black", linewidths=0.5)
            ax.text(x + max_x * 0.03, y, f" {node.name} ({state})",
                    va="center", ha="left", fontsize=9)
        else:
            # Internal node: keep "clean" UNLESS a gain/loss event occurred on
            # the branch leading here (reconstructed state differs from parent).
            # Only the few event nodes get a prominent ringed circle, colored by
            # the NEW state, labeled "old->new" (XCL1-paper style).
            parent = node.up
            pstate = (node_states.get(parent.name, "Amb")
                      if parent is not None and parent.name else "Amb")
            if parent is not None and state != pstate:
                ax.scatter(x, y, c=color, s=150, zorder=4, marker="o",
                           edgecolors="black", linewidths=1.6)
                ax.text(x, y + 0.28, f"{pstate}→{state}", ha="center",
                        va="bottom", fontsize=6, fontweight="bold", zorder=6)
            else:
                ax.scatter(x, y, c="#cfcfcf", s=8, zorder=2, marker="o",
                           edgecolors="none")

    # Legend
    legend_handles = [mpatches.Patch(color=c, label=s) for s, c in _STATE_COLORS.items()]
    ax.legend(handles=legend_handles, loc="upper left", fontsize=8, framealpha=0.8)

    # Title
    fold1, fold2 = pair_id.split("_") if "_" in pair_id else (pair_id, "")
    root_state = node_states.get(root.name if root.name else "", "?")
    ax.set_title(f"{fold1} / {fold2}\nRoot: {root_state}", fontsize=11)
    ax.set_xlabel("Branch length", fontsize=9)
    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    plt.tight_layout()
    os.makedirs(os.path.dirname(output_png), exist_ok=True)
    plt.savefig(output_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved ancestral tree: {output_png}")


# ---------------------------------------------------------------------------
# 8. ESM / MSAT concordance (post-hoc)
# ---------------------------------------------------------------------------

def _load_esm_cluster_tm(pair_id: str) -> Optional[dict]:
    """Load df_esm.csv and compute per-cluster max TM1/TM2, if available."""
    csv_path = os.path.join(DATA_DIR, pair_id, "Analysis", "df_esm.csv")
    if not os.path.isfile(csv_path):
        return None
    try:
        df = pd.read_csv(csv_path)
        if "TMscore_fold1" not in df.columns and "score_pdb1" in df.columns:
            df = df.rename(columns={"score_pdb1": "TMscore_fold1",
                                    "score_pdb2": "TMscore_fold2"})
        col = "cluster_num" if "cluster_num" in df.columns else "cluster"
        df["_tag"] = df[col].apply(_parse_cluster_tag)
        df = df.dropna(subset=["_tag"])
        df["TM1"] = pd.to_numeric(df["TMscore_fold1"], errors="coerce")
        df["TM2"] = pd.to_numeric(df["TMscore_fold2"], errors="coerce")
        result = {}
        for tag, grp in df.groupby("_tag"):
            result[tag] = (grp["TM1"].max(), grp["TM2"].max())
        return result
    except Exception:
        return None


def compute_concordance(af2_prefs: dict, other_tms: Optional[dict],
                        delta: float = 0.05) -> Optional[float]:
    """Fraction of clusters where AF2 and other method agree on fold preference."""
    if other_tms is None:
        return None
    other_prefs = assign_fold_preference(other_tms, delta)
    common = set(af2_prefs.keys()) & set(other_prefs.keys())
    if not common:
        return None
    agree = sum(1 for t in common if af2_prefs[t] == other_prefs[t])
    return agree / len(common)


# ---------------------------------------------------------------------------
# 9. Full per-pair analysis
# ---------------------------------------------------------------------------

def analyze_pair(pair_id: str, delta: float = 0.05,
                 output_dir: str = None, n_perm: int = 1000) -> dict:
    """Run full ancestral reconstruction pipeline for one pair.

    Returns a result dict with summary statistics.
    """
    if output_dir is None:
        output_dir = os.path.join(DATA_DIR, pair_id)

    # 1. Build coarse tree
    tree, tag_to_label, label_to_tag = build_coarse_tree(pair_id)

    # 2. Load AF2 TM-scores
    cluster_tms = load_cluster_af2_tm(pair_id)

    # 3. Assign fold preferences
    af2_prefs = assign_fold_preference(cluster_tms, delta)

    # Map preferences to tree leaf labels
    leaf_states = {}
    for tag, pref in af2_prefs.items():
        label = tag_to_label.get(tag)
        if label:
            leaf_states[label] = pref

    # Also assign "Amb" to any tree leaves without TM data
    for leaf in tree.iter_leaves():
        if leaf.name not in leaf_states:
            leaf_states[leaf.name] = "Amb"

    n_clusters = len(leaf_states)
    pref_counts = Counter(leaf_states.values())

    # Name internal nodes before ASR so the copy inherits them
    _name_internal_nodes(tree)

    # 4. Ancestral state reconstruction
    work_dir = os.path.join(output_dir, "Analysis", "ancestral")
    os.makedirs(work_dir, exist_ok=True)
    node_states, method, posteriors = run_asr(tree, leaf_states, work_dir)

    # 5. Count transitions
    n_trans = count_transitions(tree, node_states)
    n_fold_trans = count_fold_transitions(tree, node_states)

    # 6. Classify pattern
    pattern = classify_pattern(tree, node_states)

    # 7. Root state
    root = tree.get_tree_root()
    root_state = node_states.get(root.name, "?")

    # 8. Permutation test (only if there are multiple fold types)
    p_value = None
    if pref_counts.get("F1", 0) > 0 and pref_counts.get("F2", 0) > 0:
        p_value = run_permutation_test(tree, leaf_states, n_perm=n_perm)

    # 9. Render figure
    fig_path = os.path.join(output_dir, "output_figs",
                            f"{pair_id}_ancestral_tree.png")
    try:
        render_ancestral_tree(tree, node_states, pair_id, fig_path, posteriors)
    except Exception as e:
        print(f"  [warn] Could not render tree for {pair_id}: {e}")

    # 10. Save ancestral states CSV
    states_csv = os.path.join(work_dir, "ancestral_states.csv")
    _save_states_csv(tree, node_states, posteriors, states_csv)

    # 11. ESM concordance (post-hoc)
    esm_tms = _load_esm_cluster_tm(pair_id)
    esm_agree = compute_concordance(af2_prefs, esm_tms, delta)

    return {
        "pair_id": pair_id,
        "n_clusters": n_clusters,
        "n_f1": pref_counts.get("F1", 0),
        "n_f2": pref_counts.get("F2", 0),
        "n_amb": pref_counts.get("Amb", 0),
        "n_transitions": n_trans,
        "n_fold_transitions": n_fold_trans,
        "pattern": pattern,
        "root_state": root_state,
        "p_value": p_value,
        "method": method,
        "esm_agree_pct": round(esm_agree * 100, 1) if esm_agree is not None else None,
    }


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _name_internal_nodes(tree):
    """Assign names to unnamed internal nodes: N0, N1, ..."""
    idx = 0
    for node in tree.traverse("preorder"):
        if not node.name:
            node.name = f"N{idx}"
            idx += 1


def _save_states_csv(tree, node_states, posteriors, path):
    """Save per-node ancestral states to CSV."""
    rows = []
    for node in tree.traverse("preorder"):
        key = node.name if node.name else ""
        state = node_states.get(key, "")
        row = {"node": key, "is_leaf": node.is_leaf(), "state": state}
        if posteriors and key in posteriors:
            for s in ["F1", "F2", "Amb"]:
                row[f"posterior_{s}"] = posteriors[key].get(s, "")
        rows.append(row)
    pd.DataFrame(rows).to_csv(path, index=False)
    print(f"  Saved ancestral states: {path}")
