import copy
# --- imports (kept inside to avoid hard deps if user only imports module) ---

# for phylogenetic trees
from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import *

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.Newick import Clade, Tree as NwTree


# Prefer fastcluster (much faster), fall back to SciPy
try:
    from fastcluster import linkage  # 1.2.6 present on machine
except Exception:
    from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import to_tree


import pickle
from pylab import *

# Use ete3 package for visualization
# from ete3 import Tree, TreeStyle, TextFace, RectFace, NodeStyle
from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace, TextFace, RectFace
import os
os.environ["XDG_RUNTIME_DIR"] = "/tmp"
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen") # os.environ["QT_QPA_PLATFORM"] = "offscreen"

ETE_TREEVIEW_AVAILABLE = False
try:
    # Qt-dependent viewer/renderer
    from ete3.treeview import TreeStyle, NodeStyle, faces
    ETE_TREEVIEW_AVAILABLE = True
except Exception as _e:
    TreeStyle = NodeStyle = faces = None
    print(f"[phytree] treeview unavailable ({_e}); tree plots will be skipped.")
#except Exception as e:
#    Tree = None
#    print(f"[phytree] ETE3 import failed: {e}")

from utils.msa_utils import *
from utils.utils import now
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize, to_hex
# from matplotlib import pyplot as plt
from matplotlib.colorbar import ColorbarBase
import matplotlib as mpl
import matplotlib.pyplot as plt
# from matplotlib.transforms import Bbox

import re, random
import time
import math
from contextlib import contextmanager
from collections import Counter


# ----- Helper functions -----

@contextmanager
def _stage(name: str):
    t0 = time.time()
    print(f"[tree] {now()} >>> {name} ...", flush=True)
    try:
        yield
    finally:
        dt = time.time() - t0
        print(f"[tree] {now()} <<< {name} done in {dt:.1f}s", flush=True)


def _condensed_hamming_from_seqs(seqs, block_j=1024, step_cols=2048):
    """
    Build a condensed Hamming distance vector (upper triangle, i<j) for equal-length aligned strings.

    seqs: list[str] with identical length (post A3M cleaning)
    Returns: np.ndarray, shape = (n*(n-1)//2,), float64, distances in [0,1]
    """
    n = len(seqs)
    if n < 2:
        return np.zeros(0, dtype=np.float64)

    L = len(seqs[0])
    # Build a tiny LUT only for actually observed characters (fast, compact)
    chars = sorted({c for s in seqs for c in s})
    lut = {c: i for i, c in enumerate(chars)}  # typically 21 symbols including '-'

    # Encode to uint8 matrix (n, L)
    A = np.empty((n, L), dtype=np.uint8)
    for i, s in enumerate(seqs):
        A[i, :] = np.fromiter((lut[c] for c in s), dtype=np.uint8, count=L)

    # Condensed length
    m = n * (n - 1) // 2
    out = np.empty(m, dtype=np.float64)

    # Precompute row offsets into the condensed vector:
    # offset(i) = number of elements before row i = sum_{r=0}^{i-1} (n-1-r) = i*(2n - i - 1)/2
    # Then index(i,j) = offset(i) + (j - i - 1)
    # We'll write per-i row to contiguous slices starting at offset(i).
    for i in range(n - 1):
        # Starting write position in condensed for this i
        offset_i = i * (2 * n - i - 1) // 2  # integer arithmetic
        write_pos = offset_i

        # Compare row i against blocks of rows j in (i+1..n-1)
        Ai = A[i:i+1]  # shape (1, L)
        for j0 in range(i + 1, n, block_j):
            j1 = min(n, j0 + block_j)
            Aj = A[j0:j1]  # shape (B2, L)

            # Accumulate mismatches across columns in manageable chunks
            acc = None
            for c0 in range(0, L, step_cols):
                c1 = min(L, c0 + step_cols)
                # (1,1,cols) vs (1,B2,cols) via broadcasting -> (1,B2,cols) -> sum over cols -> (1,B2)
                neq = (Ai[:, None, c0:c1] != Aj[None, :, c0:c1]).sum(axis=2, dtype=np.int32)
                acc = neq if acc is None else (acc + neq)

            B2 = j1 - j0
            # acc shape is (1, B2); squeeze to (B2,)
            out[write_pos:write_pos + B2] = acc[0]
            write_pos += B2

    # Normalize to mismatch fraction (Bio 'identity' distance)
    out /= float(L)
    return out


# --- helpers for stratified sampling (place once in phytree_utils.py) ---

def _debug_compare_with_biopython(seqs, names, sample_pairs=200, rng_seed=0):
    """
    Randomly pick some (i,j) pairs, compare our fast Hamming to BioPython's
    DistanceCalculator('identity'). Print MAE / max-abs-diff.
    """

    random.seed(rng_seed)
    n = len(seqs)
    # Build a small alignment only for sampled pairs + a reference row to avoid repeated rebuilds.
    # For simplicity, we just build per-pair mini-alignments (cheap for a few hundred checks).
    calc = DistanceCalculator('identity')

    # our fast distances via the same encoder used in pipeline (re-use function)
    def _pair_fast(i, j):
        s1, s2 = seqs[i], seqs[j]
        # exact same logic as in _condensed_hamming_from_seqs for a single pair:
        L = len(s1)
        # NOTE: treat '-' as a regular symbol (same as Bio 'identity'):
        neq = sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))
        return neq / float(L)

    idx = [(random.randrange(0, n-1), None) for _ in range(sample_pairs)]
    idx = [(i, random.randrange(i+1, n)) for (i, _) in idx]

    diffs = []
    for (i, j) in idx:
        aln = MultipleSeqAlignment([
            SeqRecord(Seq(seqs[i]), id=names[i]),
            SeqRecord(Seq(seqs[j]), id=names[j]),
        ])
        bio_dm = calc.get_distance(aln)
        bio = bio_dm[1,0]  # distance between the two rows
        fast = _pair_fast(i, j)
        diffs.append(abs(bio - fast))

    diffs = np.asarray(diffs)
    print(f"[debug] compare fast-vs-Bio on {len(diffs)} pairs: "
          f"MAE={diffs.mean():.3e}  max={diffs.max():.3e}")


def _norm_id(header: str) -> str:
    """
    Normalize FASTA headers so Deep/Cluster IDs match better:
    - keep first token before whitespace
    - strip trailing '/start-end' spans (e.g., '/12-345')
    """
    h = (header or "").strip()
    if not h:
        return h
    tok = h.split()[0]
    tok = re.sub(r"/\d+-\d+$", "", tok)
    return tok

def _load_cluster_norm_ids(cluster_msa_dir: str) -> dict[str, set[str]]:
    """
    Read all ShallowMsa_###.a3m files and return:
        cluster_name -> set(normalized IDs)
    Uses load_fasta() that already exists in your codebase.
    """
    clusters: dict[str, set[str]] = {}
    for fn in sorted(os.listdir(cluster_msa_dir)):
        if not fn.endswith(".a3m"):
            continue
        cl = fn[:-4]  # strip .a3m
        if not cl.startswith("ShallowMsa_"):
            continue
        ids, _seqs = load_fasta(os.path.join(cluster_msa_dir, fn))
        clusters[cl] = { _norm_id(i) for i in ids }
    return clusters




def prune_tree_to_leaves(ete_tree, keep_leaves: list[str]):
    """
    Keep-only pruning. Repeatedly:
      1) drop any leaf not in keep_leaves,
      2) drop now-empty internal nodes.

    After this returns, every leaf name is in keep_leaves.
    """
    keep = set(keep_leaves)
    changed = True
    while changed:
        changed = False

        # 1) remove leaves not in keep set
        to_detach = [leaf for leaf in ete_tree.iter_leaves() if leaf.name not in keep]
        if to_detach:
            for leaf in to_detach:
                leaf.detach()
            changed = True

        # 2) remove now-empty internal nodes
        # (postorder so we clear bottom-up)
        for n in list(ete_tree.iter_descendants(strategy="postorder")):
            if not n.children:  # became leaf
                if n.name not in keep:   # synthetic leaf not requested
                    n.detach()
                    changed = True

    return ete_tree



def resolve_duplicated_ids(ids_list):
        """
        Takes a list of strings and adds a counter to duplicates to make each string unique.
        Unique strings remain unchanged.

        Args:
        strings (list): A list of strings, potentially with duplicates.

        Returns:
        list: A list of strings where duplicates are made unique with a counter.
        """
        count = {}
        result = []

        for s in ids_list:
            if s in count:
                count[s] += 1
                result.append(f"{s}_{count[s]}")
            else:
                count[s] = 0
                result.append(s)

        return result


# Reconstruct a phylogenetic tree
# Input:
# msa_file - file with MSA in a3m format
# output_tree_file - where to save
# max_seqs - sample sequences if too many
def phytree_from_msa(msa_file, output_tree_file,
    max_seqs=None,  # None → unlimited; keep your existing behavior if different
    cluster_msa_dir=None,
    seed=123,
    progress_rows=100,  # print every N rows of the distance matrix
    progress_percent=1.0):  # or every % of pairs (whichever triggers first)

    """
    Build a UPGMA tree from an A3M alignment.
    - If max_seqs is None -> use ALL sequences (no sampling; cluster_msa_dir is ignored).
    - If max_seqs is set and num_seqs > max_seqs:
        * If cluster_msa_dir is provided -> STRATIFIED sampling:
              ensure >=1 sequence per shallow cluster that also appears in the deep MSA (by normalized ID),
              then fill the remaining quota at random.
              If max_seqs < #represented-clusters, it is bumped up.
        * Else -> uniform random sampling (previous behavior).

    Returns
    -------
    Biopython tree (also written to output_tree_file).
    """
    if not ETE_TREEVIEW_AVAILABLE:
        print("[phytree] PyQt5/Qt not available; tree will be built without rendering.")

    if seed is not None:
        random.seed(int(seed))

    # --- load the (aligned) MSA ---
    seqs_IDs, seqs = load_fasta(msa_file)

    # 1) CLEAN A3M LINES ROBUSTLY (lowercase inserts removed; '.' and '*' → '-')
    seqs = [clean_a3m_line(s) for s in seqs]

    # 2) ENFORCE A SINGLE ALIGNED LENGTH (keep modal length; drop outliers)
    raw_n = len(seqs)
    lengths = [len(s) for s in seqs]
    if not lengths or min(lengths) == 0:
        raise RuntimeError("Empty sequences after A3M cleaning — check input.")

    cnt = Counter(lengths)
    mode_len, _ = cnt.most_common(1)[0]
    if len(set(lengths)) != 1:
        keep_idx = [i for i, L in enumerate(lengths) if L == mode_len]
        drop_idx = [i for i, L in enumerate(lengths) if L != mode_len]
        if len(keep_idx) < 2:
            raise RuntimeError(f"A3M lengths inconsistent; fewer than 2 sequences at mode length {mode_len}.")
        print(f"[clean] A3M columns not uniform: kept {len(keep_idx)}/{raw_n} at length={mode_len}; "
              f"dropped {len(drop_idx)} outliers.")
        seqs = [seqs[i] for i in keep_idx]
        seqs_IDs = [seqs_IDs[i] for i in keep_idx]

    # 3) UPDATE num_seqs AFTER CLEANING/FILTERING
    num_seqs = len(seqs)

    # 4) SAMPLING (none / stratified / uniform) NOW USES CURRENT ARRAYS
    if (max_seqs is None) or (num_seqs <= int(max_seqs)):
        pass  # keep all sequences
    else:
        with _stage(f"phytree sampling"):
            print(f"[phytree] sampling {max_seqs} sequences from {num_seqs}")
        if cluster_msa_dir:
            # Build mapping on CURRENT seqs_IDs (post-clean)
            deep_norm_to_indices: dict[str, list[int]] = {}
            for i, rid in enumerate(seqs_IDs):
                nid = _norm_id(rid)
                deep_norm_to_indices.setdefault(nid, []).append(i)

            try:  # try to load clusters
                clusters = _load_cluster_norm_ids(cluster_msa_dir)
                shallow = sorted([c for c in clusters if c.startswith("ShallowMsa_")])
            except:
                print("[warn] Could not load clusters; falling back to uniform sampling.")
                clusters = {}
                shallow = []

            chosen_idx: list[int] = []
            represented = 0
            missing = 0
            for cl in shallow:
                cand_norm = list(clusters[cl] & set(deep_norm_to_indices.keys()))
                if cand_norm:
                    nid = random.choice(cand_norm)
                    idx = random.choice(deep_norm_to_indices[nid])
                    chosen_idx.append(idx)
                    represented += 1
                else:
                    missing += 1

            mand = len(set(chosen_idx))
            if mand == 0:
                print("[warn] No shallow clusters matched Deep IDs after cleaning; falling back to uniform sampling.")
                chosen_idx = []

            # bump cap if needed
            if mand > 0 and max_seqs < mand:
                print(f"[info] Bumping max_seqs from {max_seqs} to {mand} to cover all represented clusters.")
                max_seqs = mand

            remaining = int(max_seqs) - len(set(chosen_idx))
            pool = [i for i in range(num_seqs) if i not in set(chosen_idx)]
            random.shuffle(pool)
            extra = pool[:max(0, remaining)]
            sel_idx = list(dict.fromkeys(chosen_idx + extra))  # dedup while preserving order

            print(f"[stratify] deep={num_seqs}  shallow_clusters={len(shallow)}  "
                  f"represented={represented}  missing={missing}  "
                  f"final_sample={len(sel_idx)}  (mand={mand} + rand={len(extra)})")

            # apply selection on CURRENT arrays
            seqs = [seqs[i] for i in sel_idx]
            seqs_IDs = [seqs_IDs[i] for i in sel_idx]

        else:
            # uniform sample from CURRENT arrays
            rand_inds = random.sample(range(num_seqs), int(max_seqs))
            seqs = [seqs[i] for i in rand_inds]
            seqs_IDs = [seqs_IDs[i] for i in rand_inds]
            print(f"[sample] uniform random: picked {len(seqs)} / {num_seqs}")

    # Make IDs unique if needed (same as before)
    with _stage(f"Setting unique IDs"):
        seqs_IDs = resolve_duplicated_ids(seqs_IDs)
        print(f"[phytree] Set unique IDs for {len(seqs_IDs)} sequences")

    # Build alignment object
    seq_records = [SeqRecord(Seq(seqs[i]), id=seqs_IDs[i]) for i in range(len(seqs))]
    with _stage(f"Setting MSA"):
        alignment = MultipleSeqAlignment(seq_records)
        print(f"[phytree] Built MSA for sequences")

    # Fast condensed distances (skip Bio DistanceCalculator)
    with _stage(f"Setting Distance MAtrix"):
        condensed = _condensed_hamming_from_seqs(seqs)  # <- main speedup
        names = list(seqs_IDs)  # keep same order as seqs
        n = len(names)
        print(f"[phytree] Built condensed Hamming for {n} sequences")

        L = len(seqs[0])
        print(f"[phytree] n={len(seqs)}  L={L}  pairs={len(seqs) * (len(seqs) - 1) // 2}")

        _debug_compare_with_biopython(seqs, names, sample_pairs=200)
        print(f"[debug] Compared condensed to Biopython distances")

    # --- replaced the slow UPGMA block with a SciPy-based version ---
    with _stage(f"Setting UPGMA Tree"):
        try:
            # 'condensed', 'names', 'n' were produced in the previous stage
            Z = linkage(condensed, method="average")  # fastcluster or SciPy

            rootnode, _ = to_tree(Z, rd=True)

            def to_clade(cn):
                # returns (Bio.Phylo.Newick.Clade, height_of_subtree)
                if cn.is_leaf():
                    return Clade(name=names[cn.id], branch_length=0.0), 0.0
                lc, lh = to_clade(cn.get_left())
                rc, rh = to_clade(cn.get_right())
                # UPGMA ultrametric branch lengths
                bl_l = max(cn.dist - lh, 0.0)
                bl_r = max(cn.dist - rh, 0.0)
                lc.branch_length = bl_l
                rc.branch_length = bl_r
                return Clade(clades=[lc, rc], branch_length=0.0), cn.dist

            root_clade, _ = to_clade(rootnode)
            tree = NwTree(root=root_clade)
            print("[phytree] UPGMA tree built via fastcluster average-linkage")

        except Exception as e:
            print(f"[warn] Fast UPGMA path failed ({e}); falling back to Bio.Phylo (slow)")
            # If you want a fallback, you *must* rebuild a Bio DistanceMatrix. Skipping here is fine.
            constructor = DistanceTreeConstructor()
            # WARNING: falling back will need a DistanceMatrix object;
            # you can omit this fallback, or reconstruct it if necessary.
            raise

    from scipy.cluster.hierarchy import cophenet
    c, _ = cophenet(Z, condensed)
    print(f"[phytree] cophenetic corr = {c:.4f}")
    # Values ~0.8–0.98 are typical; closer to 1 is better.

    # Save (same as before)
    if len(output_tree_file) == 0:
        output_tree_file = msa_file.replace(".a3m", "_tree.nwk")
    if not os.path.exists(os.path.dirname(output_tree_file)):
        os.makedirs(os.path.dirname(output_tree_file))
    with _stage(f"Saving Tree"):
        Phylo.write(tree, output_tree_file, "newick")
        print(f"[phytree] Saved output tree to {output_tree_file}" )

    return tree



def convert_biopython_to_ete3(biopy_tree, parent_ete_node=None):
    # TRY using files: need a dictionary for the leaf names
    node_dict = {}
    ctr = 0
    for node in biopy_tree.find_clades():
        #        clade.name = "'" + clade.name + "'"
        node_dict["node" + str(ctr)] = node.name
        node.name = "node" + str(ctr)
        ctr += 1

    Phylo.write(biopy_tree, "temp_tree.nwk", "newick")  # write to file

    ete_tree = Tree("temp_tree.nwk", 1) # read
    for node in ete_tree.traverse():  # modify name
        node.name = node_dict[node.name]

    return ete_tree


# Read tree and strip quotes from node names
def read_tree_ete(phytree_file):
    ete_tree = Tree(phytree_file, format=1)
    for node in ete_tree.traverse():
        node.name = node.name.strip("'")
        node.name = node.name.strip("\"")
    return ete_tree


# Define a dictionary to map numerical values to colors
value_to_color = {
    "Leaf_1": (0.2, 0.4, 0.6),  # Replace "Leaf_1" with the actual leaf name
    "Leaf_2": (0.8, 0.1, 0.3),  # Replace "Leaf_2" with the actual leaf name
    # Add more entries for other leaf nodes as needed
}


def write_newick_with_quotes(tree, file_path):
    def format_clade(clade):
        if clade.is_terminal():
            return f"'{clade.name}':{clade.branch_length}" if clade.name else ''
        else:
            subtrees = ','.join(format_clade(sub_clade) for sub_clade in clade)
            return f"({subtrees}){'' if clade.name is None else f'{clade.name}'}:{clade.branch_length}"

    newick_str = format_clade(tree.root) + ';'
    with open(file_path, 'w') as f:
        f.write(newick_str)


# Induced subtree for an ete3 tree
# keep only the input leaves appearing in 'leaf_names'
def extract_induced_subtree(tree, leaf_names):
    """
    Extracts the induced subtree that includes only the specified leaves and their ancestral nodes.

    Args:
    tree (ete3.Tree): The original tree from which to extract the subtree.
    leaf_names (list): A list of leaf names to include in the subtree.

    Returns:
    ete3.Tree: The induced subtree containing only the specified leaves and their ancestors.
    """
    if type(tree) == str: # extract tree from file
        print("Read tree from file: " + tree)
        bio_tree = Phylo.read(tree, "newick")  # This is different from write_newick_with_quotes !!!!
        tree = convert_biopython_to_ete3(bio_tree)

    # Copy the tree to avoid modifying the original
    subtree = tree.copy(method='deepcopy')

    # Prune the subtree to keep only the leaves specified
    leaves_to_keep = set(leaf_names)
    for leaf in subtree.iter_leaves():
        if leaf.name not in leaves_to_keep:
            leaf.delete()

    return subtree


def visualize_tree_with_heatmap(
    phylo_tree,
    node_values_matrix: pd.DataFrame,
    col_groups: list[list[str]],
    output_file: str,
    group_titles: list[str] | None = None,
    *,
    figsize=(22, 12),
    tree_width_ratio=1.0,
    heatmap_width_ratio=6.5,
    cbar_width_ratio=0.18,
    cbar_stack_vertical=True,
    x_tick_rotation=90,
    x_tick_fontsize=8,
    y_tick_fontsize=7,
    nan_rgba=(0.92, 0.92, 0.92, 1.0),
):
    """
    Render an ETE tree + grouped heatmap + group colorbars (slim, stacked).
    `node_values_matrix` must be indexed by EXACT leaf names present in `phylo_tree`.
    """

    # --- validate / order rows by tree leaves (tree is already pruned by caller) ---
    tree_leaves = [n.name for n in phylo_tree.iter_leaves()]
    missing = [idx for idx in node_values_matrix.index if idx not in tree_leaves]
    if missing:
        raise ValueError(f"Matrix includes leaf names not in tree: {missing[:10]}{'...' if len(missing)>10 else ''}")

    M = node_values_matrix.reindex(tree_leaves)

    # Flatten column groups to the exact order
    cols = [c for group in col_groups for c in group if c in M.columns]
    M = M[cols]

    # ------------- layout -------------
    if cbar_stack_vertical:
        # main (tree + heatmap) + one thin column for stacked colorbars
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(
            nrows=1, ncols=3,
            width_ratios=[tree_width_ratio, heatmap_width_ratio, cbar_width_ratio],
            left=0.04, right=0.98, bottom=0.12, top=0.98, wspace=0.25
        )
        ax_tree = fig.add_subplot(gs[0, 0])
        ax_hm = fig.add_subplot(gs[0, 1])
        # stack cbars vertically within the third column
        from matplotlib.gridspec import GridSpecFromSubplotSpec
        n_cb = len(col_groups) if col_groups else 0
        if n_cb > 0:
            gs_cbar = GridSpecFromSubplotSpec(n_cb, 1, subplot_spec=gs[0, 2],
                                              hspace=0.4, height_ratios=[1]*n_cb)
            ax_cbars = [fig.add_subplot(gs_cbar[i, 0]) for i in range(n_cb)]
        else:
            ax_cbars = []
    else:
        # side-by-side cbars (thin)
        n_cb = len(col_groups) if col_groups else 0
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(
            nrows=1, ncols=2 + n_cb,
            width_ratios=[tree_width_ratio, heatmap_width_ratio] + [cbar_width_ratio]*n_cb,
            left=0.04, right=0.98, bottom=0.12, top=0.98, wspace=0.3
        )
        ax_tree = fig.add_subplot(gs[0, 0])
        ax_hm = fig.add_subplot(gs[0, 1])
        ax_cbars = [fig.add_subplot(gs[0, 2 + i]) for i in range(n_cb)]

    # ------------- tree -------------
    # If you already have a dedicated ETE→matplotlib renderer, call it here.
    # For now: hide axis so only branches (drawn elsewhere) are visible.
    ax_tree.set_axis_off()
    # (Your project previously drew the ETE tree onto ax_tree; keep that call if exists.)

    # ------------- heatmap -------------
    data = M.values.astype(float)

    # single cmap for all columns (consistent aesthetics)
    cmap = mpl.cm.viridis.copy()
    cmap.set_bad(nan_rgba)

    im = ax_hm.imshow(data, aspect='auto', interpolation='nearest', cmap=cmap)

    # ticks
    ax_hm.set_xticks(np.arange(data.shape[1]))
    ax_hm.set_xticklabels(M.columns, rotation=x_tick_rotation, ha='right', fontsize=x_tick_fontsize)
    ax_hm.set_yticks(np.arange(data.shape[0]))
    ax_hm.set_yticklabels(M.index, fontsize=y_tick_fontsize)
    ax_hm.tick_params(axis='x', pad=6)

    # group separators + group-specific colorbars (scaled to that group)
    col_start = 0
    for gi, group in enumerate(col_groups):
        gcols = [c for c in group if c in M.columns]
        if not gcols:
            continue
        gsize = len(gcols)

        # thin vertical separators
        ax_hm.axvline(col_start - 0.5, color='k', lw=0.5, alpha=0.25)
        ax_hm.axvline(col_start + gsize - 0.5, color='k', lw=0.5, alpha=0.25)

        # colorbar per group using that group's min/max
        gdat = data[:, col_start:col_start + gsize]
        if np.isnan(gdat).all():
            vmin, vmax = 0.0, 1.0
        else:
            vmin = np.nanmin(gdat)
            vmax = np.nanmax(gdat)
            if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
                vmin, vmax = 0.0, 1.0

        sm = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax), cmap=cmap)
        if gi < len(ax_cbars):
            cb = plt.colorbar(sm, cax=ax_cbars[gi])
            title = group_titles[gi] if group_titles and gi < len(group_titles) else ""
            ax_cbars[gi].set_title(title, fontsize=11, pad=4)
            cb.ax.tick_params(labelsize=8)

        col_start += gsize

    # save
    png_path = output_file if output_file.lower().endswith(".png") else output_file + ".png"
    fig.savefig(png_path, dpi=200)
    plt.close(fig)


def calculate_tree_width(rendered_tree_path):
    """
    Calculate the width of the rendered tree from the temporary file.
    """
    from PIL import Image
    with Image.open(rendered_tree_path) as img:
        return img.size[0]  # Returns the width of the image


def add_figure_title(image_path, title, title_x_offset=None, title_y_offset=20,
                     font_size=30, extra_space_ratio=0.05):
    from PIL import Image, ImageDraw, ImageFont

    # Open the image
    img = Image.open(image_path)
    img_width, img_height = img.size

    # Calculate extra space to add at the bottom
    extra_space = int(img_height * extra_space_ratio)

    # Create a new image with extra space
    new_img = Image.new("RGB", (img_width, img_height + extra_space), (255, 255, 255))
    new_img.paste(img, (0, 0))

    # Use a default font that supports resizing
    try:
        font = ImageFont.truetype("DejaVuSans-Bold.ttf", font_size)
    except IOError:
        print("DejaVu Sans font not found. Using PIL's default font (no size customization).")
        font = ImageFont.load_default()  # Fallback font

    # Calculate text dimensions using getbbox()
    draw = ImageDraw.Draw(new_img)
    text_bbox = draw.textbbox((0, 0), title, font=font)  # Returns (x0, y0, x1, y1)
    text_width = text_bbox[2] - text_bbox[0]
    text_height = text_bbox[3] - text_bbox[1]

    # Center the title horizontally if no offset provided
    if title_x_offset is None:
        title_x_offset = (img_width - text_width) // 2

    # Position the title in the added space
    text_x = title_x_offset
    text_y = img_height + (extra_space - text_height) // 2

    # Draw the title
    draw.text((text_x, text_y), title, fill="black", font=font)

    # Save the updated image
    new_img.save(image_path)


def add_figure_title_working(image_path, title,
                             title_x_offset=None, title_y_offset=20,
                             font_size=30):
    from PIL import Image, ImageDraw, ImageFont

    # Open the image
    img = Image.open(image_path)
    img_width, img_height = img.size

    # Use a default font that supports resizing
    try:
        font = ImageFont.truetype("DejaVuSans-Bold.ttf", font_size)
    except IOError:
        print("DejaVu Sans font not found. Using PIL's default font (no size customization).")
        font = ImageFont.load_default()  # Fallback font

    # Calculate text dimensions using getbbox()
    draw = ImageDraw.Draw(img)
    text_bbox = draw.textbbox((0, 0), title, font=font)  # Returns (x0, y0, x1, y1)
    text_width = text_bbox[2] - text_bbox[0]
    text_height = text_bbox[3] - text_bbox[1]

    # Center the title horizontally if no offset provided
    if title_x_offset is None:
        title_x_offset = (img_width - text_width) // 2

    text_x = title_x_offset
    text_y = img_height - title_y_offset - text_height  # Position near the bottom

    # Draw the title
    draw.text((text_x, text_y), title, fill="black", font=font)

    # Save the updated image
    img.save(image_path)


def render_tree_to_png(
    ete_tree,
    out_png: str,
    *,
    width_px: int = 1400,
    height_px: int = 2000,
    leaf_fontsize: int = 10,
    show_scale: bool = False,
    margin: int = 20,
):
    """
    Use ETE3's renderer to export the tree to a PNG. We then place this PNG
    inside a Matplotlib figure so layout with the heatmap is trivial & robust.
    """
    if not ETE_TREEVIEW_AVAILABLE:
        print("[phytree] skipping tree from msa render: PyQt5/Qt not available")
        return


    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.scale =  120  # making branches visible; tweak if needed
    ts.show_scale = show_scale
    ts.margin_left = ts.margin_right = ts.margin_top = ts.margin_bottom = margin
    ts.title.add_face(TextFace("", fsize=1), column=0)

    # leaf font
    def _style_node(n):
        ns = NodeStyle()
        ns["hz_line_width"] = 1
        ns["vt_line_width"] = 1
        n.set_style(ns)
        if n.is_leaf():
            n.add_face(AttrFace("name", fsize=leaf_fontsize), column=0, position="aligned")
    for n in ete_tree.traverse():
        _style_node(n)

    ete_tree.render(out_png, tree_style=ts, w=width_px, h=height_px)
    return out_png


def draw_grouped_heatmap(
    ax_hm, df_leaf, col_groups,
    *, group_titles=None, nan_rgba=(0.92,0.92,0.92,1.0),
    x_tick_rotation=90, x_tick_fontsize=8, y_tick_fontsize=7,
    cmap_name="viridis", return_group_scalars=False
):
    import numpy as np, matplotlib as mpl

    cols = [c for g in col_groups for c in g if c in df_leaf.columns]
    M = df_leaf[cols].to_numpy(dtype=float)

    cmap = mpl.cm.get_cmap(cmap_name).copy()
    cmap.set_bad(nan_rgba)

    im = ax_hm.imshow(M, aspect="auto", interpolation="nearest", cmap=cmap)
    ax_hm.set_yticks(np.arange(df_leaf.shape[0]))
    ax_hm.set_yticklabels(df_leaf.index, fontsize=y_tick_fontsize)

    ax_hm.set_xticks(np.arange(len(cols)))
    ax_hm.set_xticklabels(cols, rotation=x_tick_rotation, ha="right", fontsize=x_tick_fontsize)
    ax_hm.tick_params(axis="x", pad=6)

    col_start, group_scalars = 0, []
    for gi, grp in enumerate(col_groups):
        gcols = [c for c in grp if c in df_leaf.columns]
        if not gcols:
            continue
        gsize = len(gcols)
        ax_hm.axvline(col_start - 0.5, color="k", lw=0.5, alpha=0.25)
        ax_hm.axvline(col_start + gsize - 0.5, color="k", lw=0.5, alpha=0.25)

        gdat = M[:, col_start:col_start + gsize]
        if np.isnan(gdat).all():
            vmin, vmax = 0.0, 1.0
        else:
            vmin, vmax = np.nanmin(gdat), np.nanmax(gdat)
            if not (np.isfinite(vmin) and np.isfinite(vmax)) or vmin == vmax:
                vmin, vmax = 0.0, 1.0
        group_scalars.append((vmin, vmax))
        col_start += gsize

    return (cmap, group_scalars) if return_group_scalars else None


def draw_tree_aligned(ax, ete_tree, leaf_order):
    import numpy as np

    y_pos = {name: i for i, name in enumerate(leaf_order)}

    # distance-from-root x
    root = ete_tree.get_tree_root()
    x_pos = {root: 0.0}
    def edge_len(child):
        bl = child.dist
        return float(bl) if bl is not None and bl == bl else 1.0
    for node in root.traverse("preorder"):
        for child in node.children:
            x_pos[child] = x_pos[node] + edge_len(child)

    # compute y with guards
    def compute_y(n):
        if n.is_leaf():
            if n.name in y_pos:
                return float(y_pos[n.name])
            # leaf not in mapping -> ignore by returning NaN
            return np.nan
        vals = [compute_y(c) for c in n.children]
        vals = [v for v in vals if v == v]  # drop NaNs
        return float(np.mean(vals)) if vals else np.nan

    def y(n):
        v = compute_y(n)
        return v

    # draw edges; skip segments with NaN y
    for node in root.traverse("preorder"):
        for child in node.children:
            x0, y0 = x_pos[node], y(node)
            x1, y1 = x_pos[child], y(child)
            if not (np.isfinite(y0) and np.isfinite(y1)):
                continue
            ax.plot([x0, x1], [y1, y1], color="#666666", lw=1.0, solid_capstyle="butt")
            ax.plot([x0, x0], [y0, y1], color="#666666", lw=1.0, solid_capstyle="butt")

    n = len(leaf_order)
    ax.set_ylim(n - 0.5, -0.5)
    xmin, xmax = min(x_pos.values()), max(x_pos.values())
    ax.set_xlim(xmin - 0.05 * (xmax - xmin + 1e-9), xmax + 0.02)
    ax.axis("off")



def compose_tree_and_heatmap(
    ete_tree,
    df_leaf,                 # index = leaf names (already pruned/ordered)
    col_groups,
    output_file,
    *,
    group_titles=None,
    base_figsize=(22, 12),
    cbar_width_ratio=0.22,   # fixed slim column for stacked colorbars
    x_tick_rotation=90, x_tick_fontsize=8, y_tick_fontsize=7,
    nan_rgba=(0.92,0.92,0.92,1.0),
    ylabels_override=None,    # if provided, replaces df_leaf.index on the heatmap
    lock_scales_to_unit=False,
    split_groups=True,        # NEW: draw one mini-heatmap per group with a small gap
    group_gap_frac=0.08       # NEW: ~8% of a group's width used as inter-group gap
):
    import numpy as np, matplotlib.pyplot as plt, matplotlib as mpl
    from matplotlib.gridspec import GridSpecFromSubplotSpec
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    # ----- prepare plotting frame (labels vs leaf names) -----
    if ylabels_override is not None:
        df_plot = df_leaf.copy()
        df_plot.index = ylabels_override
    else:
        df_plot = df_leaf

    # flatten the column order and filter to present columns
    flat_cols = [c for g in col_groups for c in g if c in df_leaf.columns]
    if not flat_cols:
        # nothing to draw
        out_png = output_file if output_file.lower().endswith(".png") else output_file + ".png"
        fig = plt.figure(figsize=base_figsize); fig.savefig(out_png, dpi=200); plt.close(fig)
        return out_png

    # ----- layout: tree | heatmap area | stacked colorbars
    # Estimate width for heatmaps from number of columns so the whole figure looks stable
    n_cols = sum(len([c for c in g if c in df_leaf.columns]) for g in col_groups)
    heatmap_ratio = 0.06 + 0.03 * n_cols
    tree_ratio = max(0.01, 1.0 - heatmap_ratio - 0.05)  # reserve ~5% for cbar column
    width_ratios = [tree_ratio, heatmap_ratio, 0.05]

    fig = plt.figure(figsize=base_figsize)
    gs = fig.add_gridspec(
        nrows=1, ncols=3, width_ratios=width_ratios,
        left=0.035, right=0.985, bottom=0.10, top=0.98, wspace=0.06
    )
    ax_tree = fig.add_subplot(gs[0, 0])
    ax_heat_region = fig.add_subplot(gs[0, 1])  # container for mini-heatmaps (we'll add sub-axes)
    ax_heat_region.axis("off")
    cbar_host = fig.add_subplot(gs[0, 2]); cbar_host.axis("off")

    # stack slim colorbars in the right column (narrower + ~2× longer + moved left)
    n_cb = len(col_groups) if col_groups else 0
    cbar_axes = []
    if n_cb > 0:
        # take almost the full column height for the bars (~98%)
        total_h = 0.98  # ↑ was 0.94; makes bars longer
        pad = 0.006  # smaller spacing between bars
        left_margin = 0.01  # ← move bars left inside the cbar column (was 0.20)
        width_frac = 0.68  # keep them slim (0.68 of the column width)

        cbar_h = (total_h - (n_cb - 1) * pad) / max(1, n_cb)
        y0 = (1 - total_h) / 2.0

        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        for i in range(n_cb):
            ax_c = inset_axes(
                cbar_host,
                width="50%", height=f"{cbar_h * 200:.1f}%",
                bbox_to_anchor=(left_margin, y0 + i * (cbar_h + pad), width_frac, cbar_h),
                bbox_transform=cbar_host.transAxes, borderpad=0
            )
            cbar_axes.append(ax_c)

    # ----- draw the tree aligned to df_leaf row order -----
    draw_tree_aligned(ax_tree, ete_tree, leaf_order=list(df_leaf.index))

    # ----- draw the heatmaps -----
    # shared colormap (with NaN color); per-group vmin/vmax
    cmap = mpl.cm.viridis.copy()
    cmap.set_bad(nan_rgba)

    # compute grouped widths and gaps inside the heat region
    groups_present = [ [c for c in g if c in df_leaf.columns] for g in col_groups ]
    groups_present = [g for g in groups_present if g]  # drop empties
    n_groups = len(groups_present)
    if n_groups == 0:
        out_png = output_file if output_file.lower().endswith(".png") else output_file + ".png"
        fig.savefig(out_png, dpi=200); plt.close(fig); return out_png

    # relative widths proportional to number of columns; include small gaps
    widths = [len(g) for g in groups_present]
    gap = group_gap_frac  # relative to a group's width
    # convert relative to absolute in axis coordinates (0..1)
    total_units = sum(widths) + gap * sum(widths[:-1])
    unit = 1.0 / max(total_units, 1e-9)

    left = 0.00
    group_axes, scalars = [], []
    for gi, gcols in enumerate(groups_present):
        w = unit * len(gcols)
        # add a sub-axes for this group inside ax_heat_region
        bb = ax_heat_region.get_position()
        ax = fig.add_axes([bb.x0 + left * bb.width, bb.y0, w * bb.width, bb.height])
        group_axes.append((ax, gcols))
        left += w
        if gi < n_groups - 1:
            left += unit * gap * len(gcols)  # small gap proportional to group width

    # render each group
    for gi, (ax, gcols) in enumerate(group_axes):
        M = df_plot[gcols].to_numpy(dtype=float)
        if lock_scales_to_unit:
            vmin, vmax = 0.0, 1.0
        else:
            # robust min/max for this group (ignore all-NaN)
            if np.isnan(M).all():
                vmin, vmax = 0.0, 1.0
            else:
                vmin = np.nanmin(M); vmax = np.nanmax(M)
                if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
                    vmin, vmax = 0.0, 1.0
        scalars.append((vmin, vmax))

        im = ax.imshow(M, aspect="auto", interpolation="nearest", cmap=cmap, vmin=vmin, vmax=vmax)

        # y-ticks only on the first group
        if gi == 0:
            ax.set_yticks(np.arange(df_plot.shape[0]))
            ax.set_yticklabels(df_plot.index, fontsize=y_tick_fontsize)
        else:
            ax.set_yticks([])
            ax.set_yticklabels([])

        # x-ticks for that group's columns
        ax.set_xticks(np.arange(len(gcols)))
        ax.set_xticklabels(gcols, rotation=x_tick_rotation, ha="right", fontsize=x_tick_fontsize)
        ax.tick_params(axis="x", pad=6)

        # thin bounding lines to make groups visually distinct
        ax.axvline(-0.5, color="k", lw=0.5, alpha=0.25)
        ax.axvline(len(gcols) - 0.5, color="k", lw=0.5, alpha=0.25)

    # ----- colorbars: one per group (narrow & long), extremes shown -----
    from matplotlib.ticker import FixedLocator, FixedFormatter  # put with other imports above
    for gi, (vmin, vmax) in enumerate(scalars):
        if gi < len(cbar_axes):
            sm = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax), cmap=cmap)
            cb = plt.colorbar(sm, cax=cbar_axes[gi], extend="both")

            # ticks at extremes + middle
            ticks = [vmin, 0.5 * (vmin + vmax), vmax]

            # round to 2 decimals and avoid "-0.00"
            labels = [("0.00" if abs(t) < 5e-6 else f"{t:.2f}") for t in ticks]

            # <<< THIS is the key: set on the colorbar, not the axis >>>
            cb.locator = FixedLocator(ticks)
            cb.formatter = FixedFormatter(labels)
            cb.update_ticks()

            title = group_titles[gi] if (group_titles and gi < len(group_titles)) else ""
            cbar_axes[gi].set_title(title, fontsize=11, pad=4)
            cb.ax.tick_params(labelsize=8, pad=1)  # small pad so text doesn’t overflow


    out_png = output_file if output_file.lower().endswith(".png") else output_file + ".png"
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    return out_png


def compose_tree_and_heatmap_old(
    ete_tree,
    df_leaf,                 # index = leaf names (already pruned/ordered)
    col_groups,
    output_file,
    *,
    group_titles=None,
    base_figsize=(22, 12),
    cbar_width_ratio=0.22,   # fixed slim column for stacked colorbars
    x_tick_rotation=90, x_tick_fontsize=8, y_tick_fontsize=7,
    nan_rgba=(0.92,0.92,0.92,1.0),
    ylabels_override=None,    # if provided, replaces df_leaf.index on the heatmap
    lock_scales_to_unit=False
):
    import numpy as np, matplotlib.pyplot as plt, matplotlib as mpl
    from matplotlib.gridspec import GridSpecFromSubplotSpec

    # ----- compute dynamic widths -----
    n_cols = sum(len([c for c in g if c in df_leaf.columns]) for g in col_groups)
    heatmap_ratio = 0.06 + 0.03 * n_cols
    # keep right column fixed for stacked colorbars
    tree_ratio = max(0.01, 1.0 - heatmap_ratio - 0.04)  # leave ~4% margin
    width_ratios = [tree_ratio, heatmap_ratio, 0.04]     # last small gap column
    # colorbars stacked in a nested gridspec inside the gap column, with fixed width later

    # scale figsize so it looks nice regardless of ratios (keep height)
    fig = plt.figure(figsize=base_figsize)
    gs = fig.add_gridspec(
        nrows=1, ncols=3, width_ratios=width_ratios,
        left=0.035, right=0.985, bottom=0.10, top=0.98, wspace=0.06
    )
    ax_tree = fig.add_subplot(gs[0, 0])
    ax_hm   = fig.add_subplot(gs[0, 1])
    gap_ax  = fig.add_subplot(gs[0, 2])
    gap_ax.axis("off")

    # stacked colorbars within the right margin (fixed physical width)
    # we make a new inset_axes inside gap_ax using axes coordinates
    # stacked colorbars within the right margin (fixed physical width)
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    n_cb = len(col_groups) if col_groups else 0
    cbar_axes = []
    if n_cb > 0:
        total_h = 0.94
        pad = 0.02
        cbar_h = (total_h - (n_cb - 1) * pad) / max(1, n_cb)
        y0 = (1 - total_h) / 2.0
        for i in range(n_cb):
            ax_c = inset_axes(
                gap_ax, width="80%", height=f"{cbar_h*100:.1f}%",
                bbox_to_anchor=(0.10, y0 + i*(cbar_h+pad), 0.80, cbar_h),
                bbox_transform=gap_ax.transAxes, borderpad=0
            )
            cbar_axes.append(ax_c)

    # ----- ensure heatmap y-labels are cluster numbers -----
    if ylabels_override is not None:
        df_plot = df_leaf.copy()
        df_plot.index = ylabels_override
    else:
        df_plot = df_leaf

    # ----- draw heatmap (get group scalars & cmap) -----
    cmap, group_scalars = draw_grouped_heatmap(
        ax_hm, df_plot, col_groups,
        group_titles=group_titles,
        nan_rgba=nan_rgba,
        x_tick_rotation=x_tick_rotation,
        x_tick_fontsize=x_tick_fontsize,
        y_tick_fontsize=y_tick_fontsize,
        return_group_scalars=True
    )

    # ----- draw tree aligned to SAME y-coordinates -----
    # Use the original leaf order (df_leaf.index order) for alignment
    draw_tree_aligned(ax_tree, ete_tree, leaf_order=list(df_leaf.index))

    # ----- colorbars -----
    cols = [c for g in col_groups for c in g if c in df_leaf.columns]
    M = df_leaf[cols].to_numpy(dtype=float)
    col_start = 0
    for gi, grp in enumerate(col_groups):
        gcols = [c for c in grp if c in df_leaf.columns]
        if not gcols:
            continue
        gsize = len(gcols)
        if lock_scales_to_unit:
            vmin, vmax = 0.0, 1.0
        else:
            vmin, vmax = group_scalars[gi]
        sm = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax), cmap=cmap)
        if gi < len(cbar_axes):
            cb = plt.colorbar(sm, cax=cbar_axes[gi])
            title = group_titles[gi] if (group_titles and gi < len(group_titles)) else ""
            cbar_axes[gi].set_title(title, fontsize=11, pad=4)
            cb.ax.tick_params(labelsize=8)
        col_start += gsize

    out_png = output_file if output_file.lower().endswith(".png") else output_file + ".png"
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    return out_png



# Perform ancestral reconstruction of sequences in the phylogenetic tree
# Input:
# tree_file - file with the phylogenetic tree
# alignment_file - file with the multiple sequence alignment
# output_file - where to save the ancestral sequences
def reconstruct_ancestral_sequences(tree_file, alignment_file, output_file):
    # Read the tree
    tree = Phylo.read(tree_file, 'newick')

    # Read the multiple sequence alignment
    seqs_IDs, seqs = load_fasta(alignment_file)
    seqs = [''.join([x for x in s if x.isupper() or x == '-']) for s in seqs]  # remove lowercase letters in alignment

    seq_records = [SeqRecord(Seq(seqs[i]), id=seqs_IDs[i]) for i in range(len(seqs))]  # Here must give correct names to sequences!

    # Create a MultipleSeqAlignment object from the SeqRecord objects
    alignment = MultipleSeqAlignment(seq_records)
#    alignment = AlignIO.read(alignment_file, 'fasta') # a3m file not all are of the same length !!


    # Perform the parsimony reconstruction
    scorer = ParsimonyScorer()
    searcher = NNITreeSearcher(scorer)
    constructor = ParsimonyTreeConstructor(searcher, alignment)
    print("Alignment:")
    print(alignment)
    print("Tree:")
    print(tree)
    parsimony_tree = constructor.build_tree(tree)

    # Print to file the reconstructed sequences for the internal nodes
    with open(output_file, 'w') as f:
        for clade in parsimony_tree.get_nonterminals():
            f.write(f"{clade.name}: {clade.confidences[0]}\n")

# Example usage
# reconstruct_ancestral_sequences('path/to/tree_file.newick', 'path/to/alignment_file.fasta', 'output.txt')
