from config import *
import re, os

# sys.path.append(os.path.join(os.path.dirname(__file__), '.'))

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
import numpy as np


# === CSV loaders & builders ===

def _pair_analysis_dir(pair_id: str) -> str:
    return os.path.join("Pipeline", pair_id, "Analysis")

def _load_pair_csvs(pair_id: str):
    """Load per-pair analysis CSVs; return (df_af, df_esm, df_cmap) or None if missing."""
    anal = _pair_analysis_dir(pair_id)
    df_af   = load_csv_or_none(os.path.join(anal, "df_af.csv"))
    df_esm  = load_csv_or_none(os.path.join(anal, "df_esm.csv"))
    df_cmap = load_csv_or_none(os.path.join(anal, "df_cmap.csv"))
    return df_af, df_esm, df_cmap

def _normalize_cluster_tag(s: str) -> str:
    """
    Normalize cluster tag into canonical labels used across the code:
    - 'DeepMsa', 'MSA_deep' -> 'DeepMsa'
    - Shallow numbers -> 'ShallowMsa_###' (zero-padded)
    """
    if not isinstance(s, str):
        s = str(s)
    t = s.strip()
    if "deep" in t.lower():
        return "DeepMsa"
    m = re.search(r"(\d+)$", t)
    if m:
        return f"ShallowMsa_{int(m.group(1)):03d}"
    # already canonical? keep
    if t.startswith("ShallowMsa_"):
        return t
    return t

def _build_tmscores_from_pair_csvs(df_af, df_esm, cluster_index: list[str]) -> pd.DataFrame:
    """
    Build columns (subset of): TM-AF1, TM-AF2, TM-ESM1, TM-ESM2
    Assumes the CSVs have: 'cluster_num', 'TMscore_fold1', 'TMscore_fold2'
    """
    import pandas as pd
    out = pd.DataFrame(index=cluster_index)

    def _pick(df, prefix):
        if df is None or df.empty:
            return
        d = df.copy()
        # accept either 'cluster_num' or 'cluster'
        if "cluster_num" not in d.columns and "cluster" in d.columns:
            d = d.rename(columns={"cluster": "cluster_num"})
        d["_tag"] = d["cluster_num"].astype(str).map(_normalize_cluster_tag)
        t1 = pd.to_numeric(d["TMscore_fold1"], errors="coerce")
        t2 = pd.to_numeric(d["TMscore_fold2"], errors="coerce")
        agg = d.assign(TM1=t1, TM2=t2).groupby("_tag")[["TM1","TM2"]].max()
        out[f"TM-{prefix}1"] = agg.reindex(out.index).TM1
        out[f"TM-{prefix}2"] = agg.reindex(out.index).TM2

    _pick(df_af,  "AF")
    _pick(df_esm, "ESM")
    keep = [c for c in ["TM-AF1","TM-AF2","TM-ESM1","TM-ESM2"] if c in out.columns]
    return out[keep]


def _build_msat_metrics_from_df_cmap(df_cmap, cluster_index: list[str]) -> pd.DataFrame:
    """
    Build the MSAT block from df_cmap.csv.

    We prefer MCC columns (common_mcc, uniq1_mcc, uniq2_mcc) and fall back to
    F1/recall/jaccard when MCC is missing. Rows are indexed by canonical cluster tags.
    """
    import pandas as pd
    out = pd.DataFrame(index=cluster_index)
    if df_cmap is None or df_cmap.empty:
        return out

    d = df_cmap.copy()

    # unify cluster id column to canonical tags
    if "cluster_num" in d.columns:
        d["_tag"] = d["cluster_num"].astype(str).map(_normalize_cluster_tag)
    elif "cluster" in d.columns:
        d["_tag"] = d["cluster"].astype(str).map(_normalize_cluster_tag)
    else:
        return out

    def pick(colnames):
        """Return the first present column as numeric series, else NaNs."""
        for c in colnames:
            if c in d.columns:
                return pd.to_numeric(d[c], errors="coerce")
        return pd.Series(index=d.index, dtype=float)

    cols = {
        "RE-MSAT-COM": pick(["common_mcc", "common_f1", "common_recall", "common_jaccard"]),
        "RE-MSAT1":    pick(["uniq1_mcc", "t1_mcc", "t1_f1", "t1_recall", "uniq1_recall", "t1_jaccard"]),
        "RE-MSAT2":    pick(["uniq2_mcc", "t2_mcc", "t2_f1", "t2_recall", "uniq2_recall", "t2_jaccard"]),
    }
    dd = pd.DataFrame(cols)
    dd["_tag"] = d["_tag"]

    # group by canonical cluster tag and take max (each cluster appears once anyway)
    agg = dd.groupby("_tag").max()

    # keep only the columns that actually exist and reindex to heatmap row order
    keep = [c for c in ["RE-MSAT-COM","RE-MSAT1","RE-MSAT2"] if c in agg.columns]
    return agg.reindex(cluster_index)[keep]


def _cluster_short_label_from_leaf(leaf_name, ete_leaves_cluster_ids):
    """Return a compact cluster label (e.g., '3' for ShallowMsa_003, 'D' for Deep)."""
    raw = ete_leaves_cluster_ids.get(leaf_name)
    if not raw or raw == 'p':
        return ""  # unknown
    tag = _normalize_cluster_tag(raw)
    if tag.lower().startswith("deep"):
        return "D"
    m = re.search(r"(\d+)$", tag)
    if m:
        return str(int(m.group(1)))  # no leading zeros
    return tag


def _cluster_tag_from_leaf(leaf_name, ete_leaves_cluster_ids):
    raw = ete_leaves_cluster_ids.get(leaf_name)
    if not raw or raw == "p":
        return None
    return _normalize_cluster_tag(raw)

def _group_leaves_by_cluster(ete_tree, ete_leaves_cluster_ids):
    """Return {cluster_tag: [leaf1, leaf2, ...]} using only leaves present in the tree."""
    groups = {}
    for leaf in ete_tree.iter_leaves():
        tag = _cluster_tag_from_leaf(leaf.name, ete_leaves_cluster_ids)
        if not tag:
            continue
        groups.setdefault(tag, []).append(leaf.name)
    return groups

def _choose_representative_leaves(cluster_to_leaves):
    """
    Choose exactly ONE leaf per cluster. For now, just pick the first;
    swap to a heuristic later (e.g., shortest branch, etc.) if you like.
    """
    reps = {tag: leaves[0] for tag, leaves in cluster_to_leaves.items() if leaves}
    return reps  # {cluster_tag: rep_leaf_name}

def _cluster_short_label_from_tag(tag):
    if not tag:
        return ""
    if tag.lower().startswith("deep"):
        return "D"
    m = re.search(r"(\d+)$", tag)
    return str(int(m.group(1))) if m else tag



def _collect_cluster_index_from_tree(ete_tree, ete_leaves_cluster_ids, cluster_node_values_keys) -> list[str]:
    """
    Build a *canonical* cluster index to index rows (ShallowMsa_### or DeepMsa)
    by scanning:
      - leaves present in the tree (preferred order),
      - and any additional keys we computed metrics for (to not drop data).
    """
    order = []
    seen = set()

    # from leaves
    for n in ete_tree.iter_leaves():
        cid = ete_leaves_cluster_ids.get(n.name)
        if not cid or cid == "p":
            continue
        tag = _normalize_cluster_tag(cid)
        if tag not in seen:
            order.append(tag); seen.add(tag)

    # include any extra computed keys (e.g. when leaf names are partial)
    for k in cluster_node_values_keys:
        tag = _normalize_cluster_tag(k)
        if tag not in seen:
            order.append(tag); seen.add(tag)
    return order


def _cluster_metrics_to_leaf_df(
    df_cluster,
    ete_tree,
    ete_leaves_cluster_ids,
    *,
    fillna_with=None,          # keep None to preserve NaNs; set 0.0 to fill
    drop_all_nan_rows=True     # drop leaves whose row is all-NaN
):
    """
    Map a cluster-indexed table -> leaf-indexed table (tree order).
    Optionally fill NaNs and/or drop leaves whose entire row is NaN.
    """
    import pandas as pd
    import numpy as np

    leaf_names = [n.name for n in ete_tree.iter_leaves()]
    rows = []
    keep = []
    for leaf in leaf_names:
        raw_cid = ete_leaves_cluster_ids.get(leaf, None)
        tag = _normalize_cluster_tag(raw_cid) if raw_cid and raw_cid != "p" else None
        if tag and tag in df_cluster.index:
            s = df_cluster.loc[tag]
        else:
            s = pd.Series(index=df_cluster.columns, dtype=float)  # all NaN
        rows.append(s)
        keep.append(leaf)

    df_leaf = pd.DataFrame(rows, index=keep)
    df_leaf = df_leaf[df_cluster.columns]  # preserve column order

    if drop_all_nan_rows:
        mask = ~df_leaf.isna().all(axis=1)
        df_leaf = df_leaf.loc[mask]

    if fillna_with is not None:
        df_leaf = df_leaf.fillna(fillna_with)

    return df_leaf


def _prepare_tree_heatmap_inputs(
    foldpair_id: str,
    ete_tree,
    ete_leaves_cluster_ids: dict[str,str],
    cluster_node_values: dict[str, tuple] | None = None
):
    """
    Build the DataFrame and grouping for visualize_tree_with_heatmap() *robustly*.
    We try CSV first; if missing, we fall back to in-memory computed metrics
    (cluster_node_values with shape (shared, foldA, foldB)).
    Returns: (df_for_heatmap, col_groups: list[list[str]], group_titles: list[str])
    """
    df_af, df_esm, df_cmap = _load_pair_csvs(foldpair_id)

    # row index (canonical cluster tags)
    cluster_index = _collect_cluster_index_from_tree(
        ete_tree, ete_leaves_cluster_ids,
        list(cluster_node_values.keys()) if cluster_node_values else []
    )

    # 1) TM blocks from AF/ESM CSVs
    tm_df = _build_tmscores_from_pair_csvs(df_af, df_esm, cluster_index)

    # 2) MSAT block from df_cmap
    ms_df = _build_msat_metrics_from_df_cmap(df_cmap, cluster_index)

    # 3) (optional) add RE-MSAT from computed contact overlap if CSV missing
    if ms_df.empty and cluster_node_values:
        # cluster_node_values: key -> (shared, foldA, foldB)
        tmp = pd.DataFrame({ _normalize_cluster_tag(k): v for k, v in cluster_node_values.items() }).T
        tmp.columns = ["RE-MSAT-COM", "RE-MSAT1", "RE-MSAT2"]
        ms_df = tmp.reindex(cluster_index)

    # Merge whatever exists
    blocks = []
    col_groups: list[list[str]] = []
    group_titles: list[str] = []

    if not tm_df.empty:
        blocks.append(tm_df)
        # derive which TM columns exist to set labels
        tm_cols = [c for c in ["TM-AF1","TM-AF2"] if c in tm_df.columns]
        if tm_cols:
            col_groups.append(tm_cols)
            group_titles.append("AF")
        tm_cols2 = [c for c in ["TM-ESM1","TM-ESM2"] if c in tm_df.columns and c.startswith("TM-ESM")]
        if tm_cols2:
            col_groups.append(tm_cols2)
            group_titles.append("ESM")

    if not ms_df.empty:
        blocks.append(ms_df)
        col_groups.append([c for c in ["RE-MSAT-COM","RE-MSAT1","RE-MSAT2"] if c in ms_df.columns])
        group_titles.append("MSAT")

    if not blocks:
        # return an empty-but-valid frame so the caller can decide to skip
        return pd.DataFrame(index=cluster_index), [], []

    df = pd.concat(blocks, axis=1)
    # keep only finite values; leave NaN where missing (visualize fills with nan_color)
    df = df.apply(pd.to_numeric, errors="coerce")
    return df, col_groups, group_titles


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


def _plot_contacts_panel(match_predicted_cmaps, match_true_cmap, fig_dir_root, foldpair_id):
    """
    Save three images into fig_dir_root:
      1) {pair}_all_clusters_cmap.png   (grid of all clusters)
      2) {pair}_best_clusters_cmap.png  (focused, picked by recall)  [done inside plot_array_contacts_and_predictions]
      3) {pair}_deep_clusters_cmap.png  (focused deep-MSA panel, if a deep key exists)
    """
    print("Plot Array Contact Map")

    # --- 1+2) 'all' grid + 'best' focused panel (this call saves both) ---
    save_root = os.path.join(fig_dir_root, f"{foldpair_id}_all_clusters_cmap")
    plot_array_contacts_and_predictions(
        match_predicted_cmaps,
        match_true_cmap,
        save_file=save_root,          # saves ..._all_clusters_cmap.png
        foldpair_id=foldpair_id       # also saves ..._best_clusters_cmap.png inside
    )

    # --- 3) Deep-only focused panel (if we have a deep-/DeepMsa key) ---
    deep_key = None
    # try a simple substring match first
    for k in match_predicted_cmaps.keys():
        if "deep" in str(k).lower():
            deep_key = k
            break
    # or normalize into canonical 'DeepMsa'
    if deep_key is None:
        try:
            deep_key = next(k for k in match_predicted_cmaps.keys()
                            if _normalize_cluster_tag(str(k)) == "DeepMsa")
        except StopIteration:
            deep_key = None

    if deep_key is None:
        print("[plot] no 'deep' cluster found -> skipping deep panel")
        return

    # Render deep panel for the two folds (matched layout to 'best')
    import matplotlib.pyplot as plt
    plt.figure(figsize=(10, 8))
    plot_foldswitch_contacts_and_predictions(
        predictions=match_predicted_cmaps[deep_key],  # tuple/list of [pred_f1, pred_f2]
        contacts=match_true_cmap,
        title="Deep-MSA",
        show_legend=True
    )
    out = os.path.join(fig_dir_root, f"{foldpair_id}_deep_clusters_cmap.png")
    plt.savefig(out, dpi=200)
    plt.close()
    print("[plot] wrote", out)


def _render_true_structures(pair_dir, pdb1, pdb2, out_dir, out_prefix):
    """Save (1) two_structures.png (unaligned) and (2) two_structures_aligned.png."""
    if not PYMOL_AVAILABLE:
        return
    import pymol
    from pymol import cmd
    pymol.finish_launching(['pymol', '-cq'])
    try:
        cmd.delete('all')
        cmd.load(os.path.join(pair_dir, pdb1), 'fold1')
        cmd.load(os.path.join(pair_dir, pdb2), 'fold2')
        cmd.color('red',  'fold1')
        cmd.color('blue', 'fold2')
        cmd.zoom('all', buffer=10)
        cmd.png(os.path.join(out_dir, f"{out_prefix}_two_structures.png"))
        cmd.align('fold2', 'fold1')
        cmd.png(os.path.join(out_dir, f"{out_prefix}_two_structures_aligned.png"))
    finally:
        cmd.quit()


def _best_cluster_from_df(pair_id: str, model_ver: str, fold_idx: int) -> str | None:
    """
    Return canonical cluster tag ('DeepMsa' or 'ShallowMsa_###') for the best TMscore
    in Pipeline/<pair>/Analysis/df_af.csv for AF{model_ver} and fold_idx in {0,1}.
    """
    import pandas as pd, numpy as np, re
    df_path = os.path.join('Pipeline', pair_id, 'Analysis', 'df_af.csv')
    if not os.path.isfile(df_path):
        return None
    d = pd.read_csv(df_path)

    # If there's a 'model' column, keep rows for AF2/AF3 accordingly
    if 'model' in d.columns:
        want = f"AF{model_ver}"
        d = d[d['model'].astype(str).str.upper() == want.upper()]
        if d.empty:
            return None

    # choose score column for fold1/fold2
    score_col = 'TMscore_fold1' if fold_idx == 0 else 'TMscore_fold2'
    if score_col not in d.columns:
        return None

    # choose cluster col
    ccol = 'cluster_num' if 'cluster_num' in d.columns else ('cluster' if 'cluster' in d.columns else None)
    if ccol is None:
        return None

    row = d.loc[d[score_col].astype(float).fillna(-np.inf).idxmax()]
    tag = str(row[ccol]).strip()
    # normalize canonical tag
    if tag.lower().startswith('deep'):
        return 'DeepMsa'
    m = re.fullmatch(r'\d+', tag)
    if m:
        return f"ShallowMsa_{int(m.group(0)):03d}"
    return tag


def _render_true_vs_best_models(pair_id: str, pdbids: list[str], pdbchains: list[str], fig_dir_root: str, model_ver='2'):
    """For each fold, overlay true structure vs best AF{ver} prediction."""
    if not PYMOL_AVAILABLE:
        return
    for idx in (0, 1):
        cluster = _best_cluster_from_df(pair_id, model_ver, idx)
        if not cluster:
            print(f"[3D] no best cluster in df_af.csv for AF{model_ver}, fold{idx+1}")
            continue
        true_pdb = os.path.join('Pipeline', pair_id, f"{pdbids[idx]}.pdb")
        pred_pdb = os.path.join('Pipeline', pair_id,
                                f"output_AF/AF{model_ver}/{cluster}__{pdbids[idx]+pdbchains[idx]}.pdb")
        if not (os.path.isfile(true_pdb) and os.path.isfile(pred_pdb)):
            print(f"[3D] missing PDB for AF{model_ver}: {true_pdb} / {pred_pdb}")
            continue
        out = os.path.join(fig_dir_root, f"{pair_id}_fold{idx+1}_vs_best_AF{model_ver}.png")
        try:
            align_and_visualize_proteins(true_pdb, pred_pdb, out, open_environment=True)
        except Exception as e:
            print(f"[3D] AF{model_ver} fold{idx+1} overlay failed: {e}")



# ----------- Plots Generation Functions -----------


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
    fasta_file_names = {  # Take just generated chain files
        pdbids[i] + pdbchains[i]:
            os.path.join(fasta_dir, foldpair_id, "fasta_chain_files", f"{pdbids[i]}{pdbchains[i]}.fasta")
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

    # ---------- Check length match between fasta seq and cmaps ----------
    for key in (pdbids[0] + pdbchains[0], pdbids[1] + pdbchains[1]):
        n_seq = len(seqs[key])
        n_map = true_cmap[key].shape[0]
        if n_seq != n_map:
            raise ValueError(
                f"[{foldpair_id}] FASTA length {n_seq} != contact-map size {n_map} for {key}. "
                "Regenerate the FASTA using CA-only residues (see get_fasta_chain_seq).")

    # ---------- Align truth/pred indices ----------
    pairwise_alignment = Align.PairwiseAligner().align(
        seqs[pdbids[0] + pdbchains[0]],
        seqs[pdbids[1] + pdbchains[1]],
    )
    print("Get matching indices: pdbids", pdbids, "pdbchains", pdbchains)
    match_true_cmap, match_predicted_cmaps = get_matching_indices_two_cmaps(
        pairwise_alignment, true_cmap, msa_transformer_pred)

    # ---------- Plot CMAPs ----------
    if plot_contacts:
        _plot_contacts_panel(match_predicted_cmaps, match_true_cmap, fig_dir_root, foldpair_id)

    # ---------- Metrics on shared/unique contacts ----------
    shared_unique_contacts, shared_unique_contacts_metrics, contacts_united = \
        match_predicted_and_true_contact_maps(match_predicted_cmaps, match_true_cmap)

    # Keep this dict in case df_cmap.csv is missing; values = (shared, fold1, fold2)
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

    # === Build inputs for the heatmap robustly (CSV-first; fallback to computed) ===
    heat_df, col_groups, group_titles = _prepare_tree_heatmap_inputs(
        foldpair_id=foldpair_id,
        ete_tree=ete_tree,
        ete_leaves_cluster_ids=ete_leaves_cluster_ids,
        cluster_node_values=cluster_node_values)

    if heat_df is not None and not heat_df.empty and col_groups:
        # build leaf-indexed matrix (may contain many leaves per cluster)
        heat_df_leaf = _cluster_metrics_to_leaf_df(
            df_cluster=heat_df,
            ete_tree=ete_tree,
            ete_leaves_cluster_ids=ete_leaves_cluster_ids,
            fillna_with=None,
            drop_all_nan_rows=True
        )

        # === NEW: one row per CLUSTER ===
        # 1) group leaves by cluster, 2) pick one representative leaf per cluster,
        # 3) reduce the matrix to those representatives (reindex),
        # 4) rename row labels to the CLUSTER TAGS (not leaf names).
        cluster_to_leaves = _group_leaves_by_cluster(ete_tree, ete_leaves_cluster_ids)
        reps = _choose_representative_leaves(cluster_to_leaves)   # {tag: rep_leaf}
        # keep only tags that exist in the cluster-indexed data (avoid orphan clusters)
        valid_tags = [t for t in reps if t in heat_df.index]

        # prune the tree to the representative leaves, then order rows by tree order
        rep_leaves = [reps[t] for t in valid_tags]
        prune_tree_to_leaves(ete_tree, keep_leaves=rep_leaves)
        leaf_order = [n.name for n in ete_tree.iter_leaves()]

        # map leaf order -> their cluster tags
        leaf_to_tag = {v: k for k, v in reps.items()}
        ordered_tags = [leaf_to_tag[l] for l in leaf_order if l in leaf_to_tag]

        # build the final matrix: one row per cluster-tag in tree order
        df_cluster_ordered = heat_df.loc[ordered_tags]
        print("[debug] columns:", df_cluster_ordered.columns.tolist())
        print("[debug] NaNs per column:\n", df_cluster_ordered.isna().sum())

        # ---- NEW: remap index back to the representative leaf names (for the tree) ----
        df_by_leaf = df_cluster_ordered.copy()
        df_by_leaf.index = leaf_order

        # override y tick labels with short cluster numbers (e.g., '3', 'D')
        ylabels_override = [_cluster_short_label_from_tag(t) for t in ordered_tags]

        out_root = os.path.join(fig_dir_root, f"{foldpair_id}_phytree_cluster")
        print("Making tree heatmap plot with column groups...", col_groups, flush=True)
        compose_tree_and_heatmap(
            ete_tree=ete_tree,
            df_leaf=df_by_leaf,  # <-- use LEAF names here
            col_groups=col_groups,
            group_titles=group_titles,
            output_file=out_root,
            base_figsize=(22, 12),
            x_tick_rotation=90, x_tick_fontsize=9, y_tick_fontsize=9,
            nan_rgba=(0.92, 0.92, 0.92, 1.0),
            ylabels_override=ylabels_override,
        )

        print(f"[ok] saved tree heatmap -> {out_root}.png")
        concat_scores = df_cluster_ordered

    else:
        # Fallback: plot a simple tree with per-leaf metrics if group layout is unavailable
        out_root = os.path.join(fig_dir_root, f"{foldpair_id}_phytree")
        # Build a minimal leaf->values frame from the MSA leaf mapping (may be empty)
        entries = {}
        for n in ete_tree.iter_leaves():
            cid = ete_leaves_cluster_ids.get(n.name)
            if not cid or cid == 'p':
                continue
            # prefer exact key; otherwise try normalized
            tag = _normalize_cluster_tag(cid)
            key = tag if tag in cluster_node_values else (cid if cid in cluster_node_values else None)
            if key is None:
                continue
            entries[n.name] = cluster_node_values[key]
        if entries:
            ete_leaves_node_values = pd.DataFrame(entries).T
            ete_leaves_node_values.columns = ["shared", pdbids[0] + pdbchains[0], pdbids[1] + pdbchains[1]]
        else:
            ete_leaves_node_values = pd.DataFrame(columns=["shared", pdbids[0] + pdbchains[0], pdbids[1] + pdbchains[1]])

        # Here we pass a single flat list of columns for the minimal figure
        visualize_tree_with_heatmap(
            phylo_tree=ete_tree,
            node_values_matrix=ete_leaves_node_values,
            col_groups=[list(ete_leaves_node_values.columns)],
            output_file=out_root,
            group_titles=["Scores"]
        )
        concat_scores = ete_leaves_node_values


    # ---------- Pair summaries ----------
    # cmap distances: only if we actually have predicted cmaps
    if match_predicted_cmaps and len(match_predicted_cmaps) > 0:
        try:
            cmap_dists_vec = compute_cmap_distances(match_predicted_cmaps)
        except Exception as e:
            print(f"[warn] compute_cmap_distances failed: {e!r}")
            cmap_dists_vec = None
    else:
        print("[warn] No predicted cmaps after alignment; skipping cmap distance calc.")
        cmap_dists_vec = None

    # seq distances: protect against empty MSA dict just in case
    try:
        seqs_dists = compute_seq_distances(msa_clusters)
        seqs_dists_vec = float(np.mean(seqs_dists)) if len(seqs_dists) else None
    except Exception as e:
        print(f"[warn] compute_seq_distances failed: {e!r}")
        seqs_dists_vec = None

    num_seqs_msa_vec = len(seqs)


    # ---------- 3D PANELS (PyMOL headless) ----------
    # --------- Optional 3D-align plot (only on non-Linux in your original code)  --------------
    if PYMOL_AVAILABLE:
        try:
            _render_true_structures(
                pair_dir=os.path.join('Pipeline', foldpair_id),
                pdb1=f"{pdbids[0]}.pdb",
                pdb2=f"{pdbids[1]}.pdb",
                out_dir=fig_dir_root,
                out_prefix=foldpair_id
            )
        except Exception as e:
            print(f"[3D] true-true render failed: {e}")

        # AF2
        try:
            _render_true_vs_best_models(foldpair_id, pdbids, pdbchains, fig_dir_root, model_ver='2')
        except Exception as e:
            print(f"[3D] AF2 overlays failed: {e}")

        # AF3
        try:
            _render_true_vs_best_models(foldpair_id, pdbids, pdbchains, fig_dir_root, model_ver='3')
        except Exception as e:
            print(f"[3D] AF3 overlays failed: {e}")
    else:
        print("[3D] PyMOL not available; skipping 3D panels.")

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
