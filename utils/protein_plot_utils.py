"""Orchestrator for per-pair fold-switch plots.

Delegates to focused submodules:
  - utils.plot_3d       : PyMOL / py3Dmol structure rendering
  - utils.plot_contacts : contact-map grid panels and evaluation
  - utils.plot_global   : cross-pair scatter plots

Keeps:
  - Tree / CSV helper functions used by make_foldswitch_all_plots
  - make_foldswitch_all_plots itself (the main per-pair orchestrator)
  - Re-exports for backward compatibility
"""

import os
import re
import shutil

import numpy as np
import pandas as pd
from Bio import Align
from glob import glob

from config import FIGURE_RES_DIR, PER_PAIR_PUBLISH_SUBDIR
from utils.utils import load_csv_or_none, pick, find_max_keys
from utils.protein_utils import read_msa, compute_precisions, load_cmap_and_idx
from utils.align_utils import align_truth_and_preds_via_msa_first2
from utils.phytree_utils import (
    Phylo,
    convert_biopython_to_ete3,
    prune_tree_to_leaves,
    compose_tree_and_heatmap,
    visualize_tree_with_heatmap,
    add_figure_title,
)
from utils.energy_utils import read_energy_tuples, align_and_compare_residues

# ── Submodule imports (used by make_foldswitch_all_plots) ──────────────
from utils.plot_3d import (
    _render_true_structures,
    _render_true_vs_best_models_generic,
    _render_true_vs_best_esm,
    _render_ddg_aligned_image,
)
from utils.plot_contacts import (
    _plot_contacts_panel,
    match_predicted_and_true_contact_maps,
    compute_cmap_distances,
    compute_seq_distances,
    seqs_ids_to_cluster_ids,
)
from utils.plot_global import global_pairs_statistics_plots

# ── Re-exports for backward compatibility ──────────────────────────────
# (run_foldswitch_pipeline.py imports make_foldswitch_all_plots and
#  global_pairs_statistics_plots from here)
from utils.plot_3d import (                                       # noqa: F401
    align_and_visualize_proteins,
    write_py3dmol_overlay_html,
)
from utils.plot_contacts import (                                 # noqa: F401
    plot_contacts_and_predictions,
    plot_foldswitch_contacts_and_predictions,
    plot_array_contacts_and_predictions,
    plot_viz_cmap,
)


# ═══════════════════════════════════════════════════════════════════════
# Tree / CSV helper functions (used only inside this file)
# ═══════════════════════════════════════════════════════════════════════

def _publish_figs_to_docs(foldpair_id: str, fig_dir_root: str) -> None:
    """Copy all figures from output_figs/ to docs/HTML/figs/<pair_id>/ for GitHub Pages."""
    docs_dir = os.path.join("docs", PER_PAIR_PUBLISH_SUBDIR, foldpair_id)
    os.makedirs(docs_dir, exist_ok=True)
    for f in os.listdir(fig_dir_root):
        if f.lower().endswith((".png", ".html", ".jpg")):
            shutil.copy2(os.path.join(fig_dir_root, f), os.path.join(docs_dir, f))
    print(f"[publish] synced figures -> {docs_dir}")


def _add_phytree_title(png_path: str, foldpair_id: str) -> None:
    """Add pair ID as title to a phytree figure (uses PIL via add_figure_title)."""
    if not os.path.isfile(png_path):
        return
    try:
        fold1, fold2 = foldpair_id.split("_", 1)
        title = f"{fold1} / {fold2}"
        add_figure_title(png_path, title, font_size=28)
    except Exception as e:
        print(f"[warn] could not add title to {png_path}: {e}")

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
    if t.startswith("ShallowMsa_"):
        return t
    return t


def _cluster_short_disp(tag: str) -> str:
    """ShallowMsa_019 -> C#19 ; DeepMsa/MSA_deep -> D; else pass-through."""
    t = _normalize_cluster_tag(tag or "")
    if "deep" in t.lower():
        return "D"
    m = re.search(r"(\d+)$", t)
    return f"C#{int(m.group(1))}" if m else t


def _pair_analysis_dir(pair_id: str) -> str:
    return os.path.join("Pipeline", "FoldPairs", pair_id, "Analysis")


def _load_pair_csvs(pair_id: str):
    """Load per-pair analysis CSVs; return (df_af, df_esm, df_cmap) or None if missing."""
    anal = _pair_analysis_dir(pair_id)
    df_af   = load_csv_or_none(os.path.join(anal, "df_af.csv"))
    df_esm  = load_csv_or_none(os.path.join(anal, "df_esm.csv"))
    df_cmap = load_csv_or_none(os.path.join(anal, "df_cmap.csv"))
    return df_af, df_esm, df_cmap


def _build_tmscores_from_pair_csvs(df_af, df_esm, cluster_index: list[str]) -> pd.DataFrame:
    """Build columns: TM-AF1, TM-AF2, TM-ESM1, TM-ESM2."""
    out = pd.DataFrame(index=cluster_index)

    def _pick(df, prefix):
        if df is None or df.empty:
            return
        d = df.copy()
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
    """Build the MSAT block from df_cmap.csv."""
    out = pd.DataFrame(index=cluster_index)
    if df_cmap is None or df_cmap.empty:
        return out

    d = df_cmap.copy()

    if "cluster_num" in d.columns:
        d["_tag"] = d["cluster_num"].astype(str).map(_normalize_cluster_tag)
    elif "cluster" in d.columns:
        d["_tag"] = d["cluster"].astype(str).map(_normalize_cluster_tag)
    else:
        return out

    cand_cols = {
        "RE-MSAT-COM": ["common_mcc", "common_f1", "common_recall", "common_jaccard"],
        "RE-MSAT1":    ["uniq1_mcc", "t1_mcc", "t1_f1", "t1_recall", "uniq1_recall", "t1_jaccard"],
        "RE-MSAT2":    ["uniq2_mcc", "t2_mcc", "t2_f1", "t2_recall", "uniq2_recall", "t2_jaccard"],
    }
    if not any(any(c in d.columns for c in cols) for cols in cand_cols.values()):
        return out

    cols = {k: pick(d, v) for k, v in cand_cols.items()}
    dd = pd.DataFrame(cols)
    dd["_tag"] = d["_tag"]

    agg = dd.groupby("_tag").max()

    keep = [c for c in ["RE-MSAT-COM","RE-MSAT1","RE-MSAT2"] if c in agg.columns]
    out_df = agg.reindex(cluster_index)[keep]
    out_df = out_df.dropna(axis=1, how="all")
    return out_df


def _cluster_short_label_from_leaf(leaf_name, ete_leaves_cluster_ids):
    """Return a compact cluster label (e.g., '3' for ShallowMsa_003, 'D' for Deep)."""
    raw = ete_leaves_cluster_ids.get(leaf_name)
    if not raw or raw == 'p':
        return ""
    tag = _normalize_cluster_tag(raw)
    if tag.lower().startswith("deep"):
        return "D"
    m = re.search(r"(\d+)$", tag)
    if m:
        return str(int(m.group(1)))
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
    """Choose exactly ONE leaf per cluster."""
    return {tag: leaves[0] for tag, leaves in cluster_to_leaves.items() if leaves}


def _cluster_short_label_from_tag(tag):
    if not tag:
        return ""
    if tag.lower().startswith("deep"):
        return "D"
    m = re.search(r"(\d+)$", tag)
    return f"S{int(m.group(1)):02d}" if m else tag


def _collect_cluster_index_from_tree(ete_tree, ete_leaves_cluster_ids, cluster_node_values_keys) -> list[str]:
    """Build a canonical cluster index by scanning tree leaves and computed keys."""
    order = []
    seen = set()

    for n in ete_tree.iter_leaves():
        cid = ete_leaves_cluster_ids.get(n.name)
        if not cid or cid == "p":
            continue
        tag = _normalize_cluster_tag(cid)
        if tag not in seen:
            order.append(tag); seen.add(tag)

    for k in cluster_node_values_keys:
        tag = _normalize_cluster_tag(k)
        if tag not in seen:
            order.append(tag); seen.add(tag)
    return order


def _cluster_metrics_to_leaf_df(df_cluster, ete_tree, ete_leaves_cluster_ids,
                                 *, fillna_with=None, drop_all_nan_rows=True):
    """Map a cluster-indexed table -> leaf-indexed table (tree order)."""
    leaf_names = [n.name for n in ete_tree.iter_leaves()]
    rows = []
    keep = []
    for leaf in leaf_names:
        raw_cid = ete_leaves_cluster_ids.get(leaf, None)
        tag = _normalize_cluster_tag(raw_cid) if raw_cid and raw_cid != "p" else None
        if tag and tag in df_cluster.index:
            s = df_cluster.loc[tag]
        else:
            s = pd.Series(index=df_cluster.columns, dtype=float)
        rows.append(s)
        keep.append(leaf)

    df_leaf = pd.DataFrame(rows, index=keep)
    df_leaf = df_leaf[df_cluster.columns]

    if drop_all_nan_rows:
        mask = ~df_leaf.isna().all(axis=1)
        df_leaf = df_leaf.loc[mask]

    if fillna_with is not None:
        df_leaf = df_leaf.fillna(fillna_with)

    return df_leaf


_UNIFIED_METHOD_ORDER = ["AF2", "AF3", "ESM", "MSAT", "CCMpred", "DDG", "Boltz2"]


def _consensus_fold_preference(df_centered: pd.DataFrame,
                                threshold: float = 0.10) -> dict:
    """Per-cluster F1/F2/Amb call from the cross-method consensus.

    Uses the unified per-method normalized centered preference matrix
    (rows = clusters, cols = methods) from `_build_unified_method_preference`.
    For each cluster takes the row mean across methods (NaN-skip) and
    thresholds:
        mean > +threshold -> F1   (cluster more F1-leaning than the pair's
                                   method-specific average across most methods)
        mean < -threshold -> F2
        otherwise         -> Amb

    Replaces the legacy AF2-only `assign_fold_preference` call which was
    redundant with the AF2 heatmap column.
    """
    if df_centered is None or df_centered.empty:
        return {}
    means = df_centered.mean(axis=1, skipna=True)
    out: dict = {}
    for tag, m in means.items():
        if not np.isfinite(m):
            out[tag] = "Amb"
        elif m > threshold:
            out[tag] = "F1"
        elif m < -threshold:
            out[tag] = "F2"
        else:
            out[tag] = "Amb"
    return out


def _load_af2_raw_tms_from_survey(foldpair_id: str) -> dict:
    """Return {cluster_tag: (TM1_max, TM2_max)} for the pair's AF2 rows.

    Used for the right-edge F1/F2/Amb fold-preference strip, which needs
    absolute TM-scores (not the normalized signed Δ used by the heatmap).
    """
    csv_path = os.path.join("docs", "fold_diversity_survey.csv")
    if not os.path.isfile(csv_path):
        return {}
    try:
        survey = pd.read_csv(csv_path)
    except Exception:
        return {}
    sub = survey[(survey["pair_id"] == foldpair_id) & (survey["method"] == "AF2")]
    out = {}
    for _, r in sub.iterrows():
        cl = str(r.get("cluster", ""))
        tag = "DeepMsa" if cl.lower().startswith("deep") else _normalize_cluster_tag(cl)
        try:
            tm1 = float(r["TM1_max"]); tm2 = float(r["TM2_max"])
        except (KeyError, ValueError, TypeError):
            continue
        if tag:
            out[tag] = (tm1, tm2)
    return out


def _build_unified_method_preference(foldpair_id: str, cluster_index: list[str]):
    """Per-method signed preference toward F1, one column per method.

    Reads docs/fold_diversity_survey.csv (aggregated per pair x cluster x method).
    Each cell = TMdiff_max for TM-based methods, t1_recall - t2_recall for contact
    methods, ddg_kcal_per_mol for DDG. All on a per-method-normalized scale so
    the diverging colormap is comparable across columns.

    Returns
    -------
    (df, deep_row) where:
      df         : DataFrame indexed by cluster tags (subset of cluster_index)
                   with one column per method; values in [-1, +1], NaN if missing.
      deep_row   : pd.Series of the same columns, holding the DeepMsa cluster
                   value per method (also normalized to [-1, +1]). May contain
                   NaN if any method has no DeepMsa row.
    """
    csv_path = os.path.join("docs", "fold_diversity_survey.csv")
    if not os.path.isfile(csv_path):
        return pd.DataFrame(index=cluster_index, columns=_UNIFIED_METHOD_ORDER), \
               pd.Series(index=_UNIFIED_METHOD_ORDER, dtype=float)
    try:
        survey = pd.read_csv(csv_path)
    except Exception:
        return pd.DataFrame(index=cluster_index, columns=_UNIFIED_METHOD_ORDER), \
               pd.Series(index=_UNIFIED_METHOD_ORDER, dtype=float)

    sub = survey[survey["pair_id"] == foldpair_id].copy()
    if sub.empty:
        return pd.DataFrame(index=cluster_index, columns=_UNIFIED_METHOD_ORDER), \
               pd.Series(index=_UNIFIED_METHOD_ORDER, dtype=float)

    sub["TMdiff_centered"] = pd.to_numeric(sub.get("TMdiff_centered"),
                                            errors="coerce")
    sub["TMdiff_max"] = pd.to_numeric(sub.get("TMdiff_max"), errors="coerce")
    sub["_tag"] = sub["cluster"].astype(str).map(
        lambda c: "DeepMsa" if str(c).lower().startswith("deep")
        else _normalize_cluster_tag(c)
    )

    # Per-cluster cells show TMdiff_centered (per-method signed deviation from
    # the pair's median TMdiff). This matches the per-pair-centered classification
    # used by the rest of the analysis (Sec. 2.5/3.4 of the paper), so a row's
    # color reflects the same signal the concordance / BH FDR pipeline uses.
    centered = sub.pivot_table(index="_tag", columns="method",
                               values="TMdiff_centered", aggfunc="first")
    # The Deep MSA "baseline" row uses ABSOLUTE TMdiff (TMdiff_max), not centered,
    # so the reader can see the un-clustered global bias on the same colorbar.
    absolute = sub.pivot_table(index="_tag", columns="method",
                               values="TMdiff_max", aggfunc="first")

    keep_cols = [m for m in _UNIFIED_METHOD_ORDER if m in centered.columns]
    centered = centered[keep_cols] if keep_cols else centered
    absolute = absolute.reindex(columns=centered.columns)

    # Per-method symmetric normalization to [-1, +1] using max-abs over the
    # whole pair, applied jointly to the centered cluster rows and the
    # absolute Deep row so they share a single colorbar.
    norm_c = centered.copy()
    norm_a = absolute.copy()
    for col in norm_c.columns:
        vc = pd.to_numeric(norm_c[col], errors="coerce")
        va = pd.to_numeric(norm_a[col], errors="coerce")
        m = float(max(vc.abs().max(skipna=True) or 0.0,
                      va.abs().max(skipna=True) or 0.0))
        if not np.isfinite(m) or m == 0.0:
            continue
        norm_c[col] = vc / m
        norm_a[col] = va / m

    deep_row = norm_a.loc["DeepMsa"] if "DeepMsa" in norm_a.index else \
               pd.Series(index=norm_c.columns, dtype=float)

    per_cluster = norm_c.drop(index="DeepMsa", errors="ignore")
    # Reindex to the tree-derived cluster order so missing clusters become NaN.
    per_cluster = per_cluster.reindex(cluster_index)
    return per_cluster, deep_row


def _prepare_tree_heatmap_inputs(foldpair_id, ete_tree, ete_leaves_cluster_ids,
                                  cluster_node_values=None):
    """Build a unified one-block heatmap input: one column per method, values are
    per-method-normalized signed preference toward F1 (positive=F1, negative=F2).
    Returns (df, col_groups, group_titles, deep_row).
    """
    cluster_index = _collect_cluster_index_from_tree(
        ete_tree, ete_leaves_cluster_ids,
        list(cluster_node_values.keys()) if cluster_node_values else []
    )

    df, deep_row = _build_unified_method_preference(foldpair_id, cluster_index)
    if df.empty or df.dropna(how="all").empty:
        return pd.DataFrame(index=cluster_index), [], [], pd.Series(dtype=float)

    df = df.apply(pd.to_numeric, errors="coerce")
    cols = list(df.columns)
    return df, [cols], ["Method preference (F2 ← → F1)"], deep_row


def _resolve_cluster_key(raw_key: str, cluster_node_values: dict) -> str | None:
    """Map ete leaves to cluster metrics; accept several cluster-key aliases."""
    cands = []
    k = str(raw_key)

    cands.append(k)

    m = re.fullmatch(r"\d+", k)
    if m:
        cands.append(f"ShallowMsa_{int(k):03d}")

    if "deep" in k.lower():
        cands.extend(["DeepMsa", "MSA_deep"])

    if k.lower().startswith("shallow") and not k.lower().startswith("shallowmsa_"):
        nn = re.sub(r"(?i)^shallow(msa)?_?", "ShallowMsa_", k)
        nn = re.sub(r"(\D)(\d{1,3})$", lambda m: f"{m.group(1)}{int(m.group(2)):03d}", nn)
        cands.append(nn)

    cands.append(k.replace("msa_t__", ""))

    for c in cands:
        if c in cluster_node_values:
            return c
    return None


# ═══════════════════════════════════════════════════════════════════════
# Main per-pair orchestrator
# ═══════════════════════════════════════════════════════════════════════

def make_foldswitch_all_plots(
    pdbids, fasta_dir, foldpair_id, pdbchains,
    plot_contacts: bool = True,
    global_plots: bool = False, plot3dformat: str = 'static'
):
    """
    Create per-pair plots (contact maps, tree heatmaps, etc.)
    Returns (cmap_dists_vec, seqs_dists_vec, num_seqs_msa_vec, concat_scores)
    """
    print("Plot for foldpair_id: " + foldpair_id)

    # ---------- Paths & ensure a single per-pair figure dir ----------
    fig_dir_root = os.path.join("Pipeline", "FoldPairs", foldpair_id, "output_figs")
    os.makedirs(fig_dir_root, exist_ok=True)

    # ---------- Inputs ----------
    fasta_file_names = {
        pdbids[i] + pdbchains[i]:
            os.path.join(fasta_dir, foldpair_id, "fasta_chain_files", f"{pdbids[i]}{pdbchains[i]}.fasta")
        for i in range(2)
    }

    msa_pred_files = glob(os.path.join(fasta_dir, foldpair_id, "output_cmaps", "msa_transformer", "*.npz"))
    msa_files = glob(os.path.join(fasta_dir, foldpair_id, "output_msa_cluster", "*.a3m"))
    msa_clusters = {os.path.basename(fp)[:-4]: read_msa(fp) for fp in msa_files}

    if len(msa_files) == 0:
        print("Shallow alignment! No MSA Clusters! Skipping family")
        return [None] * 4
    if pdbids[0] == pdbids[1]:
        print("Same PDB-ID for both folds! Might be buggy!")

    # ---------- Load FASTA sequences ----------
    seqs = {}
    for i in range(2):
        key = pdbids[i] + pdbchains[i]
        with open(fasta_file_names[key], "r") as fh:
            seqs[key] = fh.read().splitlines()[1].strip()

    if max(len(s) for s in seqs.values()) > 1024:
        print("Long sequence! Length > 1024 (not supported by this plotting path)")
        return [None] * 4

    # ---------- Load predicted CMAPs (shallow + deep) ----------
    print("msa_pred_files num=", len(msa_pred_files))
    msa_transformer_pred = {}
    for path in msa_pred_files:
        base = os.path.splitext(os.path.basename(path))[0]
        key = base[len("msa_t__"):] if base.startswith("msa_t__") else base
        try:
            pair_cmap, pair_idx, keep_cols = load_cmap_and_idx(path)
            msa_transformer_pred[key] = (pair_cmap, pair_idx)
        except Exception:
            print("[make_all_plots] Failed to load", path)
            msa_transformer_pred[key] = np.genfromtxt(path)
    print("msa_transformer_pred keys:", list(msa_transformer_pred.keys()), " num=", len(msa_transformer_pred))

    if not msa_transformer_pred:
        print(f"[plot:{foldpair_id}] No MSA-Transformer .npz files — skipping contact panels.")

    # ---------- Contact maps: load, align, plot, metrics ----------
    # This entire section is optional for tree/3D rendering.
    # If npy files are missing or alignment fails, we skip contacts gracefully.
    cluster_node_values = {}
    match_predicted_cmaps = {}
    if msa_transformer_pred:
        try:
            # Truth contacts
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

            # Check length match
            for key in (pdbids[0] + pdbchains[0], pdbids[1] + pdbchains[1]):
                n_seq = len(seqs[key])
                n_map = true_cmap[key].shape[0]
                if n_seq != n_map:
                    raise ValueError(
                        f"[{foldpair_id}] FASTA length {n_seq} != contact-map size {n_map} for {key}. "
                        "Regenerate the FASTA using CA-only residues (see get_fasta_chain_seq).")

            # Align truth/pred indices
            pairwise_alignment = Align.PairwiseAligner().align(
                seqs[pdbids[0] + pdbchains[0]],
                seqs[pdbids[1] + pdbchains[1]])
            print("Get matching indices: pdbids", pdbids, "pdbchains", pdbchains)
            print("true_cmap type and len, and len([0]):", type(true_cmap), len(true_cmap), len(true_cmap[pdbids[0] + pdbchains[0]]))
            first_key = next(iter(msa_transformer_pred.keys()), None)
            first_len = (msa_transformer_pred[first_key][0].shape[0] if first_key is not None else 0)
            print("msa_transformer_pred keys:", list(msa_transformer_pred.keys()), "first_len:", first_len)

            print(" true_cmaps sizes: ", [cmap.shape for cmap in true_cmap.values()])
            print(" pred_cmaps sizes: ", [cmap[0].shape for cmap in msa_transformer_pred.values()])
            msa_full_default = os.path.join("Pipeline", "FoldPairs", foldpair_id, "output_get_msa", "DeepMsa.a3m")
            msa_full_trim = os.path.join("Pipeline", "FoldPairs", foldpair_id, "output_get_msa", "DeepMsa_pairtrim.a3m")
            msa_full = msa_full_trim if os.path.isfile(msa_full_trim) else msa_full_default
            keyA = pdbids[0] + pdbchains[0]
            keyB = pdbids[1] + pdbchains[1]

            match_true_cmap, match_predicted_cmaps, pred_truth_slices = align_truth_and_preds_via_msa_first2(
                true_cmap=true_cmap,
                pred_cmaps=msa_transformer_pred,
                msa_path=msa_full,
                keyA=keyA,
                keyB=keyB,
                verbose=True,
            )

            print("After align match_true_cmaps sizes: ", [cmap.shape for cmap in match_true_cmap.values()])
            print("After align match_pred_cmaps sizes: ", [cmap.shape for cmap in match_predicted_cmaps.values()])

            # Plot CMAPs
            if plot_contacts:
                print("All match_true len and match_predicted_cmaps len:",  len(match_true_cmap), len(match_predicted_cmaps))
                try:
                    _plot_contacts_panel(match_predicted_cmaps, match_true_cmap, fig_dir_root, foldpair_id, pred_truth_slices=pred_truth_slices)
                except Exception as e:
                    print(f"[warn] _plot_contacts_panel failed: {e}")

            # Metrics on shared/unique contacts
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
        except Exception as e:
            print(f"[warn] contact map section failed for {foldpair_id}: {e}")
            print("       (tree heatmap will still use AF/ESM data from CSVs)")
            match_predicted_cmaps = {}

    # ---------- Tree & overlays ----------
    phytree_file = os.path.join('Pipeline', 'FoldPairs', foldpair_id, 'output_phytree', 'DeepMsa_tree.nwk')
    print("Load Biopython treefile: " + phytree_file)
    bio_tree = Phylo.read(phytree_file, "newick")
    print("Convert to ete3 tree:")
    ete_tree = convert_biopython_to_ete3(bio_tree)

    ete_leaves_cluster_ids = seqs_ids_to_cluster_ids(
        os.path.join(fasta_dir, foldpair_id, "output_msa_cluster", "*.a3m"),
        [n.name for n in ete_tree])
    print("Converted seq ids to cluster ids:")

    # === Build inputs for the heatmap ===
    heat_df, col_groups, group_titles, deep_row = _prepare_tree_heatmap_inputs(
        foldpair_id=foldpair_id,
        ete_tree=ete_tree,
        ete_leaves_cluster_ids=ete_leaves_cluster_ids,
        cluster_node_values=cluster_node_values)

    if heat_df is not None and not heat_df.empty and col_groups:
        heat_df_leaf = _cluster_metrics_to_leaf_df(
            df_cluster=heat_df,
            ete_tree=ete_tree,
            ete_leaves_cluster_ids=ete_leaves_cluster_ids,
            fillna_with=None,
            drop_all_nan_rows=True
        )

        cluster_to_leaves = _group_leaves_by_cluster(ete_tree, ete_leaves_cluster_ids)
        reps = _choose_representative_leaves(cluster_to_leaves)
        valid_tags = [t for t in reps if t in heat_df.index]

        rep_leaves = [reps[t] for t in valid_tags]
        prune_tree_to_leaves(ete_tree, keep_leaves=rep_leaves)
        leaf_order = [n.name for n in ete_tree.iter_leaves()]

        leaf_to_tag = {v: k for k, v in reps.items()}
        ordered_tags = [leaf_to_tag[l] for l in leaf_order if l in leaf_to_tag]

        df_cluster_ordered = heat_df.loc[ordered_tags]
        print("[debug] columns:", df_cluster_ordered.columns.tolist())
        print("[debug] NaNs per column:\n", df_cluster_ordered.isna().sum())

        df_by_leaf = df_cluster_ordered.copy()
        df_by_leaf.index = leaf_order

        ylabels_override = [_cluster_short_label_from_tag(t) for t in ordered_tags]

        # ----- Compute fold preferences for colored tree tips -----
        # Multi-method consensus call (mean of normalized centered preferences
        # across the 7 methods) rather than AF2-only — avoids visually repeating
        # the AF2 signal on both the tree leaves and the AF2 heatmap column.
        _FOLD_COLORS = {"F1": "#d62728", "F2": "#1f77b4", "Amb": "#999999"}
        fold_pref_per_row = None
        leaf_colors_dict = None
        prefs = _consensus_fold_preference(df_cluster_ordered, threshold=0.10)
        if prefs:
            fold_pref_per_row = [prefs.get(t, "Amb") for t in ordered_tags]
            leaf_colors_dict = {}
            for tag, leaf_name in zip(ordered_tags, leaf_order):
                pref = prefs.get(tag, "Amb")
                leaf_colors_dict[leaf_name] = _FOLD_COLORS[pref]

        out_root = os.path.join(fig_dir_root, f"{foldpair_id}_phytree_cluster")
        print("Making tree heatmap plot with column groups...", col_groups, flush=True)
        compose_tree_and_heatmap(
            ete_tree=ete_tree,
            df_leaf=df_by_leaf,
            col_groups=col_groups,
            group_titles=group_titles,
            output_file=out_root,
            base_figsize=(22, 12),
            x_tick_rotation=90, x_tick_fontsize=9, y_tick_fontsize=9,
            nan_rgba=(0.92, 0.92, 0.92, 1.0),
            ylabels_override=ylabels_override,
            leaf_colors=leaf_colors_dict,
            # Drop the redundant Fold strip: the leaf-tip rings on the tree
            # already encode the F1/F2/Amb call.
            fold_pref_per_row=None,
            unified_diverging=True,
            extra_top_row=deep_row if (deep_row is not None and not deep_row.empty) else None,
            extra_top_row_label="Deep MSA",
            label_in_leaf=True,
        )

        print(f"[ok] saved tree heatmap -> {out_root}.png")
        _add_phytree_title(out_root + ".png", foldpair_id)
        concat_scores = df_cluster_ordered

    else:
        out_root = os.path.join(fig_dir_root, f"{foldpair_id}_phytree_cluster")
        # Build per-cluster entries from contact map metrics
        entries = {}
        for n in ete_tree.iter_leaves():
            cid = ete_leaves_cluster_ids.get(n.name)
            if not cid or cid == 'p':
                continue
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

        # Prune tree to one representative leaf per cluster (same as primary path)
        cluster_to_leaves = _group_leaves_by_cluster(ete_tree, ete_leaves_cluster_ids)
        reps = _choose_representative_leaves(cluster_to_leaves)
        rep_leaves = list(reps.values())
        if rep_leaves:
            prune_tree_to_leaves(ete_tree, keep_leaves=rep_leaves)

        # Reindex to pruned leaves and set cluster labels
        leaf_order = [n.name for n in ete_tree.iter_leaves()]
        leaf_to_tag = {v: k for k, v in reps.items()}
        ylabels_fb = [_cluster_short_label_from_tag(leaf_to_tag.get(l, ""))
                      for l in leaf_order]
        pruned_vals = ete_leaves_node_values.reindex(leaf_order)

        visualize_tree_with_heatmap(
            phylo_tree=ete_tree,
            node_values_matrix=pruned_vals,
            col_groups=[list(ete_leaves_node_values.columns)],
            output_file=out_root,
            group_titles=["Scores"],
            ylabels_override=ylabels_fb,
        )
        _add_phytree_title(out_root + ".png", foldpair_id)
        concat_scores = ete_leaves_node_values


    # ---------- Pair summaries ----------
    if match_predicted_cmaps and len(match_predicted_cmaps) > 0:
        try:
            cmap_dists_vec = compute_cmap_distances(match_predicted_cmaps)
        except Exception as e:
            print(f"[warn] compute_cmap_distances failed: {e!r}")
            cmap_dists_vec = None
    else:
        print("[warn] No predicted cmaps after alignment; skipping cmap distance calc.")
        cmap_dists_vec = None

    try:
        seqs_dists = compute_seq_distances(msa_clusters)
        seqs_dists_vec = float(np.mean(seqs_dists)) if len(seqs_dists) else None
    except Exception as e:
        print(f"[warn] compute_seq_distances failed: {e!r}")
        seqs_dists_vec = None

    num_seqs_msa_vec = len(seqs)

    # ---- ΔΔG aligned profile image ----
    try:
        _render_ddg_aligned_image(foldpair_id, pdbids, pdbchains, fig_dir_root)
    except Exception as e:
        print(f"[ΔΔG] skipped for {foldpair_id}: {e}")


    # ---------- 3D PANELS ----------
    make_interactive = plot3dformat in ('interactive', 'both')
    try:
        _render_true_structures(
            pair_dir=os.path.join('Pipeline', 'FoldPairs', foldpair_id),
            pdb1=f"{pdbids[0]}.pdb",
            pdb2=f"{pdbids[1]}.pdb",
            out_dir=fig_dir_root,
            out_prefix=foldpair_id,
            make_interactive=make_interactive
        )
    except Exception as e:
        print(f"[3D] true-true render failed: {e}")

    # AF2 / AF3
    _render_true_vs_best_models_generic(foldpair_id, pdbids, pdbchains, fig_dir_root, family='AF', ver='2',
                                        make_interactive=make_interactive)
    _render_true_vs_best_models_generic(foldpair_id, pdbids, pdbchains, fig_dir_root, family='AF', ver='3',
                                        make_interactive=make_interactive)

    # ESM2
    try:
        _render_true_vs_best_esm(foldpair_id, pdbids, pdbchains, fig_dir_root, model_ver='2', make_interactive=make_interactive)
    except Exception as e:
        print(f"[3D] ESM2 overlays failed: {e}")

    # ESM3
    try:
        _render_true_vs_best_esm(foldpair_id, pdbids, pdbchains, fig_dir_root, model_ver='3', make_interactive=make_interactive)
    except Exception as e:
        print(f"[3D] ESM3 overlays failed: {e}")

    if global_plots:
        print("Make global plots!")
        global_pairs_statistics_plots(output_dir=FIGURE_RES_DIR)

    # Auto-publish figures to docs/ for GitHub Pages
    try:
        _publish_figs_to_docs(foldpair_id, fig_dir_root)
    except Exception as e:
        print(f"[publish] failed to sync figures to docs/: {e}")

    return cmap_dists_vec, seqs_dists_vec, num_seqs_msa_vec, concat_scores
