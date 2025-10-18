from PIL import Image, ImageDraw, ImageFont  # NEW

import pandas as pd
import numpy as np
import os, re, math, shutil  # make sure shutil is imported
import subprocess, shlex, tempfile
from Bio import Align
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt

# sys.path.append(os.path.join(os.path.dirname(__file__), '.'))

PYMOL_AVAILABLE = False
try:
    import pymol
    from pymol import cmd
    PYMOL_AVAILABLE = True
except Exception:
    # Keep plotting utils importable on machines without PyMOL
    PYMOL_AVAILABLE = False

from config import *
from utils.phytree_utils import *
from scripts.MSA_Clust import *
from utils.utils import *
from utils.energy_utils import *
from utils.align_utils import *
from utils.protein_utils import *

# at top, near the PYMOL_AVAILABLE detection
PYMOL_BIN = os.environ.get("PYMOL_BIN")  # path like /sci/.../pymol-venv/bin/pymol



# === CSV loaders & builders ===
def _detect_pymol_bin() -> str | None:
    # 1) explicit env
    cand = os.environ.get("PYMOL_BIN")
    if cand and os.path.isfile(cand):
        return cand
    # 2) known per-user install (adjust if you move it)
    known = PYMOL_EXE
    if os.path.isfile(known):
        return known
    # 3) in PATH
    for name in ("pymol", "pymol-open-source"):
        path = shutil.which(name)
        if path:
            return path
    return None

# cache at import time, but functions will re-check env so you can override per-run
_PYMOL_BIN_CACHED = _detect_pymol_bin()



def _pred_pdb_path(pair_id: str, family: str, ver: str, cluster: str, pdbid: str, chain: str) -> str:
    """family = 'AF' or 'ESM'; ver e.g. '2' or '3'."""
    tag = _normalize_cluster_tag(cluster)
    base = f"output_{family}/{family}{ver}/{tag}__{pdbid}{chain}.pdb"
    return os.path.join("Pipeline", pair_id, base)

def _best_cluster_from_table(pair_id: str, table: str, model_tag: str, fold_idx: int) -> tuple[str,float] | None:
    """
    Generic best-cluster picker from df_{table}.csv.
    Returns (canonical_cluster_tag, TMscore) for fold_idx {0,1}.
    """
    df_path = os.path.join('Pipeline', pair_id, 'Analysis', f'df_{table}.csv')
    if not os.path.isfile(df_path):
        return None
    try:
        d = pd.read_csv(df_path)
    except Exception:
        return None

    if 'model' in d.columns:
        d = d[d['model'].astype(str).str.upper() == model_tag.upper()]
        if d.empty:
            return None

    score_col = 'TMscore_fold1' if fold_idx == 0 else 'TMscore_fold2'
    ccol = 'cluster_num' if 'cluster_num' in d.columns else ('cluster' if 'cluster' in d.columns else None)
    if (score_col not in d.columns) or (ccol is None):
        return None

    d = d.copy()
    d['_score'] = pd.to_numeric(d[score_col], errors='coerce')
    d = d.dropna(subset=['_score'])
    if d.empty:
        return None

    row = d.loc[d['_score'].idxmax()]
    tag = _normalize_cluster_tag(str(row[ccol]).strip())
    return tag, float(row['_score'])


def _cluster_short_disp(tag: str) -> str:
    """ShallowMsa_019 -> C#19 ; DeepMsa/MSA_deep -> D; else pass-through."""
    t = _normalize_cluster_tag(tag or "")
    if "deep" in t.lower():
        return "D"
    m = re.search(r"(\d+)$", t)
    return f"C#{int(m.group(1))}" if m else t

def _add_png_legend(img_path: str, left_label: str, right_label: str, *,
                    title: str | None = None,
                    left_rgb=(255, 0, 0), right_rgb=(0, 80, 255)) -> None:
    """Append a white strip with a red/blue legend (and optional title) at the bottom of a PNG."""
    try:
        im = Image.open(img_path).convert("RGB")
    except Exception as e:
        print(f"[3D] legend skipped (open failed): {e}")
        return

    W, H = im.size
    extra = max(70, int(H * 0.12))
    out = Image.new("RGB", (W, H + extra), (255, 255, 255))
    out.paste(im, (0, 0))
    draw = ImageDraw.Draw(out)

    # font (fallback friendly)
    try:
        f_big = ImageFont.truetype("DejaVuSans-Bold.ttf", 20)
        f = ImageFont.truetype("DejaVuSans.ttf", 18)
    except Exception:
        f_big = f = ImageFont.load_default()

    y = H + (extra - 20) // 2
    pad = 18
    box = 16

    # left (red)
    x = pad
    draw.rectangle((x, y, x + box, y + box), fill=left_rgb)
    draw.text((x + box + 8, y - 2), f"Red: {left_label}", fill=(0, 0, 0), font=f)

    # right (blue)
    x2 = x + 360
    draw.rectangle((x2, y, x2 + box, y + box), fill=right_rgb)
    draw.text((x2 + box + 8, y - 2), f"Blue: {right_label}", fill=(0, 0, 0), font=f)

    # optional title (above legend, aligned left)
    if title:
        draw.text((pad, H - 26), title, fill=(0, 0, 0), font=f_big)

    try:
        out.save(img_path)
    except Exception as e:
        print(f"[3D] legend save failed: {e}")

def _pymol_cli_render(script_text: str, exe: str | None = None):
    exe = exe or os.environ.get("PYMOL_BIN") or _PYMOL_BIN_CACHED
    if not exe or not os.path.isfile(exe):
        raise RuntimeError("PyMOL CLI not found. Set PYMOL_BIN or install pymol in PATH.")
    with tempfile.NamedTemporaryFile("w", suffix=".pml", delete=False) as tf:
        tf.write(script_text)
        tf.flush()
        cmdline = f"{shlex.quote(exe)} -cq {shlex.quote(tf.name)}"
        subprocess.run(cmdline, shell=True, check=True)


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
      2) {pair}_best_clusters_cmap.png  (focused 1-panel with the best cluster per fold)
      3) {pair}_deep_cmap.png           (Deep-MSA focused panel, if Deep exists)
    """
    print("Plot Array Contact Map")
    os.makedirs(fig_dir_root, exist_ok=True)

    # --- 1+2) 'all' grid + 'best' focused panel ---
    save_root = os.path.join(fig_dir_root, f"{foldpair_id}_all_clusters_cmap")
    n_pred = len(match_predicted_cmaps or {})
    if n_pred == 0:
        print("[plot] WARN: no predicted cmaps found after alignment; skipping cmap panels.")
        return
    plot_array_contacts_and_predictions(
        match_predicted_cmaps,
        match_true_cmap,
        save_file=save_root,          # saves ..._all_clusters_cmap.png
        foldpair_id=foldpair_id       # also saves ..._best_clusters_cmap.png inside
    )

    # --- 3) Deep-only panel ---
    # Show what keys we actually got (handy for debugging)
    try:
        klist = [str(k) for k in (match_predicted_cmaps or {}).keys()]
    except Exception:
        pass

    deep_key = None
    # 1) exact common filenames
    if isinstance(match_predicted_cmaps, dict):
        if "MSA_deep" in match_predicted_cmaps:
            deep_key = "MSA_deep"

    # 2) any key that *mentions* 'deep'
    if deep_key is None:
        for k in match_predicted_cmaps.keys():
            if "deep" in str(k).lower():
                deep_key = k
                break

    # 3) normalize aliases to 'DeepMsa'
    if deep_key is None:
        for k in match_predicted_cmaps.keys():
            try:
                if _normalize_cluster_tag(str(k)) == "DeepMsa":
                    deep_key = k
                    break
            except Exception:
                continue

    if deep_key is None:
        print("[plot] no 'deep' cluster found -> skipping deep panel")
        return

    plt.figure(figsize=(10, 8))
    plot_foldswitch_contacts_and_predictions(
        predictions=match_predicted_cmaps[deep_key],  # accepts 1 or 2 maps
        contacts=match_true_cmap,
        title=f"Deep-MSA ({deep_key})",
        show_legend=True,
        cluster_names=("D", "D"),
    )
    out = os.path.join(fig_dir_root, f"{foldpair_id}_deep_cmap.png")
    plt.savefig(out, dpi=200)
    plt.close()
    print("[plot] wrote cmaps", out)


def _render_true_structures(pair_dir, pdb1, pdb2, out_dir, out_prefix, make_interactive: bool = False):
    """Save (1) two_structures.png and (2) two_structures_aligned.png, with legends."""
    out_two = os.path.join(out_dir, f"{out_prefix}_two_structures.png")
    out_aln = os.path.join(out_dir, f"{out_prefix}_two_structures_aligned.png")
    os.makedirs(out_dir, exist_ok=True)

    label_red  = os.path.splitext(pdb1)[0]  # e.g., 2qqj
    label_blue = os.path.splitext(pdb2)[0]  # e.g., 4qds

    if not PYMOL_AVAILABLE:
        exe = os.environ.get("PYMOL_BIN") or _PYMOL_BIN_CACHED
        if exe: #if PYMOL_BIN:
            p1 = os.path.join(pair_dir, pdb1)
            p2 = os.path.join(pair_dir, pdb2)
            script = f"""
reinitialize
load "{p1}", fold1
load "{p2}", fold2
color red, fold1
color blue, fold2
bg_color white
set orthoscopic, on
viewport 1200, 900
zoom all, 10
png {out_two}, dpi=220, width=1200, height=900
align fold2, fold1
png {out_aln}, dpi=220, width=1200, height=900
quit
""".lstrip()
            _pymol_cli_render(script)
            _add_png_legend(out_two,  label_red, label_blue, title=None)
            _add_png_legend(out_aln,  label_red, label_blue, title=None)
            return
        print("[PyMOL] skipped (neither import nor CLI available)")
        return

    pymol.finish_launching(['pymol', '-cq'])
    try:
        cmd.delete('all')
        cmd.load(os.path.join(pair_dir, pdb1), 'fold1')
        cmd.load(os.path.join(pair_dir, pdb2), 'fold2')
        cmd.color('red',  'fold1')
        cmd.color('blue', 'fold2')
        cmd.bg_color('white')
        cmd.set('orthoscopic', 1)
        cmd.zoom('all', buffer=10)
        cmd.png(out_two, dpi=220, width=1200, height=900)
        cmd.align('fold2', 'fold1')
        cmd.png(out_aln, dpi=220, width=1200, height=900)
    finally:
        cmd.quit()
    _add_png_legend(out_two,  label_red, label_blue, title=None)
    _add_png_legend(out_aln,  label_red, label_blue, title=None)

    # Optional interactive (3D)
    if make_interactive:
        # Unaligned interactive
        try:
            html_unaln = os.path.join(out_dir, f"{out_prefix}_two_structures.html")
            write_py3dmol_overlay_html(
                os.path.join(pair_dir, pdb1),
                os.path.join(pair_dir, pdb2),
                html_unaln,
                left_label=label_red,
                right_label=label_blue,
                title="True vs True (unaligned)"
            )
            print(f"[3D] wrote {html_unaln}")
        except Exception as e:
            print(f"[3D] true-true interactive (unaligned) failed: {e}")

        # Aligned interactive
        try:
            html_aln = os.path.join(out_dir, f"{out_prefix}_two_structures_aligned.html")
            write_py3dmol_overlay_html(
                os.path.join(pair_dir, pdb1),
                os.path.join(pair_dir, pdb2),
                html_aln,
                left_label=label_red,
                right_label=label_blue,
                title="True vs True (aligned)"
            )
            print(f"[3D] wrote {html_aln}")
        except Exception as e:
            print(f"[3D] true-true interactive (aligned) failed: {e}")




def _best_cluster_from_df(pair_id: str, model_ver: str, fold_idx: int) -> str | None:
    """
    Return canonical cluster tag ('DeepMsa' or 'ShallowMsa_###') for the best TMscore
    in Pipeline/<pair>/Analysis/df_af.csv for AF{model_ver} and fold_idx in {0,1}.
    """
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



def _best_esm_from_df(pair_id: str, model_ver: str, fold_idx: int):
    """
    Return (cluster_tag, sample_idx) for the best TMscore in df_esm.csv
    for ESM{model_ver} and fold_idx in {0,1}.
    """
    import pandas as pd, numpy as np, re
    df_path = os.path.join('Pipeline', pair_id, 'Analysis', 'df_esm.csv')
    if not os.path.isfile(df_path):
        return None
    d = pd.read_csv(df_path)

    # Keep only this ESM model if a 'model' column exists
    want = f"ESM{model_ver}".upper()
    if 'model' in d.columns:
        d = d[d['model'].astype(str).str.upper() == want]
        if d.empty:
            return None

    # Score column for fold1/fold2
    score_col = 'TMscore_fold1' if fold_idx == 0 else 'TMscore_fold2'
    if score_col not in d.columns or d[score_col].isna().all():
        return None

    row = d.loc[d[score_col].astype(float).fillna(-np.inf).idxmax()]

    # Cluster tag -> canonical form (e.g. 'ShallowMsa_019')
    if 'cluster_num' in d.columns and pd.notna(row['cluster_num']):
        tag_raw = str(row['cluster_num']).strip()
    elif 'cluster' in d.columns and pd.notna(row['cluster']):
        tag_raw = str(row['cluster']).strip()
    else:
        return None
    cluster_tag = _normalize_cluster_tag(tag_raw)

    # Sample index (expect one of these columns; default to 1 if absent)
    sample_col = next((c for c in ('sample', 'sample_id', 'sample_num') if c in d.columns), None)
    sample_idx = int(row[sample_col]) if sample_col and pd.notna(row[sample_col]) else 1

    return cluster_tag, sample_idx


def _render_true_vs_best_esm(pair_id: str, pdbids: list[str], pdbchains: list[str], fig_dir_root: str,
                             model_ver='2', make_interactive: bool=False):
    """
    For each fold, overlay true structure vs best ESM{ver} prediction.
    ESM PDBs live in: Pipeline/<pair>/output_esm_fold/esm{ver}/ShallowMsa_XXX__sample_YYY_esm{ver}.pdb
    """
    for idx in (0, 1):
        picked = _best_esm_from_df(pair_id, model_ver, idx)
        if not picked:
            print(f"[3D] no best row in df_esm.csv for ESM{model_ver}, fold{idx+1}")
            continue
        cluster_tag, sample_idx = picked

        true_pdb = os.path.join('Pipeline', pair_id, f"{pdbids[idx]}.pdb")
        pred_pdb = os.path.join(
            'Pipeline', pair_id, 'output_esm_fold', f"esm{model_ver}",
            f"{cluster_tag}__sample_{sample_idx:03d}_esm{model_ver}.pdb"
        )
        if not (os.path.isfile(true_pdb) and os.path.isfile(pred_pdb)):
            print(f"[3D] missing PDB for ESM{model_ver}: {true_pdb} / {pred_pdb}")
            continue

        out = os.path.join(fig_dir_root, f"{pair_id}_fold{idx+1}_vs_best_ESM{model_ver}.png")
        try:
            align_and_visualize_proteins(true_pdb, pred_pdb, out, open_environment=True)
            print(f"[3D] wrote {out}")
        except Exception as e:
            print(f"[3D] ESM{model_ver} fold{idx+1} overlay failed: {e}")

        if make_interactive:
            title = f"ESM{model_ver}({_cluster_short_disp(cluster_tag)}) vs. {pdbids[idx]}{pdbchains[idx]}"
            html_out = os.path.join(fig_dir_root, f"{pair_id}_fold{idx+1}_vs_best_ESM{model_ver}.html")
            try:
                write_py3dmol_overlay_html(
                    true_pdb, pred_pdb, html_out,
                    left_label=f"{pdbids[idx]}{pdbchains[idx]}",
                    right_label=f"ESM{model_ver}",
                    title=title, color_left="red", color_right="blue"
                )
                print(f"[3D] wrote {html_out}")
            except Exception as e:
                print(f"[3D] ESM{model_ver} interactive failed: {e}")


def _render_true_vs_best_models_generic(pair_id: str, pdbids: list[str], pdbchains: list[str],
                                        fig_dir_root: str, *, family: str, ver: str,  make_interactive: bool = False):
    """
    family='AF' or 'ESM'; ver='2' or '3'
    For each fold, overlay true structure vs best {family}{ver} prediction (CLI/API) and add legend/title.
    """
    exe = os.environ.get("PYMOL_BIN") or _PYMOL_BIN_CACHED
    if not PYMOL_AVAILABLE and not exe:
        print(f"[3D] Skipping {family}{ver}: no PyMOL")
        return

    os.makedirs(fig_dir_root, exist_ok=True)
    table = 'af' if family == 'AF' else 'esm'
    model_tag = f"{family}{ver}"

    for idx in (0, 1):
        pick = _best_cluster_from_table(pair_id, table=table, model_tag=model_tag, fold_idx=idx)
        if not pick:
            print(f"[3D] no best cluster in df_{table}.csv for {model_tag}, fold{idx+1}")
            continue

        cluster, tm = pick
        true_pdb = os.path.join('Pipeline', pair_id, f"{pdbids[idx]}.pdb")
        pred_pdb = _pred_pdb_path(pair_id, family, ver, cluster, pdbids[idx], pdbchains[idx])

        if not (os.path.isfile(true_pdb) and os.path.isfile(pred_pdb)):
            print(f"[3D] missing PDB for {model_tag}: {true_pdb} / {pred_pdb}")
            continue

        out = os.path.join(fig_dir_root, f"{pair_id}_fold{idx+1}_vs_best_{model_tag}.png")
        try:
            # render
            align_and_visualize_proteins(true_pdb, pred_pdb, out, open_environment=bool(PYMOL_AVAILABLE))
            # annotate
            left = f"{pdbids[idx]}{pdbchains[idx]}"
            right = f"{model_tag}"
            title = f"{model_tag}({_cluster_short_disp(cluster)}) vs. {left}: TM={tm:.2f}"
            _add_png_legend(out, left, right, title=title)
            print(f"[3D] wrote {out}")
        except Exception as e:
            print(f"[3D] {model_tag} fold{idx+1} overlay failed: {e}")

        if make_interactive:
            html_out = os.path.join(fig_dir_root, f"{pair_id}_fold{idx+1}_vs_best_{model_tag}.html")
            try:
                write_py3dmol_overlay_html(
                   true_pdb, pred_pdb, html_out,
                    left_label=f"{pdbids[idx]}{pdbchains[idx]}",
                    right_label=model_tag,
                    title=title, color_left="red", color_right="blue"
                )
                print(f"[3D] wrote {html_out}")
            except Exception as e:
                print(f"[3D] {model_tag} interactive failed: {e}")


def _render_ddg_aligned_image(pair_id: str, pdbids: list[str], pdbchains: list[str], fig_dir_root: str) -> None:
    """
    Build the aligned per-residue ΔΔG image into output_figs/,
    using cached CSVs written by the compute_deltaG step.

    Output: Pipeline/<pair>/output_figs/<pair>_ddg_aligned.png
    """
    os.makedirs(fig_dir_root, exist_ok=True)
    pair_dir = os.path.join("Pipeline", pair_id)
    anal_dir = os.path.join(pair_dir, "Analysis")
    out_png = os.path.join(fig_dir_root, f"{pair_id}_ddg_aligned.png")

    print("fig_dir_root:", fig_dir_root)
    print("out_png:", out_png)

    # Try to load per-residue energies for both folds.
    res_lists = []
    for i in (0, 1):
        # Accept both naming schemes: df_energy_2qqjA.csv or df_energy_2qqj.csv
        token_with_chain = f"{pdbids[i]}{pdbchains[i]}"
        candidates = [
            os.path.join(anal_dir, f"df_energy_{token_with_chain}.csv"),
            os.path.join(anal_dir, f"df_energy_{pdbids[i]}.csv"),
        ]
        csv_path = next((p for p in candidates if os.path.isfile(p)), None)
        if not csv_path:
            print(f"[ΔΔG] missing per-residue CSV for fold {i+1}. "
                  f"Expected one of: {', '.join(candidates)}. "
                  f"Run: --run_mode compute_deltaG for this pair first.")
            return

        d = pd.read_csv(csv_path)
        # expect columns: residue_name, residue_index, energy
        if not set(["residue_name", "residue_index", "energy"]).issubset(d.columns):
            print(f"[ΔΔG] bad columns in {csv_path}; needs residue_name,residue_index,energy")
            return

        res_lists.append(list(d[["residue_name", "residue_index", "energy"]].itertuples(index=False, name=None)))

    # Make the image (align & draw)
    try:
        align_and_compare_residues(
            res_lists[0], res_lists[1],
            f"{pdbids[0]}{pdbchains[0]}", f"{pdbids[1]}{pdbchains[1]}",
            output_file=out_png
        )
        print(f"[ΔΔG] wrote {out_png}")
    except Exception as e:
        print(f"[ΔΔG] align/plot failed for {pair_id}: {e}")


# ----------- Plots Generation Functions -----------
def write_py3dmol_overlay_html(
    pdb1: str,
    pdb2: str,
    out_html: str,
    *,
    left_label: str = "Fold1",
    right_label: str = "Fold2",
    title: str | None = None,
    color_left: str = "red",
    color_right: str = "blue",
    width: int = 700,
    height: int = 520,
):
    import py3Dmol, html as _html

    v = py3Dmol.view(width=width, height=height)
    with open(pdb1) as f:
        v.addModel(f.read(), "pdb")
    v.setStyle({"model": 0}, {"cartoon": {"color": color_left}})
    with open(pdb2) as f:
        v.addModel(f.read(), "pdb")
    v.setStyle({"model": 1}, {"cartoon": {"color": color_right}})
    v.setBackgroundColor("white")
    v.zoomTo()

    body = v._make_html()
    # Title + legend (colored squares; no “Red/Blue” words)
    title_html = f"<h3 style='font-family:sans-serif;margin:6px 0 6px 0;'>{_html.escape(title) if title else ''}</h3>"
    legend_html = f"""
<div style="display:flex;justify-content:center;gap:18px;margin:4px 0 10px 0;font-family:sans-serif;">
  <span><span style="display:inline-block;width:12px;height:12px;background:{color_left};margin-right:6px;border:1px solid #444;"></span>{_html.escape(left_label)}</span>
  <span><span style="display:inline-block;width:12px;height:12px;background:{color_right};margin-right:6px;border:1px solid #444;"></span>{_html.escape(right_label)}</span>
</div>"""

    with open(out_html, "w") as f:
        f.write(title_html + legend_html + body)




def make_foldswitch_all_plots(
    pdbids, fasta_dir, foldpair_id, pdbchains,
    plot_tree_clusters: bool = False, plot_contacts: bool = True, global_plots: bool = False, plot3dformat: str = 'static'
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

    msa_pred_files = glob(os.path.join(fasta_dir, foldpair_id, "output_cmaps", "msa_transformer", "*.npz"))
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
    print("msa_pred_files num=", len(msa_pred_files))
    msa_transformer_pred = {}
    for path in msa_pred_files:
        base = os.path.splitext(os.path.basename(path))[0]  # e.g., "msa_t__ShallowMsa_008"
        key = base[len("msa_t__"):] if base.startswith("msa_t__") else base
        try:
            pair_cmap, pair_idx, keep_cols = load_cmap_and_idx(path) #  np.load(path, allow_pickle=True)
            msa_transformer_idx[key] = (pair_cmap, pair_idx)
        except Exception:
            msa_transformer_pred[key] = np.genfromtxt(path)
    print("msa_transformer_pred keys:", list(msa_transformer_pred.keys()), " num=", len(msa_transformer_pred))

    if not msa_transformer_pred:
        print(f"[plot:{foldpair_id}] No MSA-Transformer .npz files — skipping contact panels.")
        return  # skip

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
        seqs[pdbids[1] + pdbchains[1]])
    print("Get matching indices: pdbids", pdbids, "pdbchains", pdbchains)
    print("true_cmap type and len, and len([0]):", type(true_cmap), len(true_cmap), len(true_cmap[pdbids[0] + pdbchains[0]]))
    # new (safe):
    first_key = next(iter(msa_transformer_pred.keys()), None)
    first_len = (msa_transformer_pred[first_key].shape[0] if first_key is not None else 0)
    print("msa_transformer_pred keys:", list(msa_transformer_pred.keys()), "first_len:", first_len)

    print(" true_cmaps sizes: ", [cmap.shape for cmap in true_cmap.values()])
    print(" pred_cmaps sizes: ", [cmap.shape for cmap in msa_transformer_pred.values()])
    msa_full_default = os.path.join("Pipeline", foldpair_id, "output_get_msa", "DeepMsa.a3m")
    msa_full_trim = os.path.join("Pipeline", foldpair_id, "output_get_msa", "DeepMsa_pairtrim.a3m")
    msa_full = msa_full_trim if os.path.isfile(msa_full_trim) else msa_full_default
    keyA = pdbids[0] + pdbchains[0]  # e.g. "2qqjA"
    keyB = pdbids[1] + pdbchains[1]  # e.g. "4qdsA"

    match_true_cmap, match_predicted_cmaps = align_truth_and_preds_via_msa_first2(
        true_cmap=true_cmap,
        pred_cmaps=msa_transformer_pred,
        msa_path=msa_full,
        keyA=keyA,
        keyB=keyB,
        verbose=True,
    )

    print("After align match_true_cmaps sizes: ", [cmap.shape for cmap in match_true_cmap.values()])
    print("After align match_pred_cmaps sizes: ", [cmap.shape for cmap in match_predicted_cmaps.values()])

    # ---------- Plot CMAPs ----------
    if plot_contacts:
        print("All match_true len and match_predicted_cmaps len:",  len(match_true_cmap), len(match_predicted_cmaps))
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

    # ---- NEW: ΔΔG aligned profile image into output_figs/ ----
    try:
        _render_ddg_aligned_image(foldpair_id, pdbids, pdbchains, fig_dir_root)
    except Exception as e:
        print(f"[ΔΔG] skipped for {foldpair_id}: {e}")


    # ---------- 3D PANELS (PyMOL headless) ----------
    # --------- Optional 3D-align plot (only on non-Linux in your original code)  --------------
    # 3D PANELS
    make_interactive = plot3dformat in ('interactive', 'both')
    try:
        _render_true_structures(
            pair_dir=os.path.join('Pipeline', foldpair_id),
            pdb1=f"{pdbids[0]}.pdb",
            pdb2=f"{pdbids[1]}.pdb",
            out_dir=fig_dir_root,
            out_prefix=foldpair_id,
            make_interactive=make_interactive
        )
    except Exception as e:
        print(f"[3D] true-true render failed: {e}")


    # AF2 / AF3
#    _render_true_vs_best_models_generic(foldpair_id, pdbids, pdbchains, fig_dir_root, family='AF',  ver='2')
#    _render_true_vs_best_models_generic(foldpair_id, pdbids, pdbchains, fig_dir_root, family='AF',  ver='3')

    _render_true_vs_best_models_generic(foldpair_id, pdbids, pdbchains, fig_dir_root, family='AF', ver='2',
                                        make_interactive = make_interactive)
    _render_true_vs_best_models_generic(foldpair_id, pdbids, pdbchains, fig_dir_root, family='AF', ver='3',
                                        make_interactive = make_interactive)

#    # ESM2 / ESM3  (NEW)
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

    return cmap_dists_vec, seqs_dists_vec, num_seqs_msa_vec, concat_scores


# Align two pdb structures from two pdb files
def align_and_visualize_proteins(pdb_file1, pdb_file2, output_file, open_environment=True):
    """
    Align two protein structures and save the visualization as an image.
    Uses PyMOL Python API if available; otherwise falls back to PyMOL CLI via PYMOL_BIN.
    """
    out_dir = os.path.dirname(os.path.abspath(output_file)) or "."
    os.makedirs(out_dir, exist_ok=True)

    if PYMOL_AVAILABLE:
        try:
            if open_environment:
                pymol.finish_launching(['pymol', '-cq'])
            cmd.delete('all')
            cmd.load(pdb_file1, 'protein1')
            cmd.load(pdb_file2, 'protein2')
            cmd.align('protein2', 'protein1')
            cmd.color('red', 'protein1')
            cmd.color('blue', 'protein2')
            cmd.bg_color('white')
            cmd.zoom('all', buffer=10)
            cmd.png(output_file, dpi=220, width=1200, height=900)
        finally:
            if open_environment:
                try: cmd.quit()
                except Exception: pass
        return

    # ---- CLI fallback ----
    exe = os.environ.get("PYMOL_BIN") or _PYMOL_BIN_CACHED
    if not exe:
        raise RuntimeError("PyMOL module not available and no CLI found.")

    p1 = os.path.abspath(pdb_file1)
    p2 = os.path.abspath(pdb_file2)
    out = os.path.abspath(output_file)

    # Keep quotes for load paths (safe), but **no quotes** for png target.
    script = f"""
reinitialize
load "{p1}", protein1
load "{p2}", protein2
align protein2, protein1
color red, protein1
color blue, protein2
bg_color white
viewport 1200, 900
zoom all, 10
set antialias, 2
png {out}, dpi=220, width=1200, height=900
quit
""".lstrip()

    if '_pymol_cli_render' in globals() and callable(globals()['_pymol_cli_render']):
        _pymol_cli_render(script)
    else:
        with tempfile.NamedTemporaryFile("w", suffix=".pml", delete=False) as tf:
            tf.write(script); tf.flush()
            cmdline = f"{shlex.quote(PYMOL_BIN)} -cq {shlex.quote(tf.name)}"
            subprocess.run(cmdline, shell=True, check=True)

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
    Draw:
      • Grid of ALL clusters found in `predictions` -> {save_file}.png
      • One-panel BEST clusters (one cluster chosen per fold) -> {save_file.replace('all','best')}.png
    """
    preds = dict(predictions or {})
    n_pred = len(preds)
    if n_pred == 0:
        print("[plot] no predictions to plot.")
        return

    # choose a near-square grid
    n_row = int(math.ceil(math.sqrt(n_pred)))
    n_col = n_row - 1 if n_row * (n_row - 1) >= n_pred else n_row

    fig, axes = plt.subplots(n_row, n_col, figsize=(12, 12), layout="compressed")
    axes = np.array(axes).reshape(-1)  # flatten safely for any 1D/2D case

    # ---- draw all clusters and collect per-fold recall ----
    fold_ids = list(contacts.keys())
    recall: dict[str, dict[str, float]] = {}
    for i, name in enumerate(preds.keys()):
        ax = axes[i]
        recall[name] = plot_foldswitch_contacts_and_predictions(
            preds[name], contacts, ax=ax, title=_cluster_short_disp(name), show_legend=False)


    # hide leftover cells
    for j in range(n_pred, len(axes)):
        axes[j].axis("off")

    # ---- single global legend below the mosaic (one narrow row) ----
    legend_handles = [
        Patch(facecolor='lightgray', edgecolor='#444', label='Shared Contacts'),
        Patch(facecolor='darkgray', edgecolor='#444', label='Unique Contacts'),
        Line2D([0], [0], marker='o', linestyle='None', color='red', label='False Positives', markersize=6),
        Line2D([0], [0], marker='o', linestyle='None', color='blue', label='True Shared Positives', markersize=6),
        Line2D([0], [0], marker='o', linestyle='None', color='green', label='True Unique Positives', markersize=6),
    ]
    legend_labels = [h.get_label() for h in legend_handles]

    # put legend just below the axes area
    fig.legend(handles=legend_handles, labels=legend_labels, loc='lower center', bbox_to_anchor=(0.5, -0.02),
        ncol=5, frameon=False, fontsize=8)
    # leave a bit more room for the legend
    fig.subplots_adjust(bottom=0.18)


    if save_file:
        plt.savefig(save_file + ".png", dpi=200)
        print("Save cmap fig:", save_file + ".png")
    plt.close(fig)

    if not recall:
        return

    # exclude deep from the "best" picker (as before)
    filtered = {k: v for k, v in recall.items() if "deep" not in k.lower()}
    pick_from = filtered if filtered else recall

    best = find_max_keys(pick_from)  # {'foldA': (cluster, val), 'foldB': (...)}
    print("Best recall clusters:", best)
    best_names = {fold_ids[0]: best[fold_ids[0]][0], fold_ids[1]: best[fold_ids[1]][0]}
    print("Best recall cluster names:", best_names)

    # optional ΔG overlay
    xvec = yvec = None
    if foldpair_id:
        energy_dir = os.path.join("Pipeline", foldpair_id, "output_deltaG")
        try:
            e0 = read_energy_tuples(os.path.join(energy_dir, f"deltaG_{fold_ids[0][:4]}.txt"))
            e1 = read_energy_tuples(os.path.join(energy_dir, f"deltaG_{fold_ids[1][:4]}.txt"))
            tmp_delta = align_and_compare_residues(e0, e1, fold_ids[0][:4], fold_ids[1][:4])
            delta_energies_filtered = tmp_delta[1]
            xvec = yvec = np.array(delta_energies_filtered[:len(contacts[next(iter(contacts))])])
        except Exception as e:
            print(f"[plot] NOTE: ΔG overlay skipped for {foldpair_id}: {e}")

    # render BEST single panel
    plt.figure(figsize=(10, 8))
    plot_foldswitch_contacts_and_predictions(
        predictions=(preds[best_names[fold_ids[0]]], preds[best_names[fold_ids[1]]]),
        contacts=contacts,
        title="Best clusters",
        show_legend=True,
        cluster_names=(
            _cluster_short_disp(best_names[fold_ids[0]]),
            _cluster_short_disp(best_names[fold_ids[1]]),
        ),
        x_vector=xvec, y_vector=yvec,
    )
    if save_file:
        plt.savefig(save_file.replace("all", "best") + ".png", dpi=200)
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
        ax.set_xlabel(f"{fold_ids[0]}, recall={round(recall[fold_ids[0]], 4)}", fontsize=12)
        ax.set_ylabel(f"{fold_ids[1]}, recall={round(recall[fold_ids[1]], 4)}", fontsize=12)
    else:
        cx = _cluster_short_disp(cluster_names[0])
        cy = _cluster_short_disp(cluster_names[1])
        ax.set_xlabel(f"{fold_ids[0]}, recall={round(recall[fold_ids[0]], 4)} — {cx}", fontsize=12)
        ax.set_ylabel(f"{fold_ids[1]}, recall={round(recall[fold_ids[1]], 4)} — {cy}", fontsize=12)


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

    # Normalize model tags (robust if 'model' is missing)
    if "model" not in df_det.columns:
        df_det["model"] = ""  # create an empty column with the right length
    df_det["model"] = df_det["model"].astype(str).str.upper()

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


    # === NEW: cross-method pairwise comparisons ===
    # We use the unified *summary* table because it already contains the “Deep” and “Best Cluster/Best Prediction” summaries.
    if (not os.path.isfile(SUMMARY_RESULTS_TABLE)) or os.path.getsize(SUMMARY_RESULTS_TABLE) == 0:
        print(f"[warn] No summary CSV at {SUMMARY_RESULTS_TABLE}. Skipping cross-method plots.")
        return

    df_sum = pd.read_csv(SUMMARY_RESULTS_TABLE)

    # columns we rely on (present in DEFAULT_PREFERRED_COLS in gen_html_table.py)
    # AF2/AF3 Deep:
    cols_deep = dict(
        AF2_TM1="AF2Deep_TM1", AF2_TM2="AF2Deep_TM2",
        AF3_TM1="AF3Deep_TM1", AF3_TM2="AF3Deep_TM2",
    )
    # AF2/AF3 Best Cluster (shallow only):
    cols_bestcl = dict(
        AF2_TM1="AF2Clust_TM1", AF2_TM2="AF2Clust_TM2",
        AF3_TM1="AF3Clust_TM1", AF3_TM2="AF3Clust_TM2",
    )
    # ESM2/ESM3 Best Prediction (may be different sequence per fold):
    cols_esm = dict(
        ESM2_TM1="ESM2_TM1", ESM2_TM2="ESM2_TM2",
        ESM3_TM1="ESM3_TM1", ESM3_TM2="ESM3_TM2",
    )


    def _split_pair_id(pid: str):
        # e.g., "2qqjA_4qdsA" -> ("2qqjA","4qdsA")
        s = str(pid)
        return tuple(s.split("_", 1)) if "_" in s else (s, "fold2")

    def _plot_xy_pairs(
        df: pd.DataFrame,
        x1_col: str, x2_col: str,
        y1_col: str, y2_col: str,
        title: str,
        outfile: str
    ):
        sub = df.copy()
        # numeric-ize
        for c in [x1_col, x2_col, y1_col, y2_col]:
            if c not in sub.columns:
                print(f"[skip] Missing column {c} for {title}")
                return
            sub[c] = sub[c].map(numify)
        sub = sub.dropna(subset=[x1_col, x2_col, y1_col, y2_col])
        if sub.empty:
            print(f"[info] No rows for plot: {title}")
            return

        # build points per pair: two points (fold1, fold2) and connect them
        plt.figure(figsize=(6.4, 6.0))
        # y=x
        plt.plot([0,1],[0,1], 'k-', linewidth=1, alpha=0.8)

        for _, r in sub.iterrows():
            pid = r.get("pair_id") or r.get("fold_pair") or "pair"
            f1, f2 = _split_pair_id(pid)

            x1, y1 = r[x1_col], r[y1_col]  # fold1
            x2, y2 = r[x2_col], r[y2_col]  # fold2

            # connector
            plt.plot([x1, x2], [y1, y2], linestyle='--', linewidth=0.8, alpha=0.5)

            # points
            plt.scatter(x1, y1, s=22, marker='o')
            plt.scatter(x2, y2, s=22, marker='o')

            # light labels (small)
            try:
                plt.text(x1, y1, f1, fontsize=6, ha='right', va='bottom', alpha=0.9)
                plt.text(x2, y2, f2, fontsize=6, ha='left', va='top', alpha=0.9)
            except Exception:
                pass

        plt.xlim(0, 1); plt.ylim(0, 1)
        plt.xlabel("X method TM-score")
        plt.ylabel("Y method TM-score")
        plt.title(title)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        out_path = os.path.join(output_dir, outfile)
        plt.savefig(out_path, dpi=180)
        plt.close()
        print(f"[plot] wrote {out_path}")

    # 1) AF2 (Deep) vs AF3 (Deep)
    _plot_xy_pairs(
        df_sum,
        x1_col=cols_deep["AF2_TM1"], x2_col=cols_deep["AF2_TM2"],
        y1_col=cols_deep["AF3_TM1"], y2_col=cols_deep["AF3_TM2"],
        title="AF2 (Deep) vs AF3 (Deep) — TM-scores (per fold, pairs linked)",
        outfile="compare_AF2_AF3_DEEP.png"
    )

    # 2) AF2 (Best Shallow Cluster) vs AF3 (Best Shallow Cluster)
    _plot_xy_pairs(
        df_sum,
        x1_col=cols_bestcl["AF2_TM1"], x2_col=cols_bestcl["AF2_TM2"],
        y1_col=cols_bestcl["AF3_TM1"], y2_col=cols_bestcl["AF3_TM2"],
        title="AF2 (Best Shallow Cluster) vs AF3 (Best Shallow Cluster)",
        outfile="compare_AF2_AF3_BESTCLUST.png"
    )

    # 3) ESM2 (Best Prediction) vs ESM3 (Best Prediction)
    _plot_xy_pairs(
        df_sum,
        x1_col=cols_esm["ESM2_TM1"], x2_col=cols_esm["ESM2_TM2"],
        y1_col=cols_esm["ESM3_TM1"], y2_col=cols_esm["ESM3_TM2"],
        title="ESM2 (Best Prediction) vs ESM3 (Best Prediction)",
        outfile="compare_ESM2_ESM3_BESTPRED.png"
    )



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
