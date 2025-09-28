# postprocess_unified.py
import os, re, glob
import pandas as pd
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple
import fnmatch

from tqdm import tqdm
import argparse

from config import DATA_DIR, SUMMARY_RESULTS_TABLE, DETAILED_RESULTS_TABLE, MSA_TRANS_MODEL_FILE

from utils.utils import list_protein_pairs, pair_str_to_tuple
from utils.align_utils import compute_tmscore_align
from Analysis.cmap_analysis import compute_cmap_metrics_for_pair


PAIR_DIR = Path(DATA_DIR)


# === NEW: fast builder from cached per-pair cluster CSVs to unified CSVs ===
def _deepmsa_a3m_path(pair_id: str) -> str:
    # exact path you asked for
    return f"{DATA_DIR}/{pair_id}/output_get_msa/DeepMsa.a3m"

def _count_a3m_sequences_fast(a3m_path: str) -> int:
    """Count '>' headers in a plain-text A3M file."""
    try:
        n = 0
        with open(a3m_path, "r") as fh:
            for line in fh:
                if line.startswith(">"):
                    n += 1
        return n
    except Exception:
        return 0

def _count_shallow_clusters_fast(pair_id: str) -> int:
    """Count ShallowMsa_*.a3m files directly under output_msa_cluster/."""
    base = f"{DATA_DIR}/{pair_id}/output_msa_cluster"
    try:
        names = os.listdir(base)
    except Exception:
        return 0
    return sum(1 for name in names if fnmatch.fnmatch(name, "ShallowMsa_*.a3m"))

def _safe_read_csv(path: str) -> Optional[pd.DataFrame]:
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    try:
        return pd.read_csv(path)
    except Exception:
        return None

def _norm_tm_df(df: Optional[pd.DataFrame], default_source: str) -> pd.DataFrame:
    """
    Normalize AF/ESM TM tables into a single schema:
      columns: model (af2|af3|esm2|esm3), cluster_tag, cluster_type ('clust'|'deep'|None),
               cluster_id (e.g., '7' for ShallowMsa_007), TM1, TM2
    """
    import re
    SHALLOW_RE = re.compile(r"ShallowMsa_(\d+)", re.IGNORECASE)
    if df is None or df.empty:
        return pd.DataFrame(columns=["model","cluster_tag","cluster_type","cluster_id","TM1","TM2"])

    df = df.copy()
    # Accept both new and legacy column names
    if "TMscore_fold1" not in df.columns and "score_pdb1" in df.columns:
        df = df.rename(columns={"score_pdb1": "TMscore_fold1", "score_pdb2": "TMscore_fold2"})
    # Model normalization
    if "model" not in df.columns:
        df["model"] = default_source
    df["model"] = df["model"].astype(str).str.lower()  # 'AF2' -> 'af2', etc.

    # Cluster tag
    if "cluster_num" in df.columns:
        ctag = df["cluster_num"].astype(str)
    elif "cluster" in df.columns:
        ctag = df["cluster"].astype(str)
    else:
        ctag = pd.Series([None] * len(df))
    df["_cluster_tag"] = ctag

    # Parse type/id
    def _ctype_id(s):
        if s is None or s != s:
            return None, None
        st = str(s)
        if "deep" in st.lower():
            return "deep", None
        m = SHALLOW_RE.search(st)
        if m:
            cid = m.group(1).lstrip("0") or "0"
            return "clust", cid
        if st.isdigit():
            return "clust", st.lstrip("0") or "0"
        return None, None

    types, ids = zip(*[_ctype_id(s) for s in df["_cluster_tag"]])
    out = pd.DataFrame({
        "model": df["model"].str.lower(),
        "cluster_tag": df["_cluster_tag"],
        "cluster_type": list(types),
        "cluster_id": list(ids),
        "TM1": pd.to_numeric(df.get("TMscore_fold1"), errors="coerce"),
        "TM2": pd.to_numeric(df.get("TMscore_fold2"), errors="coerce"),
    })
    return out

def _pick_best(dfall: pd.DataFrame, model_tag: str, mode: str, which_fold: int) -> str:
    """
    mode: 'clust' or 'deep'; which_fold: 1 or 2
    Returns formatted string:
      - Clust: "0.78 (7)"
      - Deep:  "0.64"
      - "-" if unavailable
    """
    if dfall is None or dfall.empty:
        return "-"
    df = dfall[dfall["model"] == model_tag]
    if mode == "clust":
        df = df[df["cluster_type"] == "clust"]
    elif mode == "deep":
        df = df[df["cluster_type"] == "deep"]
    col = "TM1" if which_fold == 1 else "TM2"
    if df.empty or df[col].isna().all():
        return "-"
    idx = df[col].idxmax()
    val = float(df.loc[idx, col])
    if mode == "clust":
        cid = df.loc[idx, "cluster_id"]
        if cid is not None and cid == cid:
            return f"{val:.2f} ({cid})"
    return f"{val:.2f}"

def _pick_best_overall(dfall: pd.DataFrame, model_tag: str, which_fold: int) -> str:
    if dfall is None or dfall.empty:
        return "-"
    df = dfall[dfall["model"] == model_tag]
    col = "TM1" if which_fold == 1 else "TM2"
    s = pd.to_numeric(df[col], errors="coerce")
    return f"{s.max():.2f}" if s.notna().any() else "-"

def _best_max(df: Optional[pd.DataFrame], col: str) -> str:
    if df is None or df.empty or col not in df.columns:
        return "-"
    s = pd.to_numeric(df[col], errors="coerce")
    return f"{s.max():.2f}" if s.notna().any() else "-"

def _pair_max_len_from_truth(pair_id: str) -> int:
    pdb1, c1, pdb2, c2 = _truth_pdbs(pair_id)
    def _len(pdb_path, chain_id):
        n, seen = 0, set()
        try:
            with open(pdb_path) as fh:
                for line in fh:
                    if not line.startswith("ATOM"): continue
                    if line[12:16].strip() != "CA": continue
                    if chain_id and (line[21].strip() != chain_id): continue
                    resnum = (line[22:27], line[26])
                    if resnum not in seen:
                        seen.add(resnum); n += 1
        except Exception:
            pass
        return n
    return max(_len(pdb1, c1), _len(pdb2, c2))

def build_unified_tables_from_cluster_dfs(pairs: Optional[List[str]] = None,
                                          write_out: bool = True) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Fast path:
      - READ cached per-pair Analysis CSVs (df_af.csv, df_esm.csv, df_cmap.csv) – no recompute
      - COMPUTE summary (one row per pair) incl. "#RES" and "MSA DEPTH (#Clusters)"
      - COMPUTE detailed (concat of df_af/df_esm with fold_pair column)
      - Optionally write SUMMARY_RESULTS_TABLE and DETAILED_RESULTS_TABLE
    """
    if not pairs:
        pairs = [f"{a}_{b}" for a, b in list_protein_pairs()]

    all_detailed = []
    summary_rows = []

    for pair_id in pairs:
        pair_dir = _pair_dir(pair_id)
        if not pair_dir.is_dir():
            continue

        df_af   = _safe_read_csv(str(_ensure_pair_analysis(pair_id) / "df_af.csv"))
        df_esm  = _safe_read_csv(str(_ensure_pair_analysis(pair_id) / "df_esm.csv"))
        df_cmap = _safe_read_csv(str(_ensure_pair_analysis(pair_id) / "df_cmap.csv"))

        # Detailed rows (cluster-level)
        det = []
        if df_af is not None and not df_af.empty:
            det.append(df_af.copy())
        if df_esm is not None and not df_esm.empty:
            det.append(df_esm.copy())
        if det:
            det = pd.concat(det, ignore_index=True)
            det["fold_pair"] = pair_id
            all_detailed.append(det)

        # Build best-per-pair row + depth & #clusters
        tm_af  = _norm_tm_df(df_af,  "af2")
        tm_esm = _norm_tm_df(df_esm, "esm2")
        tm_all = pd.concat([tm_af, tm_esm], ignore_index=True) if len(tm_af) or len(tm_esm) else pd.DataFrame()

        row = {"fold_pair": pair_id, "#RES": _pair_max_len_from_truth(pair_id)}

        deepmsa_file = _deepmsa_a3m_path(pair_id)
        msa_depth = _count_a3m_sequences_fast(deepmsa_file)
        n_clusters = _count_shallow_clusters_fast(pair_id)
        row["MSA DEPTH (#Clusters)"] = f"{msa_depth} ({n_clusters})"

        # --- TRUE TM between the two known folds (symmetric max) ---
        pdb1, c1, pdb2, c2 = _truth_pdbs(pair_id)
        # note: compute both directions; TM-score is length-normalized
        tm12 = compute_tmscore_align(pdb1, pdb2, chain1=c1, chain2=c2)
        tm21 = compute_tmscore_align(pdb2, pdb1, chain1=c2, chain2=c1)
#        row["TRUE_TM12"] = tm12
#        row["TRUE_TM21"] = tm21
        row["PAIR_TM"] = round(max(tm12, tm21), 3)  # single symmetric number to display
        print("TM12:", tm12, "TM21:", tm21, "MAX-rounded:", row["PAIR_TM"])

        # AF: Clust + Deep, per fold
        for tag, up in (("af2","AF2"), ("af3","AF3")):
            row[f"{up}Clust_TM1"] = _pick_best(tm_all, tag, "clust", 1)
            row[f"{up}Clust_TM2"] = _pick_best(tm_all, tag, "clust", 2)
            row[f"{up}Deep_TM1"]  = _pick_best(tm_all, tag, "deep",  1)
            row[f"{up}Deep_TM2"]  = _pick_best(tm_all, tag, "deep",  2)

        # ESM overall bests
        row["ESM2_TM1"] = _pick_best_overall(tm_all, "esm2", 1)
        row["ESM2_TM2"] = _pick_best_overall(tm_all, "esm2", 2)
        row["ESM3_TM1"] = _pick_best_overall(tm_all, "esm3", 1)
        row["ESM3_TM2"] = _pick_best_overall(tm_all, "esm3", 2)

        # CMAP maxima (MSA-Transformer)
        row["MSATrans_CMAP_PR1"] = _best_max(df_cmap, "t1_precision")
        row["MSATrans_CMAP_PR2"] = _best_max(df_cmap, "t2_precision")
        row["MSATrans_CMAP_RE1"] = _best_max(df_cmap, "t1_recall")
        row["MSATrans_CMAP_RE2"] = _best_max(df_cmap, "t2_recall")

        summary_rows.append(row)

    detailed_df = pd.concat(all_detailed, ignore_index=True) if all_detailed else pd.DataFrame()
    summary_df  = pd.DataFrame(summary_rows)

    if write_out:
        Path(SUMMARY_RESULTS_TABLE).parent.mkdir(parents=True, exist_ok=True)
        summary_df.to_csv(SUMMARY_RESULTS_TABLE, index=False)
        detailed_df.to_csv(DETAILED_RESULTS_TABLE, index=False)

    return summary_df, detailed_df



def _cmap_csv_path(pair_id: str) -> str:
    return f"{DATA_DIR}/{pair_id}/Analysis/df_cmap.csv"

def _read_or_compute_cmap(pair_id: str, force: bool = False) -> pd.DataFrame:
    out_csv = _cmap_csv_path(pair_id)
    pred_dir = f"{DATA_DIR}/{pair_id}/output_cmaps/msa_transformer"

    # use cached CSV when present (unless forced)
    if (not force) and os.path.isfile(out_csv):
        return pd.read_csv(out_csv)

    # need predictions to compute metrics
    if not os.path.isdir(pred_dir) or not any(f.endswith(".npy") for f in os.listdir(pred_dir)):
        print(f"[warn] No MSA-Transformer NPYs found at {pred_dir}; skipping CMAP metrics.")
        return pd.DataFrame()

    # import the callable (local import avoids import cycles at module load time)
    return compute_cmap_metrics_for_pair(pair_id, include_deep=True, thresh=0.4, sep_min=6, index_tol=0)



def _pair_dir(pair_id: str) -> Path:
    return PAIR_DIR / pair_id

def _truth_pdbs(pair_id: str) -> tuple[str, str, str, str]:
    """Return (pdb1_path, chain1, pdb2_path, chain2). Prefers chain-sliced PDBs if present."""
    a, b = pair_str_to_tuple(pair_id)
    p1, c1 = a[:-1], a[-1]
    p2, c2 = b[:-1], b[-1]
    cand1 = _pair_dir(pair_id) / "chain_pdb_files" / f"{a}.pdb"
    cand2 = _pair_dir(pair_id) / "chain_pdb_files" / f"{b}.pdb"
    pdb1 = str(cand1 if cand1.is_file() else (_pair_dir(pair_id) / f"{p1}.pdb"))
    pdb2 = str(cand2 if cand2.is_file() else (_pair_dir(pair_id) / f"{p2}.pdb"))
    return pdb1, c1, pdb2, c2

def _ensure_pair_analysis(pair_id: str) -> Path:
    out = _pair_dir(pair_id) / "Analysis"
    out.mkdir(parents=True, exist_ok=True)
    return out

def _read_or_compute_af(pair_id: str, force: bool) -> pd.DataFrame:
    """Aggregate AF2/AF3 TM-scores across all clusters/chains."""
    out_csv = _ensure_pair_analysis(pair_id) / "df_af.csv"
    if out_csv.is_file() and not force:
        return pd.read_csv(out_csv)

    pdb1, c1, pdb2, c2 = _truth_pdbs(pair_id)
    rows = []

    # ONLY canonical top-level PDBs:
    for ver in ("AF2", "AF3"):
        top = _pair_dir(pair_id) / "output_AF" / ver
        if not top.is_dir():
            continue
        for pred in sorted(top.glob("*.pdb")):  # <-- no recursion
            name = pred.name.replace(".pdb", "")
            # cluster: DeepMsa or ShallowMsa_###
            m = re.search(r"(DeepMsa|ShallowMsa_\d+)", name)
            cluster = m.group(1) if m else "UNK"

            tm1 = compute_tmscore_align(pdb1, str(pred), chain2=None)
            tm2 = compute_tmscore_align(pdb2, str(pred), chain2=None)

            rows.append({
                "fold_pair": pair_id,
                "model": ver,  # AF2 vs AF3
                "cluster_num": cluster,
                "name": name,  # short
                "TMscore_fold1": tm1,
                "TMscore_fold2": tm2
            })

    df = pd.DataFrame(rows)
    if len(df):
        df["TMdiff"] = df["TMscore_fold1"] - df["TMscore_fold2"]
        df = df[["fold_pair", "model", "cluster_num", "name", "TMscore_fold1", "TMscore_fold2", "TMdiff"]]
        df.sort_values(["model", "cluster_num", "name"], inplace=True)
        df.to_csv(out_csv, index=False)
    return df


def _read_or_compute_esm(pair_id: str, force: bool) -> pd.DataFrame:
    """Read per-sample ESM TM-scores if cached; else compute quickly by scanning outputs."""
    out_csv = _ensure_pair_analysis(pair_id) / "df_esm.csv"
    if out_csv.is_file() and not force:
        return pd.read_csv(out_csv)

    pdb1, c1, pdb2, c2 = _truth_pdbs(pair_id)
    rows = []
    for model_tag in ("esm2", "esm3"):
        mdir = _pair_dir(pair_id) / "output_esm_fold" / model_tag
        if not mdir.is_dir():
            continue
        idx = mdir / "samples_index.tsv"
        if idx.is_file():
            df_idx = pd.read_csv(idx, sep="\t")
            for _, r in df_idx.iterrows():
                pred = str(r["pdb_path"])
                name = r["name"]
                tm1 = compute_tmscore_align(pdb1, pred, chain2=None)
                tm2 = compute_tmscore_align(pdb2, pred, chain2=None)
                rows.append({"fold_pair": pair_id, "model": model_tag, "cluster_num":
                             (re.search(r"ShallowMsa_(\d+)", name).group(1)
                              if re.search(r"ShallowMsa_(\d+)", name) else "DeepMsa"),
                             "name": name, "pdb_path": pred,
                             "TMscore_fold1": tm1, "TMscore_fold2": tm2})
        else:
            for pred in sorted(mdir.glob("*.pdb")):
                base = pred.name
                name = base.replace(f"_{model_tag}.pdb","").replace(".pdb","")
                tm1 = compute_tmscore_align(pdb1, str(pred), chain2=None)
                tm2 = compute_tmscore_align(pdb2, str(pred), chain2=None)
                rows.append({"fold_pair": pair_id, "model": model_tag,
                             "cluster_num": (re.search(r"ShallowMsa_(\d+)", name).group(1)
                                             if re.search(r"ShallowMsa_(\d+)", name) else "DeepMsa"),
                             "name": name, "pdb_path": str(pred),
                             "TMscore_fold1": tm1, "TMscore_fold2": tm2})
    df = pd.DataFrame(rows)

    # normalize: cluster label and clean sample name (no paths)
    df["cluster_num"] = df["name"].str.extract(r"(ShallowMsa_\d+)", expand=False).fillna("DeepMsa")
    df["name"] = df["name"].apply(lambda s: os.path.basename(str(s)).replace(".pdb", ""))

    # helpful deltas
    df["TMdiff"] = df["TMscore_fold1"] - df["TMscore_fold2"]

    # ⬅ Save WITHOUT pdb_path column so the CSV is clean
    cols = ["fold_pair", "model", "cluster_num", "name", "TMscore_fold1", "TMscore_fold2", "TMdiff"]
    df = df[cols].sort_values(["model", "cluster_num", "name"])

    df.to_csv(out_csv, index=False)

    return df

def _read_cmap(pair_id: str) -> pd.DataFrame:
    """Read per-pair cmap metrics produced by cmap_analysis.py (don’t recompute here)."""
    csv = _pair_dir(pair_id) / "Analysis" / "df_cmap.csv"
    return pd.read_csv(csv) if csv.is_file() else pd.DataFrame()

# New: function for computing TM-score between the two true structures chains
def compute_tmscore_all_pairs():
    all_pairs = list_protein_pairs()
    num_pairs = len(all_pairs)
    tm_pairs_scores = {pair_id: [0.0, 0.0] for pair_id in all_pairs}  # Prepare empty dict
    for pair_id in all_pairs:
        pdb1, c1, pdb2, c2 = _truth_pdbs(pair_id)
        tm_pairs_scores[pair_id][0] = compute_tmscore_align(pdb1, pdb2, chain1=c1, chain2=c2)
        tm_pairs_scores[pair_id][1] = compute_tmscore_align(pdb2, pdb1, chain1=c2, chain2=c1)
        print(f"{pair_id}: {tm_pairs_scores[pair_id][0]:.3f} | {tm_pairs_scores[pair_id][1]:.3f}")

    return tm_pairs_scores


def build_pair_cluster_table(pair_id: str) -> pd.DataFrame:
    """
    Return a tidy cluster table for this pair with columns:
      cluster, TM-AF1, TM-AF2, TM-ESM1, TM-ESM2, RE-MSAT-COM, RE-MSAT1, RE-MSAT2
    Values may be NaN if a modality is missing for a cluster.
    """
    anal = _ensure_pair_analysis(pair_id)
    df_af   = _safe_read_csv(str(anal / "df_af.csv"))
    df_esm  = _safe_read_csv(str(anal / "df_esm.csv"))
    df_cmap = _safe_read_csv(str(anal / "df_cmap.csv"))

    def _norm_tag(s: str) -> str:
        if not isinstance(s, str): s = str(s)
        t = s.strip()
        if "deep" in t.lower(): return "DeepMsa"
        m = re.search(r"(\d+)$", t)
        if m: return f"ShallowMsa_{int(m.group(1)):03d}"
        if t.startswith("ShallowMsa_"): return t
        return t

    # ---- TM block ----
    tm = pd.DataFrame()
    def _pick_tm(df, prefix):
        if df is None or df.empty: return
        d = df.copy()
        if "cluster_num" not in d.columns and "cluster" in d.columns:
            d = d.rename(columns={"cluster": "cluster_num"})
        d["_tag"] = d["cluster_num"].astype(str).map(_norm_tag)
        t1 = pd.to_numeric(d.get("TMscore_fold1"), errors="coerce")
        t2 = pd.to_numeric(d.get("TMscore_fold2"), errors="coerce")
        agg = d.assign(TM1=t1, TM2=t2).groupby("_tag")[["TM1","TM2"]].max()
        for col, src in (("TM-AF1","TM1"),("TM-AF2","TM2")) if prefix=="AF" else (("TM-ESM1","TM1"),("TM-ESM2","TM2")):
            tm[col] = agg[src]
    _pick_tm(df_af, "AF")
    _pick_tm(df_esm, "ESM")
    if tm.empty:
        tm.index = []

    # ---- MSAT block ----
    ms = pd.DataFrame()
    if df_cmap is not None and not df_cmap.empty:
        d = df_cmap.copy()
        if "cluster_num" in d.columns:
            d["_tag"] = d["cluster_num"].astype(str).map(_norm_tag)
        elif "cluster" in d.columns:
            d["_tag"] = d["cluster"].astype(str).map(_norm_tag)
        else:
            d["_tag"] = None

        def pick(cands):
            for c in cands:
                if c in d.columns:
                    return pd.to_numeric(d[c], errors="coerce")
            return pd.Series(index=d.index, dtype=float)
        ms = pd.DataFrame({
            "RE-MSAT-COM": pick(["common_mcc","common_f1","common_recall","common_jaccard"]),
            "RE-MSAT1":    pick(["uniq1_mcc","t1_mcc","t1_f1","t1_recall","uniq1_recall","t1_jaccard"]),
            "RE-MSAT2":    pick(["uniq2_mcc","t2_mcc","t2_f1","t2_recall","uniq2_recall","t2_jaccard"]),
            "_tag": d["_tag"]
        }).groupby("_tag").max()

    # ---- merge blocks on index (cluster tag) ----
    idx = sorted(set(tm.index) | set(ms.index))
    out = pd.DataFrame(index=idx)
    for block in (tm, ms):
        if not block.empty:
            out = out.join(block, how="left")

    out = out.reset_index().rename(columns={"index": "cluster"})
    # prefer compact labels: "DeepMsa" -> "Deep", "ShallowMsa_007" -> "7"
    def _short(s):
        if not s: return s
        if str(s).lower().startswith("deep"): return "Deep"
        m = re.search(r"(\d+)$", str(s))
        return m.group(1) if m else s
    out["cluster"] = out["cluster"].map(_short)

    # nice rounding
    for c in ["TM-AF1","TM-AF2","TM-ESM1","TM-ESM2","RE-MSAT-COM","RE-MSAT1","RE-MSAT2"]:
        if c in out.columns:
            out[c] = out[c].round(2)
    return out[["cluster"] + [c for c in ["TM-AF1","TM-AF2","TM-ESM1","TM-ESM2","RE-MSAT-COM","RE-MSAT1","RE-MSAT2"] if c in out.columns]]


# Main post-processing function:
def post_processing_analysis(force_rerun: bool = False, pairs: Optional[List[str]] = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    - Computes (or reads cached) AF/ESM TM tables and cmap metrics.
    - Builds:
        * summary_df: one row per pair (best TM and best cmap metric)
        * detailed_df: one row per model/cluster/sample (TM & cmap when available)
    - Writes summary -> SUMMARY_RESULTS_TABLE, detailed -> DETAILED_RESULTS_TABLE
    """
    if not pairs:
        pairs = [f"{a}_{b}" for a, b in list_protein_pairs()]

    all_detailed = []
    summary_rows = []

    for pair_id in tqdm(pairs, desc="postproc"):
        pair_dir = _pair_dir(pair_id)
        if not pair_dir.is_dir():
            continue

        # AF + ESM tables
        df_af  = _read_or_compute_af(pair_id, force=force_rerun)
        df_esm = _read_or_compute_esm(pair_id, force=force_rerun)
        df_tm  = pd.concat([df_af, df_esm], ignore_index=True) if len(df_af) or len(df_esm) else pd.DataFrame()
        df_cmap = _read_or_compute_cmap(pair_id, force=force_rerun)

        # Best TM per pair (over all predictions; maximize max(TM1,TM2))
        best_tm = None
        if len(df_tm):
            df_tm["TM_best_vs_truth"] = df_tm[["TMscore_fold1","TMscore_fold2"]].max(axis=1)
            best_row = df_tm.loc[df_tm["TM_best_vs_truth"].idxmax()].to_dict()
            best_tm = {k: best_row[k] for k in ["model","cluster_num","name","TMscore_fold1","TMscore_fold2","TM_best_vs_truth"]}

        # Cmap metrics
        best_cmap = None
        if len(df_cmap):
            # define "best" by max common_f1 (use t1_f1/t2_f1 if you prefer)
            score_key = "common_f1" if "common_f1" in df_cmap.columns else ("t1_f1" if "t1_f1" in df_cmap.columns else None)
            if score_key:
                best_cmap_row = df_cmap.loc[df_cmap[score_key].idxmax()].to_dict()
                best_cmap = {"best_cmap_file": best_cmap_row.get("file"), "best_cmap_score": best_cmap_row.get(score_key)}

        # Detailed rows (cluster-level)
        if len(df_tm):
            det = df_tm.copy()
            det["fold_pair"] = pair_id
            all_detailed.append(det)

        max_len = _pair_max_len_from_truth(pair_id)
        # Summary row (pair-level)
        summary_rows.append({
            "fold_pair": pair_id,
            "#RES": max_len,
            **(best_tm or {}),
            **(best_cmap or {}),
            "n_af_preds": int(len(df_af)),
            "n_esm_preds": int(len(df_esm)),
            "n_cmap_preds": int(len(df_cmap))
        })


    detailed_df = pd.concat(all_detailed, ignore_index=True) if all_detailed else pd.DataFrame()
    summary_df  = pd.DataFrame(summary_rows)

    # New: add tm-scores for true structures pairs
#    tm_pairs_scores = compute_tmscore_all_pairs()
#    summary_df['PAIR_TM'] = summary_df.apply(lambda row: tm_pairs_scores[row['fold_pair']][0], axis=1)

    return summary_df, detailed_df

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description="Unified post-processing: build summary/detailed CSVs for the website.")
    p.add_argument("--pairs", nargs="*", help="Pair IDs like 1fzpD_2frhA. If omitted, process all pairs.")
    p.add_argument("--force_rerun", action="store_true",
                   help="FULL mode: recompute per-pair Analysis CSVs even if they exist.")
    p.add_argument("--mode", choices=["full", "unify"], default="full",
                   help="'full' recomputes per-pair Analysis CSVs then unifies; 'unify' only rebuilds unified CSVs from cached per-pair Analysis CSVs.")
    args = p.parse_args()

    if args.mode == "full":
        # 1) Recompute/refresh per-pair Analysis CSVs (df_af.csv, df_esm.csv, df_cmap.csv)
        summary_df, detailed_df = post_processing_analysis(force_rerun=args.force_rerun, pairs=args.pairs)
        print(f"[postprocess] per-pair Analysis refreshed: summary rows={len(summary_df)} | detailed rows={len(detailed_df)}")

    # 2) Build unified CSVs from those cached per-pair Analysis CSVs (single writer)
    build_unified_tables_from_cluster_dfs(pairs=args.pairs, write_out=True)
    print(f"[postprocess] unified CSVs written:\n  {SUMMARY_RESULTS_TABLE}\n  {DETAILED_RESULTS_TABLE}")
