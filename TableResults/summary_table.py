# TableResults/summary_table.py
import os, sys, re
import pandas as pd
from typing import Optional, Tuple, List

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, ROOT)

from config import DATA_DIR, DETAILED_RESULTS_TABLE, SUMMARY_RESULTS_TABLE
from utils.utils import list_protein_pairs

SHALLOW_RE = re.compile(r"ShallowMsa_(\d+)", re.IGNORECASE)

def _truth_pdbs_for_pair(pair_id: str):
    a, b = pair_id.split("_", 1)
    p1, c1 = a[:-1], a[-1]
    p2, c2 = b[:-1], b[-1]
    base = os.path.join(DATA_DIR, pair_id)
    cand1 = os.path.join(base, "chain_pdb_files", f"{a}.pdb")
    cand2 = os.path.join(base, "chain_pdb_files", f"{b}.pdb")
    pdb1 = cand1 if os.path.isfile(cand1) else os.path.join(base, f"{p1}.pdb")
    pdb2 = cand2 if os.path.isfile(cand2) else os.path.join(base, f"{p2}.pdb")
    return pdb1, c1, pdb2, c2

def _chain_len_from_pdb(pdb_path: str, chain_id: str | None) -> int:
    n = 0
    seen = set()
    try:
        with open(pdb_path) as fh:
            for line in fh:
                if not line.startswith("ATOM"): 
                    continue
                if line[12:16].strip() != "CA":
                    continue
                if chain_id and (line[21].strip() != chain_id):
                    continue
                resnum = (line[22:27], line[26])
                if resnum not in seen:
                    seen.add(resnum); n += 1
    except Exception:
        pass
    return n

def _pair_max_len(pair_id: str) -> int:
    pdb1, c1, pdb2, c2 = _truth_pdbs_for_pair(pair_id)
    return max(_chain_len_from_pdb(pdb1, c1), _chain_len_from_pdb(pdb2, c2))


def _pick_best_overall(dfall: pd.DataFrame, model_tag: str, which_fold: int) -> str:
    if dfall is None or dfall.empty:
        return "-"
    df = dfall[dfall["model"] == model_tag]
    col = "TM1" if which_fold == 1 else "TM2"
    s = pd.to_numeric(df[col], errors="coerce")
    if s.notna().any():
        return f"{s.max():.2f}"
    return "-"


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
    ctag = None
    if "cluster_num" in df.columns:
        ctag = df["cluster_num"].astype(str)
    elif "cluster" in df.columns:
        ctag = df["cluster"].astype(str)
    else:
        ctag = pd.Series([None] * len(df))
    df["_cluster_tag"] = ctag

    # Parse type/id
    def _ctype_id(s: Optional[str]) -> Tuple[Optional[str], Optional[str]]:
        if s is None or s != s:
            return None, None
        st = str(s)
        if "deep" in st.lower():
            return "deep", None
        m = SHALLOW_RE.search(st)
        if m:
            cid = m.group(1).lstrip("0") or "0"
            return "clust", cid
        # Sometimes numbers only:
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
    mode: 'clust' or 'deep'
    which_fold: 1 or 2 (TM1 vs TM2)
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

def _best_max(df: Optional[pd.DataFrame], col: str) -> str:
    if df is None or df.empty or col not in df.columns:
        return "-"
    s = pd.to_numeric(df[col], errors="coerce")
    if s.notna().any():
        return f"{s.max():.2f}"
    return "-"

def collect_summary_tables(
    DATA_DIR_in: str,
    out_united_csv: str,
    out_best_csv: str,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Writes:
      - united per-cluster table (kept broadly compatible)
      - the *new* best-per-pair summary with the requested columns
    """
    globals()["DATA_DIR"] = DATA_DIR_in

    try:
        # strings like '1xjuB_1xjtA'
        pair_names = list_protein_pairs(parsed=False, sort_result=True)
    except Exception:
        pair_names = []

    united_frames: List[pd.DataFrame] = []
    best_rows: List[dict] = []

    for pair_id in pair_names:
        pair_dir = os.path.join(DATA_DIR, pair_id)
        df_af   = _safe_read_csv(os.path.join(pair_dir, "Analysis", "df_af.csv"))
        df_esm  = _safe_read_csv(os.path.join(pair_dir, "Analysis", "df_esm.csv"))
        df_cmap = _safe_read_csv(os.path.join(pair_dir, "Analysis", "df_cmap.csv"))

        # ---------- Build the summary row ----------
        tm_af  = _norm_tm_df(df_af,  "af2")    # default if legacy lacks 'model'
        tm_esm = _norm_tm_df(df_esm, "esm2")   # default if legacy
        tm_all = pd.concat([tm_af, tm_esm], ignore_index=True) if len(tm_af) or len(tm_esm) else pd.DataFrame()

        row = {"pair_id": pair_id, "#RES": _pair_max_len(pair_id)}
        # AF keeps Clust & Deep
        for tag, up in (("af2","AF2"), ("af3","AF3")):
            row[f"{up}Clust_TM1"] = _pick_best(tm_all, tag, "clust", 1)
            row[f"{up}Clust_TM2"] = _pick_best(tm_all, tag, "clust", 2)
            row[f"{up}Deep_TM1"]  = _pick_best(tm_all, tag, "deep",  1)
            row[f"{up}Deep_TM2"]  = _pick_best(tm_all, tag, "deep",  2)

        # ESM is overall-best (no cluster split)    
        row["ESM2_TM1"] = _pick_best_overall(tm_all, "esm2", 1)
        row["ESM2_TM2"] = _pick_best_overall(tm_all, "esm2", 2)
        row["ESM3_TM1"] = _pick_best_overall(tm_all, "esm3", 1)
        row["ESM3_TM2"] = _pick_best_overall(tm_all, "esm3", 2)

        # CMAP maxima (MSA-Transformer)
        row["MSATrans_CMAP_PR1"] = _best_max(df_cmap, "t1_precision")
        row["MSATrans_CMAP_PR2"] = _best_max(df_cmap, "t2_precision")
        row["MSATrans_CMAP_RE1"] = _best_max(df_cmap, "t1_recall")
        row["MSATrans_CMAP_RE2"] = _best_max(df_cmap, "t2_recall")

        best_rows.append(row)

        # ---------- (Optional) united per-cluster table (kept for clusters page) ----------
        # Reuse AF (TM columns as 'af_tm_*') + CMAP and merge on (pair, cluster) when possible
        def _ensure_cluster(df: pd.DataFrame) -> pd.DataFrame:
            if "cluster" in df.columns:
                df["cluster"] = df["cluster"].astype(str)
                return df
            out = df.copy()
            if "cluster_num" in out.columns:
                out["cluster"] = out["cluster_num"].astype(str)
                return out
            for name in ("file","File","filename","Filename"):
                if name in out.columns:
                    m = out[name].astype(str).str.extract(r"ShallowMsa_(\d+)", expand=False)
                    if m.notna().any():
                        out["cluster"] = m.astype(str)
                        return out
            return out

        # CMAP slice
        if df_cmap is not None and not df_cmap.empty:
            cm = df_cmap.copy()
            cm["pair_id"] = pair_id
            cm = _ensure_cluster(cm)
            cm = cm[[c for c in cm.columns if c in ("pair_id","cluster","file") or c.startswith(("t1_","t2_","common_","uniq1_","uniq2_"))]]
        else:
            cm = None

        # AF slice (standardize af_tm_* for backward compat)
        if df_af is not None and not df_af.empty:
            af = df_af.copy()
            if "score_pdb1" in af.columns: af = af.rename(columns={"score_pdb1":"af_tm_fold1","score_pdb2":"af_tm_fold2"})
            if "TMscore_fold1" in af.columns: af = af.rename(columns={"TMscore_fold1":"af_tm_fold1","TMscore_fold2":"af_tm_fold2"})
            af["pair_id"] = pair_id
            af = _ensure_cluster(af)
            af = af[[c for c in ("pair_id","cluster","af_tm_fold1","af_tm_fold2") if c in af.columns]]
        else:
            af = None

        if cm is None and af is not None:
            united_frames.append(af)
        elif cm is not None and af is None:
            united_frames.append(cm)
        elif cm is not None and af is not None:
            join_keys = ["pair_id","cluster"] if ("cluster" in cm.columns and "cluster" in af.columns) else ["pair_id"]
            united_frames.append(pd.merge(cm, af, on=join_keys, how="left"))

    # Write united per-cluster
    if united_frames:
        df_united = pd.concat(united_frames, ignore_index=True)
        metric_cols = [c for c in df_united.columns if c.startswith(("t1_","t2_","common_","uniq1_","uniq2_","af_tm_"))]
        if metric_cols:
            df_united = df_united[df_united[metric_cols].notna().any(axis=1)].copy()
    else:
        df_united = pd.DataFrame(columns=["pair_id","cluster","file"])

    os.makedirs(os.path.dirname(out_united_csv), exist_ok=True)
    df_united.to_csv(out_united_csv, index=False)

    # Write best-per-pair
    df_best = pd.DataFrame(best_rows).sort_values("pair_id")
    os.makedirs(os.path.dirname(out_best_csv), exist_ok=True)
    df_best.to_csv(out_best_csv, index=False)

    print(f"[summary] united per-cluster → {out_united_csv}  (rows={len(df_united)})")
    print(f"[summary] best-per-pair   → {out_best_csv}    (rows={len(df_best)})")
    return df_united, df_best

if __name__=='__main__':
    collect_summary_tables(DATA_DIR, DETAILED_RESULTS_TABLE, SUMMARY_RESULTS_TABLE)
