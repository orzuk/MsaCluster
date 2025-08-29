import os, sys, re
import pandas as pd
from typing import Optional, Tuple, List

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, ROOT)
from config import *

# Reuse your existing PAIR_DIR_RE and list_protein_pairs() already in this file.
# PAIR_DIR_RE = re.compile(r'^[0-9A-Za-z]{5}_[0-9A-Za-z]{5}$')
SHALLOW_RE = re.compile(r"ShallowMsa_(\d+)")

def _extract_cluster_from_string(s: str) -> Optional[str]:
    if not isinstance(s, str):
        return None
    stem = os.path.splitext(os.path.basename(s))[0]
    m = SHALLOW_RE.search(stem)
    return m.group(1) if m else None

def _ensure_cluster_column(df: pd.DataFrame) -> pd.DataFrame:
    """
    Ensure a 'cluster' column if possible:
      1) keep existing 'cluster'
      2) else derive from 'file'/'File'/'filename'/'Filename'
      3) else scan any object column for 'ShallowMsa_###'
      4) else leave missing (merge will fall back to pair-only)
    """
    if "cluster" in df.columns:
        df["cluster"] = df["cluster"].astype(str)
        return df

    for name in ("file", "File", "filename", "Filename"):
        if name in df.columns:
            df["cluster"] = df[name].map(_extract_cluster_from_string)
            if df["cluster"].notna().any():
                df["cluster"] = df["cluster"].astype(str)
                return df
            df.drop(columns=["cluster"], inplace=True)

    for col in df.columns:
        if df[col].dtype == object:
            sample = df[col].dropna().astype(str).head(200)
            if sample.map(_extract_cluster_from_string).notna().any():
                df["cluster"] = df[col].astype(str).map(_extract_cluster_from_string).astype(str)
                return df

    return df

def _format_best(value, cluster) -> str:
    """
    '0.53 (6)' if cluster present, otherwise just '0.53'.
    '-' if value missing.
    """
    if value is None or pd.isna(value):
        return "-"
    s = f"{float(value):.2f}"
    if cluster is None or pd.isna(cluster):
        return s
    # numeric clusters: show without leading zeros
    try:
        cid = str(int(str(cluster)))
    except Exception:
        cid = str(cluster)
    return f"{s} ({cid})" if cid else s

def _safe_read_csv(path: str) -> Optional[pd.DataFrame]:
    if not os.path.exists(path):
        return None
    try:
        return pd.read_csv(path)
    except Exception:
        return None

def collect_summary_tables(
    DATA_DIR: str,
    out_united_csv: str,
    out_best_csv: str,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build:
      1) united per-cluster table across all pairs (CMAP + AF) → out_united_csv
      2) best-per-pair table (one row per pair; BEST of fold-1 metrics + AF TM1) → out_best_csv

    Robustness:
      - If a pair has no data at all → still produces a best-row with '-' placeholders.
      - If cluster is missing → merge falls back to pair-only; best cells show numbers without (CID).
      - Clusters with no data are skipped from the united table.

    ESM results are intentionally ignored (add later).
    """
    # Make sure list_protein_pairs() sees the same DATA_DIR
    # (it uses the module-global DATA_DIR inside its body)
#    global DATA_DIR as _DATA_DIR_GLOBAL  # type: ignore  # just to be explicit to readers
    try:
        old_data_dir = DATA_DIR
    except Exception:
        old_data_dir = None
    # If your file already defines DATA_DIR globally, this assignment keeps them in sync:
    globals()["DATA_DIR"] = DATA_DIR  # so list_protein_pairs() uses the provided path

    # --- discover pairs using your helper ---
    try:
        pair_names = list_protein_pairs(parsed=False, sort_result=True)  # e.g. ['1dzlA_5keqF', ...]
    except Exception:
        # Fallback: no helper available → empty
        pair_names = []

    united_frames: List[pd.DataFrame] = []
    best_rows: List[dict] = []

    # Metric columns we’ll report best-of (fold-1 only)
    cmap_fold1_metrics = ["t1_precision", "t1_recall", "t1_f1", "t1_jaccard", "t1_mcc"]

    # If no pairs found, still write empty CSVs and return
    if not pair_names:
        os.makedirs(os.path.dirname(out_united_csv), exist_ok=True)
        pd.DataFrame().to_csv(out_united_csv, index=False)
        pd.DataFrame({"pair_id": [], "BEST_AF_TM_FOLD1": []}).to_csv(out_best_csv, index=False)
        print(f"[summary] No pairs under {DATA_DIR}. Wrote empty CSVs.")
        return pd.DataFrame(), pd.DataFrame()

    for pair_id in pair_names:
        pair_dir = os.path.join(DATA_DIR, pair_id)
        cmap_csv = os.path.join(pair_dir, "Analysis", "df_cmap.csv")
        af_csv   = os.path.join(pair_dir, "Analysis", "df_af.csv")

        df_cmap = _safe_read_csv(cmap_csv)
        df_af   = _safe_read_csv(af_csv)

        # ---- normalize CMAP (per-cluster) ----
        if df_cmap is not None and len(df_cmap):
            df_cmap["pair_id"] = pair_id
            df_cmap = _ensure_cluster_column(df_cmap)
            keep_cmap = ["pair_id"]
            if "cluster" in df_cmap.columns: keep_cmap.append("cluster")
            if "file"    in df_cmap.columns: keep_cmap.append("file")
            keep_cmap += [c for c in df_cmap.columns if c.startswith(("t1_", "t2_", "common_", "uniq1_", "uniq2_"))]
            df_cmap = df_cmap[[c for c in keep_cmap if c in df_cmap.columns]].copy()
        else:
            df_cmap = None  # no cmap for this pair

        # ---- normalize AF (optional) ----
        if df_af is not None and len(df_af):
            # common pipeline: columns 'cluster_num', 'score_pdb1', 'score_pdb2'
            if "cluster_num" in df_af.columns and "cluster" not in df_af.columns:
                df_af = df_af[df_af["cluster_num"] != "Query"]
                df_af["cluster"] = df_af["cluster_num"]
            # rename to standard
            colmap = {}
            if "score_pdb1" in df_af.columns: colmap["score_pdb1"] = "af_tm_fold1"
            if "score_pdb2" in df_af.columns: colmap["score_pdb2"] = "af_tm_fold2"
            df_af = df_af.rename(columns=colmap)
            df_af["pair_id"] = pair_id
            df_af = _ensure_cluster_column(df_af)
            keep_af = ["pair_id"]
            if "cluster" in df_af.columns: keep_af.append("cluster")
            if "af_tm_fold1" in df_af.columns: keep_af.append("af_tm_fold1")
            if "af_tm_fold2" in df_af.columns: keep_af.append("af_tm_fold2")
            df_af = df_af[[c for c in keep_af if c in df_af.columns]].copy()
        else:
            df_af = None  # no AF for this pair

        # ---- build united per-cluster rows ----
        if df_cmap is None and df_af is None:
            # nothing to add for united table (but best table will still get a placeholder row)
            pass
        elif df_cmap is None:
            # AF-only rows (may or may not have cluster)
            united_frames.append(df_af.copy())
        elif df_af is None:
            # CMAP-only rows
            united_frames.append(df_cmap.copy())
        else:
            # Both exist → merge (prefer (pair, cluster); fallback to pair-only)
            if "cluster" in df_cmap.columns and "cluster" in df_af.columns:
                join_keys = ["pair_id", "cluster"]
            else:
                join_keys = ["pair_id"]
            df_u = pd.merge(df_cmap, df_af, on=join_keys, how="left")
            united_frames.append(df_u)

        # ---- best-per-pair row ----
        row = {"pair_id": pair_id}

        g = None
        if united_frames:
            # collect just for this pair from the latest appended chunk(s)
            # safer: filter later from df_united; for now, use direct reads:
            # build a temporary group 'g' from whatever we have
            if df_cmap is not None and df_af is not None:
                g = df_u  # merged
            elif df_cmap is not None:
                g = df_cmap
            elif df_af is not None:
                g = df_af
            else:
                g = pd.DataFrame(columns=["pair_id"])

        # AF best fold-1
        if g is not None and "af_tm_fold1" in g.columns and len(g):
            idx = g["af_tm_fold1"].astype(float).idxmax(skipna=True)
            if pd.notna(idx):
                val = g.loc[idx, "af_tm_fold1"]
                cid = g.loc[idx, "cluster"] if "cluster" in g.columns else None
                row["BEST_AF_TM_FOLD1"] = _format_best(val, cid)
            else:
                row["BEST_AF_TM_FOLD1"] = "-"
        else:
            row["BEST_AF_TM_FOLD1"] = "-"

        # CMAP fold-1 metrics
        for m in cmap_fold1_metrics:
            if g is not None and m in g.columns and len(g):
                idx = g[m].astype(float).idxmax(skipna=True)
                if pd.notna(idx):
                    val = g.loc[idx, m]
                    cid = g.loc[idx, "cluster"] if "cluster" in g.columns else None
                    row[f"BEST_CMAP_{m.upper()}"] = _format_best(val, cid)
                else:
                    row[f"BEST_CMAP_{m.upper()}"] = "-"
            else:
                row[f"BEST_CMAP_{m.upper()}"] = "-"

        best_rows.append(row)

    # ----- finalize united table -----
    if united_frames:
        df_united = pd.concat(united_frames, ignore_index=True)
        # drop rows that have neither CMAP nor AF metrics (purely empty clusters)
        metric_cols = [c for c in df_united.columns if c.startswith(("t1_","t2_","common_","uniq1_","uniq2_","af_tm_"))]
        if metric_cols:
            mask_nonempty = df_united[metric_cols].notna().any(axis=1)
            df_united = df_united[mask_nonempty].copy()
    else:
        df_united = pd.DataFrame(columns=["pair_id","cluster","file"])

    os.makedirs(os.path.dirname(out_united_csv), exist_ok=True)
    df_united.to_csv(out_united_csv, index=False)

    # ----- best-per-pair table -----
    df_best = pd.DataFrame(best_rows).sort_values("pair_id")
    os.makedirs(os.path.dirname(out_best_csv), exist_ok=True)
    df_best.to_csv(out_best_csv, index=False)

    print(f"[summary] united per-cluster → {out_united_csv}  (rows={len(df_united)})")
    print(f"[summary] best-per-pair   → {out_best_csv}    (rows={len(df_best)})")

    return df_united, df_best


if __name__=='__main__':
    print("Run make summary table")
    collect_summary_tables(DATA_DIR, DETAILED_RESULTS_TABLE, SUMMARY_RESULTS_TABLE)


