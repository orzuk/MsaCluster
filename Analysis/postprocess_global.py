# Analysis/postprocess_global.py
from __future__ import annotations
import os, glob, re, sys
import pandas as pd
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, ROOT)

from utils.cache_utils import get_from_pair_cache
from config import PIPELINE_DIR, DOCS_DIR

SUMMARY_CSV  = os.path.join(DOCS_DIR, "summary_final_res_all_pairs_df.csv")
DETAIL_CSV   = os.path.join(DOCS_DIR, "detailed_final_res_all_pairs_df.csv")

_NUM = re.compile(r"^\s*([+-]?\d+(?:\.\d+)?)")
def _numify(x):
    if pd.isna(x): return np.nan
    if isinstance(x, (int,float)): return float(x)
    m = _NUM.match(str(x));  return float(m.group(1)) if m else np.nan

def _pair_id_from_dir(d: str) -> str:
    return os.path.basename(d.rstrip("/"))

def _read_pair_tables(pair_dir: str) -> dict:
    pair_id = _pair_id_from_dir(pair_dir)
    out = {"pair_id": pair_id}

    f_af   = os.path.join(pair_dir, "Analysis", "df_af.csv")
    f_esm  = os.path.join(pair_dir, "Analysis", "df_esm.csv")
    # OPTIONAL: f_cmap = os.path.join(pair_dir, "Analysis", "df_cmap.csv")

    # ---- AF2/AF3 Deep & Best Cluster ----
    if os.path.isfile(f_af):
        df = pd.read_csv(f_af)
        def _pick(model, kind, fold):
            sub = df[(df["model"]==model) & (df["kind"]==kind) & (df["fold"]==fold)]
            return _numify(sub["TM"].iloc[0]) if not sub.empty else np.nan
        out["AF2Deep_TM1"]   = _pick("AF2", "Deep",      1)
        out["AF2Deep_TM2"]   = _pick("AF2", "Deep",      2)
        out["AF3Deep_TM1"]   = _pick("AF3", "Deep",      1)
        out["AF3Deep_TM2"]   = _pick("AF3", "Deep",      2)
        out["AF2Clust_TM1"]  = _pick("AF2", "ClustBest", 1)
        out["AF2Clust_TM2"]  = _pick("AF2", "ClustBest", 2)
        out["AF3Clust_TM1"]  = _pick("AF3", "ClustBest", 1)
        out["AF3Clust_TM2"]  = _pick("AF3", "ClustBest", 2)

    # ---- ESM2/ESM3 Best Prediction ----
    if os.path.isfile(f_esm):
        df = pd.read_csv(f_esm)
        def _pick(model, fold):
            sub = df[(df["model"]==model) & (df["fold"]==fold)]
            return _numify(sub["TM"].iloc[0]) if not sub.empty else np.nan
        out["ESM2_TM1"] = _pick("ESM2", 1)
        out["ESM2_TM2"] = _pick("ESM2", 2)
        out["ESM3_TM1"] = _pick("ESM3", 1)
        out["ESM3_TM2"] = _pick("ESM3", 2)

    # ---- true vs true TM from cache ----
    tm_true = get_from_pair_cache(pair_dir, "true_true_tm")
    out["TrueTrue_TM"] = float(tm_true) if tm_true is not None else np.nan

    # ---- ΔG from pair-local energy CSV ----
    f_energy = os.path.join(pair_dir, "Analysis", "df_energy_global.csv")
    try:
        if os.path.isfile(f_energy):
            dfe = pd.read_csv(f_energy)
            # Identify which PDB belongs to fold1/fold2 by name tokens in pair_id
            a, b = pair_id.split("_", 1)  # "2qqjA", "4qdsA"
            # Match by basename (without .pdb)
            names = {os.path.splitext(os.path.basename(x))[0]: x for x in dfe["PDB_ID"].tolist()}
            def _pick(token):
                p = f"{token}.pdb"
                row = dfe[dfe["PDB_ID"] == p]
                return _numify(row["Delta_G"].iloc[0]) if not row.empty else np.nan
            out["ΔG1"] = _pick(a)
            out["ΔG2"] = _pick(b)
    except Exception:
        pass

    return out

def build_or_load_global_tables(force: bool = False,
                                pairs_glob: str | None = None) -> tuple[str,str]:
    if (not force) and os.path.isfile(SUMMARY_CSV) and os.path.isfile(DETAIL_CSV):
        return SUMMARY_CSV, DETAIL_CSV

    pair_dirs = sorted(glob.glob(os.path.join(PIPELINE_DIR, "*_*")))
    if pairs_glob:
        pair_dirs = sorted(glob.glob(os.path.join(PIPELINE_DIR, pairs_glob)))

    rows = []
    for d in pair_dirs:
        if os.path.isdir(os.path.join(d, "Analysis")):
            rows.append(_read_pair_tables(d))

    os.makedirs(DOCS_DIR, exist_ok=True)
    if not rows:
        pd.DataFrame([]).to_csv(SUMMARY_CSV, index=False)
        pd.DataFrame([]).to_csv(DETAIL_CSV, index=False)
        return SUMMARY_CSV, DETAIL_CSV

    df = pd.DataFrame(rows)

    keep_cols = ["pair_id",
                 "AF2Deep_TM1","AF2Deep_TM2","AF3Deep_TM1","AF3Deep_TM2",
                 "AF2Clust_TM1","AF2Clust_TM2","AF3Clust_TM1","AF3Clust_TM2",
                 "ESM2_TM1","ESM2_TM2","ESM3_TM1","ESM3_TM2",
                 "TrueTrue_TM",  '#RES', 'PAIR_TM',  'ΔG1', 'ΔG2' ]

    summary  = df[[c for c in keep_cols if c in df.columns]].copy()
    detailed = df.copy()

    summary.to_csv(SUMMARY_CSV, index=False)
    detailed.to_csv(DETAIL_CSV, index=False)
    return SUMMARY_CSV, DETAIL_CSV
