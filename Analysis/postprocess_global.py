# Analysis/postprocess_global.py
from __future__ import annotations
import os, glob, re, sys
import pandas as pd
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, ROOT)

from utils.cache_utils import get_from_pair_cache
from utils.utils import numify
from utils.align_utils import seq_identity_for_pair
from config import *

def _add_seqid_column(df):
    seqids = []
    for pid in df["pair_id"]:
        try:
            seqids.append(seq_identity_for_pair(pid))
        except Exception as e:
            # keep going if MSA missing for some pairs
            print(f"[postprocess] WARN: seqid failed for {pid}: {e}")
            seqids.append(np.nan)
    df["SEQ_ID"] = seqids
    return df



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
        def _pick_mkf(model, kind, fold):
            sub = df[(df["model"]==model) & (df["kind"]==kind) & (df["fold"]==fold)]
            return numify(sub["TM"].iloc[0]) if not sub.empty else np.nan
        out["AF2Deep_TM1"]   = _pick_mkf("AF2", "Deep",      1)
        out["AF2Deep_TM2"]   = _pick_mkf("AF2", "Deep",      2)
        out["AF3Deep_TM1"]   = _pick_mkf("AF3", "Deep",      1)
        out["AF3Deep_TM2"]   = _pick_mkf("AF3", "Deep",      2)
        out["AF2Clust_TM1"]  = _pick_mkf("AF2", "ClustBest", 1)
        out["AF2Clust_TM2"]  = _pick_mkf("AF2", "ClustBest", 2)
        out["AF3Clust_TM1"]  = _pick_mkf("AF3", "ClustBest", 1)
        out["AF3Clust_TM2"]  = _pick_mkf("AF3", "ClustBest", 2)

    # ---- ESM2/ESM3 Best Prediction ----
    if os.path.isfile(f_esm):
        df = pd.read_csv(f_esm)
        def _pick_mf(model, fold):
            sub = df[(df["model"]==model) & (df["fold"]==fold)]
            return numify(sub["TM"].iloc[0]) if not sub.empty else np.nan
        out["ESM2_TM1"] = _pick_mf("ESM2", 1)
        out["ESM2_TM2"] = _pick_mf("ESM2", 2)
        out["ESM3_TM1"] = _pick_mf("ESM3", 1)
        out["ESM3_TM2"] = _pick_mf("ESM3", 2)

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
            def _pick_token(token):
                p = f"{token}.pdb"
                row = dfe[dfe["PDB_ID"] == p]
                return numify(row["Delta_G"].iloc[0]) if not row.empty else np.nan
            out["ΔG1"] = round(_pick_token(a), 2)
            out["ΔG2"] = round(_pick_token(b), 2)
    except Exception:
        pass

    return out

def build_or_load_global_tables(force: bool = False,
                                pairs_glob: str | None = None) -> tuple[str,str]:
    """
    Builds two CSVs:
      - SUMMARY_RESULTS_TABLE: one row per pair, with AF/ESM bests + PAIR_TM + ΔG + seq_id + msa stats
      - DETAILED_RESULTS_TABLE : your detailed per-pair row(s) from _read_pair_tables(d)
    Also merges cache.json (seq_id, msa_depth, msa_width) if present.
    """
    if (not force) and os.path.isfile(SUMMARY_RESULTS_TABLE) and os.path.isfile(DETAILED_RESULTS_TABLE):
        return SUMMARY_RESULTS_TABLE, DETAILED_RESULTS_TABLE

    pair_dirs = sorted(glob.glob(os.path.join(PIPELINE_DIR, "*_*")))
    if pairs_glob:
        pair_dirs = sorted(glob.glob(os.path.join(PIPELINE_DIR, pairs_glob)))

    rows = []
    cache_rows = []
    for d in pair_dirs:
        anal = os.path.join(d, "Analysis")
        if os.path.isdir(anal):
            rows.append(_read_pair_tables(d))
            # --- pull cache.json if exists ---
            cache_p = os.path.join(anal, "cache.json")
            if os.path.isfile(cache_p) and os.path.getsize(cache_p) > 0:
                try:
                    with open(cache_p, "r") as fh:
                        C = json.load(fh)
                    pair_id = os.path.basename(d)
                    cache_rows.append({
                        "pair_id": pair_id,
                        "seq_id": C.get("seq_id", None),
                        "msa_depth": C.get("msa_depth", None),
                        "msa_width": C.get("msa_width", None),
                    })
                except Exception:
                    pass

    os.makedirs(DOCS_DIR, exist_ok=True)
    if not rows:
        pd.DataFrame([]).to_csv(SUMMARY_RESULTS_TABLE, index=False)
        pd.DataFrame([]).to_csv(DETAILED_RESULTS_TABLE, index=False)
        return SUMMARY_RESULTS_TABLE, DETAILED_RESULTS_TABLE

    # base table assembled from per-pair readers
    df = pd.DataFrame(rows)

    # merge cache additions
    if cache_rows:
        df_cache = pd.DataFrame(cache_rows)
        df = df.merge(df_cache, on="pair_id", how="left")

    # order/keep columns if present
    keep_cols = [
        "pair_id",
        "AF2Deep_TM1","AF2Deep_TM2","AF3Deep_TM1","AF3Deep_TM2",
        "AF2Clust_TM1","AF2Clust_TM2","AF3Clust_TM1","AF3Clust_TM2",
        "ESM2_TM1","ESM2_TM2","ESM3_TM1","ESM3_TM2",
        "TrueTrue_TM", "#RES", "PAIR_TM", "ΔG1", "ΔG2",
        # NEW:
        "seq_id", "msa_depth", "msa_width",
    ]

    summary  = df[[c for c in keep_cols if c in df.columns]].copy()
    detailed = df.copy()

    # ... where you finish constructing df_summary / df_detailed:
    summary = _add_seqid_column(summary)
    # if you also want it in the detailed table:
    detailed = _add_seqid_column(detailed)

    # nicer formatting
    for _df in (summary, detailed):
        num_cols = _df.select_dtypes(include="number").columns
        _df[num_cols] = _df[num_cols].round(3)

    summary.to_csv(SUMMARY_RESULTS_TABLE, index=False)
    detailed.to_csv(DETAILED_RESULTS_TABLE, index=False)
    return SUMMARY_RESULTS_TABLE, DETAILED_RESULTS_TABLE

