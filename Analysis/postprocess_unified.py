# postprocess_unified.py
import os, re, glob, sys, json
import pandas as pd
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple
import fnmatch
from tqdm import tqdm
import argparse

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, ROOT)

from config import DATA_DIR, SUMMARY_RESULTS_TABLE, DETAILED_RESULTS_TABLE, MSA_TRANS_MODEL_FILE
from utils.utils import list_protein_pairs, pair_str_to_tuple, pick
from utils.align_utils import compute_tmscore_align, get_or_compute_true_tm
from Analysis.cmap_analysis import compute_cmap_metrics_for_pair

import numpy as np
import pandas as pd


PAIR_DIR = Path(DATA_DIR)


# === NEW: fast builder from cached per-pair cluster CSVs to unified CSVs ===

def _norm_tag(s: str) -> str:
    if not isinstance(s, str): s = str(s)
    t = s.strip()
    if "deep" in t.lower(): return "DeepMsa"
    m = re.search(r"(\d+)$", t)
    if m: return f"ShallowMsa_{int(m.group(1)):03d}"
    if t.startswith("ShallowMsa_"): return t
    return t


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

# --- add near the top of the file, with the other small helpers ---
def _msa_width_a3m_columns(a3m_path: str) -> int:
    """
    Return the MSA alignment width (number of columns), counting gaps.
    We read the first sequence line after a header and strip lowercase insertions.
    """
    try:
        with open(a3m_path, "r") as fh:
            for line in fh:
                if line.startswith(">"):
                    # next non-header line is the reference sequence
                    seq = next(fh).strip()
                    # in A3M, lowercase = insertions; remove them but keep '-' gaps
                    core = "".join(ch for ch in seq if (ch == '-' or ch.isupper()))
                    return len(core)
    except Exception:
        pass
    return 0


def _count_shallow_clusters_fast(pair_id: str) -> int:
    """Count ShallowMsa_*.a3m files directly under output_msa_cluster/."""
    base = f"{DATA_DIR}/{pair_id}/output_msa_cluster"
    try:
        names = os.listdir(base)
    except Exception:
        return 0
    return sum(1 for name in names if fnmatch.fnmatch(name, "ShallowMsa_*.a3m"))

def _short_cluster_label(tag: str) -> str:
    t = str(tag)
    if t.lower().startswith("deep"):
        return "Deep"
    m = re.search(r"(\d+)$", t)
    return f"{int(m.group(1)):03d}" if m else t


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
      - READ cached per-pair Analysis CSVs – no recompute
      - COMPUTE summary (one row per pair)
      - COMPUTE detailed (ONE ROW PER CLUSTER; includes per-fold AF2/AF3 + ESM list columns)
      - Optionally write SUMMARY_RESULTS_TABLE and DETAILED_RESULTS_TABLE
    """
    if not pairs:
        pairs = [f"{a}_{b}" for a, b in list_protein_pairs()]

    all_detailed = []
    summary_rows = []

    num_pairs = len(pairs)
    ctr_pairs = 0
    for pair_id in pairs:
        pair_dir = _pair_dir(pair_id)
        if not pair_dir.is_dir():
            continue

#        print("Building tables from pair:", pair_id)

        anal = _ensure_pair_analysis(pair_id)
        df_af   = _safe_read_csv(str(anal / "df_af.csv"))
        df_esm  = _safe_read_csv(str(anal / "df_esm.csv"))
        df_cmap = _safe_read_csv(str(anal / "df_cmap.csv"))

        # === DETAILED: ONE ROW PER CLUSTER ===
        try:
            per_cluster = build_pair_cluster_table(pair_id)
#            print("Add cluster sub-table: ", per_cluster)
            if per_cluster is not None and not per_cluster.empty:
                # keep both names for downstream code
                per_cluster.insert(0, "pair_id", pair_id)
                all_detailed.append(per_cluster)
        except Exception as e:
            print(f"[unified] WARN per-cluster table for {pair_id}: {e}")

        # === SUMMARY (kept as you had) ===
        tm_af  = _norm_tm_df(df_af,  "af2")
        tm_esm = _norm_tm_df(df_esm, "esm2")
        tm_all = pd.concat([tm_af, tm_esm], ignore_index=True) if len(tm_af) or len(tm_esm) else pd.DataFrame()

        row = {"pair_id": pair_id, "#RES": _pair_max_len_from_truth(pair_id)}

        deepmsa_file = _deepmsa_a3m_path(pair_id)
        msa_depth = _count_a3m_sequences_fast(deepmsa_file)
        msa_width = _msa_width_a3m_columns(deepmsa_file)
        n_clusters = _count_shallow_clusters_fast(pair_id)
        row["MSA: DEPTH; #RES; #Clusters"] = f"{int(msa_depth)}; {int(msa_width)}; {int(n_clusters)}"
#        print("Set MSA: DEPTH; #RES; #Clusters:", row["MSA: DEPTH; #RES; #Clusters"])

        pdb1, c1, pdb2, c2 = _truth_pdbs(pair_id)
        pair_dir_str = str(_pair_dir(pair_id))
        print("Compute TrueTM: pair_id:", pair_id, " ; ", ctr_pairs,  " out of ", num_pairs)
        ctr_pairs += 1
        def _tm_sym_once(pa: str, pb: str) -> float:
            t12 = compute_tmscore_align(pa, pb, chain1=c1, chain2=c2)
            t21 = compute_tmscore_align(pb, pa, chain1=c2, chain2=c1)
            return max(t12, t21)
        row["PAIR_TM"] = round(get_or_compute_true_tm(pair_dir_str, pdb1, pdb2, _tm_sym_once), 3)

        # ΔG (pair-local energy CSV) – PDB_ID may be without chain (e.g. "2qqj.pdb")
        f_energy = f"{DATA_DIR}/{pair_id}/Analysis/df_energy_global.csv"
        if os.path.isfile(f_energy):
            try:
                dfe = pd.read_csv(f_energy)
                # Normalize column names (guard against different casing)
                colmap = {c.strip(): c for c in dfe.columns}
                pid_col  = colmap.get("PDB_ID", "PDB_ID")
                dg_col   = colmap.get("Delta_G", "Delta_G")
                pair_col = colmap.get("Pair", None)

                a, b = pair_id.split("_", 1)

                def _cands(tok: str):
                    tok = tok.strip()
                    c = [f"{tok}.pdb"]
                    if tok and tok[-1].isalpha():  # drop trailing chain, e.g. 2qqjA -> 2qqj
                        c.append(f"{tok[:-1]}.pdb")
                    return c

                def _pick(tok: str):
                    dfp = dfe if pair_col is None else dfe[dfe[pair_col] == pair_id]
                    for cand in _cands(tok):
                        sub = dfp[dfp[pid_col] == cand]
                        if not sub.empty:
                            try:
                                return float(sub[dg_col].iloc[0])
                            except Exception:
                                pass
                    return float("nan")

                row["ΔG1"] = round(_pick(a), 2)
                row["ΔG2"] = round(_pick(b), 2)
            except Exception:
                row["ΔG1"] = float("nan"); row["ΔG2"] = float("nan")
        else:
            row["ΔG1"] = float("nan"); row["ΔG2"] = float("nan")

        for tag, up in (("af2","AF2"), ("af3","AF3")):
            row[f"{up}Clust_TM1"] = _pick_best(tm_all, tag, "clust", 1)
            row[f"{up}Clust_TM2"] = _pick_best(tm_all, tag, "clust", 2)
            row[f"{up}Deep_TM1"]  = _pick_best(tm_all, tag, "deep",  1)
            row[f"{up}Deep_TM2"]  = _pick_best(tm_all, tag, "deep",  2)

        row["ESM2_TM1"] = _pick_best_overall(tm_all, "esm2", 1)
        row["ESM2_TM2"] = _pick_best_overall(tm_all, "esm2", 2)
        row["ESM3_TM1"] = _pick_best_overall(tm_all, "esm3", 1)
        row["ESM3_TM2"] = _pick_best_overall(tm_all, "esm3", 2)

        row["MSATrans_CMAP_PR1"] = _best_max(df_cmap, "t1_precision")
        row["MSATrans_CMAP_PR2"] = _best_max(df_cmap, "t2_precision")
        row["MSATrans_CMAP_RE1"] = _best_max(df_cmap, "t1_recall")
        row["MSATrans_CMAP_RE2"] = _best_max(df_cmap, "t2_recall")

        summary_rows.append(row)

#    print("ALL DETAILED list:")
#    print(all_detailed)
    detailed_df = pd.concat(all_detailed, ignore_index=True) if all_detailed else pd.DataFrame()
    summary_df  = pd.DataFrame(summary_rows)


    if write_out:
        print("WRITING SUMMARY CLUSTERS TABLE:")
        print(summary_df.head())
        Path(SUMMARY_RESULTS_TABLE).parent.mkdir(parents=True, exist_ok=True)
        summary_df.to_csv(SUMMARY_RESULTS_TABLE, index=False)

        print("WRITING DETAILED CLUSTERS TABLE:")
        print(detailed_df.head())
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
    has_cmap_outputs = (os.path.isdir(pred_dir) and
            any(f.endswith((".npz", ".npy")) for f in os.listdir(pred_dir)))
    if not has_cmap_outputs:
        print(f"[warn] No MSA-Transformer outputs found at {pred_dir}; skipping CMAP metrics.")
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
    """
    Return the per-pair analysis dir that contains the cached CSVs
    (df_af.csv, df_esm.csv, df_cmap.csv). We try several standard
    locations and prefer the first one that actually has any of them.
    """
    out = Path(DATA_DIR)  / pair_id / "Analysis"
    if out.is_dir():
        if any((out / fn).exists() for fn in ("df_af.csv", "df_esm.csv", "df_cmap.csv")):
            return out
    # Fallback (still lets us read cache.json for n/neff)
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
                "pair_id": pair_id,
                "model": ver,  # AF2 vs AF3
                "cluster_num": cluster,
                "name": name,  # short
                "TMscore_fold1": tm1,
                "TMscore_fold2": tm2
            })

    df = pd.DataFrame(rows)
    if len(df):
        df["TMdiff"] = df["TMscore_fold1"] - df["TMscore_fold2"]
        df = df[["pair_id", "model", "cluster_num", "name", "TMscore_fold1", "TMscore_fold2", "TMdiff"]]
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
                rows.append({"pair_id": pair_id, "model": model_tag, "cluster_num":
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
                rows.append({"pair_id": pair_id, "model": model_tag,
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
    cols = ["pair_id", "model", "cluster_num", "name", "TMscore_fold1", "TMscore_fold2", "TMdiff"]
    df = df[cols].sort_values(["model", "cluster_num", "name"])

    df.to_csv(out_csv, index=False)
    return df

def _read_cmap(pair_id: str) -> pd.DataFrame:
    """Read per-pair cmap metrics produced by cmap_analysis.py (don’t recompute here)."""
    csv = _pair_dir(pair_id) / "Analysis" / "df_cmap.csv"
    return pd.read_csv(csv) if csv.is_file() else pd.DataFrame()

def build_pair_cluster_table(pair_id: str) -> pd.DataFrame:
    """
    ONE row per cluster with columns (when available):
      cluster, n, neff,
      AF2_TM1, AF2_TM2, AF3_TM1, AF3_TM2,
      RE-MSAT-COM, RE-MSAT1, RE-MSAT2,
      ESM_SEQIDS, ESM_TM1_LIST, ESM_TM2_LIST
    """
    anal = _ensure_pair_analysis(pair_id)
    df_af   = _safe_read_csv(str(anal / "df_af.csv"))
    df_esm  = _safe_read_csv(str(anal / "df_esm.csv"))
    df_cmap = _safe_read_csv(str(anal / "df_cmap.csv"))

    # ---- cache clusters: sizes & neff ----
    cache_stats = {}
    cpath = anal / "cache.json"
    if cpath.is_file():
        try:
            import json
            C = json.loads(cpath.read_text())
            cache_stats = C.get("clusters", {}) or {}
        except Exception:
            cache_stats = {}

    # ---- cluster tag resolver (works with cluster_num/cluster/name) ----
    SHALLOW_RE = re.compile(r"ShallowMsa[_\-]?(\d+)")
    def _tag_series(d: Optional[pd.DataFrame]) -> pd.Series:
        if d is None or d.empty:
            return pd.Series([], dtype=object)
        cols = [str(c) for c in d.columns]
        if "cluster_num" in cols:
            base = d["cluster_num"].astype(str)
        elif "cluster" in cols:
            base = d["cluster"].astype(str)
        elif "name" in cols:
            base = d["name"].astype(str)
        else:
            base = pd.Series([None] * len(d))
        def _norm(x):
            if not isinstance(x, str):
                return None
            s = x
            if "DeepMsa" in s:
                return "DeepMsa"
            m = SHALLOW_RE.search(s)
            if m:
                return f"ShallowMsa_{int(m.group(1)):03d}"
            # last numeric suffix as fallback
            m = re.search(r"(\d+)$", s)
            return f"ShallowMsa_{int(m.group(1)):03d}" if m else None
        return base.map(_norm)

    # ---- AF per-fold maxima (per model, tolerant to 'model' column) ----
    def _tm_by_model(df: Optional[pd.DataFrame], model_pat: str, prefix: str) -> pd.DataFrame:
        if df is None or df.empty:
            return pd.DataFrame(columns=[f"{prefix}_TM1", f"{prefix}_TM2"])
        d = df.copy()
        d["_tag"] = _tag_series(d)
        if "model" in d.columns:
            mm = d["model"].astype(str).str.upper()
            d = d[mm.str.contains(model_pat, na=False)]  # e.g. "AF2" or "AF3"
        t1 = pd.to_numeric(d.get("TMscore_fold1"), errors="coerce")
        t2 = pd.to_numeric(d.get("TMscore_fold2"), errors="coerce")
        g = d.assign(TM1=t1, TM2=t2).groupby("_tag", dropna=True)[["TM1","TM2"]].max()
        return g.rename(columns={"TM1": f"{prefix}_TM1", "TM2": f"{prefix}_TM2"})

    af2 = _tm_by_model(df_af,  "AF2", "AF2")
    af3 = _tm_by_model(df_af,  "AF3", "AF3")

    # ---- CMAP / MSAT block ----
    ms = pd.DataFrame()
    if df_cmap is not None and not df_cmap.empty:
        d = df_cmap.copy()
        d["_tag"] = _tag_series(d)
        ms = pd.DataFrame({
            "RE-MSAT-COM": d.get("common_f1", d.get("common_mcc", None)),
            "RE-MSAT1":    d.get("t1_f1",    d.get("t1_mcc",    None)),
            "RE-MSAT2":    d.get("t2_f1",    d.get("t2_mcc",    None)),
            "_tag": d["_tag"]
        }).groupby("_tag", dropna=True).max()

    # ---- ESM list columns (IDs cleaned; 2 decimals; top-10 by TM1) ----
    esm_lists = pd.DataFrame()
    if df_esm is not None and not df_esm.empty:
        d = df_esm.copy()
        d["_tag"] = _tag_series(d)

        # choose ID column & clean
        id_col = "name" if "name" in d.columns else ("sample_id" if "sample_id" in d.columns else None)
        if id_col is None:
            d["ID"] = ["" for _ in range(len(d))]
        else:
            def _clean_id(x):
                s = str(x)
                if "__" in s:  # strip prefix
                    s = s.split("__")[-1]
                m = re.match(r"sample_(\d+)$", s)
                return f"{int(m.group(1)):03d}" if m else s
            d["ID"] = d[id_col].map(_clean_id)

        t1 = pd.to_numeric(d.get("TMscore_fold1"), errors="coerce").round(2)
        t2 = pd.to_numeric(d.get("TMscore_fold2"), errors="coerce").round(2)

        blocks = []
        for tag, sub in d.assign(TM1=t1, TM2=t2).groupby("_tag", dropna=True):
            sub = sub.sort_values("TM1", ascending=False).head(10)
            ids  = "; ".join([str(x) for x in sub["ID"].fillna("")])
            tm1s = "; ".join([("{:.2f}".format(v) if pd.notna(v) else "") for v in sub["TM1"]])
            tm2s = "; ".join([("{:.2f}".format(v) if pd.notna(v) else "") for v in sub["TM2"]])
            blocks.append(pd.DataFrame({"ESM_SEQIDS":[ids], "ESM_TM1_LIST":[tm1s], "ESM_TM2_LIST":[tm2s]}, index=[tag]))
        if blocks:
            esm_lists = pd.concat(blocks, axis=0)

    # ---- assemble per-cluster table ----
    idx = sorted(set(af2.index) | set(af3.index) | set(ms.index) | set(esm_lists.index) | set(cache_stats.keys()))
    out = pd.DataFrame(index=idx)
    for block in (af2, af3, ms, esm_lists):
        if not block.empty:
            out = out.join(block, how="left") if not out.empty else block.copy()

    # n & neff from cache (keys like "ShallowMsa_XXX" or "DeepMsa")
    def _get_stat(tag, key):
        s = cache_stats.get(tag)
        if s is None and tag != "DeepMsa":
            m = re.search(r"(\d+)$", str(tag))
            if m:
                s = cache_stats.get(f"ShallowMsa_{int(m.group(1)):03d}")
        return (s or {}).get(key, None)

    out["n"]    = [ _get_stat(t, "n")    for t in out.index ]
    out["neff"] = [ _get_stat(t, "neff") for t in out.index ]

    # finalize
    out = out.reset_index().rename(columns={"index": "cluster"})
    def _short(s):
        if s is None or (isinstance(s, float) and pd.isna(s)): return ""
        st = str(s).strip()
        if st == "" or "unknown" in st.lower(): return ""
        if st.lower().startswith("deep"): return "Deep"
        m = re.search(r"(\d+)$", st)
        return m.group(1) if m else st
#    print("Inside build-pair finished finalized, out cluster: ", out["cluster"])
    if 'cluster' in out:
#        print("cluster name before: ", out["cluster"])
        out["cluster"] = out["cluster"].map(_short)
#        print("cluster name after: ", out["cluster"])
    if '_tag' in out:
#        print("cluster-tag name before: ", out["_tag"])
        out["_tag"] = out["_tag"].map(_short)
#        print("cluster-tag name after: ", out["_tag"])
        out.rename(columns={"_tag":"cluster"}, inplace=True)  # change to cluster

    wanted = [
        "cluster","n","neff",
        "AF2_TM1","AF2_TM2","AF3_TM1","AF3_TM2",
        "RE-MSAT-COM","RE-MSAT1","RE-MSAT2",
        "ESM_SEQIDS","ESM_TM1_LIST","ESM_TM2_LIST",
    ]
    out = out[[c for c in wanted if c in out.columns]]
    num_cols = out.select_dtypes(include="number").columns

    if len(num_cols):
        out[num_cols] = out[num_cols].round(3)
    return out


def build_pair_cluster_table_one_row(
    pair_id: str,
    max_esm: int = 10,
) -> pd.DataFrame:
    """
    Return ONE row per cluster for `pair_id` with columns:
      pair_id, cluster, n, neff,
      AF2_TM1, AF2_TM2,
      AF3_TM1, AF3_TM2,
      CMAP_RE1, CMAP_RE2,
      ESM_SEQIDS, ESM_TM1_LIST, ESM_TM2_LIST
    Rules:
      - AF2/AF3: for each cluster & fold, take the MAX TMscore across predictions (names).
      - CMAP: use Recall per fold (RE1/RE2) if found; fallback to common names.
      - ESM: keep up to `max_esm` samples per cluster; lists are “; ”-joined in same order.
      - cluster label is "Deep" or a 3-digit index "000","001",...
      - round numeric columns to 3 decimals.
    """
    # Always use the cached per-pair Analysis directory where df_af/df_esm/df_cmap live
    anal = _ensure_pair_analysis(pair_id)

    df_af   = _safe_read_csv(anal / "df_af.csv")
    df_esm  = _safe_read_csv(anal / "df_esm.csv")
    df_cmap = _safe_read_csv(anal / "df_cmap.csv")

    # cache stats: n, neff
    n_map, neff_map = {}, {}
    cpath = anal / "cache.json"
    if cpath.is_file():
        try:
            C = json.loads(cpath.read_text())
            for k, v in (C.get("clusters", {}) or {}).items():
                n_map[_norm_tag(k)]    = v.get("n")
                neff_map[_norm_tag(k)] = v.get("neff")
        except Exception:
            pass

    # collect all cluster tags we see anywhere
    tags: set[str] = set()
    def add_tags_from_df(df):
        if df is None or df.empty: return
        d = df.copy()
        if "cluster_num" not in d.columns and "cluster" in d.columns:
            d = d.rename(columns={"cluster": "cluster_num"})
        if "cluster_num" in d.columns:
            tags.update(_norm_tag(str(x)) for x in d["cluster_num"].astype(str))
    add_tags_from_df(df_af); add_tags_from_df(df_esm); add_tags_from_df(df_cmap)
    tags.update(n_map.keys())
    if not tags:
        return pd.DataFrame(columns=["pair_id","cluster","n","neff",
                                     "AF2_TM1","AF2_TM2","AF3_TM1","AF3_TM2",
                                     "CMAP_RE1","CMAP_RE2",
                                     "ESM_SEQIDS","ESM_TM1_LIST","ESM_TM2_LIST"])

    # --- AF2/AF3 max per fold ---
    def _tm_max_by_model(df: Optional[pd.DataFrame], model_name: str) -> pd.DataFrame:
        if df is None or df.empty:
            return pd.DataFrame(columns=["_tag", f"{model_name}_TM1", f"{model_name}_TM2"])
        d = df.copy()
        if "cluster_num" not in d.columns and "cluster" in d.columns:
            d = d.rename(columns={"cluster": "cluster_num"})
        d["_tag"] = d["cluster_num"].astype(str).map(_norm_tag)
        if "model" in d.columns:
            d = d[d["model"].astype(str).str.upper() == model_name.upper()]
        # TM columns
        t1 = pd.to_numeric(d.get("TMscore_fold1"), errors="coerce")
        t2 = pd.to_numeric(d.get("TMscore_fold2"), errors="coerce")
        d = d.assign(TM1=t1, TM2=t2)
        g = d.groupby("_tag")[["TM1","TM2"]].max().rename(
            columns={"TM1": f"{model_name}_TM1", "TM2": f"{model_name}_TM2"}
        )
        return g

    af2 = _tm_max_by_model(df_af,  "AF2")
    af3 = _tm_max_by_model(df_af,  "AF3")

    # --- CMAP recall per fold (robust column finder) ---
    def _pick_first_col(d: pd.DataFrame, candidates: List[str]) -> Optional[pd.Series]:
        for c in candidates:
            if c in d.columns:
                return pd.to_numeric(d[c], errors="coerce")
        return None

    cmap = pd.DataFrame()
    if df_cmap is not None and not df_cmap.empty:
        dc = df_cmap.copy()
        if "cluster_num" not in dc.columns and "cluster" in dc.columns:
            dc = dc.rename(columns={"cluster": "cluster_num"})
        dc["_tag"] = dc["cluster_num"].astype(str).map(_norm_tag)
        re1 = _pick_first_col(dc, ["MSATrans_CMAP_RE1","RE1","t1_re","t1_recall","RE-MSAT1"])
        re2 = _pick_first_col(dc, ["MSATrans_CMAP_RE2","RE2","t2_re","t2_recall","RE-MSAT2"])
        cmap = pd.DataFrame({"_tag": dc["_tag"]})
        if re1 is not None: cmap["CMAP_RE1"] = re1
        if re2 is not None: cmap["CMAP_RE2"] = re2
        cmap = cmap.groupby("_tag").max()  # max across any duplicates

    # --- ESM: keep up to max_esm samples; 3 string columns ---
    esm = pd.DataFrame(columns=["ESM_SEQIDS","ESM_TM1_LIST","ESM_TM2_LIST"])
    if df_esm is not None and not df_esm.empty:
        de = df_esm.copy()
        if "cluster_num" not in de.columns and "cluster" in de.columns:
            de = de.rename(columns={"cluster": "cluster_num"})
        de["_tag"] = de["cluster_num"].astype(str).map(_norm_tag)
        # 'name' holds a per-sample id; keep as is
        t1 = pd.to_numeric(de.get("TMscore_fold1"), errors="coerce")
        t2 = pd.to_numeric(de.get("TMscore_fold2"), errors="coerce")
        de = de.assign(TM1=t1, TM2=t2)
        blocks = []
        for tag, sub in de.groupby("_tag"):
            # sort by TM1 desc, take top max_esm
            sub = sub.sort_values("TM1", ascending=False).head(max_esm)
            ids = [str(x) for x in sub.get("name", pd.Series([None]*len(sub))).fillna("").tolist()]
            tm1 = [("" if pd.isna(x) else f"{float(x):.3f}") for x in sub["TM1"]]
            tm2 = [("" if pd.isna(x) else f"{float(x):.3f}") for x in sub["TM2"]]
            row = pd.DataFrame(
                {"ESM_SEQIDS": ["; ".join(ids)],
                 "ESM_TM1_LIST": ["; ".join(tm1)],
                 "ESM_TM2_LIST": ["; ".join(tm2)]},
                index=[tag]
            )
            blocks.append(row)
        if blocks:
            esm = pd.concat(blocks, axis=0)

    # --- assemble per-cluster table ---
    idx = sorted(set(tags))
    out = pd.DataFrame(index=idx)
    for block in (af2, af3, cmap, esm):
        if block is not None and not block.empty:
            out = out.join(block, how="left")

    out["n"]    = [n_map.get(t)    for t in out.index]
    out["neff"] = [neff_map.get(t) for t in out.index]

    out = out.reset_index().rename(columns={"index":"_tag"})
    out.insert(0, "pair_id", pair_id)
    out.insert(1, "cluster", out["_tag"].map(_short_cluster_label))
    out = out.drop(columns=["_tag"])

    # order and round
    preferred = ["pair_id","cluster","n","neff",
                 "AF2_TM1","AF2_TM2","AF3_TM1","AF3_TM2",
                 "CMAP_RE1","CMAP_RE2",
                 "ESM_SEQIDS","ESM_TM1_LIST","ESM_TM2_LIST"]
    cols = [c for c in preferred if c in out.columns] + [c for c in out.columns if c not in preferred]
    out = out[cols]
    num_cols = out.select_dtypes(include="number").columns
    if len(num_cols):
        out[num_cols] = out[num_cols].round(3)
    return out


def build_all_pairs_clusters_table(
    pairs_list_csv: str = DETAILED_RESULTS_TABLE,
    out_csv: str = SUMMARY_RESULTS_TABLE,
    max_esm: int = 10,
) -> str:
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    # get list of pairs robustly (from the summary CSV)
    if not os.path.exists(pairs_list_csv):
        raise FileNotFoundError(f"Cannot find pairs list at {pairs_list_csv}")
    df_pairs = pd.read_csv(pairs_list_csv)
    pair_col = "fold_pair" if "fold_pair" in df_pairs.columns else ("pair_id" if "pair_id" in df_pairs.columns else None)
    if pair_col is None:
        raise KeyError("Expected 'fold_pair' (or 'pair_id') in pairs list CSV.")
    frames = []
    for pid in df_pairs[pair_col].dropna().astype(str).unique():
        try:
            frames.append(build_pair_cluster_table_one_row(pid, max_esm=max_esm))
        except Exception as e:
            print(f"[clusters] WARN: {pid}: {e}")
    if not frames:
        pd.DataFrame(columns=["pair_id","cluster"]).to_csv(out_csv, index=False)
        return out_csv
    out = pd.concat(frames, axis=0, ignore_index=True)
    out.to_csv(out_csv, index=False)
    print(f"[clusters] wrote: {out_csv} (rows={len(out)})")
    return out_csv


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

        print("Post-processing pair:", pair_id)
        # AF/ESM tables
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
            det["pair_id"] = pair_id
            all_detailed.append(det)

        max_len = _pair_max_len_from_truth(pair_id)
        # Summary row (pair-level)
        summary_rows.append({
            "pair_id": pair_id,
            "#RES": max_len,
            **(best_tm or {}),
            **(best_cmap or {}),
            "n_af_preds": int(len(df_af)),
            "n_esm_preds": int(len(df_esm)),
            "n_cmap_preds": int(len(df_cmap))
        })

    detailed_df = pd.concat(all_detailed, ignore_index=True) if all_detailed else pd.DataFrame()
    summary_df  = pd.DataFrame(summary_rows)

    return summary_df, detailed_df


if __name__ == "__main__":
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

# 3) Build global per-cluster table (one row per cluster)
    build_all_pairs_clusters_table(pairs=args.pairs, write_out=True)
    print(f"[postprocess] clusters table written:\n  {Path(SUMMARY_RESULTS_TABLE).parent / 'clusters_global_table.csv'}")