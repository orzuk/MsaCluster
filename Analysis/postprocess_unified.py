# postprocess_unified.py
import os, re, glob
import pandas as pd
from pathlib import Path
from typing import Optional, Dict, Any, List
from tqdm import tqdm
import argparse

from config import DATA_DIR, SUMMARY_RESULTS_TABLE, DETAILED_RESULTS_TABLE, MSA_TRANS_MODEL_FILE

from utils.utils import list_protein_pairs, pair_str_to_tuple
from utils.align_utils import compute_tmscore_align
from Analysis.cmap_analysis import compute_cmap_metrics_for_pair

PAIR_DIR = Path(DATA_DIR)


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

        # Summary row (pair-level)
        summary_rows.append({
            "fold_pair": pair_id,
            **(best_tm or {}),
            **(best_cmap or {}),
            "n_af_preds": int(len(df_af)),
            "n_esm_preds": int(len(df_esm)),
            "n_cmap_preds": int(len(df_cmap))
        })

    detailed_df = pd.concat(all_detailed, ignore_index=True) if all_detailed else pd.DataFrame()
    summary_df  = pd.DataFrame(summary_rows)

    # Write global tables
    Path(SUMMARY_RESULTS_TABLE).parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(SUMMARY_RESULTS_TABLE, index=False)
    detailed_df.to_csv(DETAILED_RESULTS_TABLE, index=False)
    return summary_df, detailed_df


if __name__ == "__main__":

    p = argparse.ArgumentParser(description="Unified post-processing: TM-scores (AF/ESM), CMAP metrics, summary tables.")
    p.add_argument("--pairs", nargs="*", help="Pair IDs like 1fzpD_2frhA. If omitted, process all pairs.")
    p.add_argument("--force_rerun", action="store_true", help="Recompute per-pair CSVs even if they exist.")
    args = p.parse_args()

    summary_df, detailed_df = post_processing_analysis(force_rerun=args.force_rerun, pairs=args.pairs)
    print(f"[postprocess] wrote:\n  {SUMMARY_RESULTS_TABLE}\n  {DETAILED_RESULTS_TABLE}")
    print(f"[postprocess] summary rows={len(summary_df)} | detailed rows={len(detailed_df)}")
