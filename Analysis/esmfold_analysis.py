import os, sys
import pandas as pd
from pathlib import Path

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, ROOT)

from utils.align_utils import compute_tmscore_align
from utils.utils import list_protein_pairs
from config import DATA_DIR
from tqdm import tqdm

def compute_df_esm(fold_pair: str, write: bool = True) -> pd.DataFrame:
    """
    Build df_esm.csv for a single pair:
      - reads Pipeline/<pair>/output_esm_fold/{esm2,esm3}/
      - uses samples_index.tsv if present, else scans *.pdb
      - computes TM vs truth (chain-sliced if available)
      - returns a DataFrame and optionally writes Analysis/df_esm.csv
    """
    pair_dir = Path(DATA_DIR) / fold_pair
    (pair_dir / "Analysis").mkdir(exist_ok=True)
    folds = fold_pair.split("_")
    chains = [folds[0][-1], folds[1][-1]]

    # prefer chain-sliced truth
    pdb1_root, pdb2_root = folds[0][:-1], folds[1][:-1]
    cand1 = pair_dir / "chain_pdb_files" / f"{folds[0]}.pdb"
    cand2 = pair_dir / "chain_pdb_files" / f"{folds[1]}.pdb"
    truth1 = cand1 if cand1.is_file() else pair_dir / f"{pdb1_root}.pdb"
    truth2 = cand2 if cand2.is_file() else pair_dir / f"{pdb2_root}.pdb"

    rows = []
    for model_tag in ("esm2", "esm3"):
        mdir = pair_dir / "output_esm_fold" / model_tag
        if not mdir.is_dir():
            continue
        tsv = mdir / "samples_index.tsv"
        if tsv.is_file():
            idx = pd.read_csv(tsv, sep="\t")  # name,pdb_path
            it = idx[["name","pdb_path"]].itertuples(index=False, name=None)
        else:
            pdbs = sorted(p for p in mdir.glob("*.pdb"))
            it = []
            for p in pdbs:
                base = p.name
                name = base.replace(f"_{model_tag}.pdb","").replace(".pdb","")
                it.append((name, str(p)))

        for name, pred_pdb in it:
            # TM vs each truth (respect chain IDs of the *truth* structure)
            TM1 = compute_tmscore_align(str(truth1), pred_pdb, chain1=chains[0])
            TM2 = compute_tmscore_align(str(truth2), pred_pdb, chain1=chains[1])
            rows.append({
                "fold_pair": fold_pair,
                "model": model_tag,
                "name": name,
                "pdb_path": pred_pdb,
                "TMscore_fold1": TM1,
                "TMscore_fold2": TM2,
            })

    if not rows:
        # Return an empty DF with the expected columns so callers are robust
        cols = ["fold_pair","model","name","pdb_path","TMscore_fold1","TMscore_fold2",
                "cluster_num","TMdiff","TM_mean_cluster_pdb1","TM_mean_cluster_pdb2","sample_count","class"]
        return pd.DataFrame(columns=cols)

    df = pd.DataFrame(rows)
    df["cluster_num"] = df["name"].str.extract(r"ShallowMsa_(\d+)", expand=False).fillna("DeepMsa")
    df["TMdiff"] = df["TMscore_fold1"] - df["TMscore_fold2"]
    grp = ["fold_pair","model","cluster_num"]
    df["TM_mean_cluster_pdb1"] = df.groupby(grp)["TMscore_fold1"].transform("mean")
    df["TM_mean_cluster_pdb2"] = df.groupby(grp)["TMscore_fold2"].transform("mean")
    df["sample_count"]         = df.groupby(grp)["TMscore_fold1"].transform("count")

    # optional coarse class
    tau, delta = 0.50, 0.05
    df["class"] = "neither"
    df.loc[(df.TMscore_fold1 >= tau) & (df.TMscore_fold1 - df.TMscore_fold2 >= delta), "class"] = "toward_fold1"
    df.loc[(df.TMscore_fold2 >= tau) & (df.TMscore_fold2 - df.TMscore_fold1 >= delta), "class"] = "toward_fold2"
    df.loc[(df[["TMscore_fold1","TMscore_fold2"]].min(axis=1) >= tau) & (df["TMdiff"].abs() < delta), "class"] = "both"

    df.sort_values(by=["fold_pair","model","cluster_num","name"], inplace=True)

    if write:
        out_csv = pair_dir / "Analysis" / "df_esm.csv"
        df.to_csv(out_csv, index=False)

    return df

if __name__ == "__main__":
    # keep the old batch behavior but call the new API
    fold_pairs = list_protein_pairs()
    for fp in fold_pairs:
        # original script expected subdir format "<pdb1>_<pdb2>"
        fp_sub = f"{fp[0]}_{fp[1]}" if isinstance(fp, (list, tuple)) else fp
        out_dir = Path(DATA_DIR) / fp_sub / "output_esm_fold"
        if not out_dir.is_dir() or not any(out_dir.iterdir()):
            print(f"[skip] No ESMFold outputs for {fp_sub} at {out_dir}")
            continue
        compute_df_esm(fp_sub, write=True)
    print("Finish all ESMFold TMscore computations!")

