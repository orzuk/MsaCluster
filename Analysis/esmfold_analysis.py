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
      - scans Pipeline/<pair>/output_esm_fold/{esm2,esm3}/*.pdb
      - if samples_index.tsv exists AND has rows, it is used; otherwise glob the directory
      - computes TM vs truth (chain-sliced if available)
    """
    pair_dir = Path(DATA_DIR) / fold_pair
    (pair_dir / "Analysis").mkdir(exist_ok=True)

    # Parse truth pdbs (prefer chain-sliced, fall back to whole PDBs)
    left, right = fold_pair.split("_", 1)
    pdb1_root, ch1 = left[:-1], left[-1]
    pdb2_root, ch2 = right[:-1], right[-1]
    cand1 = pair_dir / "chain_pdb_files" / f"{left}.pdb"
    cand2 = pair_dir / "chain_pdb_files" / f"{right}.pdb"
    truth1 = cand1 if cand1.is_file() else pair_dir / f"{pdb1_root}.pdb"
    truth2 = cand2 if cand2.is_file() else pair_dir / f"{pdb2_root}.pdb"

    def _iter_preds(model_tag: str):
        """Yield (name, path) for all predicted PDBs of a given model."""
        mdir = pair_dir / "output_esm_fold" / model_tag
        if not mdir.is_dir():
            return
        tsv = mdir / "samples_index.tsv"
        used_index = False
        if tsv.is_file():
            try:
                idx = pd.read_csv(tsv, sep="\t")
                if {"name","pdb_path"}.issubset(idx.columns) and len(idx) > 0:
                    used_index = True
                    for name, pth in idx[["name","pdb_path"]].itertuples(index=False, name=None):
                        yield name, str(pth)
            except Exception:
                pass  # fall through to glob

        if not used_index:
            # Glob all per-sample PDBs that end with _<model>.pdb
            for p in sorted(mdir.glob(f"*_{model_tag}.pdb")):
                # Strip the "_esmX.pdb" suffix to recover the sample name
                base = p.name
                name = base[: -(len(model_tag) + 5)]  # remove f"_{model_tag}.pdb"
                # safety: handle plain ".pdb" (shouldn't happen with our writer)
                if name.endswith(".pdb"):
                    name = name[:-4]
                yield name, str(p)

    rows = []
    for model_tag in ("esm2", "esm3"):
        for name, pred_pdb in _iter_preds(model_tag):
            # TM vs each truth (chain1 refers to TRUTH chain we compare to)
            print("[esm_tm] ", fold_pair, model_tag, name)
            TM1 = compute_tmscore_align(str(truth1), pred_pdb, chain1=ch1)
            TM2 = compute_tmscore_align(str(truth2), pred_pdb, chain1=ch2)
            rows.append({
                "fold_pair": fold_pair,
                "model": model_tag,
                "name": name,
                "pdb_path": pred_pdb,
                "TMscore_fold1": TM1,
                "TMscore_fold2": TM2,
            })

    # If nothing found, return an empty, schema-compatible DF
    if not rows:
        cols = ["fold_pair","model","name","pdb_path","TMscore_fold1","TMscore_fold2",
                "cluster_num","TMdiff","TM_mean_cluster_pdb1","TM_mean_cluster_pdb2","sample_count","class"]
        return pd.DataFrame(columns=cols)

    df = pd.DataFrame(rows)

    # cluster id: ShallowMsa_### or DeepMsa
    cl = df["name"].str.extract(r"(ShallowMsa_(\d+)|DeepMsa)", expand=True)
    df["cluster_num"] = cl[1].fillna("DeepMsa")

    df["TMdiff"] = df["TMscore_fold1"] - df["TMscore_fold2"]
    grp = ["fold_pair","model","cluster_num"]
    df["TM_mean_cluster_pdb1"] = df.groupby(grp)["TMscore_fold1"].transform("mean")
    df["TM_mean_cluster_pdb2"] = df.groupby(grp)["TMscore_fold2"].transform("mean")
    df["sample_count"]         = df.groupby(grp)["TMscore_fold1"].transform("count")

    # Simple coarse class
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

