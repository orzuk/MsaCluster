import pandas as pd
import os, sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, ROOT)

from utils.align_utils import *
from utils.utils import list_protein_pairs
from config import *
from tqdm import tqdm

'''
To run this script on your local computer you have to specify where the TM align software is installed on your computer and adapt the function
compute_tmscore_align changing your path.
And to set in pipeline_folder your path where the results are on your machine
'''

# Compute TMscores of AF predictions for all clusters in a fold pair
def compute_esmfold_pred_tmscores(fold_pair):
    final_res       = []
    print(f'Fold_pair: {fold_pair}')
    path = f'{DATA_DIR}/{fold_pair}'
    print(f'Path: {path}')
    if not os.path.exists(f'{path}/Analysis'):
        print("Mkdir: " + f'{path}/Analysis')
        os.mkdir(f'{path}/Analysis')

    folds = fold_pair.split("_")
    chains = [folds[0][-1],folds[1][-1]]

    # Prefer chain-sliced PDBs if present; else fall back to root PDBs
    pdb1, pdb2 = folds[0][:-1], folds[1][:-1]
    cand1 = f"{DATA_DIR}/{fold_pair}/chain_pdb_files/{folds[0]}.pdb"
    cand2 = f"{DATA_DIR}/{fold_pair}/chain_pdb_files/{folds[1]}.pdb"
    truth1_pdb = cand1 if os.path.isfile(cand1) else f"{DATA_DIR}/{fold_pair}/{pdb1}.pdb"
    truth2_pdb = cand2 if os.path.isfile(cand2) else f"{DATA_DIR}/{fold_pair}/{pdb2}.pdb"

    # Iterate normalized per-model outputs via samples_index.tsv (preferred)
    for model_tag in ("esm2", "esm3"):
        mdir = f"{DATA_DIR}/{fold_pair}/output_esm_fold/{model_tag}"
        tsv = os.path.join(mdir, "samples_index.tsv")

        if os.path.isfile(tsv):
            idx = pd.read_csv(tsv, sep="\t")  # columns: name, pdb_path
            for _, r in idx.iterrows():
                pred_name = r["name"]
                pred_pdb  = r["pdb_path"]
                TMscore_fold1 = compute_tmscore_align(truth1_pdb, pred_pdb)
                TMscore_fold2 = compute_tmscore_align(truth2_pdb, pred_pdb)
                final_res.append({
                    "fold_pair": fold_pair,
                    "model": model_tag,
                    "name": pred_name,
                    "pdb_path": pred_pdb,
                    "TMscore_fold1": TMscore_fold1,
                    "TMscore_fold2": TMscore_fold2,
                })
        else:
            # Fallback: scan *.pdb if TSV missing
            if not os.path.isdir(mdir):
                continue
            pdb_files = sorted([p for p in os.listdir(mdir) if p.endswith(".pdb")])
            for pred_file in pdb_files:
                pred_path = os.path.join(mdir, pred_file)
                # infer a clean sample name (works with ShallowMsa_003__sample_007_esm2.pdb)
                base = os.path.basename(pred_file)
                pred_name = base.replace(f"_{model_tag}.pdb", "").replace(".pdb", "")
                TMscore_fold1 = compute_tmscore_align(truth1_pdb, pred_path)
                TMscore_fold2 = compute_tmscore_align(truth2_pdb, pred_path)
                final_res.append({
                    "fold_pair": fold_pair,
                    "model": model_tag,
                    "name": pred_name,
                    "pdb_path": pred_path,
                    "TMscore_fold1": TMscore_fold1,
                    "TMscore_fold2": TMscore_fold2,
                })

    if not final_res:
        print("[warn] No ESM predictions found for", fold_pair)
        return

    df = pd.DataFrame(final_res)
    print("df columns:", df.columns.tolist())

    # Cluster id from the standardized sample name; DeepMsa if present
    df["cluster_num"] = df["name"].str.extract(r"ShallowMsa_(\d+)", expand=False).fillna("DeepMsa")

    # Helpful differential signal and per-model, per-cluster aggregates
    df["TMdiff"] = df["TMscore_fold1"] - df["TMscore_fold2"]
    grp = ["fold_pair", "model", "cluster_num"]
    df["TM_mean_cluster_pdb1"] = df.groupby(grp)["TMscore_fold1"].transform("mean")
    df["TM_mean_cluster_pdb2"] = df.groupby(grp)["TMscore_fold2"].transform("mean")
    df["sample_count"]         = df.groupby(grp)["TMscore_fold1"].transform("count")

    # Optional: quick classification
    tau, delta = 0.50, 0.05
    df["class"] = "neither"
    df.loc[(df.TMscore_fold1 >= tau) & (df.TMscore_fold1 - df.TMscore_fold2 >= delta), "class"] = "toward_fold1"
    df.loc[(df.TMscore_fold2 >= tau) & (df.TMscore_fold2 - df.TMscore_fold1 >= delta), "class"] = "toward_fold2"
    df.loc[(df[["TMscore_fold1","TMscore_fold2"]].min(axis=1) >= tau) & (df["TMdiff"].abs() < delta), "class"] = "both"

    df.sort_values(by=["fold_pair","model","cluster_num","name"], inplace=True)
    df.to_csv(f"{path}/Analysis/df_af.csv", index=False)


if __name__=='__main__':
    fold_pairs = list_protein_pairs()
    print("All pairs to compute ESMFold TMscores: ", fold_pairs)

    for fold_pair in fold_pairs:
        fold_pair_subdir = fold_pair[0] + "_" + fold_pair[1]
        if 'sh' in fold_pair_subdir:
            continue
        out_dir = f"{DATA_DIR}/{fold_pair_subdir}/output_esm_fold"
        if not os.path.isdir(out_dir) or len(os.listdir(out_dir)) == 0:
            print(f"[skip] No ESMFold outputs for {fold_pair_subdir} at {out_dir}")
            continue

        compute_esmfold_pred_tmscores(fold_pair_subdir)

    print("Finish all ESMFold TMscores computations!")


