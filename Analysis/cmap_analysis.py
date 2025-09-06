import os, sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, ROOT)

from tqdm import tqdm
from Analysis.PlotUtils import *
from config import *
import pandas as pd
from utils.utils import *
from utils.align_utils import *
import argparse


def align_and_resize_contact_maps(cmap1, cmap2, window_size=10, step_size=1):
    """
    Align two contact maps and resize them to match the smaller map's dimensions.
    :param cmap1: First contact map (2D symmetric numpy array)
    :param cmap2: Second contact map (2D symmetric numpy array)
    :param window_size: Size of the sliding window for comparison
    :param step_size: Step size for sliding the window
    :return: Tuple of (aligned_cmap1, aligned_cmap2), both with dimensions of the smaller input map
    """

    # Determine which map is larger
    if cmap1.shape[0] * cmap1.shape[1] >= cmap2.shape[0] * cmap2.shape[1]:
        larger_cmap, smaller_cmap = cmap1, cmap2
    else:
        larger_cmap, smaller_cmap = cmap2, cmap1

    best_score = float('-inf')
    best_offset = (0, 0)

    # Find the best alignment
    for i in range(0, larger_cmap.shape[0] - smaller_cmap.shape[0] + 1, step_size):
        for j in range(0, larger_cmap.shape[1] - smaller_cmap.shape[1] + 1, step_size):
            window = larger_cmap[i:i + smaller_cmap.shape[0], j:j + smaller_cmap.shape[1]]
            score = np.sum(window * smaller_cmap)  # Simple dot product for similarity

            if score > best_score:
                best_score = score
                best_offset = (i, j)

    # Extract the aligned portion from the larger map
    aligned_larger = larger_cmap[best_offset[0]:best_offset[0] + smaller_cmap.shape[0],
                     best_offset[1]:best_offset[1] + smaller_cmap.shape[1]]

    # Ensure we return cmap1 and cmap2 in the correct order
    if cmap1.shape[0] * cmap1.shape[1] >= cmap2.shape[0] * cmap2.shape[1]:
        return aligned_larger, smaller_cmap
    else:
        return smaller_cmap, aligned_larger


def load_cmap(path):
    return np.load(path)

def get_only_cmaps(cmap1,cmap2):
    diff_folds = cmap1 - cmap2
    only_fold1 = diff_folds.copy()
    only_fold1[only_fold1 == -1] = 0
    only_fold2 = diff_folds.copy()
    only_fold2[only_fold2 == 1] = 0
    only_fold2[only_fold2 == -1] = 1
    return only_fold1,only_fold2



import numpy as np


# Evaluate a predicted contact map against 2 ground truth contact maps
def evaluate_pred_cmap(
    pred_map: np.ndarray,
    truth1_bin: np.ndarray,
    truth2_bin: np.ndarray,
    *,
    thresh: float = 0.4,   # probability threshold if pred_map is probabilistic
    sep_min: int = 6,      # ignore |i-j| < sep_min in sequence positional distance
    index_tol: int = 0,    # allow ±k residue index tolerance in sequence positional distance (0 = strict)
    symmetrize: bool = True,
) -> dict:
    """Flat metrics dict for one predicted map vs 2 truths (+ common/unique sets)."""

    # 1) binarize & symmetrize
    if pred_map.dtype.kind in "fc":
        pred_bin = (pred_map > thresh).astype(np.uint8)
    else:
        pred_bin = (pred_map > 0).astype(np.uint8)

    t1 = (truth1_bin > 0).astype(np.uint8)
    t2 = (truth2_bin > 0).astype(np.uint8)

    if symmetrize:
        pred_bin = ((pred_bin | pred_bin.T) > 0).astype(np.uint8)
        t1 = ((t1 | t1.T) > 0).astype(np.uint8)
        t2 = ((t2 | t2.T) > 0).astype(np.uint8)

    n = pred_bin.shape[0]
    assert t1.shape == (n, n) and t2.shape == (n, n), "All maps must be same NxN"

    # 2) optional index tolerance (pure NumPy, small & clear)
    pred_eff = pred_bin.copy()
    tol = int(index_tol)
    if tol > 0:
        for di in range(-tol, tol+1):
            for dj in range(-tol, tol+1):
                if di == 0 and dj == 0:
                    continue
                sh = np.zeros_like(pred_bin)
                si0, si1 = max(0, -di), min(n, n-di)
                sj0, sj1 = max(0, -dj), min(n, n-dj)
                di0, dj0 = max(0, di), max(0, dj)
                if si1 > si0 and sj1 > sj0:
                    sh[di0:di0+(si1-si0), dj0:dj0+(sj1-sj0)] = pred_bin[si0:si1, sj0:sj1]
                pred_eff |= sh
        pred_eff = pred_eff.astype(np.uint8)

    # 3) score only upper triangle with sequence-separation filter
    ksep = max(1, int(sep_min))
    iu = np.triu_indices(n, k=ksep)

    def compute_metrics(truth_bin: np.ndarray) -> dict:
        t = truth_bin[iu].astype(np.uint8)
        p = pred_eff[iu].astype(np.uint8)
        tp = int(np.sum((t == 1) & (p == 1)))
        fp = int(np.sum((t == 0) & (p == 1)))
        fn = int(np.sum((t == 1) & (p == 0)))
        tn = int(np.sum((t == 0) & (p == 0)))
        prec = tp / (tp + fp) if (tp + fp) else 0.0
        rec  = tp / (tp + fn) if (tp + fn) else 0.0
        f1   = (2*prec*rec) / (prec + rec) if (prec + rec) else 0.0
        jac  = tp / (tp + fp + fn) if (tp + fp + fn) else 0.0
        denom = (tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)
        mcc  = ((tp*tn - fp*fn) / np.sqrt(denom)) if denom else 0.0
        return {
            "tp": tp, "fp": fp, "fn": fn, "tn": tn,
            "precision": round(prec, 4),
            "recall":    round(rec,  4),
            "f1":        round(f1,   4),
            "jaccard":   round(jac,  4),
            "mcc":       round(float(mcc), 4),
            "contacts_truth": int(tp + fn),
            "contacts_pred":  int(tp + fp),
        }

    # sets
    common = (t1 & t2).astype(np.uint8)
    uniq1  = (t1 & (1 - t2)).astype(np.uint8)
    uniq2  = (t2 & (1 - t1)).astype(np.uint8)

    # 4) assemble flat dict with clear prefixes (no dict-comprehension gymnastics)
    out = {
        "n": n,
        "sep_min": ksep,
        "index_tol": tol,
        "thresh": float(thresh),
    }

    for label, truth in (("t1", t1), ("t2", t2), ("common", common), ("uniq1", uniq1), ("uniq2", uniq2)):
        metrics = compute_metrics(truth)
        for mname, val in metrics.items():
            out[f"{label}_{mname}"] = val

    return out


# --- replace your existing main with this ---
if __name__ == "__main__":
    # Parse optional CLI pair
    parser = argparse.ArgumentParser(description="CMAP analysis for one pair or all pairs")
    parser.add_argument("pair", nargs="?", help="Optional pair like 1dzlA_5keqF")
    args = parser.parse_args()

    DEFAULT_PYCHARM_PAIR = "1dzlA_5keqF"

    if args.pair:  # explicit CLI pair
        fold_pairs = [pair_str_to_tuple(args.pair)]
        mode = "single-arg"
    elif is_pycharm():  # running from PyCharm, no CLI arg
        fold_pairs = [pair_str_to_tuple(DEFAULT_PYCHARM_PAIR)]
        mode = "pycharm-default"
    else:  # no arg, non-PyCharm → run all pairs
        fold_pairs = list_protein_pairs()
        mode = "all"

    print(f"[cmap_analysis] Mode={mode} | total pairs={len(fold_pairs)}")

    fold_pair_errors = []
    for fold_pair in fold_pairs:
        chains = (fold_pair[0][-1], fold_pair[1][-1])
        # ensure tuple form like ('1dzlA','5keqF')
        if isinstance(fold_pair, str):
            fold_pair = pair_str_to_tuple(fold_pair)

        fold_pair_subdir = f"{fold_pair[0]}_{fold_pair[1]}"
        res = []
        try:
            print("Set plot_tool for: ", fold_pair_subdir)
            plot_tool = PlotTool(folder=DATA_DIR, fold_pair=fold_pair_subdir)

            print("Extract sequences for: ", fold_pair_subdir)
            seq1 = extract_protein_sequence(f"{DATA_DIR}/{fold_pair_subdir}/{fold_pair[0][:-1]}.pdb",
                                            chain=chains[0], ca_only=True)
            seq2 = extract_protein_sequence(f"{DATA_DIR}/{fold_pair_subdir}/{fold_pair[1][:-1]}.pdb",
                                            chain=chains[1], ca_only=True)  # Add chain!

            print("Compute contact maps: ", fold_pair_subdir)
            cmap_pdb1, _ = pdb_to_contact_map(f"{DATA_DIR}/{fold_pair_subdir}/{fold_pair[0][:-1]}.pdb",
                                              chain=chains[0])
            cmap_pdb2, _ = pdb_to_contact_map(f"{DATA_DIR}/{fold_pair_subdir}/{fold_pair[1][:-1]}.pdb",
                                              chain=chains[1])

            print(f"[debug] cmap sizes: n1={cmap_pdb1.shape[0]}, n2={cmap_pdb2.shape[0]}")
            print(f"[debug] seq lengths: L1={len(seq1)}, L2={len(seq2)}")

#            print("DEBUG SEQ ALIGNMENT: ", fold_pair_subdir)
#            pair_alignment = get_align_indexes(seq1, seq2, mode="standard", debug=True)
#            print("DEBUG SEQ ALIGNMENT: BLOSUM ", fold_pair_subdir)
#            pair_alignment_blosum = get_align_indexes(seq1, seq2, mode="blosum", debug=True)

            print("NEW!!! Align cmaps for: ", fold_pair_subdir)
#            cmap_aligned_pdb1, cmap_aligned_pdb2 = align_and_resize_contact_maps(
#                cmap_pdb1, cmap_pdb2, window_size=1, step_size=1)  # OLD!!!
            cmap_aligned_pdb1, cmap_aligned_pdb2, (idx1, idx2) = align_cmaps_by_sequence(
                cmap_pdb1, seq1, cmap_pdb2, seq2, mode="standard") #"standard") # 'blosum')

            print("Get only cmaps for: ", fold_pair_subdir)
            cmap_only_pdb1, cmap_only_pdb2 = get_only_cmaps(cmap_aligned_pdb1, cmap_aligned_pdb2)

            # === NEW: metrics only, no visualization files ===
            pred_dir = f"{DATA_DIR}/{fold_pair_subdir}/output_cmaps/msa_transformer"
            cmaps = sorted(os.listdir(pred_dir))
            print("Loop on cmaps for: ", fold_pair_subdir)
            print("List of cmaps: ", cmaps)

            # Aligned truths (same NxN already)
            truth1_aln = cmap_aligned_pdb1.astype(np.uint8)
            truth2_aln = cmap_aligned_pdb2.astype(np.uint8)
            n_common = truth1_aln.shape[0]

            # For mapping predictions onto the same common frame:
            # idx1 are positions in seq1 that survived alignment (map to common frame)
            # idx2 are positions in seq2 that survived alignment (map to common frame)
            idx1 = np.asarray(idx1, dtype=int)
            idx2 = np.asarray(idx2, dtype=int)
            L1 = len(seq1)
            L2 = len(seq2)

            for fname in tqdm(cmaps):
                try:
                    # Only shallow predicted cmaps
                    if not (fname.endswith(".npy") and "ShallowMsa" in fname):
                        continue
                    if "visualization_map" in fname or "VizCmaps" in fname or "deep" in fname:
                        continue

                    stem = os.path.splitext(os.path.basename(fname))[0]  # e.g., 'msa_t__ShallowMsa_000'
                    m = re.search(r"ShallowMsa_(\d+)", stem)
                    cluster = m.group(1) if m else "UNK"

                    pred_path = os.path.join(pred_dir, fname)
                    pred_full = np.load(pred_path)  # float or binary
                    if pred_full.ndim != 2 or pred_full.shape[0] != pred_full.shape[1]:
                        print("[warn]", fname, "-> prediction is not square, skipping")
                        continue

                    Np = pred_full.shape[0]
                    # Map prediction to the same common aligned frame used for truth1_aln/truth2_aln
                    if Np == L1:
                        # prediction is on seq1 indexing -> take common positions idx1
                        pred_common = pred_full[np.ix_(idx1, idx1)]
                        src_seq = "seq1"
                    elif Np == L2:
                        # prediction is on seq2 indexing -> take common positions idx2
                        pred_common = pred_full[np.ix_(idx2, idx2)]
                        src_seq = "seq2"
                    elif Np == n_common:
                        # already on the common frame (unlikely, but handle it)
                        pred_common = pred_full
                        src_seq = "common"
                    else:
                        print("[warn]", fname, f"-> size {Np} does not match len(seq1)={L1}, len(seq2)={L2}, "
                                               f"or common={n_common}; skipping")
                        continue

                    if pred_common.shape[0] != n_common:
                        print("[warn]", fname, "-> mapped prediction still not NxN common; skipping")
                        continue

                    # Evaluate on the common frame against the two aligned truths
                    # Evaluate
                    metrics = evaluate_pred_cmap(
                        pred_map=pred_common,
                        truth1_bin=truth1_aln,
                        truth2_bin=truth2_aln,
                        thresh=0.4,  # used only if pred is probabilistic
                        sep_min=6,  # ignore |i-j| < 6 (diagonal/near-diagonal)
                        index_tol=0,  # 0 strict; 1–2 if you want “nearby” matches allowed
                        symmetrize=True
                    )

                    # Keep metrics only (drop pred_source_seq, n, sep_min, index_tol, thresh)
                    metric_keys = [k for k in metrics.keys()
                                   if k.startswith(("t1_", "t2_", "common_", "uniq1_", "uniq2_"))]

                    row = {"file": fname}  # keep just the filename for reference (optional)
                    row.update({k: metrics[k] for k in metric_keys})
                    res.append(row)


                except Exception as e:
                    print("[warn]", fname, "->", e)
                    continue

            # Save data-frame for each protein pair
            os.makedirs(f"{DATA_DIR}/{fold_pair_subdir}/Analysis", exist_ok=True)
            df_cmap = pd.DataFrame(res)
            df_cmap.to_csv(f"{DATA_DIR}/{fold_pair_subdir}/Analysis/df_cmap.csv", index=False)
            print("Saved metrics:", f"{DATA_DIR}/{fold_pair_subdir}/Analysis/df_cmap.csv", "| rows:", len(df_cmap))

        except Exception as err:
            print("[error pair]", fold_pair_subdir, "->", err)
    print("Finish all MSA-Transformer CMAP Similarity computations!")