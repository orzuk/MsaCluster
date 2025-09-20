import os, sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, ROOT)

from tqdm import tqdm
from config import *
from utils.utils import *
from utils.align_utils import *
import argparse


import re
import numpy as np

CONTACT_CUTOFF = 8.0  # Å; adjust if you use a different truth definition

def _read_two_seeds_from_a3m(a3m_path: str):
    """
    Return (seed1, seed2) strings from an A3M; if only 1 seed is present,
    seed2 may be None. Lines are concatenated and kept exactly as in A3M
    (uppercase = aligned residue, lowercase = insertion, '-' = gap).
    """
    if not os.path.isfile(a3m_path):
        return None, None
    seqs = []
    take = False
    cur = []
    with open(a3m_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur)); cur = []
                take = True
                continue
            if take:
                s = line.strip()
                if not s:
                    continue
                cur.append(s)
        if cur:
            seqs.append("".join(cur))
    # top two
    s1 = seqs[0] if len(seqs) >= 1 else None
    s2 = seqs[1] if len(seqs) >= 2 else None
    return s1, s2


def _uppercase_cols(seed: str):
    """Return list of alignment-column indices to KEEP (seed has uppercase AAs)."""
    if not seed:
        return []
    return [i for i, ch in enumerate(seed) if ch.isalpha() and ch.isupper()]


def _common_cols(seed1: str, seed2: str):
    """Columns where BOTH seeds have uppercase residues (shared-by-residue frame)."""
    if not seed1 or not seed2:
        return []
    keep1 = _uppercase_cols(seed1)
    keep2 = _uppercase_cols(seed2)
    K1 = set(keep1)
    K2 = set(keep2)
    return sorted(list(K1.intersection(K2)))


def _coords_from_pdb_chain(pdb_path: str, chain_id: str | None):
    """Extract CA coords (N×3). Falls back to CB if needed."""
    if not os.path.isfile(pdb_path):
        return None
    ca = []
    cb = []
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            if chain_id and line[21].strip() != chain_id:
                continue
            atom = line[12:16].strip()
            if atom not in ("CA", "CB"):
                continue
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except Exception:
                continue
            if atom == "CA":
                ca.append((x,y,z))
            elif atom == "CB":
                cb.append((x,y,z))
    arr = np.array(ca if ca else cb, dtype=np.float32)
    return arr if len(arr) else None


def _truth_contacts(coords: np.ndarray, cutoff=CONTACT_CUTOFF, sep_min=6):
    """Return boolean NxN contact mask from coordinates."""
    if coords is None or len(coords) == 0:
        return None
    X = coords
    d2 = np.sum((X[:,None,:] - X[None,:,:])**2, axis=2)
    mask = d2 <= (cutoff*cutoff)
    # remove diagonal & short-range
    n = mask.shape[0]
    for i in range(n):
        lo = max(0, i - (sep_min-1))
        hi = min(n, i + (sep_min))
        mask[i, lo:hi] = False
    np.fill_diagonal(mask, False)
    # make symmetric
    mask = np.triu(mask, 1) | np.tril(mask, -1)
    return mask


def _bin_pred(mat: np.ndarray, thresh=0.4, sep_min=6, index_tol=0):
    """Threshold a probability CMAP; apply sep_min and optional +/- index tolerance."""
    n = mat.shape[0]
    pred = (mat >= thresh)
    # remove diagonal & short-range
    for i in range(n):
        lo = max(0, i - (sep_min-1))
        hi = min(n, i + (sep_min))
        pred[i, lo:hi] = False
    np.fill_diagonal(pred, False)
    pred = np.triu(pred, 1) | np.tril(pred, -1)

    if index_tol > 0:
        # simple dilation on the binary grid within a Chebyshev radius 'index_tol'
        pred_eff = pred.copy()
        for di in range(-index_tol, index_tol+1):
            for dj in range(-index_tol, index_tol+1):
                if di == 0 and dj == 0:
                    continue
                sh = np.zeros_like(pred)
                si0, si1 = max(0, -di), min(n, n-di)
                sj0, sj1 = max(0, -dj), min(n, n-dj)
                di0, dj0 = max(0, di), max(0, dj)
                if si1 > si0 and sj1 > sj0:
                    sh[di0:di0+(si1-si0), dj0:dj0+(sj1-sj0)] = pred[si0:si1, sj0:sj1]
                pred_eff |= sh
        pred = pred_eff
    return pred


def _metrics(pred_bin: np.ndarray, truth_bin: np.ndarray):
    """Return precision, recall, F1, Jaccard, MCC on upper-tri entries."""
    if pred_bin is None or truth_bin is None:
        return dict(precision=np.nan, recall=np.nan, f1=np.nan, jaccard=np.nan, mcc=np.nan)
    iu = np.triu_indices_from(truth_bin, 1)
    p = pred_bin[iu].astype(bool)
    t = truth_bin[iu].astype(bool)
    tp = int(np.sum(p & t))
    fp = int(np.sum(p & ~t))
    fn = int(np.sum(~p & t))
    tn = int(np.sum(~p & ~t))
    prec = tp / (tp + fp) if (tp+fp) else np.nan
    rec  = tp / (tp + fn) if (tp+fn) else np.nan
    f1   = 2*prec*rec / (prec+rec) if prec>0 and rec>0 else (1.0 if (prec==1 and rec==1) else (0.0 if (prec==0 or rec==0) else np.nan))
    jac  = tp / (tp + fp + fn) if (tp+fp+fn) else np.nan
    # MCC
    denom = (tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)
    mcc = ((tp*tn - fp*fn) / np.sqrt(denom)) if denom else np.nan
    return dict(precision=prec, recall=rec, f1=f1, jaccard=jac, mcc=mcc)


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


def compute_cmap_metrics_for_pair(
    subdir: str,
    include_deep: bool = True,
    thresh: float = 0.4,
    sep_min: int = 6,
    index_tol: int = 0,
):
    """
    Write Pipeline/<pair>/Analysis/df_cmap.csv with per-file metrics:
      t1_* / t2_* (per fold), and common_* / uniq1_* / uniq2_* (shared & contact-unique).
    Symmetric: truth2 does not depend on truth1 — both are aligned via _seed_both.a3m.
    """
    pair_dir = os.path.join(DATA_DIR, subdir)
    out_dir  = os.path.join(pair_dir, "Analysis")
    os.makedirs(out_dir, exist_ok=True)
    out_csv  = os.path.join(out_dir, "df_cmap.csv")

    # ---- 1) paths & truths
    a, b = subdir.split("_", 1)
    p1, c1 = a[:-1], a[-1]
    p2, c2 = b[:-1], b[-1]
    # prefer chain-specific PDBs if present
    cand1 = os.path.join(pair_dir, "chain_pdb_files", f"{a}.pdb")
    cand2 = os.path.join(pair_dir, "chain_pdb_files", f"{b}.pdb")
    pdb1  = cand1 if os.path.isfile(cand1) else os.path.join(pair_dir, f"{p1}.pdb")
    pdb2  = cand2 if os.path.isfile(cand2) else os.path.join(pair_dir, f"{p2}.pdb")

    coords1 = _coords_from_pdb_chain(pdb1, c1)
    coords2 = _coords_from_pdb_chain(pdb2, c2)
    if coords1 is None or coords2 is None:
        raise RuntimeError("Could not load truth coordinates for one or both chains")

    T1 = _truth_contacts(coords1, cutoff=CONTACT_CUTOFF, sep_min=sep_min)
    T2 = _truth_contacts(coords2, cutoff=CONTACT_CUTOFF, sep_min=sep_min)
    L1, L2 = T1.shape[0], T2.shape[0]

    # ---- 2) alignment frames from _seed_both.a3m (symmetric mapping)
    both_a3m = os.path.join(pair_dir, "_seed_both.a3m")
    s1, s2 = _read_two_seeds_from_a3m(both_a3m)
    keep1 = _uppercase_cols(s1) if s1 else []
    keep2 = _uppercase_cols(s2) if s2 else []
    common_cols = _common_cols(s1, s2)

    # map 'common' columns to residue indices in each fold
    # For each alignment column with uppercase in both seeds, we know its residue index
    # within each seed's ungapped sequence.
    def _ungapped_positions(seed):
        idxs = []
        resi = 0
        for i, ch in enumerate(seed or ""):
            if ch.isalpha():
                if ch.isupper():
                    resi += 1
                    idxs.append(resi-1)  # 0-based residue index in that seed
                # lowercase: insertion (doesn't count for residue index)
            else:
                # '-' gap
                pass
        return idxs

    # Build position maps from alignment column -> residue index
    def _col_to_res(seed):
        col2res = {}
        resi = 0
        for i, ch in enumerate(seed or ""):
            if ch.isalpha():
                if ch.isupper():
                    col2res[i] = resi  # map to current residue index (0-based)
                    resi += 1
            # lowercase or '-' → no mapping
        return col2res

    col2res1 = _col_to_res(s1) if s1 else {}
    col2res2 = _col_to_res(s2) if s2 else {}

    # For the common frame, we’ll need arrays that map common positions -> residue indices in each fold
    c_to_r1 = [col2res1[c] for c in common_cols] if common_cols else []
    c_to_r2 = [col2res2[c] for c in common_cols] if common_cols else []

    # Sanity: lengths expected after collapsing
    exp1 = len(keep1)
    exp2 = len(keep2)
    expC = len(common_cols)

    # ---- 3) iterate predicted maps
    pred_dir = os.path.join(pair_dir, "output_cmaps", "msa_transformer")
    files = []
    if os.path.isdir(pred_dir):
        for fn in sorted(os.listdir(pred_dir)):
            if fn.startswith("msa_t__") and fn.endswith(".npy"):
                if (not include_deep) and ("MSA_deep" in fn):
                    continue
                files.append(os.path.join(pred_dir, fn))

    rows = []
    for f in files:
        P = np.load(f)  # NxN in alignment columns
        N = P.shape[0]
        if P.shape[0] != P.shape[1]:
            # skip weird shapes
            continue

        # collapse to fold1 frame
        pred1 = None
        if exp1 and max(keep1) < N:
            pred1 = P[np.ix_(keep1, keep1)]
            if pred1.shape[0] != L1:
                # If unusual mismatch (e.g., missing residues), skip t1_* but still compute common
                pred1 = None

        # collapse to fold2 frame
        pred2 = None
        if exp2 and max(keep2) < N:
            pred2 = P[np.ix_(keep2, keep2)]
            if pred2.shape[0] != L2:
                pred2 = None

        # common frame (symmetric; does not depend on choosing fold1 or fold2)
        predC = None
        if expC and max(common_cols) < N:
            predC = P[np.ix_(common_cols, common_cols)]
            # Build T1/T2 in that same common frame by remapping residue indices:
            if expC and len(c_to_r1) == expC and len(c_to_r2) == expC:
                T1C = T1[np.ix_(c_to_r1, c_to_r1)]
                T2C = T2[np.ix_(c_to_r2, c_to_r2)]
            else:
                T1C = None; T2C = None
        else:
            T1C = None; T2C = None

        # --- metrics ---
        # Fold1
        t1_prec = t1_rec = t1_f1 = t1_jac = t1_mcc = np.nan
        if pred1 is not None:
            pb = _bin_pred(pred1, thresh=thresh, sep_min=sep_min, index_tol=index_tol)
            m  = _metrics(pb, T1)
            t1_prec, t1_rec, t1_f1, t1_jac, t1_mcc = m["precision"], m["recall"], m["f1"], m["jaccard"], m["mcc"]

        # Fold2
        t2_prec = t2_rec = t2_f1 = t2_jac = t2_mcc = np.nan
        if pred2 is not None:
            pb = _bin_pred(pred2, thresh=thresh, sep_min=sep_min, index_tol=index_tol)
            m  = _metrics(pb, T2)
            t2_prec, t2_rec, t2_f1, t2_jac, t2_mcc = m["precision"], m["recall"], m["f1"], m["jaccard"], m["mcc"]

        # Common / uniq-by-contact categories (symmetric)
        com_prec = com_rec = com_f1 = com_jac = com_mcc = np.nan
        u1_prec = u1_rec = u1_f1 = u1_jac = u1_mcc = np.nan
        u2_prec = u2_rec = u2_f1 = u2_jac = u2_mcc = np.nan

        if predC is not None and T1C is not None and T2C is not None:
            pb = _bin_pred(predC, thresh=thresh, sep_min=sep_min, index_tol=index_tol)

            common_mask = (T1C & T2C)
            uniq1_mask  = (T1C & ~T2C)
            uniq2_mask  = (T2C & ~T1C)

            # evaluate each category by masking the truth and reusing the same pb grid
            def _masked_metrics(pb, truth_mask):
                return _metrics(pb, truth_mask)

            mcom = _masked_metrics(pb, common_mask)
            mu1  = _masked_metrics(pb, uniq1_mask)
            mu2  = _masked_metrics(pb, uniq2_mask)

            com_prec, com_rec, com_f1, com_jac, com_mcc = mcom["precision"], mcom["recall"], mcom["f1"], mcom["jaccard"], mcom["mcc"]
            u1_prec,  u1_rec,  u1_f1,  u1_jac,  u1_mcc  = mu1["precision"],  mu1["recall"],  mu1["f1"],  mu1["jaccard"],  mu1["mcc"]
            u2_prec,  u2_rec,  u2_f1,  u2_jac,  u2_mcc  = mu2["precision"],  mu2["recall"],  mu2["f1"],  mu2["jaccard"],  mu2["mcc"]

        # File/meta
        base = os.path.basename(f)
        if "MSA_deep" in base:
            clus = "DeepMsa"
        else:
            m = re.search(r"(ShallowMsa_\d+)", base)
            clus = m.group(1) if m else "Unknown"

        rows.append({
            "pair_id": subdir,
            "file": base,
            "cluster": clus,
            "thresh": thresh,
            "sep_min": sep_min,
            # per fold
            "t1_precision": t1_prec, "t1_recall": t1_rec, "t1_f1": t1_f1, "t1_jaccard": t1_jac, "t1_mcc": t1_mcc,
            "t2_precision": t2_prec, "t2_recall": t2_rec, "t2_f1": t2_f1, "t2_jaccard": t2_jac, "t2_mcc": t2_mcc,
            # symmetric categories
            "common_precision": com_prec, "common_recall": com_rec, "common_f1": com_f1, "common_jaccard": com_jac, "common_mcc": com_mcc,
            "uniq1_precision": u1_prec, "uniq1_recall": u1_rec, "uniq1_f1": u1_f1, "uniq1_jaccard": u1_jac, "uniq1_mcc": u1_mcc,
            "uniq2_precision": u2_prec, "uniq2_recall": u2_rec, "uniq2_f1": u2_f1, "uniq2_jaccard": u2_jac, "uniq2_mcc": u2_mcc,
        })

    # ---- 4) write CSV
    import pandas as pd
    df = pd.DataFrame(rows)
    df.to_csv(out_csv, index=False)
    return df


# Compute and save contact maps metrics
if __name__ == "__main__":
    ap = argparse.ArgumentParser("Compute MSA-Transformer CMAP metrics for pair(s)")
    ap.add_argument("pairs", nargs="*", help="pair IDs like 1dzlA_5keqF; if empty, process all")
    ap.add_argument("--include_deep", action="store_true")
    ap.add_argument("--thresh", type=float, default=0.4)
    ap.add_argument("--sep_min", type=int, default=6)
    ap.add_argument("--index_tol", type=int, default=0)
    args = ap.parse_args()

    if args.pairs:
        todo = args.pairs
    else:
        from utils.utils import list_protein_pairs
        todo = list_protein_pairs(parsed=False, sort_result=True)

    for pid in todo:
        try:
            compute_cmap_metrics_for_pair(
                pid,
                include_deep=args.include_deep,
                thresh=args.thresh,
                sep_min=args.sep_min,
                index_tol=args.index_tol,
            )
            print("[cmap] wrote", os.path.join(DATA_DIR, pid, "Analysis", "df_cmap.csv"))
        except Exception as e:
            print("[cmap] ERROR", pid, "→", e)

