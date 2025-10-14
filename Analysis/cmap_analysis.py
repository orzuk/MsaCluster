import os, sys
import re
import numpy as np
import pandas as pd
from tqdm import tqdm
import argparse
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, ROOT)

from config import *
from utils.utils import *
from utils.align_utils import *


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


def load_cmap_and_idx(path):
    """
    Returns: (cmap: np.ndarray[L,L], idx: np.ndarray[L_idx], keep_cols or None)
    - .npz (preferred): expects keys 'cmap', 'idx' (and optional 'keep_cols')
    - legacy .npy dict: {'cmap':..., 'idx':...}
    - legacy plain .npy: cmap only + sidecar .idx.npy (if present)
    """
    if path.endswith(".npz"):
        z = np.load(path)
        return z["cmap"], z["idx"], (z["keep_cols"] if "keep_cols" in z.files else None)
    arr = np.load(path, allow_pickle=True)
    if isinstance(arr, np.ndarray) and arr.dtype == object:
        d = arr.item()
        return d["cmap"], d["idx"], d.get("keep_cols", None)
    # plain array case:
    cmap = arr
    idx_path = path.replace(".npy", ".idx.npy")
    if os.path.isfile(idx_path):
        idx = np.load(idx_path)
        return cmap, idx, None
    raise FileNotFoundError(f"{path}: legacy npy without idx sidecar; re-run MSAT to produce .npz.")


def seed_residue_indices(seed_a3m_path: str, tag: str):
    # tag is "4dxrA" or "4dxtA"
    lines = [ln.rstrip("\n") for ln in open(seed_a3m_path) if ln.strip()]
    seqs, heads, cur = [], [], []
    for ln in lines:
        if ln.startswith(">"):
            if cur: seqs.append("".join(cur)); cur=[]
            heads.append(ln[1:].split()[0])
        else:
            cur.append(ln)
    if cur: seqs.append("".join(cur))
    if heads[0].startswith(tag):
        s = seqs[0]
    elif heads[1].startswith(tag):
        s = seqs[1]
    else:
        raise RuntimeError(f"{tag} not found in {seed_a3m_path}")
    cols = [c for c in s if c.isupper() or c == '-']
    return np.asarray([i for i,c in enumerate(cols) if c.isupper()], dtype=np.int32)


def get_only_cmaps(cmap1,cmap2):
    diff_folds = cmap1 - cmap2
    only_fold1 = diff_folds.copy()
    only_fold1[only_fold1 == -1] = 0
    only_fold2 = diff_folds.copy()
    only_fold2[only_fold2 == 1] = 0
    only_fold2[only_fold2 == -1] = 1
    return only_fold1,only_fold2


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

    # ---- 2) alignment frames from _seed_both.a3m (symmetric mapping)
    both_a3m = os.path.join(pair_dir, "_seed_both.a3m")
    pairtrim_seed = os.path.join(pair_dir, "_seed_both_pairtrim.a3m")
    if os.path.isfile(pairtrim_seed):
        both_a3m = pairtrim_seed
    s1, s2 = _read_two_seeds_from_a3m(both_a3m)
    keep1 = _uppercase_cols(s1) if s1 else []
    keep2 = _uppercase_cols(s2) if s2 else []
    common_cols = _common_cols(s1, s2)


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


    # ---- 3) iterate predicted maps (npz with cmap+idx; also supports legacy npy)
    pred_dir = os.path.join(pair_dir, "output_cmaps", "msa_transformer")
    files = []
    if os.path.isdir(pred_dir):
        for fn in sorted(os.listdir(pred_dir)):
            if not fn.startswith("msa_t_"):  # keep all msa_t_* variants
                continue
            if not (fn.endswith(".npz") or fn.endswith(".npy")):
                continue
            if (not include_deep) and ("DeepMsa" in fn or "MSA_deep" in fn):
                continue
            files.append(os.path.join(pred_dir, fn))

    # Build residue-index mappings from seed columns -> residue indices for each fold
    # (col2res1/2 map a SEED column index -> residue index in that fold's truth matrix)
    col2res1 = _col_to_res(s1) if s1 else {}
    col2res2 = _col_to_res(s2) if s2 else {}

    rows = []
    for f in files:
        # Load prediction and the seed-frame idx it uses
        try:
            P, idx, keep_cols = load_cmap_and_idx(f)   # P: (L,L) ; idx: seed-frame indices row-1 covers
        except Exception as e:
            print("[cmap] WARN skipping", f, "->", e)
            continue

        if P.ndim != 2 or P.shape[0] != P.shape[1]:
            print("[cmap] WARN skipping non-square", f, P.shape)
            continue

        # Map SEED indices -> P-axis positions
        pos = {int(c): i for i, c in enumerate(idx)}

        # --- Build evaluation frames (SEED set ∩ pred idx), then map to P's axis
        # Seed sets (as columns): keep1 = non-gaps in 4dxrA; keep2 = non-gaps in 4dxtA; common_cols = both non-gaps
        I1_cols = sorted(set(keep1).intersection(idx))             # fold1 columns present in prediction
        I2_cols = sorted(set(keep2).intersection(idx))             # fold2 columns present in prediction
        IC_cols = sorted(set(common_cols).intersection(idx))       # common columns present in prediction
        IU1_cols = sorted(set(keep1).difference(keep2).intersection(idx))  # uniq1 columns present
        IU2_cols = sorted(set(keep2).difference(keep1).intersection(idx))  # uniq2 columns present

        # Convert SEED columns -> P-axis indices
        M1  = np.fromiter((pos[c] for c in I1_cols), dtype=int)
        M2  = np.fromiter((pos[c] for c in I2_cols), dtype=int)
        MC  = np.fromiter((pos[c] for c in IC_cols), dtype=int)
        MU1 = np.fromiter((pos[c] for c in IU1_cols), dtype=int)
        MU2 = np.fromiter((pos[c] for c in IU2_cols), dtype=int)

        # Slice predictions into each frame
        P1  = P[np.ix_(M1,  M1)]
        P2  = P[np.ix_(M2,  M2)]
        PC  = P[np.ix_(MC,  MC)]
        PU1 = P[np.ix_(MU1, MU1)]
        PU2 = P[np.ix_(MU2, MU2)]

        # --- Build truth frames by mapping SEED columns -> residue indices
        R1  = np.fromiter((col2res1[c] for c in I1_cols), dtype=int)  # length == P1.shape[0]
        R2  = np.fromiter((col2res2[c] for c in I2_cols), dtype=int)
        RC1 = np.fromiter((col2res1[c] for c in IC_cols), dtype=int)
        RC2 = np.fromiter((col2res2[c] for c in IC_cols), dtype=int)
        RU1 = np.fromiter((col2res1[c] for c in IU1_cols), dtype=int)
        RU2 = np.fromiter((col2res2[c] for c in IU2_cols), dtype=int)

        # Fold-specific truths
        T1S = T1[np.ix_(R1, R1)] if len(R1) else None
        T2S = T2[np.ix_(R2, R2)] if len(R2) else None

        # Common/uniq truths (in the SAME per-frame sizes)
        T1C = T1[np.ix_(RC1, RC1)] if len(RC1) else None
        T2C = T2[np.ix_(RC2, RC2)] if len(RC2) else None
        TU1 = T1[np.ix_(RU1, RU1)] if len(RU1) else None
        TU2 = T2[np.ix_(RU2, RU2)] if len(RU2) else None

        # --- Metrics
        def safe_metrics(Psub, Tsub):
            if Psub is None or Tsub is None or Psub.shape != Tsub.shape or Psub.size == 0:
                return dict(precision=np.nan, recall=np.nan, f1=np.nan, jaccard=np.nan, mcc=np.nan)
            pb = _bin_pred(Psub, thresh=thresh, sep_min=sep_min, index_tol=index_tol)
            return _metrics(pb, Tsub)

        m_t1 = safe_metrics(P1,  T1S)
        m_t2 = safe_metrics(P2,  T2S)

        # For common/uniq, build truth masks by combining T1/T2 in the SAME sub-frame sizes
        def comb_metrics(Psub, Tsub1, Tsub2, mode):
            if Psub is None or Tsub1 is None or Tsub2 is None or Psub.shape != Tsub1.shape or Psub.shape != Tsub2.shape or Psub.size == 0:
                return dict(precision=np.nan, recall=np.nan, f1=np.nan, jaccard=np.nan, mcc=np.nan)
            if mode == "common":
                truth = (Tsub1.astype(bool) & Tsub2.astype(bool))
            elif mode == "uniq1":
                truth = (Tsub1.astype(bool) & ~Tsub2.astype(bool))
            else:
                truth = (Tsub2.astype(bool) & ~Tsub1.astype(bool))
            pb = _bin_pred(Psub, thresh=thresh, sep_min=sep_min, index_tol=index_tol)
            return _metrics(pb, truth.astype(np.uint8))

        m_com = comb_metrics(PC,  T1C, T2C, "common")
        m_u1  = comb_metrics(PU1, T1C, T2C, "uniq1") if (len(RC1)==len(RU1)) else comb_metrics(PU1, TU1, TU2, "uniq1")
        m_u2  = comb_metrics(PU2, T1C, T2C, "uniq2") if (len(RC2)==len(RU2)) else comb_metrics(PU2, TU1, TU2, "uniq2")

        base = os.path.basename(f)
        if "MSA_deep" in base or "DeepMsa" in base:
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
            "t1_precision": m_t1["precision"], "t1_recall": m_t1["recall"], "t1_f1": m_t1["f1"], "t1_jaccard": m_t1["jaccard"], "t1_mcc": m_t1["mcc"],
            "t2_precision": m_t2["precision"], "t2_recall": m_t2["recall"], "t2_f1": m_t2["f1"], "t2_jaccard": m_t2["jaccard"], "t2_mcc": m_t2["mcc"],
            # symmetric categories
            "common_precision": m_com["precision"], "common_recall": m_com["recall"], "common_f1": m_com["f1"], "common_jaccard": m_com["jaccard"], "common_mcc": m_com["mcc"],
            "uniq1_precision": m_u1["precision"],  "uniq1_recall":  m_u1["recall"],  "uniq1_f1":  m_u1["f1"],  "uniq1_jaccard":  m_u1["jaccard"],  "uniq1_mcc":  m_u1["mcc"],
            "uniq2_precision": m_u2["precision"],  "uniq2_recall":  m_u2["recall"],  "uniq2_f1":  m_u2["f1"],  "uniq2_jaccard":  m_u2["jaccard"],  "uniq2_mcc":  m_u2["mcc"],
        })


    # ---- 4) write CSV
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

