#!/usr/bin/env python3
"""Fold-free per-cluster sequence composition features.

These feed the concordance false-positive calibration: pseudo-methods built
from random projections of these features read the same MSA the real
predictors read, but carry no structural information, so the number of
families they flag measures the false-positive rate of the concordance test.

Only composition is used. Nothing here can know about either fold.

Usage: cluster_composition.py OUT_CSV [PIPELINE_ROOT]
"""
import os, sys
import numpy as np

AA20 = "ACDEFGHIKLMNPQRSTVWY"
AA_IDX = {a: i for i, a in enumerate(AA20)}
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5,
      "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
      "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
      "Y": -1.3, "V": 4.2}
CHARGE = {"K": 1.0, "R": 1.0, "D": -1.0, "E": -1.0}
KD_VEC = np.array([KD[a] for a in AA20])
CH_VEC = np.array([CHARGE.get(a, 0.0) for a in AA20])

DEFAULT_ROOT = "/sci/labs/orzuk/orzuk/github/MsaCluster/Pipeline/FoldPairs"
COLS = (["pair_id", "cluster", "n_members", "gap_frac", "log_size",
         "net_charge", "hydrophobicity"] + [f"f_{a}" for a in AA20])


def read_a3m(path):
    seqs, cur = [], []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                    cur = []
            else:
                cur.append(line)
    if cur:
        seqs.append("".join(cur))
    return seqs


def features(members):
    """(20 aa frequencies, gap_frac, log_size, net_charge, hydrophobicity).

    Counts run over match-column residues only: upper case is a match-state
    residue, '-' a deletion at a match column, lower case an insertion that
    carries no column and is ignored, exactly as the trait code does.
    """
    C = np.zeros((len(members), 20))
    gaps = np.zeros(len(members))
    for i, s in enumerate(members):
        ng = 0
        for c in s:
            j = AA_IDX.get(c)
            if j is not None:
                C[i, j] += 1.0
            elif c == "-":
                ng += 1
        gaps[i] = ng / (C[i].sum() + ng) if (C[i].sum() + ng) else np.nan
    tot = C.sum(axis=1)
    ok = tot > 0
    if not ok.any():
        return None
    F = C[ok] / tot[ok, None]
    freq = F.mean(axis=0)
    return dict(
        n_members=len(members),
        gap_frac=float(np.nanmean(gaps)) if np.isfinite(gaps).any() else np.nan,
        log_size=float(np.log(len(members))),
        net_charge=float((C[ok] @ CH_VEC).mean()),
        hydrophobicity=float((F @ KD_VEC).mean()),
        **{f"f_{a}": float(freq[i]) for i, a in enumerate(AA20)},
    )


def write_table(out_path, root):
    """Walk every pair's clusters and write one composition row each."""
    n_ok = n_skip = 0
    with open(out_path, "w") as out:
        out.write(",".join(COLS) + "\n")
        for pair in sorted(os.listdir(root)):
            cdir = os.path.join(root, pair, "output_msa_cluster")
            if not os.path.isdir(cdir):
                continue
            for fn in sorted(os.listdir(cdir)):
                if not (fn.startswith("ShallowMsa_") and fn.endswith(".a3m")):
                    continue
                seqs = read_a3m(os.path.join(cdir, fn))
                if len(seqs) <= 1:          # row 0 is the query/medoid
                    n_skip += 1
                    continue
                d = features(seqs[1:])
                if d is None:
                    n_skip += 1
                    continue
                d["pair_id"] = pair
                d["cluster"] = fn[:-4]
                out.write(",".join(str(d[c]) for c in COLS) + "\n")
                n_ok += 1
            print(f"{pair} done ({n_ok} clusters)", flush=True)
    print(f"TOTAL clusters={n_ok} skipped={n_skip}", flush=True)
    return n_ok, n_skip


if __name__ == "__main__":
    write_table(sys.argv[1],
                sys.argv[2] if len(sys.argv) > 2 else DEFAULT_ROOT)
