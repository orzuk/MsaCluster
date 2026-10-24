#!/usr/bin/env python3
"""Divergence baseline for the XCL1 ancestor calibration.

The calibration reports that predictors rank four reconstructed ancestors and
human XCL1 consistently with their NMR-measured alternative-fold populations.
A rank correlation that high is only evidence about structure prediction if a
trivial sequence property cannot do the same. This scores three such
baselines against the same measured populations, with the same exact
permutation test over all 5! = 120 orderings:

  evo_order     position in the reconstructed lineage (time alone)
  seq_len       length of the sequence
  pid_to_XCL1   gap-aware identity to human XCL1 in the aligned fasta

Anything a baseline explains is not evidence that the predictor read fold.
"""
import itertools
import os
import sys

import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV = os.path.join(REPO, "docs", "ancestor_calibration_gpu.csv")
FASTA = os.path.join(REPO, "data", "ancestors", "dishman_ancestors_aligned.fasta")


def read_fasta(path):
    seqs, name, cur = {}, None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(cur)
            name, cur = line[1:].split()[0], []
        elif line:
            cur.append(line)
    if name:
        seqs[name] = "".join(cur)
    return seqs


def rankdata(a):
    a = np.asarray(a, float)
    order = np.argsort(a, kind="mergesort")
    r = np.empty(len(a), float)
    r[order] = np.arange(1, len(a) + 1)
    # average ties, as scipy does
    for v in np.unique(a):
        m = a == v
        if m.sum() > 1:
            r[m] = r[m].mean()
    return r


def rho(x, y):
    rx, ry = rankdata(x), rankdata(y)
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    d = np.sqrt((rx ** 2).sum() * (ry ** 2).sum())
    return float((rx * ry).sum() / d) if d > 0 else np.nan


def exact_p(x, y):
    """One-sided exact p over all 120 orderings of the five labels."""
    obs = rho(x, y)
    n = 0
    perms = list(itertools.permutations(range(len(x))))
    for p in perms:
        if rho(np.asarray(x)[list(p)], y) >= obs - 1e-12:
            n += 1
    return obs, n / len(perms), len(perms)


def gapaware_pid(a, b):
    num = den = 0
    for ca, cb in zip(a, b):
        if ca == "-" and cb == "-":
            continue
        den += 1
        if ca == cb:
            num += 1
    return num / den if den else np.nan


def main(out_csv=None):
    df = pd.read_csv(CSV)
    anc = (df[df.kind == "ancestor"]
           .drop_duplicates(subset=["label"])
           .sort_values("evo_order"))
    print("ancestors:", list(anc.label))
    pct = anc.pct_alt_measured.to_numpy(float)
    print("measured %alt:", pct)

    seqs = read_fasta(FASTA)
    print("fasta entries:", list(seqs))

    # Map each label to its fasta record.
    keymap = {}
    for lab in anc.label:
        hit = [k for k in seqs if k.replace("_", ".").lower() == lab.lower()]
        if not hit:
            hit = [k for k in seqs if lab.lower().replace(".", "") in
                   k.lower().replace("_", "").replace(".", "")]
        keymap[lab] = hit[0] if hit else None
    print("label -> fasta:", keymap)

    base = {"evo_order": anc.evo_order.to_numpy(float),
            "seq_len": anc.seq_len.to_numpy(float)}

    if all(keymap.values()):
        ref = seqs[keymap["XCL1"]] if "XCL1" in keymap else None
        if ref is not None:
            base["pid_to_XCL1"] = np.array(
                [gapaware_pid(seqs[keymap[l]], ref) for l in anc.label])

    rows = []
    for name, v in base.items():
        r, p, n = exact_p(v, pct)
        rows.append(dict(baseline=name, rho=r, p_exact_greater=p, n_perms=n,
                         values=";".join(f"{x:.4g}" for x in v)))
    out = pd.DataFrame(rows)
    print()
    print(out.to_string(index=False))
    if out_csv:
        out.to_csv(out_csv, index=False)
    return out


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
