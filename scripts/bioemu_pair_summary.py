#!/usr/bin/env python3
"""Per-pair BioEmu fold-preference summary from each df_bioemu.csv.

Per pair: n clusters scored, mean frac_toward_f1 (population bias to fold1),
how many clusters lean fold1 vs fold2, mean best-TM to each fold (does the
ensemble actually reach the folds), and the spread of frac across clusters
(high spread = MSA-cluster-dependent fold preference = the fold-switch signal).

Usage:  python3 scripts/bioemu_pair_summary.py
        python3 scripts/bioemu_pair_summary.py --sort spread   # most bimodal first
"""
import csv
import glob
import os
import statistics as st
import sys


def load(path):
    fr, b1, b2 = [], [], []
    with open(path) as f:
        for row in csv.DictReader(f):
            try:
                fr.append(float(row["frac_toward_f1"]))
                b1.append(float(row["best_tm_fold1"]))
                b2.append(float(row["best_tm_fold2"]))
            except (KeyError, ValueError):
                pass
    return fr, b1, b2


def main():
    by = "spread" if ("--sort" in sys.argv and "spread" in sys.argv) else "frac"
    rows = []
    for c in glob.glob("Pipeline/FoldPairs/*/Analysis/df_bioemu.csv"):
        name = os.path.basename(os.path.dirname(os.path.dirname(c)))
        fr, b1, b2 = load(c)
        if not fr:
            continue
        spread = st.pstdev(fr) if len(fr) > 1 else 0.0
        rows.append((name, len(fr), st.mean(fr), spread,
                     sum(x > 0.5 for x in fr), sum(x < 0.5 for x in fr),
                     st.mean(b1), st.mean(b2)))
    rows.sort(key=lambda r: -r[3] if by == "spread" else -r[2])
    print(f"{'pair':18s} {'n':>3} {'meanFrac':>8} {'spread':>6} "
          f"{'F1':>3}/{'F2':>3} {'bTM1':>5} {'bTM2':>5}")
    for name, n, mf, sp, n1, n2, mb1, mb2 in rows:
        print(f"{name:18s} {n:3d} {mf:8.2f} {sp:6.2f} {n1:3d}/{n2:3d} "
              f"{mb1:5.2f} {mb2:5.2f}")
    print(f"\n{len(rows)} pairs. meanFrac~1=fold1, ~0=fold2. "
          f"High spread = clusters disagree = candidate fold-switch signal. "
          f"bTM1/bTM2 = mean best TM-score reached to each fold "
          f"(both low => ensemble doesn't reach either fold).")


if __name__ == "__main__":
    main()
