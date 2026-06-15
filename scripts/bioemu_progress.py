#!/usr/bin/env python3
"""Per-pair BioEmu progress.

For every fold pair, print the reference sequence length (residues) and
clusters_done/total, where clusters_done = rows in Analysis/df_bioemu.csv and
total = number of ShallowMsa_*.a3m clusters. Sorted least-complete first so
gaps are obvious; a final line gives the overall total.

Usage (from the repo root):
    python3 scripts/bioemu_progress.py
    python3 scripts/bioemu_progress.py --sort len   # sort by length desc instead
"""
import glob
import os
import sys

DATA = "Pipeline/FoldPairs"


def ref_len(pair_dir):
    fs = sorted(glob.glob(os.path.join(pair_dir, "fasta_chain_files", "*.fasta")))
    if not fs:
        return 0
    return sum(len(l.strip()) for l in open(fs[0]) if not l.startswith(">"))


def csv_rows(pair_dir):
    c = os.path.join(pair_dir, "Analysis", "df_bioemu.csv")
    if not os.path.isfile(c):
        return 0
    return max(0, sum(1 for _ in open(c)) - 1)


def n_clusters(pair_dir):
    return len(glob.glob(os.path.join(pair_dir, "output_msa_cluster", "ShallowMsa_*.a3m")))


def main():
    by_len = "--sort" in sys.argv and "len" in sys.argv
    rows = []
    for pair_dir in sorted(glob.glob(os.path.join(DATA, "*"))):
        if not os.path.isdir(pair_dir):
            continue
        tot = n_clusters(pair_dir)
        if tot == 0:
            continue  # not a runnable pair
        rows.append((os.path.basename(pair_dir), ref_len(pair_dir),
                     csv_rows(pair_dir), tot))

    if by_len:
        rows.sort(key=lambda r: -r[1])
    else:  # least-complete first
        rows.sort(key=lambda r: (r[2] / r[3], -r[1]))

    for name, L, done, tot in rows:
        flag = "" if done >= tot else "  <-- incomplete"
        print(f"{name:18s}  L={L:4d}  {done:3d}/{tot:3d}{flag}")

    td = sum(r[2] for r in rows)
    tt = sum(r[3] for r in rows)
    nfull = sum(1 for r in rows if r[2] >= r[3])
    print(f"\nTOTAL: {td}/{tt} clusters scored across {len(rows)} pairs "
          f"({nfull} pairs complete).")


if __name__ == "__main__":
    main()
