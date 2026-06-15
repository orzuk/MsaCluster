#!/usr/bin/env python3
"""Why is each incomplete BioEmu pair missing clusters?

For every pair with done < total, break the shortfall into:
  skip  = clusters skipped for length (medoid > max_len)   ['SKIP query length']
  fail  = clusters that errored                             ['failed for']  (X = 'non-valid protein character' subset)
  UNEXPLAINED = missing - skip - fail
    -> clusters that are neither skipped nor failed in the log: either the job
       was still running / killed before reaching them, or they completed but
       didn't make it into the CSV. THESE are the ones that 'should have more'.
Sorted by UNEXPLAINED (then short L first) so real gaps surface at the top.

Usage:  python3 scripts/bioemu_gap_reasons.py
"""
import glob
import os

DATA = "Pipeline/FoldPairs"


def count(path, needle):
    if not os.path.isfile(path):
        return 0
    n = 0
    with open(path, errors="ignore") as f:
        for line in f:
            if needle in line:
                n += 1
    return n


def ref_len(pd):
    fa = sorted(glob.glob(os.path.join(pd, "fasta_chain_files", "*.fasta")))
    if not fa:
        return 0
    return sum(len(l.strip()) for l in open(fa[0]) if not l.startswith(">"))


rows = []
for pd in sorted(glob.glob(os.path.join(DATA, "*"))):
    if not os.path.isdir(pd):
        continue
    name = os.path.basename(pd)
    tot = len(glob.glob(os.path.join(pd, "output_msa_cluster", "ShallowMsa_*.a3m")))
    if tot == 0:
        continue
    csv = os.path.join(pd, "Analysis", "df_bioemu.csv")
    done = max(0, sum(1 for _ in open(csv)) - 1) if os.path.isfile(csv) else 0
    if done >= tot:
        continue
    log = os.path.join(pd, "jobs", f"run_bioemu_{name}.out")
    skip = count(log, "SKIP query length")
    fail = count(log, "failed for")
    res = count(log, "non-valid protein character")
    L = ref_len(pd)
    miss = tot - done
    unexp = miss - skip - fail
    rows.append((name, L, done, tot, miss, skip, fail, res, unexp))

rows.sort(key=lambda r: (-r[8], r[1]))   # UNEXPLAINED desc, then short L first
for name, L, done, tot, miss, skip, fail, res, unexp in rows:
    flag = "  <== should have more" if unexp > 0 else ""
    print(f"{name:18s} L={L:4d}  {done:3d}/{tot:3d}  miss={miss:3d}  "
          f"skip>500={skip:3d}  fail={fail:3d}(X={res})  UNEXPLAINED={unexp:3d}{flag}")

print(f"\n{len(rows)} incomplete pairs; "
      f"{sum(1 for r in rows if r[8] > 0)} have UNEXPLAINED gaps "
      f"(still-running jobs also show up here until they finish).")
