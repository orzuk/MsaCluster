#!/usr/bin/env python3
"""Structural spread of the per-cluster ESMFold2 predictions.

For every target (a fold-switch pair or a single-fold control) this computes
TM-align over all pairs of per-cluster predicted structures and summarises how
far apart the clusters land. Fold switchers should spread further than controls
if the per-cluster predictions carry fold information rather than noise.

The full pairwise matrix is saved per target, so that statistics sensitive to
how many clusters a target has (above all the minimum, an extreme-value
statistic whose expectation falls as the number of compared pairs grows) can be
recomputed later on a fixed number of clusters without rerunning TM-align.

Usage: esmfold2_diversity.py OUT_CSV SHARD NSHARD [MATRIX_DIR]
"""
import os, re, subprocess, sys, itertools
import numpy as np

TMALIGN = "/sci/labs/orzuk/orzuk/software/tmalign/TMalign"
ROOT = "/sci/labs/orzuk/orzuk/github/MsaCluster/Pipeline"
OUT, SHARD, NSHARD = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
MATDIR = sys.argv[4] if len(sys.argv) > 4 else \
    "/sci/labs/orzuk/orzuk/tmp/esm2div/tm_matrix"
os.makedirs(MATDIR, exist_ok=True)

targets = []
for grp in sorted(os.listdir(os.path.join(ROOT, "Controls"))):
    gdir = os.path.join(ROOT, "Controls", grp)
    if not os.path.isdir(gdir):
        continue
    for name in sorted(os.listdir(gdir)):
        d = os.path.join(gdir, name, "output_esmfold2")
        if os.path.isdir(d):
            targets.append((name, grp, d))
for name in sorted(os.listdir(os.path.join(ROOT, "FoldPairs"))):
    d = os.path.join(ROOT, "FoldPairs", name, "output_esmfold2")
    if os.path.isdir(d):
        targets.append((name, "FoldSwitch", d))

TM_RE = re.compile(r"TM-score=\s*([0-9.]+)")

def tm(a, b):
    try:
        r = subprocess.run([TMALIGN, a, b], capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        return None
    v = TM_RE.findall(r.stdout)
    return max(float(x) for x in v[:2]) if len(v) >= 2 else None

with open(OUT, "w") as fh:
    fh.write("target,group,n_clusters,n_pairs,mean_tm,median_tm,min_tm,p10_tm,frac_below_0.5\n")
    for i, (name, grp, d) in enumerate(targets):
        if i % NSHARD != SHARD:
            continue
        cifs = sorted(os.path.join(d, f) for f in os.listdir(d) if f.endswith(".cif"))
        if len(cifs) < 2:
            continue
        n = len(cifs)
        # NaN marks a pair TM-align could not score; it is skipped everywhere.
        mat = np.full((n, n), np.nan, dtype=np.float32)
        np.fill_diagonal(mat, 1.0)
        vals = []
        for a, b in itertools.combinations(range(n), 2):
            v = tm(cifs[a], cifs[b])
            if v is None:
                continue
            mat[a, b] = mat[b, a] = v
            vals.append(v)
        if not vals:
            continue
        np.savez_compressed(
            os.path.join(MATDIR, f"{grp}__{name}.npz"),
            tm=mat, group=grp, target=name,
            files=np.array([os.path.basename(c) for c in cifs]))
        v = np.array(vals)
        fh.write(f"{name},{grp},{n},{len(v)},{v.mean():.4f},"
                 f"{np.median(v):.4f},{v.min():.4f},{np.quantile(v, 0.1):.4f},"
                 f"{float((v < 0.5).mean()):.4f}\n")
        fh.flush()
        print(f"[{SHARD}] {name} ({grp}) {n} clusters mean TM {v.mean():.3f}", flush=True)
