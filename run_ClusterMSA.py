#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Cluster an A3M into ShallowMsa_###.a3m using one of:
  - HDBSCAN   (--cluster_alg hdbscan)  [existing behavior]
  - TREE      (--cluster_alg tree)     [requires --tree_path *.nwk]
  - AHC       (--cluster_alg ahc)      [default]

All algorithms end with the same post-processing:
  * enforce min_output_size
  * keep at most max_clusters (largest first)
  * compute Neff (identity threshold) and drop clusters with neff < min_neff
  * write ShallowMsa_000.a3m, ShallowMsa_001.a3m, ...

AHC and TREE ensure all sequences are assigned to some cluster.
HDBSCAN keeps your existing sample→cluster→assign heuristic (and can be set to reassign noise).
"""

from __future__ import annotations
import argparse, sys, os, math, random
from scipy.cluster.hierarchy import linkage, fcluster
from pathlib import Path
from typing import List, Tuple, Optional, Dict
import numpy as np
import pandas as pd
import time

# ------------- tiny I/O utils (compatible with your project) -------------
def _now() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")

class StageTimer:
    def __init__(self, name, verbose=True):
        self.name = name
        self.t0 = time.perf_counter()
        self.verbose = verbose
        if self.verbose:
            print(f"[{_now()}] >>> {self.name} ...", flush=True)
    def done(self, note: str = ""):
        t = time.perf_counter() - self.t0
        if self.verbose:
            extra = f" ({note})" if note else ""
            print(f"[{_now()}] <<< {self.name} done in {t:.2f}s{extra}", flush=True)

def _log(msg: str, verbose=True):
    if verbose:
        print(f"[{_now()}] {msg}", flush=True)


def ensure_dir(p: str) -> None:
    Path(p).mkdir(parents=True, exist_ok=True)

def load_a3m(a3m_path: str) -> Tuple[List[str], List[str]]:
    """Read A3M, strip lowercase inserts, keep '-' and uppercase AA."""
    H, S = [], []
    h, cur = None, []
    with open(a3m_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line: continue
            if line.startswith(">"):
                if h is not None:
                    s = "".join(cur)
                    s = "".join(ch for ch in s if ch == '-' or ch.isupper())
                    H.append(h); S.append(s)
                h = line[1:].strip()
                cur = []
            else:
                cur.append(line.strip())
        if h is not None:
            s = "".join(cur)
            s = "".join(ch for ch in s if ch == '-' or ch.isupper())
            H.append(h); S.append(s)
    if not S:
        raise ValueError("Empty A3M.")
    L = len(S[0])
    if any(len(s) != L for s in S):
        raise ValueError("Not all sequences have the same aligned length after cleaning.")
    return H, S

def write_fasta(ids: List[str], seqs: List[str], out_path: str) -> None:
    with open(out_path, "w") as f:
        for h, s in zip(ids, seqs):
            f.write(f">{h}\n{s}\n")

# ------------- distances / medoids / Neff -------------------
def gapaware_pid(a: str, b: str) -> float:
    """PID over positions where both are residues (not '-')."""
    both = [(x, y) for x, y in zip(a, b) if x != '-' and y != '-']
    if not both: return 0.0
    eq = sum(1 for x, y in both if x == y)
    return eq / float(len(both))

def gapaware_hamming(a: str, b: str) -> float:
    return 1.0 - gapaware_pid(a, b)

def condensed_from_square(D: np.ndarray) -> np.ndarray:
    """Square (n,n) -> condensed (n*(n-1)/2,) for SciPy linkage."""
    n = D.shape[0]
    out = np.empty(n*(n-1)//2, dtype=float)
    k = 0
    for i in range(n-1):
        for j in range(i+1, n):
            out[k] = D[i, j]; k += 1
    return out

def medoid_index(D: np.ndarray, idxs: np.ndarray) -> int:
    """Return index (from idxs) of the medoid (min avg distance)."""
    subD = D[np.ix_(idxs, idxs)]
    s = subD.sum(axis=1)
    rel = int(idxs[int(np.argmin(s))])
    return rel

def compute_neff(seqs: List[str], id_thresh: float = 0.8) -> float:
    """Meff with identity threshold id_thresh (simple O(n^2) for cluster sizes ~1e3 ok)."""
    n = len(seqs)
    if n == 0: return 0.0
    weights = np.zeros(n, dtype=float)
    for i in range(n):
        cnt = 0
        for j in range(n):
            if gapaware_pid(seqs[i], seqs[j]) >= id_thresh:
                cnt += 1
        weights[i] = 1.0 / float(cnt)
    return float(weights.sum())

# ------------- HDBSCAN path (reuse your existing behavior) -------------
def cluster_hdbscan(headers: List[str], seqs: List[str], args) -> Dict[int, List[int]]:
    """
    Keep your current sample→cluster→assign approach.
    We build a small distance matrix on a sample, run HDBSCAN, make medoids, and assign all.
    If HDBSCAN yields 0 clusters, we can (optionally) fall back to K=1.
    """
    try:
        import hdbscan  # noqa
    except Exception as e:
        raise RuntimeError(f"HDBSCAN not available: {e}")

    rng = np.random.default_rng(int(args.sample_seed))
    n = len(seqs)
    # Filter by gaps
    L = len(seqs[0])
    keep = np.array([(s.count('-') / L) < args.frac_gaps_cutoff for s in seqs], dtype=bool)
    if keep.sum() < 2:
        raise RuntimeError("Too few sequences after gap filtering.")
    idx_all = np.where(keep)[0]
    headers = [headers[i] for i in idx_all]
    seqs    = [seqs[i]    for i in idx_all]
    n = len(seqs)
    print(f"[INFO] A3M: {len(idx_all)} kept after gap filter (cutoff={args.frac_gaps_cutoff}). Aln len={L}")

    # Sample
    if n > args.sample_cap:
        sample_idx = np.array(sorted(rng.choice(n, size=int(args.sample_cap), replace=False)))
        print(f"[sample] {len(sample_idx)}/{n}")
    else:
        sample_idx = np.arange(n)

    # Build D on sample
    m = len(sample_idx)
    D = np.zeros((m, m), dtype=float)
    for ii, i in enumerate(sample_idx):
        si = seqs[i]
        for jj, j in enumerate(sample_idx[ii+1:], start=ii+1):
            D[ii, jj] = D[jj, ii] = gapaware_hamming(si, seqs[j])

    # HDBSCAN
    import hdbscan
    min_samples = None if (str(args.min_samples).lower() == "auto") else int(args.min_samples)
    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=int(args.min_cluster_size),
        min_samples=min_samples,
        metric='precomputed',
        cluster_selection_method=args.cluster_selection
    ).fit(D)

    labels = clusterer.labels_  # -1 for noise
    uniq = sorted(list(set(int(x) for x in labels if x >= 0)))
    if not uniq:
        print("[WARN] HDBSCAN found 0 clusters on the sample; writing empty metadata (all noise).")
        return {}  # handled upstream (no clusters)

    # medoids on sample clusters
    medoid_sample_rows = []
    for c in uniq:
        rows = np.where(labels == c)[0]
        med = medoid_index(D, rows)
        medoid_sample_rows.append(med)
    medoid_global_idx = sample_idx[np.array(medoid_sample_rows, dtype=int)]

    # assign all sequences to nearest medoid
    C = len(medoid_global_idx)
    centers = [seqs[i] for i in medoid_global_idx]
    assigns: Dict[int, List[int]] = {ci: [] for ci in range(C)}
    for i in range(n):
        si = seqs[i]
        dmin, cmin = 1e9, 0
        for ci, cs in enumerate(centers):
            d = gapaware_hamming(si, cs)
            if d < dmin:
                dmin, cmin = d, ci
        assigns[cmin].append(i)

    # rename cluster IDs 0..C-1 and return membership in local indices
    return assigns

# ------------- AHC path (sample → linkage → medoids → assign all) -------------
def cluster_ahc(headers: List[str], seqs: List[str], args) -> Dict[int, List[int]]:
    rng = np.random.default_rng(int(args.sample_seed))
    verbose = TRUE # bool(int(args.verbose))
    every = max(1, int(args.progress_every))

    st = StageTimer("AHC: gap-filter", verbose)
    L = len(seqs[0])
    keep = np.array([(s.count('-') / L) < args.frac_gaps_cutoff for s in seqs], dtype=bool)
    if keep.sum() < 2:
        raise RuntimeError("Too few sequences after gap filtering.")
    headers = [h for h, k in zip(headers, keep) if k]
    seqs    = [s for s, k in zip(seqs, keep) if k]
    n = len(seqs)
    st.done(f"kept={n}, L={L}, cutoff={args.frac_gaps_cutoff}")

    # sampling
    st = StageTimer("AHC: sampling", verbose)
    if n > args.sample_cap:
        sample_idx = np.array(sorted(rng.choice(n, size=int(args.sample_cap), replace=False)))
        _log(f"[sample] {len(sample_idx)}/{n}", verbose)
    else:
        sample_idx = np.arange(n)
        _log(f"[sample] using all {n} sequences", verbose)
    m = len(sample_idx)
    # quick memory estimate for dense sample distances
    bytes_est = m*m*8
    _log(f"[estimate] distance matrix {m}x{m} ~ {bytes_est/1e9:.2f} GB (float64)", verbose)
    st.done()

    # distance on sample with progress
    st = StageTimer("AHC: compute distance matrix (gap-aware Hamming)", verbose)
    D = np.zeros((m, m), dtype=float)
    for ii, i in enumerate(sample_idx):
        si = seqs[i]
        if verbose and (ii % every == 0 or ii == m-1):
            _log(f"[dist] row {ii+1}/{m}", verbose)
        for jj, j in enumerate(sample_idx[ii+1:], start=ii+1):
            D[ii, jj] = D[jj, ii] = gapaware_hamming(si, seqs[j])
    st.done()

    # Build condensed distances directly for speed/memory
    cond = np.empty(m * (m - 1) // 2, dtype=float)
    k = 0
    for ii in range(m - 1):
        si = seqs[sample_idx[ii]]
        if bool(int(args.verbose)) and (ii % max(1, int(args.progress_every)) == 0 or ii == m - 2):
            print(f"[{_now()}] [dist] row {ii + 1}/{m}", flush=True)
        for jj in range(ii + 1, m):
            sj = seqs[sample_idx[jj]]
            cond[k] = gapaware_hamming(si, sj);
            k += 1

    # Linkage
    method = args.ahc_linkage
    if '_HAS_FASTCLUSTER' in globals() and _HAS_FASTCLUSTER:
        Z = _fastcluster.linkage_vector(cond, method=method)
    else:
        from scipy.cluster.hierarchy import linkage
        Z = linkage(cond, method=method)

    # Flat cut
    from scipy.cluster.hierarchy import fcluster
    if args.ahc_cut_mode == "maxclust":
        lab = fcluster(Z, t=int(args.max_clusters), criterion="maxclust") - 1
    else:
        thr = float(args.ahc_distance_threshold)
        if thr <= 0.0:
            # auto-guess a conservative threshold from the condensed dist quantiles
            q50, q75 = np.quantile(cond, [0.50, 0.75])
            thr = (q50 + q75) / 2.0
            print(f"[AHC] Using auto distance threshold thr={thr:.4f}", flush=True)
        lab = fcluster(Z, t=thr, criterion="distance") - 1
    C = int(lab.max() + 1) if m > 0 else 0
    print(f"[cut] clusters on sample: {C}", flush=True)

    # Merge small clusters to satisfy min_output_size on the sample
    min_size = int(args.min_output_size)
    st = StageTimer(f"AHC: merge sample clusters < {min_size}", verbose)

    def cluster_rows(c): return np.where(lab == c)[0]

    def recompute_medoids(active: List[int]) -> Dict[int, int]:
        med = {}
        for c in active:
            rows = cluster_rows(c)
            if len(rows) == 0: continue
            med[c] = medoid_index(D, rows)
        return med

    active = list(range(C))
    med = recompute_medoids(active)

    def nearest_cluster(c_from, cands):
        dmin, cbest = 1e9, None
        for c in cands:
            if c == c_from: continue
            d = D[med[c_from], med[c]]
            if d < dmin:
                dmin, cbest = d, c
        return cbest

    merges = 0
    while True:
        sizes = {c: int((lab == c).sum()) for c in active}
        small = [c for c in active if sizes[c] < min_size]
        if not small:
            break
        for c in small:
            other = [x for x in active if x != c]
            if not other:
                continue
            cn = nearest_cluster(c, other)
            if cn is None:
                continue
            lab[lab == c] = cn
            active.remove(c)
            med = recompute_medoids(active)
            merges += 1
        _log(f"[merge] merges so far: {merges}, active clusters: {len(active)}", verbose)
    st.done(f"final sample clusters={len(active)}")

    # medoids of final sample clusters
    st = StageTimer("AHC: medoids (sample clusters)", verbose)
    uniq = sorted(set(int(x) for x in lab))
    med_sample_rows = [medoid_index(D, np.where(lab == c)[0]) for c in uniq]
    medoid_global_idx = sample_idx[np.array(med_sample_rows, dtype=int)]
    st.done(f"centers={len(medoid_global_idx)}")

    # assign all sequences
    st = StageTimer("AHC: assign all sequences to nearest medoid", verbose)
    centers = [seqs[i] for i in medoid_global_idx]
    assigns: Dict[int, List[int]] = {ci: [] for ci in range(len(centers))}
    for i in range(n):
        if verbose and (i % (every*10) == 0 or i == n-1):
            _log(f"[assign] seq {i+1}/{n}", verbose)
        si = seqs[i]
        dmin, cmin = 1e9, 0
        for ci, cs in enumerate(centers):
            d = gapaware_hamming(si, cs)
            if d < dmin:
                dmin, cmin = d, ci
        assigns[cmin].append(i)
    st.done()

    # normalize cluster ids
    mapping = {old: new for new, old in enumerate(sorted(assigns.keys()))}
    assigns = {mapping[k]: v for k, v in assigns.items()}
    _log(f"[assign] produced {len(assigns)} clusters (pre output filters)", verbose)
    return assigns


# ------------- TREE path (monophyletic flat cut with size/Neff constraints) -------------
def cluster_tree(headers: List[str], seqs: List[str], args) -> Dict[int, List[int]]:
    """
    Cut an input phylogenetic tree (Newick) into monophyletic clusters such that
    each cluster has >= min_output_size; small clades are merged upward.
    If tree is missing, abort with a clear message.
    """
    if not args.tree_path or not os.path.isfile(args.tree_path):
        raise SystemExit("Missing tree! Run tree reconstruction first to use --cluster_alg tree, or choose hdbscan/ahc.")
    from Bio import Phylo

    name_to_idx: Dict[str, int] = {}
    for i, h in enumerate(headers):
        # normalize like your code: take first token before space; strip trailing /start-end
        tok = (h or "").split()[0]
        tok = tok.split("/")[0] if "/" in tok else tok
        name_to_idx[tok] = i

    tree = Phylo.read(args.tree_path, "newick")

    # get leaves and map to indices; we only keep leaves that we can map
    def clade_indices(clade) -> List[int]:
        ids = []
        for leaf in clade.get_terminals():
            nm = (leaf.name or "").strip().strip("'").strip('"')
            nm = nm.split()[0]
            nm = nm.split("/")[0] if "/" in nm else nm
            if nm in name_to_idx:
                ids.append(name_to_idx[nm])
        return ids

    # postorder traversal; accept clade if it is big enough; else merge upward
    accepted: List[List[int]] = []

    def walk(clade) -> List[int]:
        if clade.is_terminal():
            idxs = clade_indices(clade)
            return idxs
        child_sets = [walk(c) for c in clade.clades]
        merged = [i for sub in child_sets for i in sub]
        # accept this clade if large enough; else bubble up
        if len(merged) >= int(args.min_output_size):
            accepted.append(merged)
            return []  # consumed here
        else:
            return merged

    leftover = walk(tree.root)
    if leftover:
        # attach leftover to the largest accepted cluster (or create one)
        if accepted:
            j = int(np.argmax([len(x) for x in accepted]))
            accepted[j].extend(leftover)
        else:
            accepted = [leftover]

    # if too many clusters, keep largest max_clusters and merge the rest into nearest by medoid
    maxC = int(args.max_clusters)
    if len(accepted) > maxC:
        # compute medoid for each, then greedily merge smallest into nearest until count == maxC
        # build a quick distance cache
        Dfull = None  # build lazily when needed
        def dist(i, j):
            nonlocal Dfull
            if Dfull is None:
                n = len(seqs)
                Dfull = np.zeros((n, n), dtype=float)
                for a in range(n-1):
                    sa = seqs[a]
                    for b in range(a+1, n):
                        Dfull[a, b] = Dfull[b, a] = gapaware_hamming(sa, seqs[b])
            return Dfull[i, j]

        def medoid_of(idxs):
            M = np.zeros((len(idxs), len(idxs)))
            for ii, gi in enumerate(idxs):
                for jj, gj in enumerate(idxs[ii+1:], start=ii+1):
                    M[ii, jj] = M[jj, ii] = dist(gi, gj)
            s = M.sum(axis=1)
            return idxs[int(np.argmin(s))]

        groups = [list(sorted(set(g))) for g in accepted]
        while len(groups) > maxC:
            sizes = [len(g) for g in groups]
            smin = int(np.argmin(sizes))
            m_from = medoid_of(groups[smin])
            # find nearest target
            best, bestd = None, 1e9
            for k, g in enumerate(groups):
                if k == smin: continue
                m_to = medoid_of(g)
                d = dist(m_from, m_to)
                if d < bestd:
                    best, bestd = k, d
            groups[best].extend(groups[smin])
            groups.pop(smin)
        accepted = groups

    # return membership dict over local indices
    assigns: Dict[int, List[int]] = {i: sorted(set(g)) for i, g in enumerate(accepted)}
    return assigns

# ------------- common post-processing and writing ------------------------
def postprocess_and_write(assigns: Dict[int, List[int]],
                          headers: List[str], seqs: List[str], args) -> pd.DataFrame:
    """
    Enforce min_output_size, cap at max_clusters, compute Neff and drop low-Neff,
    then write ShallowMsa_###.a3m. Return a per-cluster dataframe.
    """
    verbose = TRUE #  bool(int(args.verbose))
    st = StageTimer("Postprocess & write clusters", verbose)

    if not assigns:
        ensure_dir(args.outdir)
        meta = pd.DataFrame(columns=["cluster", "n", "neff", "path", "kept"], dtype=object)
        meta.to_csv(Path(args.outdir, f"{args.keyword}_clusters.csv"), index=False)
        _log("[WARN] No clusters produced; wrote empty metadata.", verbose)
        st.done()
        return meta

    # drop tiny clusters (< min_output_size)
    min_size = int(args.min_output_size)
    assigns = {c: ids for c, ids in assigns.items() if len(ids) >= min_size}
    if not assigns:
        print("[WARN] All clusters were < min_output_size; nothing to write.")
        return pd.DataFrame(columns=["cluster","n","neff","path","kept"])

    # cap at max_clusters (largest first)
    maxC = int(args.max_clusters)
    order = sorted(assigns.keys(), key=lambda c: len(assigns[c]), reverse=True)
    keep_keys = order[:maxC]
    assigns = {i_new: assigns[k] for i_new, k in enumerate(keep_keys)}

    # write & compute Neff
    ensure_dir(args.outdir)
    rows = []
    for c, idxs in assigns.items():
        ids = [headers[i] for i in idxs]
        ss  = [seqs[i]    for i in idxs]
        neff = compute_neff(ss, id_thresh=float(args.neff_id_thresh))
        outp = Path(args.outdir, f"{args.keyword}_{c:03d}.a3m")
        write_fasta(ids, ss, str(outp))
        rows.append(dict(cluster=c, n=len(idxs), neff=neff, path=str(outp), kept=True))

    df = pd.DataFrame(rows).sort_values(["kept","n"], ascending=[False, False])

    # drop low-Neff if requested
    min_neff = int(args.min_neff)
    if min_neff > 0:
        df["kept"] = df["neff"] >= float(min_neff)
        if (df["kept"] == False).any():
            # remove files of dropped clusters
            for p in df.loc[~df["kept"], "path"].tolist():
                try: os.remove(p)
                except Exception: pass

    # write metadata CSV
    meta_path = Path(args.outdir, f"{args.keyword}_clusters.csv")
    df.to_csv(meta_path, index=False)
    st.done(f"kept={(df['kept']==True).sum() if not df.empty else 0}")
    return df

# ------------- main dispatch --------------------------------------------
def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Cluster an A3M into ShallowMsa_###.a3m using HDBSCAN / TREE / AHC.")
    p.add_argument("--a3m", required=True, help="Input A3M path.")
    p.add_argument("-o", "--outdir", required=True, help="Output directory.")
    p.add_argument("--keyword", default="ShallowMsa", help="Prefix for cluster files.")
    p.add_argument("--cluster_alg", choices=["hdbscan","tree","ahc"], default="ahc",
                   help="Clustering algorithm to use (default: ahc).")
    p.add_argument("--tree_path", default=None, help="Newick file for --cluster_alg tree.")
    p.add_argument("--metric", choices=["hamming"], default="hamming", help="Distance metric (currently: hamming).")
    # shared thresholds
    p.add_argument("--max_clusters", type=int, default=100, help="Maximum clusters to output.")
    p.add_argument("--min_output_size", type=int, default=200, help="Minimum sequences per cluster (size filter).")
    p.add_argument("--min_neff", type=int, default=50, help="Minimum Neff per cluster (drop if below).")
    p.add_argument("--neff_id_thresh", type=float, default=0.8, help="Identity threshold for Neff.")
    p.add_argument("--frac_gaps_cutoff", type=float, default=0.6, help="Drop seqs with gap fraction >= this.")
    # sampling / seeds (HDBSCAN & AHC)
    p.add_argument("--sample_cap", type=int, default=5000, help="Sample size for building cluster cores.")
    p.add_argument("--sample_seed", type=int, default=12345, help="Random seed for sampling.")
    # HDBSCAN knobs (kept for compatibility)
    p.add_argument("--min_cluster_size", type=int, default=200, help="HDBSCAN min_cluster_size.")
    p.add_argument("--min_samples", default="auto", help="'auto' or integer for HDBSCAN density strictness.")
    p.add_argument("--cluster_selection", choices=["eom","leaf"], default="eom",
                   help="HDBSCAN cluster selection method; eom merges leaves to stable parents.")

    # AHC knobs
    p.add_argument("--ahc_linkage", choices=["average", "complete", "single"], default="average",
                   help="Linkage for AHC (default: average).")
    p.add_argument("--ahc_cut_mode", choices=["maxclust", "distance"], default="maxclust",
                   help="How to cut the AHC tree: maxclust or distance.")
    p.add_argument("--ahc_distance_threshold", type=float, default=0.0,
                   help="If --ahc_cut_mode=distance, cut at this distance threshold (0.0 means auto-guess).")

    # verbosity printing
    p.add_argument("--verbose", type=int, default=1, help="1=log stages (default), 0=quiet")
    p.add_argument("--progress_every", type=int, default=500,
                   help="Print progress every N rows while building distances.")

    return p

def main():
    args = build_argparser().parse_args()

    headers, seqs = load_a3m(args.a3m)

    if args.cluster_alg == "hdbscan":
        assigns = cluster_hdbscan(headers, seqs, args)
    elif args.cluster_alg == "ahc":
        assigns = cluster_ahc(headers, seqs, args)
    elif args.cluster_alg == "tree":
        assigns = cluster_tree(headers, seqs, args)
    else:
        raise SystemExit(f"Unknown --cluster_alg {args.cluster_alg}")

    df = postprocess_and_write(assigns, headers, seqs, args)
    kept = df[df.kept] if not df.empty else df
    dropped = df[~df.kept] if not df.empty else df
    print(f"[SUMMARY] Kept {0 if df.empty else len(kept)} clusters, Dropped {0 if df.empty else len(dropped)} (by Neff).")
    if not df.empty and not kept.empty:
        print(kept[["cluster","n","neff","path"]].to_string(index=False))

if __name__ == "__main__":
    main()
