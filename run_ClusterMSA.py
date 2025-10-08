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

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, ROOT)

from utils.utils import ensure_dir, now
from utils.msa_utils import clean_a3m_line, write_fasta, compute_neff_from_a3m

# ------------- tiny I/O utils (compatible with your project) -------------
def now() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")

class StageTimer:
    def __init__(self, name, verbose=True):
        self.name = name
        self.t0 = time.perf_counter()
        self.verbose = verbose
        if self.verbose:
            print(f"[{now()}] >>> {self.name} ...", flush=True)
    def done(self, note: str = ""):
        t = time.perf_counter() - self.t0
        if self.verbose:
            extra = f" ({note})" if note else ""
            print(f"[{now()}] <<< {self.name} done in {t:.2f}s{extra}", flush=True)

def _log(msg: str, verbose=True):
    if verbose:
        print(f"[{now()}] {msg}", flush=True)


def _cond_idx(i,j,m):
    if i==j: return None
    if i>j: i,j = j,i
    # SciPy condensed index
    return i*m - i*(i+1)//2 + (j - i - 1)

def medoid_of_rows(rows: np.ndarray, seqs: list[str], idx_map: np.ndarray | None = None) -> int:
    if len(rows) == 1: return int(rows[0])
    r = rows.tolist(); R = len(r)
    Msum = np.zeros(R, dtype=float)
    for a in range(R-1):
        sa = seqs[sample_idx[r[a]]]
        for b in range(a+1, R):
            sb = seqs[sample_idx[r[b]]]
            d = gapaware_hamming(sa, sb)
            Msum[a] += d; Msum[b] += d
    return int(r[int(np.argmin(Msum))])


def medoid_of_rows_condensed(rows: np.ndarray, m: int, cond_vec: np.ndarray) -> int:
    r = rows.astype(int)
    R = len(r)
    if R == 1: return int(r[0])
    sums = np.zeros(R, dtype=float)
    for a in range(R-1):
        ia = r[a]
        for b in range(a+1, R):
            ib = r[b]
            k = _cond_idx(ia, ib, m)
            d = cond_vec[k]
            sums[a] += d; sums[b] += d
    return int(r[int(np.argmin(sums))])


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
    verbose = bool(int(getattr(args, "verbose", 1)))
    every = max(1, int(getattr(args, "progress_every", 500)))

    # --- Gap filter ---
    st = StageTimer("AHC: gap-filter", verbose)
    L = len(seqs[0])
    keep = np.array([(s.count('-') / L) < args.frac_gaps_cutoff for s in seqs], dtype=bool)
    if keep.sum() < 2: raise RuntimeError("Too few sequences after gap filtering.")
    seqs    = [s for s,k in zip(seqs, keep) if k]
    n = len(seqs)
    st.done(f"kept={n}, L={L}, cutoff={args.frac_gaps_cutoff}")

    # --- Sample ---
    st = StageTimer("AHC: sampling", verbose)
    if n > args.sample_cap:
        sample_idx = np.array(sorted(rng.choice(n, size=int(args.sample_cap), replace=False)))
        _log(f"[sample] {len(sample_idx)}/{n}", verbose)
    else:
        sample_idx = np.arange(n)
        _log(f"[sample] using all {n} sequences", verbose)
    m = len(sample_idx)
    _log(f"[estimate] condensed vector size ~ {m*(m-1)//2:,} distances", verbose)
    st.done()

    # --- Condensed distances ---
    st = StageTimer("AHC: compute condensed distances (gap-aware Hamming)", verbose)
    cond = np.empty(m*(m-1)//2, dtype=float)
    k = 0
    for ii in range(m-1):
        si = seqs[sample_idx[ii]]
        if verbose and (ii % every == 0 or ii == m-2):
            _log(f"[dist] row {ii+1}/{m}", verbose)
        for jj in range(ii+1, m):
            sj = seqs[sample_idx[jj]]
            cond[k] = gapaware_hamming(si, sj); k += 1
    st.done()

    # --- Linkage ---
    method = getattr(args, "ahc_linkage", "average")
    st = StageTimer(f"AHC: linkage ({method})", verbose)
    if '_HAS_FASTCLUSTER' in globals() and _HAS_FASTCLUSTER:
        Z = _fastcluster.linkage_vector(cond, method=method)
    else:
        from scipy.cluster.hierarchy import linkage
        Z = linkage(cond, method=method)
    st.done()

    # --- Flat cut ---
    from scipy.cluster.hierarchy import fcluster
    st = StageTimer("AHC: flat cut", verbose)
    cut_mode = getattr(args, "ahc_cut_mode", "maxclust")
    if cut_mode == "maxclust":
        lab = fcluster(Z, t=int(args.max_clusters), criterion="maxclust") - 1
    else:
        thr = float(getattr(args, "ahc_distance_threshold", 0.0))
        if thr <= 0.0:
            q50, q75 = np.quantile(cond, [0.50, 0.75]); thr = (q50 + q75) / 2.0
            _log(f"[AHC] Using auto distance threshold thr={thr:.4f}", verbose)
        lab = fcluster(Z, t=thr, criterion="distance") - 1
    C = int(lab.max() + 1) if m > 0 else 0
    _log(f"[cut] clusters on sample: {C}", verbose)
    st.done()

    # --- Helpers (no full D) ---
    def rows_of(c): return np.where(lab == c)[0]
    def center_seq_of_c(c:int)->str:
        mr = medoid_of_rows(rows_per[c]); return seqs[sample_idx[mr]]
    def dist_c(c1, c2):
        return gapaware_hamming(center_seq_of_c(c1), center_seq_of_c(c2))

    # --- Build per-cluster bookkeeping on sample ---
    active = list(range(C))
    rows_per = {c: rows_of(c) for c in active}
    size_of  = {c: int(rows_per[c].size) for c in active}

    # --- (1) Remove *micro* centers before assignment (fast O(K^2)) ---
    min_center = int(getattr(args, "ahc_min_center_size", 0))
    if min_center > 0:
        st = StageTimer(f"AHC: merge sample clusters < {min_center} (pre-assign, cached)", verbose)
        active = sorted(set(int(x) for x in lab))
        rows_per = {c: np.where(lab == c)[0] for c in active}
        size_of = {c: int(rows_per[c].size) for c in active}

        # one medoid per cluster (on sample rows)
        centers_seq = {}
        i = 0
        for c in active:
            if i % 10 == 0:
                _log(f"[compute-kxk-cluster-centers] Cluster i={i} out of {len(active)}; Cluster label={c}", verbose)
            i += 1
            r = rows_per[c]
            centers_seq[c] = seqs[sample_idx[medoid_of_rows(r)]]

        # KxK center distance matrix
        act_list = list(active)
        idx_of = {c: i for i, c in enumerate(act_list)}
        K = len(act_list)
        CD = np.zeros((K, K), dtype=float)
        for i, c1 in enumerate(act_list):
            if i%10 == 0:
                _log(f"[compute-kxk-distance-metric] Cluster i={i} out of {len(act_list)}; Cluster label={c1}", verbose)
            s1 = centers_seq[c1]
            for j in range(i + 1, K):
                c2 = act_list[j]
                CD[i, j] = CD[j, i] = gapaware_hamming(s1, centers_seq[c2])

        # freeze clusters once they reach >= min_center
        frozen = {c for c in active if size_of[c] >= min_center}
        merges = 0
        for c in list(active):
            if c in frozen or size_of[c] >= min_center or c not in active:
                continue
            i = idx_of[c]
            # prefer non-frozen targets first
            cand = [x for x in active if x != c and x not in frozen]
            if not cand:
                cand = [x for x in active if x != c]
            if not cand: continue
            # nearest by cached distances
            cn = min(cand, key=lambda t: CD[i, idx_of[t]])
            lab[lab == c] = cn
            size_of[cn] += size_of[c]
            size_of[c] = 0
            active.remove(c)
            merges += 1
            if size_of[cn] >= min_center:
                frozen.add(cn)
            if verbose and merges % 10 == 0:
                _log(f"[pre-assign merge] merges={merges}; active={len(active)}", verbose)

        _log(f"[pre-assign] centers after micro-merge: {len(active)} (merges={merges})", verbose)
        st.done()
    else:
        _log("[pre-assign] micro-merge disabled (ahc_min_center_size=0).", verbose)

    # --- (2) Optional: sample-merge to min_output_size (small↔small first) ---
    min_size = int(args.min_output_size)
    if int(getattr(args, "ahc_merge_on_sample", 0)) == 1:
        st = StageTimer(f"AHC: merge sample clusters < {min_size} (small↔small first)", verbose)
        frozen = set()  # clusters that reached ≥ min_size
        changed = True
        while changed:
            changed = False
            small = [c for c in active if c not in frozen and size_of[c] < min_size]
            if not small: break
            for c in small:
                if c not in active or c in frozen: continue
                # prefer small neighbors not yet frozen; else any neighbor
                cand = [x for x in small if x != c and x in active]
                if not cand: cand = [x for x in active if x != c]
                if not cand: continue
                cn = min(cand, key=lambda t: dist_c(c, t))
                lab[lab == c] = cn
                rows_per[cn] = rows_of(cn)
                size_of[cn]  = int(rows_per[cn].size)
                active.remove(c); rows_per.pop(c, None); size_of.pop(c, None)
                # freeze target if big enough now
                if size_of[cn] >= min_size: frozen.add(cn)
                changed = True
        _log(f"[sample-merge] centers after ≥{min_size}: {len(active)}", verbose)
        st.done()
    else:
        _log("[AHC] Skipping sample-merge to min_output_size; will enforce after assignment.", verbose)

    # --- (3) Cap number of centers used in assignment ---
    cap = int(getattr(args, "ahc_center_cap", 500))
    if cap > 0 and len(active) > cap:
        st = StageTimer(f"AHC: center cap → merge smallest into nearest until ≤ {cap}", verbose)
        # min-heap by size for efficiency
        import heapq
        heap = [(size_of[c], c) for c in active]
        heapq.heapify(heap)
        while len(active) > cap:
            while True:
                if not heap: break
                sz, c = heapq.heappop(heap)
                if c in active and size_of[c] == sz: break
            if not heap: break
            # find nearest neighbor among active\{c}
            cand = [x for x in active if x != c]
            if not cand: break
            cn = min(cand, key=lambda t: dist_c(c, t))
            lab[lab == c] = cn
            rows_per[cn] = rows_of(cn)
            size_of[cn]  = int(rows_per[cn].size)
            active.remove(c); rows_per.pop(c, None); size_of.pop(c, None)
            heapq.heappush(heap, (size_of[cn], cn))
        _log(f"[cap] centers after cap: {len(active)}", verbose)
        st.done()
    else:
        _log(f"[cap] centers used for assignment: {len(active)}", verbose)

    # --- Final centers (medoids) on sample clusters ---
    st = StageTimer("AHC: medoids (sample clusters)", verbose)
    uniq = sorted(set(int(x) for x in lab))
    med_sample_rows = [medoid_of_rows_condensed(np.where(lab == c)[0], m, cond) for c in uniq]

    medoid_global_idx = sample_idx[np.array(med_sample_rows, dtype=int)]
    st.done(f"centers={len(medoid_global_idx)}")

    # --- Assign ALL sequences to nearest medoid ---
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

    # --- Normalize ids ---
    mapping = {old: new for new, old in enumerate(sorted(assigns.keys()))}
    assigns = {mapping[k]: v for k, v in assigns.items()}
    _log(f"[assign] produced {len(assigns)} clusters (pre output filters)", verbose)
    return assigns


# ------------- TREE path (monophyletic flat cut with size/Neff constraints) -------------
def cluster_tree(headers: List[str], seqs: List[str], args) -> Dict[int, List[int]]:
    """
    Cut an input phylogenetic tree (Newick) into monophyletic clusters (subtrees),
    meeting size/branch constraints and respecting long internal branches.
    """
    if not args.tree_path or not os.path.isfile(args.tree_path):
        raise SystemExit("Missing tree! Run with --tree_path *.nwk or choose --cluster_alg ahc/hdbscan.")

    from Bio import Phylo
    verbose = bool(int(getattr(args, "verbose", 1)))

    # Map header names -> sequence indices (normalize like your code)
    name_to_idx: Dict[str, int] = {}
    for i, h in enumerate(headers):
        tok = (h or "").split()[0]
        tok = tok.split("/")[0] if "/" in tok else tok
        name_to_idx[tok] = i

    tree = Phylo.read(args.tree_path, "newick")

    min_size = int(args.min_output_size)
    max_clusters = int(args.max_clusters)
    min_tot_br = float(getattr(args, "tree_min_total_branch", 0.0))
    max_size = int(getattr(args, "tree_max_size", 0))
    max_tot_br = float(getattr(args, "tree_max_total_branch", 0.0))
    split_min_ib = float(getattr(args, "tree_split_min_internal_branch", 0.0))

    # --- Relative size cap: aim for at least K clusters (K = tree_min_num_clusters) ---
    K = int(getattr(args, "tree_min_num_clusters", 5))
    # Count how many leaves in the tree map to our headers
    all_idx = []
    for leaf in tree.get_terminals():
        nm = (leaf.name or "").strip().strip("'").strip('"')
        nm = nm.split()[0]
        nm = nm.split("/")[0] if "/" in nm else nm
        if nm in name_to_idx:
            all_idx.append(name_to_idx[nm])
    total_n = len(all_idx)

    # Desired relative max size by total size / K
    rel_max_size = (total_n // K) if K > 0 else 0

    # Final effective max size rule:
    # - If rel_max_size < min_size, disable the size cap (so we can return <K clusters if needed)
    # - Otherwise use rel_max_size; if an absolute max_size>0 exists, take the MIN of the two.
    if rel_max_size > 0 and rel_max_size >= min_size:
        eff_max_size = rel_max_size
        if max_size > 0:
            eff_max_size = min(eff_max_size, max_size)
    else:
        eff_max_size = 0  # disable size cap (allows fewer clusters than K when data is small)

    if verbose:
        print(f"[TREE] total_n={total_n}, K={K}, rel_max_size={rel_max_size}, "
              f"min_output_size={min_size}, eff_max_size={eff_max_size} (0=disabled)")

    # Utilities to compute subtree data fast in a single postorder pass
    def leaf_idxs(clade) -> List[int]:
        ids = []
        for leaf in clade.get_terminals():
            nm = (leaf.name or "").strip().strip("'").strip('"')
            nm = nm.split()[0]
            nm = nm.split("/")[0] if "/" in nm else nm
            if nm in name_to_idx:
                ids.append(name_to_idx[nm])
        return ids

    accepted: List[List[int]] = []

    def postorder(clade):
        """
        Returns: (idxs, total_branch_len, longest_internal_branch_here, accepted_here)
        total_branch_len counts all descendant branch lengths within this clade.
        longest_internal_branch_here = max(branch_length of immediate children)
        """
        if clade.is_terminal():
            idxs = leaf_idxs(clade)
            # terminal: no internal branches beneath
            return idxs, 0.0, 0.0

        child_data = [postorder(c) for c in clade.clades]
        child_idxs = [x[0] for x in child_data]
        child_totL = [x[1] for x in child_data]
        # immediate children branch lengths (internal edges emanating from this node)
        child_bl = [(c.branch_length or 0.0) for c in clade.clades]

        merged = []
        for ids in child_idxs:
            merged.extend(ids)

        total_len = float(sum(child_totL) + sum(child_bl))
        longest_ib = float(max(child_bl) if child_bl else 0.0)

        # Decide acceptance at this node
        n = len(merged)
        too_large = (eff_max_size > 0 and n > eff_max_size) or (max_tot_br > 0.0 and total_len > max_tot_br)
        has_min_branch = (min_tot_br <= 0.0) or (total_len >= min_tot_br)

        # Prefer splitting at long internal branches when both children are sizeable
        children_big_enough = all(len(ids) >= min_size for ids in child_idxs) if child_idxs else False
        prefer_split = (split_min_ib > 0.0 and longest_ib >= split_min_ib and children_big_enough)

        if (n >= min_size) and has_min_branch and not too_large and not prefer_split:
            # accept this clade and do not let ancestors reuse these leaves
            accepted.append(merged)
            return [], 0.0, 0.0  # consumed here
        else:
            # bubble up leaves; children may still be accepted lower down
            return merged, total_len, longest_ib

    leftovers, _, _ = postorder(tree.root)

    # If anything bubbled to the top, attach to the largest accepted cluster or keep as its own if big enough
    if leftovers:
        if len(leftovers) >= min_size and (min_tot_br <= 0.0):
            accepted.append(leftovers)
        elif accepted:
            j = int(np.argmax([len(x) for x in accepted]))
            accepted[j].extend(leftovers)
        else:
            accepted = [leftovers]

    # Map cluster index -> a leaf name present in the tree (for tree.distance)
    def any_leaf_name(idxs):
        # find a header (name) that exists in the tree; we used normalized names for name_to_idx
        for i in idxs:
            # original header token (normalized)
            nm = headers[i].split()[0]
            nm = nm.split("/")[0] if "/" in nm else nm
            if nm in (leaf.name for leaf in tree.get_terminals()):
                return nm
        # as fallback, return the first header token (even if not in tree)
        return headers[idxs[0]].split()[0]

    rep_name_of = {i: any_leaf_name(g) for i, g in enumerate(accepted)}

    print(f"[counts] pre-postprocess clusters: {len(accepted)}  total_assigned={sum(len(g) for g in accepted)}")

    # If too many clusters, merge smallest into nearest by representative PID,
    # but avoid merging clusters whose MRCA separation is "long" (>= split_min_ib)
    if len(accepted) > max_clusters:
        # Build representatives (medoids) to define "nearest"
        rng = np.random.default_rng(int(getattr(args, "sample_seed", 0)))
        REP_CAP = 200

        def rep_of(idxs: List[int]) -> str:
            if len(idxs) <= 1:
                return seqs[idxs[0]] if idxs else ""
            ids = idxs if len(idxs) <= REP_CAP else rng.choice(idxs, size=REP_CAP, replace=False).tolist()
            S = [seqs[i] for i in ids]
            m = len(S)
            if m == 1: return S[0]
            Msum = np.zeros(m, dtype=float)
            for a in range(m - 1):
                sa = S[a]
                for b in range(a + 1, m):
                    d = gapaware_hamming(sa, S[b])
                    Msum[a] += d; Msum[b] += d
            return S[int(np.argmin(Msum))]

        reps = [rep_of(g) for g in accepted]

        # Use tree distances between cluster representatives to block merges across long separations.
        # We compare the path length between representative leaves to a threshold derived from split_min_ib.
        # (Rule of thumb: if path >= 2 * split_min_ib, the MRCA separation implies long internal edges → block.)
        def allow_merge(i, j):
            if split_min_ib <= 0.0:
                return True
            # pick any mapped leaf name for each cluster representative
            rep_i = rep_name_of[i]  # we’ll set rep_name_of below
            rep_j = rep_name_of[j]
            try:
                d_tree = float(tree.distance(rep_i, rep_j))
            except Exception:
                d_tree = 0.0
            return not (d_tree >= 2.0 * split_min_ib)

        while len(accepted) > max_clusters:
            sizes = [len(g) for g in accepted]
            smin = int(np.argmin(sizes))
            # find nearest target by rep PID
            best, bestd = None, 1e9
            for k in range(len(accepted)):
                if k == smin: continue
                d = gapaware_hamming(reps[smin], reps[k])
                if d < bestd and allow_merge(smin, k):
                    best, bestd = k, d
            if best is None:
                break
            accepted[best].extend(accepted[smin])
            if verbose:
                print(f"[TREE] accept clade: n={n}, total_branch={total_len:.2f}, longest_ib={longest_ib:.2f}")

            # refresh rep of the target
            reps[best] = rep_of(accepted[best])
            # remove source
            accepted.pop(smin); reps.pop(smin)

    print(f"[counts] after merging clusters: {len(accepted)}  total_assigned={sum(len(g) for g in accepted)}")


    # ---- Assign sequences missing from the tree to nearest cluster representative ----
    if int(getattr(args, "tree_assign_missing", 1)) == 1:
        present = set()
        for g in accepted:
            present.update(g)
        missing = [i for i in range(len(seqs)) if i not in present]
        if missing:
            rng = np.random.default_rng(int(getattr(args, "sample_seed", 0)))
            REP_CAP = int(getattr(args, "tree_rep_cap", 200))

            # representatives (approx medoids) for each accepted cluster (reuse reps list if already built)
            if 'reps' not in locals():
                def rep_of(idxs: List[int]) -> str:
                    if len(idxs) <= 1:
                        return seqs[idxs[0]] if idxs else ""
                    ids = idxs if len(idxs) <= REP_CAP else rng.choice(idxs, size=REP_CAP, replace=False).tolist()
                    S = [seqs[i] for i in ids]
                    m = len(S)
                    if m == 1: return S[0]
                    Msum = np.zeros(m, dtype=float)
                    for a in range(m - 1):
                        sa = S[a]
                        for b in range(a + 1, m):
                            d = gapaware_hamming(sa, S[b])
                            Msum[a] += d;
                            Msum[b] += d
                    return S[int(np.argmin(Msum))]

                reps = [rep_of(g) for g in accepted]

            # Build assigns dict and add missing to nearest rep
            assigns_tmp: Dict[int, List[int]] = {ci: list(set(g)) for ci, g in enumerate(accepted)}
            for i in missing:
                s = seqs[i]
                # nearest representative by gap-aware Hamming
                dmin, cmin = 1e9, 0
                for ci, rseq in enumerate(reps):
                    d = gapaware_hamming(s, rseq)
                    if d < dmin:
                        dmin, cmin = d, ci
                assigns_tmp[cmin].append(i)
            accepted = [assigns_tmp[k] for k in sorted(assigns_tmp.keys())]

    assigns: Dict[int, List[int]] = {i: sorted(set(g)) for i, g in enumerate(accepted)}

    if verbose:
        print(f"[TREE] accepted {len(assigns)} clusters (pre post-process).")
    return assigns


# ------------- common post-processing and writing ------------------------
def postprocess_and_write(assigns: Dict[int, List[int]],
                          headers: List[str], seqs: List[str], args) -> pd.DataFrame:
    """
    Enforce min_output_size, cap at max_clusters, compute Neff and drop low-Neff,
    then write ShallowMsa_###.a3m. Return a per-cluster dataframe.
    """
    verbose = bool(int(args.verbose))
    st = StageTimer("Postprocess & write clusters", verbose)

    if not assigns:
        ensure_dir(args.outdir)
        meta = pd.DataFrame(columns=["cluster", "n", "neff", "path", "kept"], dtype=object)
        meta.to_csv(Path(args.outdir, f"{args.keyword}_clusters.csv"), index=False)
        _log("[WARN] No clusters produced; wrote empty metadata.", verbose)
        st.done()
        return meta

    # drop tiny clusters (< min_output_size)
    # --- Greedy tiny-merge: repeatedly merge the smallest cluster (< min_size) into its nearest cluster ---
    min_size = int(args.min_output_size)

    st = StageTimer("Greedy tiny-merge to reach min_output_size", verbose)

    rng = np.random.default_rng(int(getattr(args, "sample_seed", 0)))
    REP_CAP = 200  # compute cluster representative on at most this many sequences

    # Work on a local dict we can mutate
    clusters: Dict[int, List[int]] = {int(c): list(ids) for c, ids in assigns.items()}

    # Cache a representative (approx medoid) sequence per cluster, recompute lazily when cluster changes
    rep_cache: Dict[int, str] = {}

    def _rep_seq(cluster_id: int) -> str:
        """Return a representative sequence (approx medoid) for a cluster using up to REP_CAP members."""
        if cluster_id in rep_cache:
            return rep_cache[cluster_id]
        ids = clusters[cluster_id]
        if not ids:
            # shouldn't happen, but guard
            rep_cache[cluster_id] = ""
            return rep_cache[cluster_id]
        # sample up to REP_CAP members to make medoid affordable
        if len(ids) > REP_CAP:
            sample_ids = rng.choice(ids, size=REP_CAP, replace=False).tolist()
        else:
            sample_ids = ids
        S = [seqs[i] for i in sample_ids]
        m = len(S)
        if m == 1:
            rep_cache[cluster_id] = S[0]
            return rep_cache[cluster_id]
        # compute medoid on the sample
        Msum = np.zeros(m, dtype=float)
        for a in range(m - 1):
            sa = S[a]
            for b in range(a + 1, m):
                d = gapaware_hamming(sa, S[b])
                Msum[a] += d
                Msum[b] += d
        rep_cache[cluster_id] = S[int(np.argmin(Msum))]
        return rep_cache[cluster_id]

    # replace the whole _closest_cluster with this version
    def _closest_cluster(c_from: int) -> int | None:
        """Return id of nearest *other* cluster to c_from; None if no candidate."""
        # work on the local 'clusters' dict we mutate
        cand = [c for c in clusters.keys() if c != c_from and len(clusters[c]) > 0]
        if not cand:
            return None
        s_from = _rep_seq(c_from)  # <-- reuse the cached representative
        best, best_d = None, float("inf")
        for cj in cand:
            d = gapaware_hamming(s_from, _rep_seq(cj))
            if d < best_d:
                best_d = d
                best = cj
        return best

        # distance by tree representatives if available, else by gap-aware Hamming of medoids
        def _rep_seq(ci):
            # reuse any representative you cache; else pick a quick medoid
            ids = assigns[ci]
            return seqs[ids[0]] if len(ids) == 1 else seqs[medoid_of_rows(np.array(ids, dtype=int))]

        s_from = _rep_seq(c_from)
        for cj in cand:
            d = gapaware_hamming(s_from, _rep_seq(cj))
            if d < best_d:
                best_d = d
                best = cj
        return best

    print("[2025-.. ..:..:..] >>> Greedy tiny-merge to reach min_output_size ...")
    eff_min = int(getattr(args, "min_output_size", 0))
    merges = 0

    # No merging if disabled or nothing to merge
    if eff_min <= 0:
        print(f"[tiny-merge] disabled (min_output_size={eff_min})")
    elif len(assigns) <= 1:
        # Nothing to merge into. Accept the single cluster as-is.
        print(f"[tiny-merge] skipped: only {len(assigns)} cluster(s); "
              f"min_output_size={eff_min} > total_n={sum(len(v) for v in assigns.values())}")
    else:
        # Iteratively merge clusters smaller than eff_min into their nearest neighbor
        # Stop when all clusters are >= eff_min or when no valid target exists.
        while True:
            small = [c for c in sorted(clusters.keys()) if len(clusters[c]) < eff_min]

            if not small:
                break
            progressed = False
            for cmin in small:
                target = _closest_cluster(cmin)
                if target is None or target == cmin:
                    # No valid merge target; skip this one
                    continue
                clusters[target].extend(clusters[cmin])
                clusters[cmin] = []
                merges += 1
                progressed = True
            if not progressed:
                # no more merges possible
                break


    st.done(f"Greedy tiny-merge: merges={merges}; clusters_now={len(clusters)}")

    # If everything was tiny and merged to nothing (shouldn't happen), guard:
    if not clusters:
        print("[WARN] All clusters vanished during tiny-merge; nothing to write.")
        return pd.DataFrame(columns=["cluster","n","neff","path","kept"])

    # Reindex cluster ids to 0..C-1 for downstream
    assigns = {i_new: clusters[k] for i_new, k in enumerate(sorted(clusters.keys()))}


    # cap at max_clusters (largest first)
    maxC = int(args.max_clusters)
    order = sorted(assigns.keys(), key=lambda c: len(assigns[c]), reverse=True)
    keep_keys = order[:maxC]
    assigns = {i_new: assigns[k] for i_new, k in enumerate(keep_keys)}


    st = StageTimer("Capped max clusters", verbose)
    st.done(f"kept_centers={len(assigns)}")

    # Announce Neff phase
    _log(f"Finished clustering (clusters={len(assigns)}). "
         f"Computing Neff for clusters (mode={args.neff_mode}) ...", verbose)

    # Write clusters & compute Neff (A3M → msa_utils.compute_neff_from_a3m)
    ensure_dir(args.outdir)
    rows = []
    total = len(assigns)

    for k, (c, idxs) in enumerate(assigns.items(), start=1):
        ids = [headers[i] for i in idxs]
        ss = [seqs[i] for i in idxs]
        outp = Path(args.outdir, f"{args.keyword}_{c:03d}.a3m")
        write_fasta(ids, ss, str(outp))

        _log(f"[neff] cluster {k}/{total} (size={len(idxs)}) using {args.neff_mode}", verbose)

        if args.neff_mode == "approx":
            neff = compute_neff_from_a3m(str(outp),
                        id_thresh = float(args.neff_id_thresh),
                        use_query_columns = bool(int(args.neff_use_query_columns)),
                        strip_inserts = bool(int(args.neff_strip_inserts)),
                        mode = "approx",
                        approx_n_hashes = int(args.neff_approx_n_hashes),
                        approx_sig_len = int(args.neff_approx_sig_len),
                        approx_bucket_cap = int(args.neff_approx_bucket_cap),
                        approx_seed = int(args.sample_seed),
                        approx_subsample_cap = int(args.neff_approx_subsample_cap))
        else:
            neff = compute_neff_from_a3m(str(outp),
                        id_thresh = float(args.neff_id_thresh),
                        use_query_columns = bool(int(args.neff_use_query_columns)),
                        strip_inserts = bool(int(args.neff_strip_inserts)),
                        mode = "exact")

        rows.append(dict(cluster=c, n=len(idxs), neff=neff, path=str(outp), kept=True))

    st = StageTimer("Computed Neff for clusters", verbose)

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

    st = StageTimer("Dropped low-Neff", verbose)

    final_n = sum(len(idxs) for idxs in assigns.values())
    print(f"[counts] final clusters: {len(assigns)}  total_assigned={final_n}")
#    print(f"[counts] total_unassigned={n_after_gaps - final_n}")

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

    p.add_argument("--metric", choices=["hamming"], default="hamming", help="Distance metric (currently: hamming).")

    # Tree clustering knobs
    p.add_argument("--tree_path", default=None, help="Newick file for --cluster_alg tree.")
    p.add_argument("--tree_min_total_branch", type=float, default=0.0,
                               help = "Require subtree total branch length >= this to accept a clade (0=ignore).")
    p.add_argument("--tree_max_size", type=int, default=0,
                               help = "If >0, clades with size > this are considered too large and will be split (0=ignore).")
    p.add_argument("--tree_max_total_branch", type=float, default=0.0,
                               help = "If >0, clades with total branch length > this are considered too large and will be split (0=ignore).")
    p.add_argument("--tree_split_min_internal_branch", type=float, default=0.0,
                               help = "Prefer splitting a clade if the longest internal branch at this node is >= threshold (0=ignore).")
    p.add_argument("--tree_min_num_clusters", type=int, default=5,
                   help="Target minimum number of clusters; sets a relative max cluster size total_n / K. "
                        "If that would give clusters smaller than min_output_size, size cap is disabled.")
    p.add_argument("--tree_assign_missing", type=int, default=1,
                   help="If 1, assign sequences not present in the tree to nearest cluster representative (default: 1).")
    p.add_argument("--tree_rep_cap", type=int, default=200,
                   help="Max sequences to use when computing a cluster representative/medoid (default: 200).")

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

      # Neff mode & approximation knobs
    p.add_argument("--neff_mode", choices=["exact", "approx"], default="approx",
                                       help = "How to compute Neff per cluster (default: approx).")
    p.add_argument("--neff_use_query_columns", type=int, default=1,
                                       help = "1=restrict Neff PID to columns where the first row has residues (default: 1).")
    p.add_argument("--neff_strip_inserts", type=int, default=1,
                                       help = "1=strip lowercase inserts from A3M before Neff (default: 1).")
    p.add_argument("--neff_approx_n_hashes", type=int, default=3,
                                       help = "Approx mode: number of random signatures for blocking (default: 3).")
    p.add_argument("--neff_approx_sig_len", type=int, default=32,
                                       help = "Approx mode: signature length in columns (default: 32).")
    p.add_argument("--neff_approx_bucket_cap", type=int, default=5000,
                                       help = "Approx mode: if a bucket would induce >cap pairs, sample inside (default: 5000).")
    p.add_argument("--neff_approx_subsample_cap", type=int, default=0,
                                       help = "Approx mode: if n>cap, subsample sequences before blocking (0=disabled).")

# AHC knobs
    p.add_argument("--ahc_linkage", choices=["average", "complete", "single"], default="average",
                   help="Linkage for AHC (default: average).")
    p.add_argument("--ahc_cut_mode", choices=["maxclust", "distance"], default="maxclust",
                   help="How to cut the AHC tree: maxclust or distance.")
    p.add_argument("--ahc_distance_threshold", type=float, default=0.18,
                   help="If --ahc_cut_mode=distance, cut at this distance threshold (0.0 means auto-guess).")
    p.add_argument("--ahc_merge_on_sample", type=int, default=0,
                   help="If 1, merge sample clusters smaller than min_output_size before assigning all sequences. "
                        "Default 0 (skip merge on sample; enforce size after assignment).")
    p.add_argument("--ahc_min_center_size", type=int, default=10,
                   help="On the sample, clusters smaller than this are merged into nearest before assignment (0=disable).")
    p.add_argument("--ahc_center_cap", type=int, default=500,
                   help="Cap the number of centers used for assignment by merging the smallest clusters until this count (0=disable).")

    # verbosity printing
    p.add_argument("--verbose", type=int, default=1, help="1=log stages (default), 0=quiet")
    p.add_argument("--progress_every", type=int, default=500,
                   help="Print progress every N rows while building distances.")

    return p

def main():
    args = build_argparser().parse_args()

    headers, seqs = load_a3m(args.a3m)

    raw_n = len(headers)  # or however you store MSA lines
    L = len(seqs[0]) if seqs else 0
    print(f"[counts] MSA total={raw_n}  L={L}")

    # Build gap filter (works for all algorithms)
    frac = np.array([s.count('-') / len(s) for s in seqs], dtype=float)
    keep_mask = (frac <= float(args.frac_gaps_cutoff))
    if keep_mask.mean() < 1.0:
        print(f"[filter] dropping {int((~keep_mask).sum())} sequences by gap fraction > {args.frac_gaps_cutoff}")
    # Apply mask
    headers = [h for h, k in zip(headers, keep_mask) if k]
    seqs = [s for s, k in zip(seqs, keep_mask) if k]

    # Suppose you build keep_mask by frac_gaps_cutoff (see next patch)
    n_after_gaps = int(keep_mask.sum())
    print(f"[counts] after gap-filter (frac_gaps <= {args.frac_gaps_cutoff}): {n_after_gaps} kept, {raw_n - n_after_gaps} dropped")

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
