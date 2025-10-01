#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Cluster MSA (A3M) with HDBSCAN, enforce limits on number/size of clusters,
and filter by N_effective (Meff). Writes A3M clusters (query on top) suitable
for AlphaFold / MSA-Transformer.

Design:
- Prefer coarse clusters: HDBSCAN(cluster_selection_method="eom")
- Avoid micro-clusters via larger min_cluster_size + post-filters
- Cap total clusters with --max_clusters (keep largest K by size)
- Require both raw size (--min_output_size) and information content (--min_neff)

This file is self-contained, but if the project provides msa_utils.parse_a3m_records
and msa_utils.write_fasta, they will be used automatically.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from glob import glob
from typing import List, Tuple, Optional
from hdbscan import HDBSCAN

# Optional: faster edit distances if present
try:
    from polyleven import levenshtein as _levenshtein
    HAS_POLYLEVEN = True
except Exception:
    HAS_POLYLEVEN = False

# ---------- Try to use project utilities if available ----------
# We try to import parse_a3m_records and write_fasta; if not present,
# fallbacks defined below will be used.
_parse_a3m_records = None
_write_fasta = None
try:
    from utils.msa_utils import parse_a3m_records as _parse_a3m_records  # type: ignore
except Exception:
    pass
try:
    from utils.msa_utils import write_fasta as _write_fasta  # type: ignore
except Exception:
    pass


# ---------- Fallback utilities (used only if project utils are missing) ----------
def _fallback_parse_a3m_records(lines: List[str],
                                strip_inserts: bool = True) -> Tuple[List[str], List[str]]:
    """
    Minimal A3M reader.
    Returns (headers, sequences) with lowercase insertions stripped if requested.
    """
    headers, seqs = [], []
    cur_h = None
    cur_s = []
    for ln in lines:
        ln = ln.rstrip("\n")
        if not ln:
            continue
        if ln.startswith(">"):
            if cur_h is not None:
                s = "".join(cur_s)
                if strip_inserts:
                    # keep uppercase letters and '-' gaps only
                    s = "".join(ch for ch in s if (ch == '-' or ch.isupper()))
                headers.append(cur_h)
                seqs.append(s)
            cur_h = ln
            cur_s = []
        else:
            cur_s.append(ln.strip())
    if cur_h is not None:
        s = "".join(cur_s)
        if strip_inserts:
            s = "".join(ch for ch in s if (ch == '-' or ch.isupper()))
        headers.append(cur_h)
        seqs.append(s)
    return headers, seqs


def _fallback_write_fasta(headers: List[str], seqs: List[str], outfile: str) -> None:
    with open(outfile, "w") as f:
        for h, s in zip(headers, seqs):
            if not h.startswith(">"):
                h = ">" + h
            f.write(h + "\n")
            # wrap to 80 chars
            for i in range(0, len(s), 80):
                f.write(s[i:i+80] + "\n")


# bind utils (prefer project versions if present)
parse_a3m_records = _parse_a3m_records or _fallback_parse_a3m_records
write_fasta = _write_fasta or _fallback_write_fasta


# ---------- Core: Meff ----------
def compute_neff_from_seqs(seqs: List[str],
                           id_thresh: float = 0.8,
                           use_query_columns: bool = True) -> float:
    """
    Compute N_effective (Meff) from in-memory aligned sequences.
    - Identity computed over positions with residues in BOTH sequences.
    - If use_query_columns=True, restrict to columns where the QUERY (seqs[0]) has residues.
    Return a float Meff (sum of 1 / neighborhood counts).
    """
    n = len(seqs)
    if n == 0:
        return 0.0
    L = len(seqs[0])
    if any(len(s) != L for s in seqs):
        raise ValueError("All sequences must be the same aligned length.")

    # build array (n, L), characters
    arr = np.array([list(s) for s in seqs], dtype='<U1')

    # restrict to query residue columns if requested
    if use_query_columns:
        keep = (arr[0] != '-')
        if keep.sum() == 0:
            return 0.0
        arr = arr[:, keep]
        L = arr.shape[1]

    # identity(i,j) = matches / non_gap_positions, with both residues present
    def pid_row(i, sl):
        A = arr[i]                 # (L,)
        B = arr[sl]                # (B, L)
        both = (A != '-') & (B != '-')
        denom = both.sum(axis=1)   # (B,)
        matches = (A == B) & both
        num = matches.sum(axis=1)
        with np.errstate(divide='ignore', invalid='ignore'):
            pid = np.where(denom > 0, num / denom, 0.0)
        return pid

    counts = np.zeros(n, dtype=np.int32)
    block = 1024
    for i in range(n):
        c = 0
        for start in range(0, n, block):
            end = min(n, start + block)
            pid = pid_row(i, slice(start, end))
            c += int((pid >= id_thresh).sum())
        counts[i] = max(c, 1)
    weights = 1.0 / counts.astype(np.float64)
    return float(weights.sum())


def compute_neff_from_a3m(a3m_path: str,
                          id_thresh: float = 0.8,
                          use_query_columns: bool = True,
                          strip_inserts: bool = True) -> int:
    with open(a3m_path, "r") as f:
        lines = f.readlines()
    headers, seqs = parse_a3m_records(lines, strip_inserts=strip_inserts)
    meff = compute_neff_from_seqs(seqs, id_thresh=id_thresh, use_query_columns=use_query_columns)
    return int(round(meff))


# ---------- Distance matrix (precomputed) ----------
def normalized_hamming_distance_matrix(seqs: List[str],
                                       use_query_columns: bool = True) -> np.ndarray:
    """
    1 - pid over positions where both have residues (like AF conventions).
    If use_query_columns=True: only columns where query (seqs[0]) has residues.
    Returns (n,n) symmetric matrix with zeros on diag.
    """
    n = len(seqs)
    L = len(seqs[0])
    arr = np.array([list(s) for s in seqs], dtype='<U1')
    if use_query_columns:
        keep = (arr[0] != '-')
        arr = arr[:, keep]
    # compute pairwise distances in blocks
    D = np.zeros((n, n), dtype=np.float64)
    block = 256
    for i in range(0, n, block):
        A = arr[i:i+block]
        for j in range(i, n, block):
            B = arr[j:j+block]
            # both residues mask: (A_b,1,L) & (1,B_b,L)
            both = (A[:, None, :] != '-') & (B[None, :, :] != '-')
            denom = both.sum(axis=2).astype(np.float64) + 1e-9
            matches = (A[:, None, :] == B[None, :, :]) & both
            pid = matches.sum(axis=2).astype(np.float64) / denom
            dist = 1.0 - pid
            D[i:i+A.shape[0], j:j+B.shape[0]] = dist
            if j != i:
                D[j:j+B.shape[0], i:i+A.shape[0]] = dist.T
    np.fill_diagonal(D, 0.0)
    return D


def normalized_levenshtein_distance_matrix(seqs: List[str]) -> np.ndarray:
    """
    Normalized Levenshtein distance (edit distance / length).
    Requires polyleven; falls back to simple Python if unavailable (slower).
    """
    n = len(seqs)
    L = len(seqs[0])
    D = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        si = seqs[i]
        for j in range(i+1, n):
            sj = seqs[j]
            if HAS_POLYLEVEN:
                d = _levenshtein(si, sj)
            else:
                # slower fallback
                d = _plain_lev(si, sj)
            dist = d / float(L)
            D[i, j] = D[j, i] = dist
    return D


def _plain_lev(a: str, b: str) -> int:
    # very simple Levenshtein for fallback (only used if polyleven missing)
    la, lb = len(a), len(b)
    dp = list(range(lb + 1))
    for i in range(1, la + 1):
        prev, dp[0] = dp[0], i
        for j in range(1, lb + 1):
            cur = dp[j]
            cost = 0 if a[i-1] == b[j-1] else 1
            dp[j] = min(dp[j] + 1, dp[j-1] + 1, prev + cost)
            prev = cur
    return dp[lb]


# ---------- IO helpers ----------
def read_a3m(a3m_path: str,
             strip_inserts: bool = True) -> Tuple[List[str], List[str]]:
    with open(a3m_path, "r") as f:
        lines = f.readlines()
    headers, seqs = parse_a3m_records(lines, strip_inserts=strip_inserts)
    # sanity
    if len(seqs) == 0:
        raise ValueError("Empty A3M.")
    L = len(seqs[0])
    if any(len(s) != L for s in seqs):
        raise ValueError("Sequences are not aligned to the same length.")
    return headers, seqs


def ensure_dir(p: str):
    os.makedirs(p, exist_ok=True)


# ---------- Main clustering ----------
def run_clustering(a3m_path: str,
                   outdir: str,
                   keyword: str,
                   metric: str = "hamming",
                   min_cluster_size: int = 200,
                   min_samples: Optional[int] = None,
                   cluster_selection: str = "eom",
                   max_clusters: int = 100,
                   min_output_size: int = 200,
                   min_neff: int = 100,
                   neff_id_thresh: float = 0.8,
                   frac_gaps_cutoff: float = 0.5,
                   use_query_columns: bool = True,
                   max_n_for_precomputed: int = 6000) -> pd.DataFrame:
    """
    Returns a DataFrame with cluster metadata. Writes cluster A3Ms to outdir.
    """
    ensure_dir(outdir)
    headers, seqs = read_a3m(a3m_path, strip_inserts=True)
    query_h = headers[0]
    query_s = seqs[0]
    L = len(query_s)
    n_all = len(seqs)

    # filter by fraction of gaps
    def frac_gaps(s): return s.count('-') / L
    keep_mask = np.array([frac_gaps(s) < frac_gaps_cutoff for s in seqs], dtype=bool)
    if keep_mask.sum() < 2:
        raise ValueError("Too few sequences after gap filtering.")
    headers = [h for h, k in zip(headers, keep_mask) if k]
    seqs = [s for s, k in zip(seqs, keep_mask) if k]
    # ensure query is first
    if headers[0] != query_h:
        # find old query and move it to front
        qi = headers.index(query_h) if query_h in headers else 0
        if qi != 0:
            headers = [headers[qi]] + headers[:qi] + headers[qi+1:]
            seqs = [seqs[qi]] + seqs[:qi] + seqs[qi+1:]

    n = len(seqs)
    print(f"[INFO] A3M: {n_all} -> {n} after gap filter (cutoff={frac_gaps_cutoff}). Aln len={L}")

    # ---------- choose clustering path ----------
    # Path A: small/medium n -> full precomputed, gap-aware PID (your current precise logic)
    # Path B: large n        -> feature-mode with metric='hamming' (no n×n matrix)
    USE_PRECOMPUTED = (metric == "hamming" and n <= max_n_for_precomputed) or (metric == "lev" and n <= 4000)
#    USE_PRECOMPUTED = (not args.force_feature_mode) and \
#                      ((metric == "hamming" and n <= max_n_for_precomputed) or (metric == "lev" and n <= 4000))

    if metric == "hamming" and USE_PRECOMPUTED:
        # ===== Path A: precomputed, gap-aware PID =====
        D = normalized_hamming_distance_matrix(seqs, use_query_columns=use_query_columns)
        hdb = HDBSCAN(min_cluster_size=min_cluster_size,
                      min_samples=min_samples,
                      cluster_selection_method=cluster_selection,
                      metric='precomputed',
                      core_dist_n_jobs=0)
        D = np.asarray(D, dtype=np.float64, order="C")
        labels = hdb.fit_predict(D)

    elif metric == "hamming" and not USE_PRECOMPUTED:
        # ===== Path B: large-n fallback: feature-mode Hamming =====
        # Build an int8 matrix of shape (n, Lq), where Lq are query residue columns only
        # Map residues to small integers; treat gaps as 255 so gap vs residue != match.
        arr = np.array([list(s) for s in seqs], dtype='<U1')
        if use_query_columns:
            keep = (arr[0] != '-')
            arr = arr[:, keep]
        # compact code per character (gap -> 255, residues -> 0..24)
        alphabet = np.array(list("ACDEFGHIKLMNPQRSTVWYBJZX-"), dtype='<U1')  # include common ambiguous + gap '-'
        codebook = {ch: i for i, ch in enumerate(alphabet)}
        GAP_CODE = 255
        X = np.empty(arr.shape, dtype=np.uint8)
        # vectorized map
        flat = arr.ravel()
        mapped = np.fromiter((codebook.get(ch, GAP_CODE) for ch in flat), count=flat.size, dtype=np.uint16)
        X[:, :] = mapped.reshape(arr.shape).astype(np.uint8)
        # Fit HDBSCAN directly on features with Hamming; no precomputed matrix is built.
        hdb = HDBSCAN(min_cluster_size=min_cluster_size,
                      min_samples=min_samples,
                      cluster_selection_method=cluster_selection,
                      metric='hamming',
                      algorithm='best',
                      approx_min_span_tree=True,
                      core_dist_n_jobs=0)
        labels = hdb.fit_predict(X)

    elif metric == "lev":
        if n > 4000:
            raise RuntimeError("Levenshtein with precomputed distances is O(n^2); "
                               "for n>4000 switch to --metric hamming (feature-mode).")
        D = normalized_levenshtein_distance_matrix(seqs)
        hdb = HDBSCAN(min_cluster_size=min_cluster_size,
                      min_samples=min_samples,
                      cluster_selection_method=cluster_selection,
                      metric='precomputed',
                      core_dist_n_jobs=0)
        D = np.asarray(D, dtype=np.float64, order="C")
        labels = hdb.fit_predict(D)
    else:
        raise ValueError("metric must be 'hamming' or 'lev'.")

    labels = np.asarray(labels, dtype=int)
    # post-filter tiny clusters
    valid = labels >= 0
    cluster_ids, counts = np.unique(labels[valid], return_counts=True)
    small = set(cluster_ids[counts < min_output_size])
    labels = np.array([-1 if (lab in small) else lab for lab in labels], dtype=int)

    # cap number of clusters by size
    cluster_ids, counts = np.unique(labels[labels >= 0], return_counts=True)
    kept_clusters = set()
    if len(cluster_ids) > max_clusters:
        order = np.argsort(counts)[::-1]
        kept_clusters = set(cluster_ids[order[:max_clusters]])
        labels = np.array([lab if lab in kept_clusters else -1 for lab in labels], dtype=int)

    # reindex cluster labels to 0..C-1
    uniq = sorted(list(set(labels[labels >= 0])))
    remap = {old: i for i, old in enumerate(uniq)}
    labels = np.array([remap[lab] if lab in remap else -1 for lab in labels], dtype=int)

    # reindex cluster labels to 0..C-1
    uniq = sorted(list(set(labels[labels >= 0])))
    remap = {old: i for i, old in enumerate(uniq)}
    labels = np.array([remap[lab] if lab in remap else -1 for lab in labels], dtype=int)

    # If no clusters at all, write an empty metadata CSV and return cleanly
    if len(uniq) == 0:
        meta_csv = os.path.join(outdir, f"{keyword}_clusters.csv")
        pd.DataFrame(columns=["cluster", "kept", "n", "neff", "path"]).to_csv(meta_csv, index=False)
        print(f"[INFO] No clusters found; all sequences labeled noise. Wrote empty {meta_csv}. "
              f"Hint: lower --min_cluster_size/--min_output_size, set a smaller --min_samples, "
              f"or try --cluster_selection leaf.")
        return pd.DataFrame(columns=["cluster", "kept", "n", "neff", "path"])

    # clean previous outputs
    for old in glob(os.path.join(outdir, f"{keyword}_*.a3m")):
        try:
            os.remove(old)
        except Exception:
            pass

    # write clusters (query on top), and compute Meff
    rows = []
    n_noise = int((labels < 0).sum())
    print(f"[INFO] clusters found (before Meff): {len(uniq)} ; noise={n_noise}")

    for c in sorted(set(labels) - {-1}):
        idx = np.where(labels == c)[0]
        # ensure query is included in each cluster file by prepending it
        # (common practice for AF/MSA-transformer)
        cl_headers = [headers[0]] + [headers[i] for i in idx if i != 0]
        cl_seqs = [seqs[0]] + [seqs[i] for i in idx if i != 0]

        out_path = os.path.join(outdir, f"{keyword}_{c:03d}.a3m")
        write_fasta(cl_headers, cl_seqs, outfile=out_path)

        meff = compute_neff_from_seqs(cl_seqs, id_thresh=neff_id_thresh,
                                      use_query_columns=True)
        n_raw = len(cl_seqs)

        if (min_neff is not None and min_neff > 0) and (meff < min_neff):
            # drop if not enough information content
            try:
                os.remove(out_path)
            except Exception:
                pass
            rows.append(dict(cluster=c, kept=False, n=n_raw, neff=meff, path=out_path))
            continue

        rows.append(dict(cluster=c, kept=True, n=n_raw, neff=meff, path=out_path))

    df_meta = pd.DataFrame(rows)
    if not df_meta.empty:
        df_meta = df_meta.sort_values(["kept", "n"], ascending=[False, False])
    meta_csv = os.path.join(outdir, f"{keyword}_clusters.csv")
    df_meta.to_csv(meta_csv, index=False)
    print(f"[INFO] Wrote cluster metadata: {meta_csv}")
    return df_meta


# ---------- CLI ----------
def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Cluster an A3M with HDBSCAN and filter by size/Neff.")
    p.add_argument("--a3m", required=True, help="Input A3M path.")
    p.add_argument("-o", "--outdir", required=True, help="Output directory.")
    p.add_argument("--keyword", default="ShallowMsa", help="Prefix for cluster files.")
    p.add_argument("--metric", choices=["hamming", "lev"], default="hamming",
                   help="Distance metric (precomputed). 'hamming' respects aligned columns and gaps.")
    p.add_argument("--min_cluster_size", type=int, default=200,
                   help="HDBSCAN min_cluster_size (coarser -> fewer clusters).")
    p.add_argument("--min_samples", type=str, default="auto",
                   help="'auto' or integer. Larger -> fewer clusters.")
    p.add_argument("--cluster_selection", choices=["eom", "leaf"], default="eom",
                   help="Use 'eom' for fewer, larger clusters.")
    p.add_argument("--max_clusters", type=int, default=100,
                   help="Keep only largest K clusters; rest -> noise.")
    p.add_argument("--min_output_size", type=int, default=200,
                   help="Drop clusters with raw size < this after HDBSCAN.")
    p.add_argument("--min_neff", type=int, default=100,
                   help="Drop clusters whose Meff < this threshold.")
    p.add_argument("--neff_id_thresh", type=float, default=0.8,
                   help="Identity threshold for Meff (e.g., 0.8).")
    p.add_argument("--frac_gaps_cutoff", type=float, default=0.5,
                   help="Filter out sequences with fraction of gaps >= cutoff.")
    p.add_argument("--use_query_columns", type=int, default=1,
                   help="1: compute distances/Meff on columns where query has residues; 0: all columns.")
    p.add_argument("--max_n_for_precomputed", type=int, default=6000,
                   help="Safety limit for O(n^2) distance computations (hamming).")
    p.add_argument("--force_feature_mode", type=int, default=0,
                   help="Set to 1 to always use feature-mode Hamming (no precomputed matrix).")

    return p


def main():
    args = build_argparser().parse_args()
    min_samples = None if args.min_samples == "auto" else int(args.min_samples)
    df = run_clustering(
        a3m_path=args.a3m,
        outdir=args.outdir,
        keyword=args.keyword,
        metric=args.metric,
        min_cluster_size=int(args.min_cluster_size),
        min_samples=min_samples,
        cluster_selection=args.cluster_selection,
        max_clusters=int(args.max_clusters),
        min_output_size=int(args.min_output_size),
        min_neff=int(args.min_neff),
        neff_id_thresh=float(args.neff_id_thresh),
        frac_gaps_cutoff=float(args.frac_gaps_cutoff),
        use_query_columns=bool(int(args.use_query_columns)),
        max_n_for_precomputed=int(args.max_n_for_precomputed),
    )
    # brief console summary
    kept = df[df.kept]
    dropped = df[~df.kept]
    print(f"[SUMMARY] Kept {len(kept)} clusters, Dropped {len(dropped)} (by Neff).")
    if not kept.empty:
        print(kept[["cluster", "n", "neff", "path"]].to_string(index=False))


if __name__ == "__main__":
    main()
