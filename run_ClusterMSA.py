#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Cluster MSA (A3M) with a scalable sample→cluster→assign heuristic:

1) Randomly sample up to --sample_cap sequences.
2) Build full precomputed gap-aware distance on the sample and run HDBSCAN.
3) Compute cluster centers (medoids) on the sample clusters.
4) Assign all non-sampled sequences to the nearest center (optionally only if <= --max_center_dist).
5) Enforce limits (min size, max clusters) and Neff threshold, then write A3Ms.

Notes:
- Gap-aware Hamming: distance = 1 - PID over positions where BOTH sequences have residues.
- Restrict to query-residue columns if --use_query_columns=1 (default).
- Centers are MEDOIDS (actual sequences), so assignment uses the same distance definition.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from glob import glob
from typing import List, Tuple, Optional, Dict
from pathlib import Path
from hdbscan import HDBSCAN

# ---------- Optional fast edit distances ----------
try:
    from polyleven import levenshtein as _levenshtein
    HAS_POLYLEVEN = True
except Exception:
    HAS_POLYLEVEN = False

# ---------- Prefer project utilities if available ----------
_parse_a3m_records = None
_write_fasta = None
try:
    # Try typical project layout (adjusted to your repo earlier)
    from msa_utils import parse_a3m_records as _parse_a3m_records  # type: ignore
except Exception:
    try:
        from utils.msa_utils import parse_a3m_records as _parse_a3m_records  # type: ignore
    except Exception:
        pass

try:
    from msa_utils import write_fasta as _write_fasta  # type: ignore
except Exception:
    try:
        from utils.msa_utils import write_fasta as _write_fasta  # type: ignore
    except Exception:
        pass


# ---------- Fallback utilities (only if project utils are missing) ----------
def _fallback_parse_a3m_records(lines: List[str],
                                strip_inserts: bool = True) -> Tuple[List[str], List[str]]:
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
            for i in range(0, len(s), 80):
                f.write(s[i:i+80] + "\n")


# bind utils
parse_a3m_records = _parse_a3m_records or _fallback_parse_a3m_records
write_fasta = _write_fasta or _fallback_write_fasta


# ---------- IO ----------
def read_a3m(a3m_path: str,
             strip_inserts: bool = True) -> Tuple[List[str], List[str]]:
    lines = Path(a3m_path).read_text().splitlines(True)
    headers, seqs = parse_a3m_records(lines, strip_inserts=strip_inserts)
    if not seqs:
        raise ValueError("Empty A3M.")
    L = len(seqs[0])
    if any(len(s) != L for s in seqs):
        raise ValueError("Sequences are not aligned to the same length.")
    return headers, seqs


def ensure_dir(p: str) -> None:
    Path(p).mkdir(parents=True, exist_ok=True)


# ---------- Meff ----------
def compute_neff_from_seqs(seqs: List[str],
                           id_thresh: float = 0.8,
                           use_query_columns: bool = True) -> float:
    n = len(seqs)
    if n == 0:
        return 0.0
    arr = np.array([list(s) for s in seqs], dtype='<U1')
    if use_query_columns:
        keep = (arr[0] != '-')
        if keep.sum() == 0:
            return 0.0
        arr = arr[:, keep]

    def pid_row(i, sl):
        A = arr[i]
        B = arr[sl]
        both = (A != '-') & (B != '-')
        denom = both.sum(axis=1)
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
    return float((1.0 / counts.astype(np.float64)).sum())


# ---------- Distances (gap-aware, query-column option) ----------
def _project_cols(arr: np.ndarray, use_query_columns: bool) -> np.ndarray:
    if not use_query_columns:
        return arr
    keep = (arr[0] != '-')
    return arr[:, keep]


def dist_pid_gapaware_matrix(seqs: List[str],
                             use_query_columns: bool = True) -> np.ndarray:
    """Full (n,n) matrix of gap-aware distances (1-PID over both-residue positions)."""
    arr = np.array([list(s) for s in seqs], dtype='<U1')
    arr = _project_cols(arr, use_query_columns)
    n = arr.shape[0]
    D = np.zeros((n, n), dtype=np.float64)
    block = 256
    for i in range(0, n, block):
        A = arr[i:i+block]
        for j in range(i, n, block):
            B = arr[j:j+block]
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


def dist_one_vs_many_gapaware(one: str, many: List[str], use_query_columns: bool = True) -> np.ndarray:
    """Distances from ONE sequence to MANY (1-PID, gap-aware)."""
    arr_many = np.array([list(s) for s in many], dtype='<U1')
    arr_one = np.array([list(one)], dtype='<U1')
    if use_query_columns:
        keep = (arr_many[0] != '-')
        arr_many = arr_many[:, keep]
        arr_one = arr_one[:, keep]
    A = arr_one[0]         # (L,)
    B = arr_many           # (m,L)
    both = (A != '-') & (B != '-')
    denom = both.sum(axis=1).astype(np.float64) + 1e-9
    matches = (A == B) & both
    pid = matches.sum(axis=1).astype(np.float64) / denom
    return 1.0 - pid


# ---------- Core: sample → cluster → assign ----------
def _cluster_sample_with_hdbscan(sample_seqs: List[str],
                                 min_cluster_size: int,
                                 min_samples: Optional[int],
                                 cluster_selection: str,
                                 use_query_columns: bool) -> np.ndarray:
    """Run HDBSCAN on the SAMPLE set with precomputed gap-aware distances."""
    D = dist_pid_gapaware_matrix(sample_seqs, use_query_columns=use_query_columns)
    D = np.asarray(D, dtype=np.float64, order="C")
    hdb = HDBSCAN(min_cluster_size=min_cluster_size,
                  min_samples=min_samples,
                  cluster_selection_method=cluster_selection,
                  metric='precomputed',
                  core_dist_n_jobs=0)
    return hdb.fit_predict(D)


def _medoids_from_labels(D: np.ndarray, labels: np.ndarray) -> Dict[int, int]:
    """Return {cluster_id: medoid_index} for clusters in the SAMPLE; indices are within the SAMPLE."""
    medoids = {}
    for c in sorted(set(labels) - {-1}):
        idx = np.where(labels == c)[0]
        if idx.size == 0:
            continue
        sub = D[np.ix_(idx, idx)]
        # medoid = argmin row-sum of distances
        medoid_local = int(sub.sum(axis=1).argmin())
        medoids[c] = int(idx[medoid_local])
    return medoids


def _assign_rest_to_medoids(rest_seqs: List[str],
                            centers: List[str],
                            use_query_columns: bool,
                            max_center_dist: Optional[float]) -> np.ndarray:
    """
    For each sequence in 'rest', assign to nearest medoid index (0..C-1 in 'centers').
    If max_center_dist is not None, assign only if distance <= max_center_dist; else label -1 (noise).
    Returns an array of length len(rest_seqs) with center indices or -1.
    """
    if not rest_seqs:
        return np.zeros(0, dtype=int)
    C = len(centers)
    # Compute distances to all centers (C is small, <= max_clusters)
    # Do this center by center to keep memory small.
    n_rest = len(rest_seqs)
    best = np.full(n_rest, np.inf, dtype=np.float64)
    best_c = np.full(n_rest, -1, dtype=int)
    for c, center_seq in enumerate(centers):
        # vectorize: compute 1-vs-many (center vs rest)
        d = np.array([dist_one_vs_many_gapaware(s, [center_seq], use_query_columns=use_query_columns)[0]
                      for s in rest_seqs], dtype=np.float64)
        improve = d < best
        best[improve] = d[improve]
        best_c[improve] = c
    if max_center_dist is not None:
        best_c = np.where(best <= float(max_center_dist), best_c, -1)
    return best_c


# ---------- Main driver ----------
def run_clustering(a3m_path: str,
                   outdir: str,
                   keyword: str,
                   metric: str = "hamming",
                   min_cluster_size: int = 200,
                   min_samples: Optional[int] = None,
                   cluster_selection: str = "leaf",
                   max_clusters: int = 100,
                   min_output_size: int = 200,
                   min_neff: int = 100,
                   neff_id_thresh: float = 0.8,
                   frac_gaps_cutoff: float = 0.5,
                   use_query_columns: bool = True,
                   # sampling heuristic
                   sample_cap: int = 10000,
                   sample_seed: int = 12345,
                   max_center_dist: Optional[float] = None,
                   # legacy/alt:
                   max_n_for_precomputed: int = 6000,
                   force_feature_mode: int = 0) -> pd.DataFrame:
    """
    Returns a DataFrame with cluster metadata. Writes cluster A3Ms to outdir.

    Heuristic path (default):
      If n > sample_cap and metric='hamming', use sample→cluster→assign.
      Else, fall back to original behavior (precomputed for small n; feature-mode for huge n if forced).

    max_center_dist: if set (e.g. 0.6), a rest sequence is assigned to nearest center only if dist <= threshold;
                     otherwise left as noise (-1). Default None means "assign all".
    """
    ensure_dir(outdir)
    headers, seqs = read_a3m(a3m_path, strip_inserts=True)
    query_h, query_s = headers[0], seqs[0]
    L = len(query_s)
    n_all = len(seqs)

    # Filter by gaps
    def frac_gaps(s): return s.count('-') / L
    keep = np.array([frac_gaps(s) < frac_gaps_cutoff for s in seqs], dtype=bool)
    if keep.sum() < 2:
        raise ValueError("Too few sequences after gap filtering.")
    headers = [h for h, k in zip(headers, keep) if k]
    seqs = [s for s, k in zip(seqs, keep) if k]

    # Ensure query is first
    if headers[0] != query_h:
        qi = headers.index(query_h) if query_h in headers else 0
        if qi != 0:
            headers = [headers[qi]] + headers[:qi] + headers[qi+1:]
            seqs = [seqs[qi]] + seqs[:qi] + seqs[qi+1:]

    n = len(seqs)
    print(f"[INFO] A3M: {n_all} -> {n} after gap filter (cutoff={frac_gaps_cutoff}). Aln len={L}")

    # ------------------ Heuristic switch ------------------
    use_heuristic = (metric == "hamming") and (n > sample_cap) and (not force_feature_mode)
    labels = None

    if use_heuristic:
        # 1) sample
        rng = np.random.default_rng(sample_seed)
        sample_idx = np.arange(n)
        rng.shuffle(sample_idx)
        sample_idx = sample_idx[:sample_cap]
        sample_idx.sort()
        rest_idx = np.setdiff1d(np.arange(n), sample_idx, assume_unique=True)

        sample_headers = [headers[i] for i in sample_idx]
        sample_seqs = [seqs[i] for i in sample_idx]

        # 2) cluster sample (precomputed gap-aware)
        D_sample = dist_pid_gapaware_matrix(sample_seqs, use_query_columns=use_query_columns)
        sample_labels = _cluster_sample_with_hdbscan(
            sample_seqs=sample_seqs,
            min_cluster_size=min_cluster_size,
            min_samples=min_samples,
            cluster_selection=cluster_selection,
            use_query_columns=use_query_columns,
        )

        # reindex sample cluster ids to 0..C-1
        uniq = sorted(list(set(sample_labels) - {-1}))
        remap = {old: i for i, old in enumerate(uniq)}
        sample_labels = np.array([remap[lab] if lab in remap else -1 for lab in sample_labels], dtype=int)

        if len(uniq) == 0:
            # No structure in sample; fallback to putting everything as noise (or one big cluster?)
            print("[WARN] No clusters found in sample; writing empty metadata (all noise).")
            meta_csv = os.path.join(outdir, f"{keyword}_clusters.csv")
            pd.DataFrame(columns=["cluster", "kept", "n", "neff", "path"]).to_csv(meta_csv, index=False)
            return pd.DataFrame(columns=["cluster", "kept", "n", "neff", "path"])

        # 3) medoids on sample clusters
        medoid_local = _medoids_from_labels(D_sample, sample_labels)  # {c: idx_in_sample}
        center_sample_idxs = [medoid_local[c] for c in sorted(medoid_local.keys())]
        center_global_idxs = [int(sample_idx[i]) for i in center_sample_idxs]
        centers = [seqs[i] for i in center_global_idxs]
        # 4) assign rest to nearest centers
        rest_seqs = [seqs[i] for i in rest_idx.tolist()]
        assigned_centers = _assign_rest_to_medoids(
            rest_seqs=rest_seqs,
            centers=centers,
            use_query_columns=use_query_columns,
            max_center_dist=max_center_dist
        )
        # assemble global labels: initialize all -1, fill sample labels and rest assignments
        labels = np.full(n, -1, dtype=int)
        # sample labels already 0..C-1
        labels[sample_idx] = sample_labels
        # rest: assigned_centers indexes are 0..C-1 (order of 'centers')
        if rest_idx.size > 0:
            labels[rest_idx] = np.where(assigned_centers >= 0, assigned_centers, -1)

    else:
        # ---------- original behavior ----------
        if metric == "hamming":
            USE_PRECOMPUTED = (n <= max_n_for_precomputed) and (not force_feature_mode)
            if USE_PRECOMPUTED:
                D = dist_pid_gapaware_matrix(seqs, use_query_columns=use_query_columns)
                D = np.asarray(D, dtype=np.float64, order="C")
                hdb = HDBSCAN(min_cluster_size=min_cluster_size,
                              min_samples=min_samples,
                              cluster_selection_method=cluster_selection,
                              metric='precomputed',
                              core_dist_n_jobs=0)
                labels = hdb.fit_predict(D)
            else:
                # feature-mode: encode residues to ints, gaps distinct
                arr = np.array([list(s) for s in seqs], dtype='<U1')
                if use_query_columns:
                    keep_cols = (arr[0] != '-')
                    arr = arr[:, keep_cols]
                alphabet = np.array(list("ACDEFGHIKLMNPQRSTVWYBJZX-"), dtype='<U1')
                codebook = {ch: i for i, ch in enumerate(alphabet)}
                GAP_CODE = 255
                X = np.empty(arr.shape, dtype=np.uint8)
                flat = arr.ravel()
                mapped = np.fromiter((codebook.get(ch, GAP_CODE) for ch in flat), count=flat.size, dtype=np.uint16)
                X[:, :] = mapped.reshape(arr.shape).astype(np.uint8)
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
                                   "for n>4000 switch to --metric hamming.")
            # simple normalized Levenshtein
            L = len(seqs[0])
            D = np.zeros((n, n), dtype=np.float64)
            for i in range(n):
                si = seqs[i]
                for j in range(i + 1, n):
                    sj = seqs[j]
                    if HAS_POLYLEVEN:
                        d = _levenshtein(si, sj)
                    else:
                        d = _plain_lev(si, sj)
                    D[i, j] = D[j, i] = d / float(L)
            hdb = HDBSCAN(min_cluster_size=min_cluster_size,
                          min_samples=min_samples,
                          cluster_selection_method=cluster_selection,
                          metric='precomputed',
                          core_dist_n_jobs=0)
            labels = hdb.fit_predict(D)
        else:
            raise ValueError("metric must be 'hamming' or 'lev'.")

    labels = np.asarray(labels, dtype=int)

    # -------- Post-filters: drop small clusters, cap K, reindex 0..C-1 --------
    valid = labels >= 0
    cluster_ids, counts = np.unique(labels[valid], return_counts=True)
    small = set(cluster_ids[counts < min_output_size])
    if small:
        labels = np.array([-1 if (lab in small) else lab for lab in labels], dtype=int)

    cluster_ids, counts = np.unique(labels[labels >= 0], return_counts=True)
    if len(cluster_ids) > max_clusters:
        order = np.argsort(counts)[::-1]
        keep = set(cluster_ids[order[:max_clusters]])
        labels = np.array([lab if lab in keep else -1 for lab in labels], dtype=int)

    uniq = sorted(list(set(labels) - {-1}))
    remap = {old: i for i, old in enumerate(uniq)}
    labels = np.array([remap.get(lab, -1) for lab in labels], dtype=int)

    # -------- Write clusters & compute Meff --------
    if len(uniq) == 0:
        meta_csv = os.path.join(outdir, f"{keyword}_clusters.csv")
        pd.DataFrame(columns=["cluster", "kept", "n", "neff", "path"]).to_csv(meta_csv, index=False)
        print(f"[INFO] No clusters found; wrote empty {meta_csv}")
        return pd.DataFrame(columns=["cluster", "kept", "n", "neff", "path"])

    # clean previous outputs
    for old in glob(os.path.join(outdir, f"{keyword}_*.a3m")):
        try: os.remove(old)
        except: pass

    rows = []
    n_noise = int((labels < 0).sum())
    print(f"[INFO] clusters found (before Meff): {len(uniq)} ; noise={n_noise}")

    for c in sorted(set(labels) - {-1}):
        idx = np.where(labels == c)[0]
        # write with query first
        cl_headers = [headers[0]] + [headers[i] for i in idx if i != 0]
        cl_seqs = [seqs[0]] + [seqs[i] for i in idx if i != 0]
        out_path = os.path.join(outdir, f"{keyword}_{c:03d}.a3m")
        write_fasta(cl_headers, cl_seqs, outfile=out_path)

        meff = compute_neff_from_seqs(cl_seqs, id_thresh=neff_id_thresh,
                                      use_query_columns=True)
        n_raw = len(cl_seqs)
        kept = True
        if (min_neff is not None and min_neff > 0) and (meff < min_neff):
            kept = False
            try: os.remove(out_path)
            except: pass
        rows.append(dict(cluster=c, kept=kept, n=n_raw, neff=meff, path=out_path))

    df_meta = pd.DataFrame(rows)
    if not df_meta.empty:
        df_meta = df_meta.sort_values(["kept", "n"], ascending=[False, False])
    meta_csv = os.path.join(outdir, f"{keyword}_clusters.csv")
    df_meta.to_csv(meta_csv, index=False)
    print(f"[INFO] Wrote cluster metadata: {meta_csv}")
    return df_meta


# ---------- Fallback Levenshtein ----------
def _plain_lev(a: str, b: str) -> int:
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


# ---------- CLI ----------
def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Cluster an A3M with HDBSCAN (sample→cluster→assign heuristic).")
    p.add_argument("--a3m", required=True, help="Input A3M path.")
    p.add_argument("-o", "--outdir", required=True, help="Output directory.")
    p.add_argument("--keyword", default="ShallowMsa", help="Prefix for cluster files.")
    p.add_argument("--metric", choices=["hamming", "lev"], default="hamming",
                   help="Distance metric for clustering (precomputed for sample; gap-aware Hamming recommended).")

    # HDBSCAN core params
    p.add_argument("--min_cluster_size", type=int, default=200,
                   help="HDBSCAN min_cluster_size.")
    p.add_argument("--min_samples", type=str, default="auto",
                   help="'auto' or integer (density strictness).")
    p.add_argument("--cluster_selection", choices=["eom","leaf"], default="leaf",
                   help="Cluster selection method; 'leaf' yields finer clusters.")

    # Output controls
    p.add_argument("--max_clusters", type=int, default=100,
                   help="Keep only the largest K clusters; rest -> noise.")
    p.add_argument("--min_output_size", type=int, default=200,
                   help="Drop clusters with raw size < this after assignment.")
    p.add_argument("--min_neff", type=int, default=100,
                   help="Drop clusters whose Meff < this threshold.")
    p.add_argument("--neff_id_thresh", type=float, default=0.8,
                   help="Identity threshold for Meff (e.g. 0.8).")

    # MSA processing
    p.add_argument("--frac_gaps_cutoff", type=float, default=0.5,
                   help="Filter out sequences with fraction of gaps >= cutoff (on full aligned length).")
    p.add_argument("--use_query_columns", type=int, default=1,
                   help="1: compute distances/Meff on columns where query has residues; 0: all columns.")

    # Heuristic controls
    p.add_argument("--sample_cap", type=int, default=10000,
                   help="Max number of sequences to cluster exactly (precomputed matrix) in the sample.")
    p.add_argument("--sample_seed", type=int, default=12345,
                   help="Random seed for sampling.")
    p.add_argument("--max_center_dist", type=float, default=None,
                   help="If set (0..1), assign rest to nearest center only if distance <= this; else leave as noise.")

    # Legacy/alt switches
    p.add_argument("--max_n_for_precomputed", type=int, default=6000,
                   help="If not using heuristic, limit for O(n^2) precomputed hamming.")
    p.add_argument("--force_feature_mode", type=int, default=0,
                   help="Set 1 to force feature-mode Hamming for the non-heuristic branch.")
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
        sample_cap=int(args.sample_cap),
        sample_seed=int(args.sample_seed),
        max_center_dist=(None if args.max_center_dist is None else float(args.max_center_dist)),
        max_n_for_precomputed=int(args.max_n_for_precomputed),
        force_feature_mode=int(args.force_feature_mode),
    )

    kept = df[df.kept] if not df.empty else df
    dropped = df[~df.kept] if not df.empty else df
    print(f"[SUMMARY] Kept {0 if df.empty else len(kept)} clusters,"
          f" Dropped {0 if df.empty else len(dropped)} (by Neff).")
    if not df.empty and not kept.empty:
        print(kept[["cluster","n","neff","path"]].to_string(index=False))


if __name__ == "__main__":
    main()
