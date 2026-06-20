#!/usr/bin/env python3
"""Test whether per-cluster fold preference (F1c/F2c/Amc) is phylogenetically
structured for each fold-switching pair, on the COARSE cluster-level tree.

For every pair p and every method m in {AF2, ESM, MSAT, DDG}:

  1. Build a COARSE cluster-level tree with one leaf per ShallowMsa_NNN
     cluster (plus DeepMsa, when present), via
     ``utils.ancestral_utils.build_coarse_tree``. That helper
       a. loads the per-sequence Newick at
          Pipeline/FoldPairs/<p>/output_phytree/DeepMsa_tree.nwk,
       b. assigns each leaf to its cluster via the per-cluster a3m files
          at Pipeline/FoldPairs/<p>/output_msa_cluster/ShallowMsa_*.a3m,
       c. picks one representative leaf per cluster and prunes the tree,
       d. renames each surviving leaf 'S00', 'S01', ..., 'D' (DeepMsa).
  2. Read each (pair, method, cluster)'s centered preference from
     docs/fold_diversity_survey.csv (column ``pref_centered``, with a
     fallback to recomputing it on the fly if the CSV lacks that column).
  3. Compute three phylogenetic-structure statistics on the F1c vs F2c
     binary partition of the COARSE tree (Amc leaves treated as missing):
       a. Fitch parsimony score (count of state changes along edges)
       b. k=1 / k=3 nearest-neighbour concordance via cophenetic distance
       c. Phylogenetic D statistic (Fritz & Purvis 2010), if computable
     Each statistic gets an empirical permutation-null p-value. The null
     shuffles the F1c/F2c labels among the cluster-level leaves; the
     coarse-tree shuffle is more conservative than a sequence-level shuffle
     because the number of "free" labels equals the number of clusters
     (typically 8-15), not the number of sequences (hundreds).
  4. Output one row per (pair, method) to docs/phylo_placement.csv.

Per-pair figures show the COARSE tree (10-20 nodes) with each leaf
coloured by F1c / F2c / Amc.

WHY COARSE?
-----------
The fold-preference label is defined PER CLUSTER, so the natural unit of
analysis is the cluster, not the sequence. Running the same statistics on
a fine per-sequence tree
   * makes within-cluster leaves trivially share a label, inflating
     NN-concordance toward 1.0;
   * leaves Fitch parsimony invariant (a constant subtree contributes 0),
     but makes the permutation null too easy to reject because most
     "shuffles" of sequence-level labels still preserve the cluster
     structure that the survey already imposed.
The coarse tree fixes both issues and produces a much cleaner figure.

Usage
-----
    python3 scripts/phylo_placement.py
    python3 scripts/phylo_placement.py --pairs 2namA_1uxmK,2n54B_2hdmA
    python3 scripts/phylo_placement.py --pairs ALL --make-figures all
"""

from __future__ import annotations

import os
# Headless rendering env vars MUST be set BEFORE importing matplotlib / ete3.
# Without these, on a headless compute node ete3.TreeStyle hangs waiting on Qt
# and matplotlib fails to find a writable config dir.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-cache-orzuk")
os.environ.setdefault("PYTHONUNBUFFERED", "1")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)

import argparse
import random
import re
import signal
import sys
import time
from collections import Counter
from contextlib import contextmanager
from typing import Dict, List, Optional, Sequence, Tuple


@contextmanager
def _per_pair_timeout(seconds: int):
    """SIGALRM-based per-pair timeout. Skips a single hung pair without
    killing the whole run. Linux-only; uses SIGALRM."""
    def _handler(signum, frame):
        raise TimeoutError(f"per-pair timeout after {seconds}s")
    old = signal.signal(signal.SIGALRM, _handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)

import numpy as np
import pandas as pd

# Allow running from repo root and from anywhere
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

print(f"[phylo] importing modules...", flush=True)
_t0 = time.time()
from config import DATA_DIR, PAIR_DIR_RE, TABLES_RES  # noqa: E402
from utils.ancestral_utils import build_coarse_tree, _name_internal_nodes  # noqa: E402
print(f"[phylo] imports done in {time.time()-_t0:.1f}s", flush=True)


# Methods we score. Methods missing from the survey CSV are silently
# skipped (n_pairs=0). Add new methods here when fold_diversity_survey.py
# starts emitting them.
SCORED_METHODS: Tuple[str, ...] = (
    "AF2", "AF3", "Boltz2",          # MSA-conditioned structure predictors
    "BioEmu",                         # generative equilibrium-ensemble (9th method)
    "ESM",                            # single-sequence ESMFold (v1) or ESMFold2 (v2)
    "MSAT", "CCMpred",                # contact-map predictors
    "DDG",                            # Rosetta ΔΔG per cluster
    "S4PRED",                         # secondary-structure-based fold preference
)


# ---------------------------------------------------------------------------
# Cluster-tag normalization (matches scripts/fold_diversity_survey.py)
# ---------------------------------------------------------------------------

_MSA_TAG_RE = re.compile(r"(?:Shallow)?Msa_(\d+)", re.IGNORECASE)


def normalize_cluster(s: str) -> str:
    """Normalize a cluster tag to either 'ShallowMsa_NNN' or 'DeepMsa'."""
    s = str(s).strip()
    s = re.sub(r"^(msa_t_)?(clusters_|cmap_)?", "", s)
    if s.startswith("DeepMsa") or s.lower() in (
        "deep", "deepmsa", "msa_deep", "deep_msa", "query"
    ):
        return "DeepMsa"
    m = _MSA_TAG_RE.fullmatch(s)
    if m:
        return f"ShallowMsa_{int(m.group(1)):03d}"
    m2 = re.search(r"ShallowMsa_(\d+)", s)
    if m2:
        return f"ShallowMsa_{int(m2.group(1)):03d}"
    if s.isdigit():
        return f"ShallowMsa_{int(s):03d}"
    return s


# ---------------------------------------------------------------------------
# Fold-preference ingestion (from fold_diversity_survey.csv)
# ---------------------------------------------------------------------------

def bh_adjust(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (NaN-safe, monotone)."""
    p = np.asarray(pvals, dtype=float)
    mask = np.isfinite(p)
    out = np.full_like(p, np.nan, dtype=float)
    if not mask.any():
        return out
    pv = p[mask]
    n = len(pv)
    order = np.argsort(pv)
    ranks = np.arange(1, n + 1)
    adj_sorted = pv[order] * n / ranks
    adj_sorted = np.minimum.accumulate(adj_sorted[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.clip(adj_sorted, 0.0, 1.0)
    out[mask] = adj
    return out


def load_corrected_pref(survey_csv: str) -> Dict[Tuple[str, str, str], str]:
    """Build (pair, method, cluster) -> residual-corrected preference label.

    Reads the 'pref_corrected' column produced by
    scripts/seq_divergence_correction.py from
    docs/fold_diversity_survey_corrected.csv.
    Used when --labels corrected is selected; corresponds to F1c/F2c/Amc
    labels computed from the regression residuals (PID-to-query stripped).
    """
    df = pd.read_csv(survey_csv)
    if "pref_corrected" not in df.columns:
        raise SystemExit(
            f"Expected column 'pref_corrected' in {survey_csv}. "
            "Run scripts/seq_divergence_correction.py first."
        )
    cluster_col = "cluster_norm" if "cluster_norm" in df.columns else "cluster"
    df[cluster_col] = df[cluster_col].astype(str).apply(normalize_cluster)

    out: Dict[Tuple[str, str, str], str] = {}
    for _, r in df.iterrows():
        label = r.get("pref_corrected")
        if not isinstance(label, str) or label not in ("F1", "F2", "Amb"):
            label = "Amb"
        key = (str(r["pair_id"]), str(r["method"]), str(r[cluster_col]))
        out[key] = label
    return out


def load_fold_value(survey_csv: str, value_col: str
                    ) -> Dict[Tuple[str, str, str], Optional[float]]:
    """Build (pair, method, cluster) -> CONTINUOUS signed fold preference.

    Used by the Moran's-I phylogenetic-signal test (no thresholding). The
    value is the method's native signed preference; Moran's I is invariant
    to per-method affine rescaling, so no normalization is needed.

      raw labels       -> value_col = 'TMdiff_centered' (median-centred)
      corrected labels -> value_col = 'TMdiff_residual' (seq-divergence
                          regression residual; produced by
                          scripts/seq_divergence_correction.py)
    """
    df = pd.read_csv(survey_csv)
    if value_col not in df.columns:
        raise SystemExit(
            f"Expected continuous column '{value_col}' in {survey_csv}."
        )
    cluster_col = "cluster_norm" if "cluster_norm" in df.columns else "cluster"
    df[cluster_col] = df[cluster_col].astype(str).apply(normalize_cluster)
    out: Dict[Tuple[str, str, str], Optional[float]] = {}
    for _, r in df.iterrows():
        v = pd.to_numeric(r.get(value_col), errors="coerce")
        key = (str(r["pair_id"]), str(r["method"]), str(r[cluster_col]))
        out[key] = float(v) if pd.notna(v) else None
    return out


def load_fold_pref(survey_csv: str, delta_tm: float = 0.05,
                   delta_ddg: float = 1.0
                   ) -> Dict[Tuple[str, str, str], str]:
    """Build (pair, method, cluster) -> centered preference (F1/F2/Amb).

    If the survey CSV already has the centered column ``pref_centered``,
    use it directly. Otherwise we recompute centering on the fly so older
    survey CSVs still work.
    """
    df = pd.read_csv(survey_csv)
    df["cluster"] = df["cluster"].apply(normalize_cluster)

    if "pref_centered" not in df.columns:
        df["pref_centered"] = "Amb"
        for (pair_id, method), grp in df.groupby(["pair_id", "method"]):
            shallow = grp[grp["cluster"] != "DeepMsa"]
            tdvals = pd.to_numeric(shallow.get("TMdiff_max"), errors="coerce")
            tdvals = tdvals.dropna()
            if len(tdvals) >= 2:
                med = float(np.median(tdvals))
            else:
                med = 0.0
            delta = delta_ddg if method == "DDG" else delta_tm
            for idx, r in grp.iterrows():
                tdm = pd.to_numeric(r.get("TMdiff_max"), errors="coerce")
                if pd.isna(tdm):
                    continue
                centered = float(tdm) - (med if r["cluster"] != "DeepMsa" else 0.0)
                if centered > delta:
                    df.at[idx, "pref_centered"] = "F1"
                elif centered < -delta:
                    df.at[idx, "pref_centered"] = "F2"
                else:
                    df.at[idx, "pref_centered"] = "Amb"

    out: Dict[Tuple[str, str, str], str] = {}
    for _, r in df.iterrows():
        key = (str(r["pair_id"]), str(r["method"]), str(r["cluster"]))
        pc = r.get("pref_centered")
        if not isinstance(pc, str) or pc not in ("F1", "F2", "Amb"):
            pc = "Amb"
        out[key] = pc
    return out


# ---------------------------------------------------------------------------
# Phylogenetic-structure statistics on the coarse tree
# ---------------------------------------------------------------------------

def _fitch_parsimony_binary(tree, leaf_states: Dict[str, str]) -> int:
    """Fitch small-parsimony score for a binary trait {'F1', 'F2'}.

    Leaves with state == 'Amb' or absent from the dict are treated as
    missing data (empty character set, no constraint on internal nodes).

    Score = number of state changes on the bottom-up Fitch pass, computed
    as (k-1) for each internal node where the children's character sets
    fall into k disjoint groups (typically 1 or 2 for a binary trait).
    """
    score = 0
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            s = leaf_states.get(node.name)
            node._fitch_set = {s} if s in ("F1", "F2") else set()
            continue
        present = [c._fitch_set for c in node.children if c._fitch_set]
        if not present:
            node._fitch_set = set()
            continue
        inter = set.intersection(*present)
        if inter:
            node._fitch_set = inter
        else:
            node._fitch_set = set.union(*present)
            groups: List[set] = []
            for s in present:
                placed = False
                for g in groups:
                    if g & s:
                        g |= s
                        placed = True
                        break
                if not placed:
                    groups.append(set(s))
            score += max(0, len(groups) - 1)
    return score


def _cophenetic_neighbours(tree, leaves_with_state: Sequence[str]
                           ) -> Dict[str, List[Tuple[float, str]]]:
    """For each leaf, sorted (distance, other-leaf) list to other leaves
    in ``leaves_with_state``. Distance uses branch lengths; if the tree
    has no branch lengths, every edge counts as length 1.
    """
    name_to_node: dict = {}
    depth: dict = {}

    def _walk(node, d):
        name_to_node[node.name] = node
        depth[node.name] = d
        for c in node.children:
            cd = c.dist if (c.dist and c.dist > 0) else 1.0
            _walk(c, d + cd)

    _walk(tree, 0.0)

    ancestors: Dict[str, List[str]] = {}
    for name in leaves_with_state:
        node = name_to_node.get(name)
        chain: List[str] = []
        while node is not None:
            chain.append(node.name)
            node = node.up
        ancestors[name] = chain

    out: Dict[str, List[Tuple[float, str]]] = {}
    leaves_list = list(leaves_with_state)
    for i, a in enumerate(leaves_list):
        a_anc_set = {n: k for k, n in enumerate(ancestors[a])}
        a_depth = depth.get(a, 0.0)
        nbrs: List[Tuple[float, str]] = []
        for j, b in enumerate(leaves_list):
            if i == j:
                continue
            lca = None
            for n in ancestors[b]:
                if n in a_anc_set:
                    lca = n
                    break
            if lca is None:
                continue
            d_ab = (a_depth - depth.get(lca, 0.0)) + (depth.get(b, 0.0) - depth.get(lca, 0.0))
            nbrs.append((d_ab, b))
        nbrs.sort(key=lambda x: (x[0], x[1]))
        out[a] = nbrs
    return out


def nn_concordance(tree, leaf_states: Dict[str, str], k: int = 1
                   ) -> Tuple[float, int]:
    """Mean fraction of a leaf's k tree-nearest neighbours sharing its
    fold preference. Only F1/F2 leaves participate (Amc treated as missing).
    """
    queryable = [n for n, s in leaf_states.items() if s in ("F1", "F2")]
    if len(queryable) < 2:
        return float("nan"), len(queryable)

    nbrs = _cophenetic_neighbours(tree, queryable)
    fracs: List[float] = []
    for leaf, sorted_nbrs in nbrs.items():
        own = leaf_states[leaf]
        considered: List[str] = []
        for _, other in sorted_nbrs:
            os_state = leaf_states.get(other)
            if os_state in ("F1", "F2"):
                considered.append(os_state)
                if len(considered) >= k:
                    break
        if not considered:
            continue
        fracs.append(sum(1 for s in considered if s == own) / len(considered))
    if not fracs:
        return float("nan"), len(queryable)
    return float(np.mean(fracs)), len(queryable)


def d_statistic(tree, leaf_states: Dict[str, str],
                n_perm: int = 200, seed: int = 0
                ) -> Optional[float]:
    """Approximate Fritz & Purvis (2010) D statistic for a binary trait.

    D = 0  -> tightly clumped (Brownian-motion-like)
    D = 1  -> phylogenetically random
    D < 0  -> MORE clumped than Brownian
    D > 1  -> overdispersed

    Implemented as a Monte-Carlo two-null comparison (random shuffle vs
    a tree-correlated Brownian-motion proxy).
    """
    binary = {n: (1 if s == "F1" else 0) for n, s in leaf_states.items()
              if s in ("F1", "F2")}
    if len(binary) < 4:
        return None
    n_f1 = sum(binary.values())
    n_f2 = len(binary) - n_f1
    if n_f1 == 0 or n_f2 == 0:
        return None

    def _sum_sister_changes(states: Dict[str, int]) -> float:
        subtree_sum: Dict[int, float] = {}
        subtree_n: Dict[int, int] = {}
        for node in tree.traverse("postorder"):
            if node.is_leaf():
                v = states.get(node.name)
                if v is None:
                    subtree_sum[id(node)] = 0.0
                    subtree_n[id(node)] = 0
                else:
                    subtree_sum[id(node)] = float(v)
                    subtree_n[id(node)] = 1
            else:
                s = sum(subtree_sum[id(c)] for c in node.children)
                n = sum(subtree_n[id(c)] for c in node.children)
                subtree_sum[id(node)] = s
                subtree_n[id(node)] = n
        total = 0.0
        for node in tree.traverse("preorder"):
            if node.is_root():
                continue
            n_self = subtree_n[id(node)]
            n_parent = subtree_n[id(node.up)]
            if n_self == 0 or n_parent == 0:
                continue
            mean_self = subtree_sum[id(node)] / n_self
            mean_parent = subtree_sum[id(node.up)] / n_parent
            total += abs(mean_self - mean_parent)
        return total

    obs = _sum_sister_changes(binary)
    rng = random.Random(seed)
    leaf_names = list(binary.keys())
    leaf_values = list(binary.values())

    rand_vals: List[float] = []
    for _ in range(n_perm):
        v = leaf_values.copy()
        rng.shuffle(v)
        m = dict(zip(leaf_names, v))
        rand_vals.append(_sum_sister_changes(m))
    mean_rand = float(np.mean(rand_vals))

    bm_vals: List[float] = []
    for _ in range(n_perm):
        bm_at: Dict[str, float] = {}

        def _walk(node, val):
            for c in node.children:
                step = rng.gauss(0.0, max(1.0, c.dist if c.dist else 1.0) ** 0.5)
                v = val + step
                if c.is_leaf():
                    bm_at[c.name] = v
                else:
                    _walk(c, v)
        _walk(tree, 0.0)

        named_vals = sorted(bm_at.items(), key=lambda x: -x[1])
        m = {}
        for k, (nm, _) in enumerate(named_vals):
            if nm in binary:
                m[nm] = 1 if k < n_f1 else 0
        bm_vals.append(_sum_sister_changes(m))
    mean_bm = float(np.mean(bm_vals))

    if mean_rand - mean_bm == 0:
        return None
    D = (obs - mean_bm) / (mean_rand - mean_bm)
    return float(D)


# ---------------------------------------------------------------------------
# Permutation null
# ---------------------------------------------------------------------------

def permutation_p(tree, leaf_states: Dict[str, str],
                  stat_fn, observed: float,
                  n_perm: int = 1000, seed: int = 12345,
                  smaller_is_more_structured: bool = True
                  ) -> float:
    """Empirical p-value via leaf-label shuffle on the coarse tree.

    Shuffles the F1/F2 labels among the leaves that originally had F1/F2
    (Amb/missing leaves stay missing). Returns
    ``(#perms with stat at least as extreme + 1) / (n_perm + 1)``.
    """
    rng = random.Random(seed)
    names = [n for n, s in leaf_states.items() if s in ("F1", "F2")]
    states = [leaf_states[n] for n in names]
    if len(names) < 2:
        return float("nan")

    n_extreme = 0
    for _ in range(n_perm):
        v = states.copy()
        rng.shuffle(v)
        perm_states = dict(leaf_states)
        for nm, st in zip(names, v):
            perm_states[nm] = st
        s = stat_fn(tree, perm_states)
        if isinstance(s, tuple):
            s = s[0]
        if smaller_is_more_structured:
            if s <= observed:
                n_extreme += 1
        else:
            if not (isinstance(s, float) and np.isnan(s)) and s >= observed:
                n_extreme += 1
    return (n_extreme + 1) / (n_perm + 1)


# ---------------------------------------------------------------------------
# Continuous phylogenetic signal: Moran's I + permutation null
# ---------------------------------------------------------------------------

def _patristic_matrix(tree, names: Sequence[str]) -> np.ndarray:
    """Pairwise patristic (branch-length) distance matrix over ``names``.

    Reuses the same root-depth + LCA logic as ``_cophenetic_neighbours``.
    Edges with no/zero length count as 1.0 (topological distance fallback).
    """
    name_to_node: dict = {}
    depth: dict = {}

    def _walk(node, d):
        name_to_node[node.name] = node
        depth[node.name] = d
        for c in node.children:
            cd = c.dist if (c.dist and c.dist > 0) else 1.0
            _walk(c, d + cd)

    _walk(tree, 0.0)
    ancestors: Dict[str, List[str]] = {}
    for nm in names:
        node = name_to_node.get(nm)
        chain: List[str] = []
        while node is not None:
            chain.append(node.name)
            node = node.up
        ancestors[nm] = chain

    n = len(names)
    D = np.zeros((n, n), dtype=float)
    for i, a in enumerate(names):
        a_anc = {nm: k for k, nm in enumerate(ancestors[a])}
        da = depth.get(a, 0.0)
        for j in range(i + 1, n):
            b = names[j]
            lca = next((nm for nm in ancestors[b] if nm in a_anc), None)
            if lca is None:
                d_ab = da + depth.get(b, 0.0)            # disconnected fallback
            else:
                d_ab = (da - depth[lca]) + (depth.get(b, 0.0) - depth[lca])
            D[i, j] = D[j, i] = d_ab
    return D


def morans_I_test(tree, leaf_values: Dict[str, Optional[float]],
                  n_perm: int = 1000, seed: int = 12345,
                  eps: float = 1e-9) -> dict:
    """Moran's I phylogenetic autocorrelation of a CONTINUOUS leaf trait,
    with a one-sided permutation p-value (testing for positive signal =
    clade-structured fold preference).

    - Weights: row-standardized inverse patristic distance (1/d_ij).
    - I > 0: phylogenetically nearby clusters have similar preference.
    - p = (1 + #{I_perm >= I_obs}) / (1 + n_perm); floor 1/(n_perm+1).
    - Invariant to affine rescaling of the trait (so methods on different
      scales are directly comparable without normalization).
    - Returns NaN with a note when n<3 or the trait has zero variance
      (e.g. CCMpred): a constant trait has no signal to measure.
    """
    names = [n for n, v in leaf_values.items()
             if v is not None and np.isfinite(v)]
    res = {"morans_I": np.nan, "morans_p": np.nan, "morans_z": np.nan,
           "morans_lambda_max": np.nan, "morans_max_perm": np.nan,
           "morans_I_norm": np.nan,
           "morans_n": len(names), "morans_note": ""}
    if len(names) < 3:
        res["morans_note"] = "n<3"
        return res
    x = np.asarray([leaf_values[m] for m in names], dtype=float)
    xc = x - x.mean()
    den = float(xc @ xc)
    if den < eps:
        res["morans_note"] = "zero_variance"
        return res

    W = 1.0 / (_patristic_matrix(tree, names) + eps)
    np.fill_diagonal(W, 0.0)
    rs = W.sum(axis=1, keepdims=True)
    rs[rs == 0] = 1.0
    W = W / rs                                   # row-standardized -> S0 = n

    I_obs = float((xc @ (W @ xc)) / den)
    rng = np.random.default_rng(seed)
    # Vectorized permutation null: each row is a permutation of xc.
    P = np.empty((n_perm, len(xc)), dtype=float)
    for k in range(n_perm):
        P[k] = rng.permutation(xc)
    I_perm = np.einsum("ij,ij->i", P @ W, P) / den
    p = (1 + int((I_perm >= I_obs).sum())) / (n_perm + 1)
    # Standardized Moran's I: deviation from THIS pair's own permutation null,
    # in null-SD units. Comparable across pairs of different n (unlike raw I,
    # whose null E[I]=-1/(n-1) shifts with n) and does not floor like p.
    mu, sd = float(I_perm.mean()), float(I_perm.std())
    z = (I_obs - mu) / sd if sd > 1e-12 else float("nan")
    res.update({"morans_I": round(I_obs, 4), "morans_p": round(float(p), 4),
                "morans_z": round(z, 3) if np.isfinite(z) else np.nan})

    # --- Achievable-max ceiling + normalized effect size ---
    # Moran's I is bounded by the tree's weight matrix, NOT by 1. Report I as a
    # fraction of the max achievable for THIS value multiset on THIS tree:
    #   lambda_max  = eigenvalue ceiling (exact upper bound over centered vectors)
    #   max_perm    = data-specific achievable max (sort-and-match on the top
    #                 eigenvector + capped vectorized 2-opt)
    #   I_norm      = (I_obs - E[I]) / (max_perm - E[I]),  E[I] = -1/(n-1)
    try:
        nn = len(xc)
        M = 0.5 * (W + W.T)                                  # symmetric quadratic form
        Pc = np.eye(nn) - np.ones((nn, nn)) / nn             # centering projector
        Mc = Pc @ M @ Pc
        evals, evecs = np.linalg.eigh(Mc)
        lam_max = float(evals[-1])

        def _two_opt(seed, maximize):
            """Capped steepest-ascent/-descent 2-opt on the quadratic form
            arr @ M @ arr; returns the optimized Moran's-I value."""
            arr = np.array(seed, dtype=float)
            for _ in range(60):
                f = M @ arr
                DZ = arr[None, :] - arr[:, None]
                DF = f[:, None] - f[None, :]
                dN = 2.0 * DZ * DF - 2.0 * (DZ ** 2) * M    # delta of form on swap (a,b)
                if maximize:
                    np.fill_diagonal(dN, -np.inf)
                    a, b = np.unravel_index(int(np.argmax(dN)), dN.shape)
                    if dN[a, b] <= 1e-12:
                        break
                else:
                    np.fill_diagonal(dN, np.inf)
                    a, b = np.unravel_index(int(np.argmin(dN)), dN.shape)
                    if dN[a, b] >= -1e-12:
                        break
                arr[a], arr[b] = arr[b], arr[a]
            return float(arr @ (M @ arr)) / den

        # Seed the search from BOTH the spectral guess (top/bottom eigenvector)
        # AND the actual data arrangement. Ascending from the data is guaranteed
        # to reach I >= I_obs (and descending <= I_obs), so the data-seeded run
        # alone makes max_perm >= I_obs >= min_perm by construction -- no clamp,
        # no |I_norm|>1 on tiny/degenerate pairs (few clusters: MSAT/CCMpred).
        seed_top = np.empty(nn); seed_top[np.argsort(evecs[:, -1])] = np.sort(xc)
        seed_bot = np.empty(nn); seed_bot[np.argsort(evecs[:, 0])]  = np.sort(xc)
        max_perm_I = max(_two_opt(seed_top, True),  _two_opt(xc, True))
        min_perm_I = min(_two_opt(seed_bot, False), _two_opt(xc, False))

        # Normalize by the achievable range on the side I_obs falls -- positive
        # signal against the reachable max, negative against the reachable min --
        # so I_norm is a true fraction of achievable clade signal in [-1, 1].
        Emean = -1.0 / (nn - 1)
        if I_obs >= Emean:
            denom = max_perm_I - Emean
        else:
            denom = Emean - min_perm_I
        I_norm = (I_obs - Emean) / denom if denom > 1e-9 else float("nan")
        if np.isfinite(I_norm):
            I_norm = float(np.clip(I_norm, -1.0, 1.0))   # belt-and-suspenders
        res.update({
            "morans_lambda_max": round(lam_max, 4),
            "morans_max_perm": round(max_perm_I, 4),
            "morans_I_norm": round(float(I_norm), 4) if np.isfinite(I_norm) else np.nan,
        })
    except Exception:
        pass
    return res


# ---------------------------------------------------------------------------
# Per-pair, per-method analysis on the COARSE tree
# ---------------------------------------------------------------------------

def analyze_pair_method(pair_id: str, method: str,
                        tree, tag_to_label: Dict[str, str],
                        cluster_pref: Dict[str, str],
                        n_perm: int, seed: int,
                        do_d_stat: bool = True,
                        cluster_value: Optional[Dict[str, Optional[float]]] = None,
                        do_discrete_perms: bool = False
                        ) -> Optional[dict]:
    """Compute phylo-structure stats for one (pair, method) on the COARSE tree.

    ``tag_to_label`` maps each cluster tag (ShallowMsa_NNN / DeepMsa) to
    the short leaf label used in the coarse tree (e.g. 'S03', 'D').
    ``cluster_pref`` is {cluster_tag: 'F1'|'F2'|'Amb'}.
    """
    # Map coarse-tree leaf name -> centered preference
    leaf_states: Dict[str, str] = {}
    for tag, label in tag_to_label.items():
        leaf_states[label] = cluster_pref.get(tag, "Amb")
    for leaf in tree.iter_leaves():
        leaf_states.setdefault(leaf.name, "Amb")

    counts = Counter(leaf_states.values())
    n_f1 = counts.get("F1", 0)
    n_f2 = counts.get("F2", 0)
    n_amb = counts.get("Amb", 0)
    n_total = sum(counts.values())

    notes: List[str] = []
    if n_f1 == 0 or n_f2 == 0:
        notes.append("no_F1F2_split")

    out = {
        "pair_id": pair_id,
        "method": method,
        "n_leaves": n_total,                 # = n_clusters in coarse tree
        "n_F1c_leaves": n_f1,
        "n_F2c_leaves": n_f2,
        "n_Amc_leaves": n_amb,
        "parsimony_score": np.nan,
        "parsimony_p": np.nan,
        "nn_concordance": np.nan,
        "nn_concordance_p": np.nan,
        "nn3_concordance": np.nan,
        "nn3_concordance_p": np.nan,
        "D_statistic": np.nan,
        "morans_I": np.nan,
        "morans_p": np.nan,
        "morans_z": np.nan,
        "morans_lambda_max": np.nan,
        "morans_max_perm": np.nan,
        "morans_I_norm": np.nan,
        "morans_n": 0,
        "morans_note": "",
        "notes": "",
    }

    # ---- PRIMARY: continuous Moran's I (no thresholding) ----
    # Uses the raw/corrected continuous preference per cluster; every cluster
    # participates (no F1/F2/Amb gate), so coverage is method-independent.
    if cluster_value is not None:
        leaf_values: Dict[str, Optional[float]] = {}
        for tag, label in tag_to_label.items():
            leaf_values[label] = cluster_value.get(tag)
        mt = morans_I_test(tree, leaf_values, n_perm=min(n_perm, 1000),
                           seed=seed + 7)
        out["morans_I"] = mt["morans_I"]
        out["morans_p"] = mt["morans_p"]
        out["morans_z"] = mt.get("morans_z", np.nan)
        out["morans_lambda_max"] = mt.get("morans_lambda_max", np.nan)
        out["morans_max_perm"] = mt.get("morans_max_perm", np.nan)
        out["morans_I_norm"] = mt.get("morans_I_norm", np.nan)
        out["morans_n"] = mt["morans_n"]
        out["morans_note"] = mt["morans_note"]

    # ---- SECONDARY (discrete; thresholded labels) ----
    # OFF by default: Moran's I is the primary (and only) phylo test now.
    # The discrete Fitch/NN permutation tests are slow (1000 tree traversals
    # per stat per pair-method) and superseded. Code retained; enable with
    # --discrete-perms only if you want the legacy parsimony/NN/D columns.
    if do_discrete_perms and n_f1 > 0 and n_f2 > 0:
        pscore = _fitch_parsimony_binary(tree, leaf_states)
        out["parsimony_score"] = int(pscore)
        out["parsimony_p"] = round(
            permutation_p(tree, leaf_states, _fitch_parsimony_binary,
                          observed=pscore, n_perm=n_perm, seed=seed,
                          smaller_is_more_structured=True), 4
        )
        nn1, _ = nn_concordance(tree, leaf_states, k=1)
        out["nn_concordance"] = round(nn1, 4) if not np.isnan(nn1) else np.nan
        if not np.isnan(nn1):
            out["nn_concordance_p"] = round(
                permutation_p(tree, leaf_states,
                              lambda t, s: nn_concordance(t, s, k=1),
                              observed=nn1, n_perm=min(n_perm, 1000),
                              seed=seed + 1,
                              smaller_is_more_structured=False), 4
            )
        # k=3 only meaningful with at least 4 F1+F2 leaves on the coarse tree
        if (n_f1 + n_f2) >= 4:
            nn3, _ = nn_concordance(tree, leaf_states, k=3)
            out["nn3_concordance"] = round(nn3, 4) if not np.isnan(nn3) else np.nan
            if not np.isnan(nn3):
                out["nn3_concordance_p"] = round(
                    permutation_p(tree, leaf_states,
                                  lambda t, s: nn_concordance(t, s, k=3),
                                  observed=nn3,
                                  n_perm=min(n_perm, 1000),
                                  seed=seed + 2,
                                  smaller_is_more_structured=False), 4
                )

        if do_d_stat:
            try:
                d = d_statistic(tree, leaf_states, n_perm=200, seed=seed + 3)
                if d is not None:
                    out["D_statistic"] = round(d, 4)
            except Exception as e:
                notes.append(f"D_failed:{type(e).__name__}")
    out["notes"] = ";".join(notes)
    out["_leaf_states"] = leaf_states  # stash for figure pass; popped before CSV
    return out


# ---------------------------------------------------------------------------
# Coarse-tree figure (per pair, per method)
# ---------------------------------------------------------------------------

_STATE_COLORS = {"F1": "#d62728", "F2": "#1f77b4", "Amb": "#bbbbbb"}
_STATE_DISPLAY = {"F1": "F1c", "F2": "F2c", "Amb": "Amc"}


def render_phylo_placement(pair_id: str, method: str,
                           tree, leaf_states: Dict[str, str],
                           output_png: str,
                           parsimony_score: Optional[int] = None,
                           parsimony_p: Optional[float] = None,
                           nn_concordance_val: Optional[float] = None,
                           ) -> None:
    """Save a PNG of the COARSE tree, with leaves coloured by F1c/F2c/Amc
    and labelled by short cluster tag (S03, D, ...).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    leaves = list(tree.iter_leaves())
    n_leaves = len(leaves)
    if n_leaves == 0:
        return

    node_x: Dict[int, float] = {}
    node_y: Dict[int, float] = {}

    def _depth(n, d):
        node_x[id(n)] = d
        for c in n.children:
            cd = c.dist if (c.dist and c.dist > 0) else 0.1
            _depth(c, d + cd)
    _depth(tree, 0.0)

    leaf_idx = 0

    def _yset(n):
        nonlocal leaf_idx
        if n.is_leaf():
            node_y[id(n)] = leaf_idx
            leaf_idx += 1
        else:
            ys = []
            for c in n.children:
                _yset(c)
                ys.append(node_y[id(c)])
            node_y[id(n)] = float(np.mean(ys))
    _yset(tree)

    max_x = max(node_x.values()) if node_x else 1.0
    if max_x <= 0:
        max_x = 1.0

    fig_h = max(3.0, 0.45 * n_leaves + 1.0)
    fig_w = 7.0
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    for node in tree.traverse("preorder"):
        if node.is_root():
            continue
        x0, y0 = node_x[id(node.up)], node_y[id(node.up)]
        x1, y1 = node_x[id(node)], node_y[id(node)]
        ax.plot([x0, x0], [y0, y1], color="#444", linewidth=1.1, zorder=1)
        ax.plot([x0, x1], [y1, y1], color="#444", linewidth=1.1, zorder=1)

    for leaf in leaves:
        s = leaf_states.get(leaf.name, "Amb")
        col = _STATE_COLORS.get(s, _STATE_COLORS["Amb"])
        x, y = node_x[id(leaf)], node_y[id(leaf)]
        ax.scatter(x, y, c=col, s=140, zorder=3,
                   edgecolors="black", linewidths=0.5)
        ax.text(x + max_x * 0.03, y,
                f" {leaf.name} ({_STATE_DISPLAY.get(s, s)})",
                va="center", ha="left", fontsize=9)

    # Internal-node markers (small grey squares)
    for node in tree.traverse("preorder"):
        if node.is_leaf():
            continue
        ax.scatter(node_x[id(node)], node_y[id(node)],
                   c="#cccccc", s=28, zorder=2, marker="s",
                   edgecolors="black", linewidths=0.4)

    title = f"{pair_id} | {method}  (coarse cluster-level tree)"
    if parsimony_score is not None and not (isinstance(parsimony_score, float)
                                            and np.isnan(parsimony_score)):
        title += f"\nparsimony={int(parsimony_score)}"
    if parsimony_p is not None and not (isinstance(parsimony_p, float)
                                        and np.isnan(parsimony_p)):
        title += f"  p={parsimony_p:.3f}"
    if nn_concordance_val is not None and not (isinstance(nn_concordance_val, float)
                                               and np.isnan(nn_concordance_val)):
        title += f"  NN={nn_concordance_val:.2f}"
    ax.set_title(title, fontsize=10)

    legend_handles = [
        mpatches.Patch(color=_STATE_COLORS["F1"], label="F1c"),
        mpatches.Patch(color=_STATE_COLORS["F2"], label="F2c"),
        mpatches.Patch(color=_STATE_COLORS["Amb"], label="Amc"),
    ]
    ax.legend(handles=legend_handles, loc="upper right", fontsize=8,
              framealpha=0.85)
    ax.set_yticks([])
    ax.set_xlabel("Branch length (cumulative)", fontsize=8)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)

    plt.tight_layout()
    os.makedirs(os.path.dirname(output_png), exist_ok=True)
    plt.savefig(output_png, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# XCL1 (2n54B_2hdmA) Dishman cross-check at the cluster level
# ---------------------------------------------------------------------------

def xcl1_cross_check(pair_id: str,
                     tag_to_label: Dict[str, str],
                     cluster_pref_per_method: Dict[str, Dict[str, str]],
                     out_dir: str) -> Optional[dict]:
    """For 2n54B_2hdmA: print cluster-level F1c/F2c/Amc by method.

    The Dishman 2021 hypothesis places fold-switching gain on the XCL1
    lineage. At the coarse-tree level the per-cluster centred preferences
    are the right object to inspect against that hypothesis.
    """
    if pair_id != "2n54B_2hdmA":
        return None

    tags = sorted(tag_to_label.keys())
    print(f"\n[XCL1 cross-check] {pair_id} (coarse tree, {len(tags)} clusters)")
    header = f"  {'cluster':<18s} {'leaf':<6s}"
    for m in SCORED_METHODS:
        header += f" {m:>8s}"
    print(header)
    rows = []
    for tag in tags:
        label = tag_to_label[tag]
        prefs = []
        for m in SCORED_METHODS:
            p = cluster_pref_per_method.get(m, {}).get(tag, "Amb")
            prefs.append(_STATE_DISPLAY.get(p, p))
        line = f"  {tag:<18s} {label:<6s}"
        for p in prefs:
            line += f" {p:>8s}"
        print(line)
        rows.append({"cluster": tag, "leaf": label,
                     **{f"pref_{m}": prefs[i]
                        for i, m in enumerate(SCORED_METHODS)}})

    out_csv = os.path.join(out_dir, "phylo_placement_xcl1_clusters.csv")
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    print(f"  Wrote {out_csv}")
    return {"pair_id": pair_id, "n_clusters": len(tags)}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def discover_pairs() -> List[str]:
    if not os.path.isdir(DATA_DIR):
        return []
    return sorted(d for d in os.listdir(DATA_DIR)
                  if PAIR_DIR_RE.match(d)
                  and os.path.isdir(os.path.join(DATA_DIR, d)))


def _submit_self_sbatch(args, script_path: str) -> None:
    """Re-invoke this script under sbatch with --run_job_mode local.

    Mirrors the run_foldswitch_pipeline.py pattern: --run_job_mode sbatch
    submits a single sbatch job that runs the actual computation inline.
    """
    import shlex
    import subprocess
    jobs_dir = "Pipeline/FoldPairs/jobs"
    os.makedirs(jobs_dir, exist_ok=True)
    log = os.path.join(jobs_dir, f"phylo_placement_{args.labels}.out")
    # Build the inner command — same args minus --run_job_mode
    inner_args = []
    for k, v in vars(args).items():
        if k == "run_job_mode":
            continue
        if v is None or v is False:
            continue
        flag = "--" + k.replace("_", "-")
        if v is True:
            inner_args.append(flag)
        else:
            inner_args.extend([flag, str(v)])
    inner = " ".join([shlex.quote(sys.executable), script_path,
                      "--run-job-mode", "local"] + inner_args)
    sbatch_cmd = [
        "sbatch", "--parsable",
        f"--job-name=phylo_placement_{args.labels}",
        f"--cpus-per-task={args.sbatch_cpus}",
        f"--mem={args.sbatch_mem}",
        f"--time={args.sbatch_time}",
        f"--output={log}",
        "--wrap", shlex.quote(inner),
    ]
    if args.sbatch_partition:
        sbatch_cmd.insert(1, f"--partition={args.sbatch_partition}")
    cmd_str = " ".join(sbatch_cmd)
    print(f"[sbatch] {cmd_str}")
    jid = subprocess.check_output(cmd_str, shell=True, text=True).strip()
    print(f"[sbatch] submitted phylo_placement ({args.labels}) job {jid}")
    print(f"[sbatch] log: {log}")


def main():
    p = argparse.ArgumentParser(
        description="Phylogenetic-structure test for cluster fold-preference "
                    "on the COARSE cluster-level tree")
    p.add_argument("--pairs", default="ALL",
                   help="Comma-separated pair IDs, or ALL")
    p.add_argument("--labels", choices=["centered", "corrected"],
                   default="centered",
                   help="centered = use pref_centered from survey CSV (default); "
                        "corrected = use pref_corrected from "
                        "fold_diversity_survey_corrected.csv (residual-corrected, "
                        "PID-to-query stripped via seq_divergence_correction.py)")
    p.add_argument("--survey-csv", default=None,
                   help="Path to survey CSV (default depends on --labels)")
    p.add_argument("--output", default=None,
                   help="Output CSV path (default depends on --labels)")
    p.add_argument("--make-figures", default="0",
                   help="'all' = render every (pair, method) figure; "
                        "an integer N = top-N by parsimony_p per method; "
                        "0 = no figures (default).")
    p.add_argument("--figure-dir", default=None,
                   help="Where to write figures "
                        "(default: <docs>/figures/phylo_placement)")
    p.add_argument("--n-perm", type=int, default=1000)
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--delta-tm", type=float, default=0.05)
    p.add_argument("--delta-ddg", type=float, default=1.0)
    p.add_argument("--no-d-statistic", action="store_true",
                   help="Skip D statistic (slowest stat)")
    p.add_argument("--discrete-perms", action="store_true",
                   help="Also run the legacy discrete Fitch/NN permutation "
                        "tests (parsimony_p, nn_concordance_p, nn3, D). OFF by "
                        "default - Moran's I is the primary test and these are "
                        "slow. Discrete F1/F2/Amb label counts are always kept "
                        "(for figures); only their permutation p-values are gated.")
    # sbatch self-submission (matches run_foldswitch_pipeline.py pattern)
    p.add_argument("--run-job-mode", choices=["local", "sbatch"], default="local",
                   help="local = run inline (default); sbatch = self-submit "
                        "one sbatch job that runs the full analysis inline")
    p.add_argument("--sbatch-partition", default=None)
    p.add_argument("--sbatch-cpus", default=2)
    p.add_argument("--sbatch-mem", default="8G")
    p.add_argument("--sbatch-time", default="04:00:00")
    args = p.parse_args()

    # --- sbatch self-submission ---
    if args.run_job_mode == "sbatch":
        _submit_self_sbatch(args, os.path.abspath(__file__))
        return

    # --- mode-aware paths ---
    if args.labels == "corrected":
        default_survey = os.path.join(TABLES_RES, "fold_diversity_survey_corrected.csv")
        default_out = os.path.join(TABLES_RES, "phylo_placement_corrected.csv")
    else:
        default_survey = os.path.join(TABLES_RES, "fold_diversity_survey.csv")
        default_out = os.path.join(TABLES_RES, "phylo_placement.csv")
    survey_csv = args.survey_csv or default_survey
    out_csv = args.output or default_out
    figure_dir = args.figure_dir or os.path.join(
        TABLES_RES, "figures", "phylo_placement"
    )

    if not os.path.isfile(survey_csv):
        print(f"[error] survey CSV not found: {survey_csv}")
        sys.exit(2)

    if args.pairs.upper() == "ALL":
        pair_ids = discover_pairs()
    else:
        pair_ids = [s.strip() for s in args.pairs.split(",") if s.strip()]
    if not pair_ids:
        print("[error] no pairs to process")
        sys.exit(2)

    if args.make_figures.lower() == "all":
        fig_mode = "all"
        fig_topn = 0
    else:
        try:
            fig_topn = int(args.make_figures)
        except ValueError:
            print(f"[error] --make-figures must be 'all' or an integer, "
                  f"got {args.make_figures!r}")
            sys.exit(2)
        fig_mode = "topn" if fig_topn > 0 else "none"

    print(f"Loading fold preferences ({args.labels}): {survey_csv}")
    if args.labels == "corrected":
        pref_lookup = load_corrected_pref(survey_csv)
        value_col = "TMdiff_residual"
    else:
        pref_lookup = load_fold_pref(survey_csv,
                                     delta_tm=args.delta_tm,
                                     delta_ddg=args.delta_ddg)
        value_col = "TMdiff_centered"
    # Continuous per-cluster preference for the Moran's-I test (primary stat).
    value_lookup = load_fold_value(survey_csv, value_col)
    pairs_in_survey = set(k[0] for k in pref_lookup)
    print(f"  {len(pairs_in_survey)} pairs in survey")

    cached_states: Dict[Tuple[str, str], Tuple[object, Dict[str, str], dict]] = {}

    rows: List[dict] = []
    t_total = time.time()
    per_pair_timeout_s = int(os.environ.get("PHYLO_PAIR_TIMEOUT", "120"))
    for i, pair_id in enumerate(sorted(pair_ids), 1):
        t_pair = time.time()
        print(f"\n[{i}/{len(pair_ids)}] {pair_id}: building coarse tree "
              f"(timeout {per_pair_timeout_s}s)...", flush=True)
        try:
            with _per_pair_timeout(per_pair_timeout_s):
                tree, tag_to_label, _ = build_coarse_tree(pair_id)
                _name_internal_nodes(tree)
        except TimeoutError as e:
            print(f"  TIMEOUT after {per_pair_timeout_s}s on coarse-tree build; "
                  f"skipping {pair_id}", flush=True)
            continue
        except FileNotFoundError as e:
            print(f"  skipped: {e}", flush=True)
            continue
        except Exception as e:
            print(f"  coarse-tree build failed: {e}", flush=True)
            continue

        n_clusters = len(tag_to_label)
        if n_clusters < 2:
            print(f"  only {n_clusters} cluster(s) -- skip ({time.time()-t_pair:.1f}s)", flush=True)
            continue
        print(f"  coarse tree built with {n_clusters} cluster-leaves "
              f"({time.time()-t_pair:.1f}s)", flush=True)

        cluster_pref_per_method: Dict[str, Dict[str, str]] = {}
        for method in SCORED_METHODS:
            t_method = time.time()
            cp: Dict[str, str] = {
                tag: pref_lookup.get((pair_id, method, tag), "Amb")
                for tag in tag_to_label
            }
            cluster_pref_per_method[method] = cp
            cv: Dict[str, Optional[float]] = {
                tag: value_lookup.get((pair_id, method, tag))
                for tag in tag_to_label
            }

            try:
                with _per_pair_timeout(per_pair_timeout_s):
                    row = analyze_pair_method(
                        pair_id, method, tree, tag_to_label, cp,
                        n_perm=args.n_perm, seed=args.seed,
                        do_d_stat=(not args.no_d_statistic),
                        cluster_value=cv,
                        do_discrete_perms=args.discrete_perms,
                    )
            except TimeoutError:
                print(f"    {method}: TIMEOUT after {per_pair_timeout_s}s; skipped", flush=True)
                continue
            except Exception as e:
                print(f"    {method}: failed: {e}", flush=True)
                continue
            print(f"    {method}: {time.time()-t_method:.1f}s", flush=True)
            if row is None:
                continue
            row["n_clusters"] = n_clusters
            leaf_states = row.pop("_leaf_states")
            cached_states[(pair_id, method)] = (tree, leaf_states, row)
            rows.append(row)
        print(f"  pair total: {time.time()-t_pair:.1f}s  (running {time.time()-t_total:.0f}s)", flush=True)

        if pair_id == "2n54B_2hdmA":
            try:
                xcl1_cross_check(pair_id, tag_to_label,
                                 cluster_pref_per_method, TABLES_RES)
            except Exception as e:
                print(f"  [warn] XCL1 cross-check failed: {e}")

    df = pd.DataFrame(rows)
    if df.empty:
        print("[warn] no rows produced")
        return

    # BH FDR across the screen for each p-value column (always-on)
    for pcol in ("morans_p", "parsimony_p", "nn_concordance_p", "nn3_concordance_p"):
        if pcol in df.columns:
            # BH within each method, across pairs (the per-method screen).
            df[pcol + "_bh"] = np.nan
            for _m in df["method"].unique():
                mask = df["method"] == _m
                df.loc[mask, pcol + "_bh"] = bh_adjust(
                    df.loc[mask, pcol].to_numpy(dtype=float))

    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"\nWrote {out_csv} ({len(df)} rows; labels={args.labels})")

    # Per-method BH-significant counts on the PRIMARY stat (Moran's I).
    print("\n=== Moran's I phylogenetic signal: BH-significant pairs per method ===")
    for method in SCORED_METHODS:
        sub = df[(df.method == method) & df.morans_p.notna()].copy()
        if sub.empty:
            continue
        nsig = int((sub["morans_p_bh"] < 0.05).sum()) if "morans_p_bh" in sub else 0
        print(f"  {method:9} tested={len(sub):>3}  "
              f"raw_p<.05={int((sub.morans_p<0.05).sum()):>3}  BH<.05={nsig:>3}  "
              f"median_I={sub.morans_I.median():.3f}")

    for method in SCORED_METHODS:
        sub = df[(df.method == method) & df.morans_p.notna()].copy()
        if sub.empty:
            continue
        # p floors at 1/(n_perm+1); break ties by SIGNED I (clustering first)
        sub = sub.sort_values(["morans_p", "morans_I"], ascending=[True, False])
        print(f"\nTop 10 phylogenetically-structured pairs ({method}, by p then signed I):")
        cols = ["pair_id", "morans_n", "morans_I", "morans_p", "morans_p_bh",
                "n_F1c_leaves", "n_F2c_leaves", "nn_concordance_p"]
        cols = [c for c in cols if c in sub.columns]
        print(sub[cols].head(10).to_string(index=False))

    if fig_mode != "none":
        os.makedirs(figure_dir, exist_ok=True)
        if fig_mode == "all":
            print(f"\nRendering ALL coarse-tree figures -> {figure_dir}")
            targets = [(k, v) for k, v in cached_states.items()]
        else:
            print(f"\nRendering top {fig_topn} coarse-tree figures per method "
                  f"-> {figure_dir}")
            targets = []
            for method in SCORED_METHODS:
                sub = df[(df.method == method) & df.parsimony_p.notna()].copy()
                if sub.empty:
                    continue
                sub = sub.sort_values("parsimony_p")
                for _, r in sub.head(fig_topn).iterrows():
                    key = (r.pair_id, method)
                    if key in cached_states:
                        targets.append((key, cached_states[key]))
        for (pair_id, method), (tree, leaf_states, row) in targets:
            out_png = os.path.join(
                figure_dir,
                f"phylo_placement_{pair_id}_{method}.png"
            )
            try:
                render_phylo_placement(
                    pair_id, method, tree, leaf_states, out_png,
                    parsimony_score=row.get("parsimony_score"),
                    parsimony_p=row.get("parsimony_p"),
                    nn_concordance_val=row.get("nn_concordance"),
                )
            except Exception as e:
                print(f"  [warn] figure failed for {pair_id}/{method}: {e}")

    print("\nDone.")


if __name__ == "__main__":
    main()
