"""Phylogenetic tree tests for per-cluster fold preference (test registry).

Used by ``scripts/phylo_placement.py`` on the COARSE cluster tree (one leaf
per MSA cluster). Every test is registered in ``TESTS`` behind ONE common
signature::

    run(tree, leaf_values, n_perm, rng, **opts) -> dict

* ``tree``        ete3-like rooted tree (needs ``name``, ``dist``,
                  ``children``, ``up``, ``traverse``, ``is_leaf``,
                  ``is_root``, ``iter_leaves``). ete3 is NOT imported here.
* ``leaf_values`` {leaf_label: value}. Continuous tests take floats (None/NaN
                  = missing). Discrete tests take 'F1'/'F2'/'Amb' labels.
                  The pair-level test ``clade_shift_combined`` takes
                  {method: {leaf_label: float}}.
* ``n_perm``      number of permutations.
* ``rng``         ``np.random.Generator`` for numpy-based permutation nulls.
                  The legacy discrete tests (fitch / nn / dstat) use Python's
                  ``random.Random`` and read their seed from ``opts['seed']``
                  so their numbers stay identical to the pre-registry code.
* returns         a flat dict of prefixed output columns (``morans_*``,
                  ``parsimony_*``, ``nn_*``, ``D_statistic``, ``clade_*``).

Registered tests
----------------
  morans               continuous Moran's I (inverse patristic weights)
  fitch                Fitch parsimony on F1/F2 labels + permutation p
  nn                   k=1 and k=3 nearest-neighbour concordance + perm p
  dstat                Fritz & Purvis D statistic (no p-value)
  clade_shift          max-|z| Mann-Whitney scan over clades (two-sided)
  clade_shift_combined same scan, per-method z's combined (Stouffer) at each
                       split; one joint leaf permutation across methods
"""

from __future__ import annotations

import random
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np


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


def morans_weights(tree, names: Sequence[str], eps: float = 1e-9) -> np.ndarray:
    """Row-standardized inverse-patristic Moran weights over ``names``
    (zero diagonal; S0 = n). Shared by ``morans_I_test`` and the fast paths."""
    W = 1.0 / (_patristic_matrix(tree, names) + eps)
    np.fill_diagonal(W, 0.0)
    rs = W.sum(axis=1, keepdims=True)
    rs[rs == 0] = 1.0
    return W / rs                                # row-standardized -> S0 = n


def morans_spectral(W: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(M, evals, evecs): symmetric form M = (W + W^T)/2 and the eigen-
    decomposition of its doubly-centered version (ascending eigenvalues).
    Depends only on the tree / leaf set, so it can be computed once and
    reused across many traits (``morans_I_norm``)."""
    nn = W.shape[0]
    M = 0.5 * (W + W.T)                                  # symmetric quadratic form
    Pc = np.eye(nn) - np.ones((nn, nn)) / nn             # centering projector
    evals, evecs = np.linalg.eigh(Pc @ M @ Pc)
    return M, evals, evecs


def morans_I_norm(xc: np.ndarray, den: float, I_obs: float, M: np.ndarray,
                  evecs: np.ndarray) -> Tuple[float, float]:
    """(max_perm_I, I_norm) for centered values ``xc`` (den = xc @ xc).

    Moran's I is bounded by the tree's weight matrix, NOT by 1. I_norm is I
    as a fraction of the max achievable for THIS value multiset on THIS tree:
      max_perm = data-specific achievable max (sort-and-match on the top
                 eigenvector + capped vectorized 2-opt)
      I_norm   = (I_obs - E[I]) / (max_perm - E[I]),  E[I] = -1/(n-1)
    (negative I_obs is normalised against the reachable minimum instead).
    """
    nn = len(xc)

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
    Emean = -1.0 / (nn - 1)
    seed_top = np.empty(nn); seed_top[np.argsort(evecs[:, -1])] = np.sort(xc)
    max_perm_I = max(_two_opt(seed_top, True),  _two_opt(xc, True))
    # Normalize by the achievable range on the side I_obs falls -- positive
    # signal against the reachable max, negative against the reachable min --
    # so I_norm is a true fraction of achievable clade signal in [-1, 1].
    if I_obs >= Emean:
        denom = max_perm_I - Emean
    else:
        seed_bot = np.empty(nn); seed_bot[np.argsort(evecs[:, 0])] = np.sort(xc)
        min_perm_I = min(_two_opt(seed_bot, False), _two_opt(xc, False))
        denom = Emean - min_perm_I
    I_norm = (I_obs - Emean) / denom if denom > 1e-9 else float("nan")
    if np.isfinite(I_norm):
        I_norm = float(np.clip(I_norm, -1.0, 1.0))   # belt-and-suspenders
    return max_perm_I, I_norm


def morans_I_test(tree, leaf_values: Dict[str, Optional[float]],
                  n_perm: int = 1000, seed: int = 12345,
                  eps: float = 1e-9,
                  rng: Optional[np.random.Generator] = None) -> dict:
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
    - ``rng``: permutation generator; default ``np.random.default_rng(seed)``.
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

    W = morans_weights(tree, names, eps)

    I_obs = float((xc @ (W @ xc)) / den)
    if rng is None:
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
    #   lambda_max  = eigenvalue ceiling (exact upper bound over centered vectors)
    #   max_perm / I_norm: see ``morans_I_norm``
    try:
        M, evals, evecs = morans_spectral(W)
        max_perm_I, I_norm = morans_I_norm(xc, den, I_obs, M, evecs)
        res.update({
            "morans_lambda_max": round(float(evals[-1]), 4),
            "morans_max_perm": round(max_perm_I, 4),
            "morans_I_norm": round(float(I_norm), 4) if np.isfinite(I_norm) else np.nan,
        })
    except Exception:
        pass
    return res


def morans_analytic_z(W: np.ndarray, xc: np.ndarray, I_obs: float) -> float:
    """Permutation-free z of Moran's I: (I - E[I]) / sd[I] with the EXACT
    mean and variance of I under random permutation of the values over the
    leaves (Cliff & Ord 1981 "randomization" moments; valid for asymmetric W).
    It estimates the same quantity as ``morans_I_test``'s ``morans_z`` (which
    uses the mean / SD of a finite permutation sample) without permuting.
    Requires n >= 4; returns NaN otherwise."""
    n = float(len(xc))
    if n < 4:
        return float("nan")
    S0 = float(W.sum())
    Ws = W + W.T
    S1 = 0.5 * float((Ws ** 2).sum())
    S2 = float(((W.sum(axis=1) + W.sum(axis=0)) ** 2).sum())
    m2 = float((xc ** 2).sum())
    b2 = n * float((xc ** 4).sum()) / (m2 * m2)
    EI = -1.0 / (n - 1.0)
    EI2 = ((n * ((n * n - 3 * n + 3) * S1 - n * S2 + 3 * S0 * S0)
            - b2 * ((n * n - n) * S1 - 2 * n * S2 + 6 * S0 * S0))
           / ((n - 1) * (n - 2) * (n - 3) * S0 * S0))
    var = EI2 - EI * EI
    return float((I_obs - EI) / np.sqrt(var)) if var > 1e-15 else float("nan")


def morans_fast_stats(W: np.ndarray, M: np.ndarray, evecs: np.ndarray,
                      x: np.ndarray, eps: float = 1e-9,
                      with_norm: bool = True) -> Tuple[float, float, float]:
    """(I, I_norm, z_analytic) for values ``x`` (all finite, aligned to W's
    leaf order) with the tree quantities precomputed by ``morans_weights`` /
    ``morans_spectral``. I and I_norm are numerically identical to
    ``morans_I_test`` (before its 4-decimal rounding); z is the analytic
    randomization z (``morans_analytic_z``). Zero variance -> all NaN."""
    xc = np.asarray(x, dtype=float) - float(np.mean(x))
    den = float(xc @ xc)
    if len(xc) < 3 or den < eps:
        return float("nan"), float("nan"), float("nan")
    I_obs = float((xc @ (W @ xc)) / den)
    I_norm = morans_I_norm(xc, den, I_obs, M, evecs)[1] if with_norm else float("nan")
    return I_obs, I_norm, morans_analytic_z(W, xc, I_obs)




# ---------------------------------------------------------------------------
# Clade-shift scan: max-|z| Mann-Whitney over all clades (two-sided)
# ---------------------------------------------------------------------------
#
# Moran's I measures DIFFUSE autocorrelation over the whole tree and is ~0
# when a single small clade (e.g. 4 of 100 clusters) gains / loses the
# alternative fold. The clade-shift scan instead asks: is there ANY clade
# whose values are shifted relative to the rest of the tree?
#
#   * Splits: every non-root node of the tree defines clade-vs-rest; keep
#     splits with min_clade <= |clade| and |rest| >= min_clade (with the
#     default min_clade=2 this is exactly the internal non-root nodes). The
#     root (clade = everything) is skipped; identical and complementary
#     splits (e.g. the root's two children) are de-duplicated, keeping the
#     smaller clade as the reported one.
#   * Per split: standardized Mann-Whitney U of clade vs rest, average ranks
#     for ties and tie-corrected variance, no continuity correction; z > 0
#     means the clade is HIGHER.
#   * T = max |z| over splits; p from permuting the values over the leaves
#     of the same tree: p = (1 + #{T_perm >= T_obs}) / (1 + n_perm).
#   * Vectorized: a boolean (n_splits x n_leaves) membership matrix is built
#     once; ranks are computed once (a permutation of the values permutes
#     the rank vector and leaves the tie term unchanged), so rank sums for
#     all splits and a whole chunk of permutations are one matrix product.
#
# Missing values: a single-method scan uses only the leaves with a value
# (ranks within those leaves; clade sizes counted over them). The combined
# scan permutes ALL leaves jointly, so each method's missing-value mask moves
# with its values; per method, z at a split is computed over that method's
# available leaves only (rank sums of available leaves, clade size = number
# of available leaves in the clade) and the split is skipped for that method
# when either side has < min_clade available leaves.

_CLADE_COLS = ("clade_T", "clade_p", "clade_n_leaves", "clade_n_splits",
               "clade_best_size", "clade_best_members",
               "clade_best_direction", "clade_best_effect",
               "clade_best_rbc", "clade_note")
_PERM_CHUNK = 2000


def _empty_clade_result(n_leaves: int, note: str) -> dict:
    res = {c: np.nan for c in _CLADE_COLS}
    res.update({"clade_n_leaves": n_leaves, "clade_n_splits": 0,
                "clade_best_members": "", "clade_note": note})
    return res


def clade_split_matrix(tree, leaf_names: Sequence[str], min_clade: int = 2
                       ) -> Tuple[np.ndarray, List[str]]:
    """Boolean clade-membership matrix over ``leaf_names``.

    Returns ``(M, node_names)`` with ``M[s, j]`` True iff leaf j is in the
    clade of split s. Leaves of ``tree`` not in ``leaf_names`` are ignored
    (clade sizes count only the listed leaves). See the block comment above
    for which splits are kept.
    """
    idx = {nm: j for j, nm in enumerate(leaf_names)}
    n = len(leaf_names)
    members: Dict[int, np.ndarray] = {}
    rows: List[np.ndarray] = []
    names: List[str] = []
    seen: Dict[bytes, int] = {}          # canonical split key -> row index
    for node in tree.traverse("postorder"):
        v = np.zeros(n, dtype=bool)
        if node.is_leaf():
            j = idx.get(node.name)
            if j is not None:
                v[j] = True
        else:
            for c in node.children:
                v |= members[id(c)]
        members[id(node)] = v
        if node.is_root():
            continue
        k = int(v.sum())
        if k < min_clade or n - k < min_clade:
            continue
        canon = v if not v[0] else ~v
        key = canon.tobytes()
        if key in seen:
            s = seen[key]
            if k < int(rows[s].sum()):   # keep the smaller side as "the clade"
                rows[s], names[s] = v, node.name
            continue
        seen[key] = len(rows)
        rows.append(v)
        names.append(node.name)
    if not rows:
        return np.zeros((0, n), dtype=bool), []
    return np.vstack(rows), names


def average_ranks(x: np.ndarray) -> Tuple[np.ndarray, float]:
    """Average (mid) ranks 1..n of ``x`` and the tie term sum(t^3 - t)."""
    _, inv, counts = np.unique(x, return_inverse=True, return_counts=True)
    avg = np.cumsum(counts) - (counts - 1) / 2.0
    return avg[inv].astype(float), float((counts.astype(float) ** 3 - counts).sum())


def mann_whitney_z(R1: np.ndarray, n1: np.ndarray, n: float, tie_term: float
                   ) -> np.ndarray:
    """Standardized Mann-Whitney U (clade vs rest) from clade rank sums.

    U1 = R1 - n1(n1+1)/2,  E[U1] = n1 n2 / 2,
    Var[U1] = n1 n2 / 12 * ((n + 1) - tie_term / (n (n - 1))).
    Returns NaN where the variance is 0. Broadcasts over R1 / n1.
    """
    n1 = np.asarray(n1, dtype=float)
    n2 = n - n1
    U1 = R1 - n1 * (n1 + 1) / 2.0
    var = n1 * n2 / 12.0 * ((n + 1) - tie_term / (n * (n - 1)))
    with np.errstate(invalid="ignore", divide="ignore"):
        z = (U1 - n1 * n2 / 2.0) / np.sqrt(var)
    return np.where(var > 1e-12, z, np.nan)


def _scan_core(M: np.ndarray, ranks: np.ndarray, avail: np.ndarray,
               ties: np.ndarray, min_clade: int, perm: Optional[np.ndarray]
               ) -> Tuple[np.ndarray, np.ndarray]:
    """Combined split z's for a batch of leaf permutations.

    M      (S, n) bool membership; ranks / avail (K, n) per method (ranks 0
    where missing); ties (K,); perm (B, n) leaf-index permutations, or None
    for the observed (identity) arrangement.
    Returns (zc, cnt): (B, S) combined z = sum_k z_k / sqrt(cnt) and (B, S)
    number of methods contributing at each split.
    """
    Mt = M.T.astype(float)                               # (n, S)
    B = 1 if perm is None else perm.shape[0]
    zsum = np.zeros((B, M.shape[0]))
    cnt = np.zeros((B, M.shape[0]))
    for r, a, tt in zip(ranks, avail, ties):
        n_k = float(a.sum())
        R = r[None, :] if perm is None else r[perm]      # (B, n)
        A = a[None, :] if perm is None else a[perm]
        R1 = R @ Mt                                      # (B, S) rank sums
        N1 = A.astype(float) @ Mt                        # (B, S) clade sizes
        valid = (N1 >= min_clade) & (n_k - N1 >= min_clade)
        z = mann_whitney_z(R1, N1, n_k, tt)
        valid &= np.isfinite(z)
        zsum += np.where(valid, z, 0.0)
        cnt += valid
    with np.errstate(invalid="ignore", divide="ignore"):
        zc = np.where(cnt > 0, zsum / np.sqrt(cnt), np.nan)
    return zc, cnt


def _max_abs(zc: np.ndarray) -> np.ndarray:
    a = np.abs(zc)
    a[~np.isfinite(a)] = -np.inf
    t = a.max(axis=1)
    return np.where(np.isfinite(t), t, np.nan)


def _scan_with_null(M, ranks, avail, ties, min_clade, n_perm, rng):
    """Observed combined-z vector and permutation p of T = max |z|."""
    zc_obs, cnt_obs = _scan_core(M, ranks, avail, ties, min_clade, None)
    zc_obs, cnt_obs = zc_obs[0], cnt_obs[0]
    T_obs = float(_max_abs(zc_obs[None, :])[0])
    if not np.isfinite(T_obs):
        return zc_obs, cnt_obs, np.nan, np.nan
    n = M.shape[1]
    base = np.arange(n)
    n_ge = 0
    done = 0
    while done < n_perm:
        b = min(_PERM_CHUNK, n_perm - done)
        perm = rng.permuted(np.tile(base, (b, 1)), axis=1)
        T_perm = _max_abs(_scan_core(M, ranks, avail, ties, min_clade, perm)[0])
        n_ge += int((T_perm >= T_obs - 1e-9).sum())     # NaN compares False
        done += b
    return zc_obs, cnt_obs, T_obs, (1 + n_ge) / (1 + n_perm)


def clade_T_stat(M: np.ndarray, x: np.ndarray, min_clade: int = 2) -> float:
    """Observed clade-shift statistic T = max |z| over the splits of ``M``
    (as ``run_clade_shift``'s ``clade_T`` for leaves with no missing values),
    with no permutation null. ``x`` is aligned to M's columns. NaN when the
    trait is constant or no split is valid."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    if M.shape[0] == 0 or n < 2 * min_clade:
        return float("nan")
    r, tt = average_ranks(x)
    if tt >= n ** 3 - n - 1e-9:
        return float("nan")
    zc, _ = _scan_core(M, r[None, :], np.ones((1, n), bool), np.array([tt]),
                       min_clade, None)
    return float(_max_abs(zc)[0])


def _members_str(M_row: np.ndarray, leaf_names: Sequence[str],
                 label_to_tag: Optional[Dict[str, str]]) -> str:
    tags = [(label_to_tag or {}).get(nm, nm)
            for nm, inside in zip(leaf_names, M_row) if inside]
    return ";".join(sorted(tags))


def _effects(x: np.ndarray, inside: np.ndarray) -> Tuple[float, float]:
    """(mean clade - mean rest, rank-biserial 2 U1 / (n1 n2) - 1)."""
    r, _ = average_ranks(x)
    n1, n2 = int(inside.sum()), int((~inside).sum())
    diff = float(x[inside].mean() - x[~inside].mean())
    U1 = r[inside].sum() - n1 * (n1 + 1) / 2.0
    return diff, float(2.0 * U1 / (n1 * n2) - 1.0)


def run_clade_shift(tree, leaf_values: Dict[str, Optional[float]],
                    n_perm: int, rng: np.random.Generator,
                    min_clade: int = 2,
                    label_to_tag: Optional[Dict[str, str]] = None,
                    **_opts) -> dict:
    """Single-method clade-shift scan (two-sided). Leaves without a finite
    value are dropped. ``label_to_tag`` maps leaf labels ('S05') back to
    cluster tags ('ShallowMsa_005') for ``clade_best_members``."""
    names = [nm for nm, v in leaf_values.items()
             if v is not None and np.isfinite(v)]
    n = len(names)
    if n < 2 * min_clade:
        return _empty_clade_result(n, f"n<{2 * min_clade}")
    x = np.asarray([leaf_values[nm] for nm in names], dtype=float)
    r, tt = average_ranks(x)
    if tt >= n ** 3 - n - 1e-9:
        return _empty_clade_result(n, "zero_variance")
    M, _ = clade_split_matrix(tree, names, min_clade)
    if M.shape[0] == 0:
        return _empty_clade_result(n, "no_valid_splits")
    zc, _, T, p = _scan_with_null(M, r[None, :], np.ones((1, n), bool),
                                  np.array([tt]), min_clade, n_perm, rng)
    if not np.isfinite(T):
        return _empty_clade_result(n, "no_valid_splits")
    s = int(np.nanargmax(np.abs(zc)))
    diff, rbc = _effects(x, M[s])
    res = _empty_clade_result(n, "")
    res.update({
        "clade_T": round(T, 4), "clade_p": round(float(p), 6),
        "clade_n_splits": int(M.shape[0]),
        "clade_best_size": int(M[s].sum()),
        "clade_best_members": _members_str(M[s], names, label_to_tag),
        "clade_best_direction": int(np.sign(zc[s])),
        "clade_best_effect": round(diff, 4),
        "clade_best_rbc": round(rbc, 4),
    })
    return res


def run_clade_shift_combined(tree,
                             leaf_values: Dict[str, Dict[str, Optional[float]]],
                             n_perm: int, rng: np.random.Generator,
                             min_clade: int = 2,
                             label_to_tag: Optional[Dict[str, str]] = None,
                             **_opts) -> dict:
    """Multi-method clade-shift scan (one row per pair, method='COMBINED').

    ``leaf_values`` is {method: {leaf_label: value}}. At each split the
    per-method z's (signs aligned: every method's value is a signed
    F1-vs-F2 preference, so z > 0 = clade leans to F1 in every method) are
    combined as sum_k z_k / sqrt(#methods present at that split); methods
    with zero variance or < 2*min_clade values are skipped. T = max |z|.
    The null applies ONE permutation of the leaves to all methods' value
    vectors (and missing-value masks), which keeps the cross-method
    agreement within each cluster while breaking any clade association.

    Extra columns: clade_n_methods, clade_methods, clade_best_method_z
    ('AF2:+2.31;ESM:+1.90' at the best split). clade_best_effect is the
    mean over contributing methods of (mean clade - mean rest) / SD(method),
    and clade_best_rbc the mean rank-biserial.
    """
    leaf_names = [l.name for l in tree.iter_leaves()]
    idx = {nm: j for j, nm in enumerate(leaf_names)}
    n = len(leaf_names)
    used, skipped = [], []
    ranks, avail, ties, xs = [], [], [], []
    for meth, vals in leaf_values.items():
        a = np.zeros(n, dtype=bool)
        x = np.full(n, np.nan)
        for nm, v in vals.items():
            j = idx.get(nm)
            if j is not None and v is not None and np.isfinite(v):
                a[j], x[j] = True, float(v)
        k = int(a.sum())
        if k == 0:                      # method not scored for this pair
            continue
        if k < 2 * min_clade:
            skipped.append(f"{meth}(n<{2 * min_clade})")
            continue
        rk, tt = average_ranks(x[a])
        if tt >= k ** 3 - k - 1e-9:
            skipped.append(f"{meth}(zero_variance)")
            continue
        r = np.zeros(n)
        r[a] = rk
        used.append(meth)
        ranks.append(r)
        avail.append(a)
        ties.append(tt)
        xs.append(x)
    note = "skipped=" + ",".join(skipped) if skipped else ""
    extra = {"clade_n_methods": len(used), "clade_methods": ",".join(used),
             "clade_best_method_z": ""}
    if not used:
        res = _empty_clade_result(n, (note + ";" if note else "") + "no_methods")
        res.update(extra)
        return res
    M, _ = clade_split_matrix(tree, leaf_names, min_clade)
    if M.shape[0] == 0:
        res = _empty_clade_result(n, (note + ";" if note else "") + "no_valid_splits")
        res.update(extra)
        return res
    zc, cnt, T, p = _scan_with_null(M, np.vstack(ranks), np.vstack(avail),
                                    np.asarray(ties), min_clade, n_perm, rng)
    res = _empty_clade_result(n, note)
    res.update(extra)
    if not np.isfinite(T):
        res["clade_note"] = (note + ";" if note else "") + "no_valid_splits"
        return res
    s = int(np.nanargmax(np.abs(zc)))
    inside = M[s]
    effs, rbcs, mz = [], [], []
    for meth, r, a, tt, x in zip(used, ranks, avail, ties, xs):
        n1 = float((a & inside).sum())
        nk = float(a.sum())
        if n1 < min_clade or nk - n1 < min_clade:
            continue
        z = float(mann_whitney_z(np.array([r[inside].sum()]), np.array([n1]),
                                 nk, tt)[0])
        if not np.isfinite(z):
            continue
        diff, rbc = _effects(x[a], inside[a])
        sd = float(np.std(x[a]))
        effs.append(diff / sd if sd > 0 else np.nan)
        rbcs.append(rbc)
        mz.append(f"{meth}:{z:+.2f}")
    res.update({
        "clade_T": round(T, 4), "clade_p": round(float(p), 6),
        "clade_n_splits": int(M.shape[0]),
        "clade_best_size": int(inside.sum()),
        "clade_best_members": _members_str(inside, leaf_names, label_to_tag),
        "clade_best_direction": int(np.sign(zc[s])),
        "clade_best_effect": round(float(np.nanmean(effs)), 4) if effs else np.nan,
        "clade_best_rbc": round(float(np.mean(rbcs)), 4) if rbcs else np.nan,
        "clade_best_method_z": ";".join(mz),
    })
    return res


# ---------------------------------------------------------------------------
# Registry wrappers for the pre-existing tests (numbers unchanged)
# ---------------------------------------------------------------------------

def _has_f1_f2(leaf_states: Dict[str, str]) -> bool:
    vals = set(leaf_states.values())
    return "F1" in vals and "F2" in vals


def run_morans(tree, leaf_values, n_perm, rng, **_opts) -> dict:
    """Moran's I (see ``morans_I_test``); ``rng`` drives the permutations."""
    return morans_I_test(tree, leaf_values, n_perm=n_perm, rng=rng)


def run_fitch(tree, leaf_values, n_perm, rng, seed: int = 12345,
              **_opts) -> dict:
    """Fitch parsimony + label-shuffle p (python ``random.Random(seed)``)."""
    out = {"parsimony_score": np.nan, "parsimony_p": np.nan}
    if not _has_f1_f2(leaf_values):
        return out
    pscore = _fitch_parsimony_binary(tree, leaf_values)
    out["parsimony_score"] = int(pscore)
    out["parsimony_p"] = round(
        permutation_p(tree, leaf_values, _fitch_parsimony_binary,
                      observed=pscore, n_perm=n_perm, seed=seed,
                      smaller_is_more_structured=True), 4
    )
    return out


def run_nn(tree, leaf_values, n_perm, rng, seed: int = 12345,
           **_opts) -> dict:
    """k=1 (seed) and k=3 (seed+1) nearest-neighbour concordance + perm p."""
    out = {"nn_concordance": np.nan, "nn_concordance_p": np.nan,
           "nn3_concordance": np.nan, "nn3_concordance_p": np.nan}
    if not _has_f1_f2(leaf_values):
        return out
    nn1, _ = nn_concordance(tree, leaf_values, k=1)
    out["nn_concordance"] = round(nn1, 4) if not np.isnan(nn1) else np.nan
    if not np.isnan(nn1):
        out["nn_concordance_p"] = round(
            permutation_p(tree, leaf_values,
                          lambda t, s: nn_concordance(t, s, k=1),
                          observed=nn1, n_perm=n_perm, seed=seed,
                          smaller_is_more_structured=False), 4
        )
    # k=3 only meaningful with at least 4 F1+F2 leaves on the coarse tree
    if sum(1 for s in leaf_values.values() if s in ("F1", "F2")) >= 4:
        nn3, _ = nn_concordance(tree, leaf_values, k=3)
        out["nn3_concordance"] = round(nn3, 4) if not np.isnan(nn3) else np.nan
        if not np.isnan(nn3):
            out["nn3_concordance_p"] = round(
                permutation_p(tree, leaf_values,
                              lambda t, s: nn_concordance(t, s, k=3),
                              observed=nn3, n_perm=n_perm, seed=seed + 1,
                              smaller_is_more_structured=False), 4
            )
    return out


def run_dstat(tree, leaf_values, n_perm, rng, seed: int = 12345,
              **_opts) -> dict:
    """Fritz & Purvis D (Monte-Carlo with ``n_perm`` draws per null)."""
    out = {"D_statistic": np.nan}
    if not _has_f1_f2(leaf_values):
        return out
    d = d_statistic(tree, leaf_values, n_perm=n_perm, seed=seed)
    if d is not None:
        out["D_statistic"] = round(d, 4)
    return out


def run_morans_combined(tree, leaf_values, n_perm, rng, eps: float = 1e-9,
                        **_opts) -> dict:
    """Multi-method Moran's I, pooled like ``run_clade_shift_combined``
    (one row per pair, method='COMBINED').

    ``leaf_values`` is {method: {leaf_label: value}}. For each method, z_m is
    the Cliff--Ord randomization z of Moran's I over the leaves that method
    scores (inverse-patristic weights, row-standardized on that subset); the
    pooled statistic is Z = sum_m z_m / sqrt(#methods). The null applies ONE
    permutation of the leaves to all methods' value vectors and missing-value
    masks, which keeps the cross-method agreement within each cluster while
    breaking any tree association. One-sided (positive autocorrelation):
    p = (1 + #{Z_perm >= Z_obs}) / (1 + n_perm). Methods with fewer than four
    values or zero variance are skipped.
    """
    names = [l.name for l in tree.iter_leaves()]
    idx = {nm: j for j, nm in enumerate(names)}
    n = len(names)
    Winv = 1.0 / (_patristic_matrix(tree, names) + eps)
    np.fill_diagonal(Winv, 0.0)
    used, skipped, X, A = [], [], [], []
    for meth, vals in leaf_values.items():
        a = np.zeros(n, dtype=bool)
        x = np.full(n, np.nan)
        for nm, v in vals.items():
            j = idx.get(nm)
            if j is not None and v is not None and np.isfinite(v):
                a[j], x[j] = True, float(v)
        k = int(a.sum())
        if k == 0:
            continue
        if k < 4 or float(np.nanvar(x[a])) < eps:
            skipped.append(f"{meth}(n<4_or_zero_variance)")
            continue
        used.append(meth); X.append(x); A.append(a)
    out = {"morans_I": np.nan, "morans_z": np.nan, "morans_p": np.nan,
           "morans_n": n, "morans_note": "", "morans_n_methods": len(used),
           "morans_methods": ";".join(used)}
    if not used:
        out["morans_note"] = "no_method"
        return out

    def _w(a):
        sub = np.where(a)[0]
        W = Winv[np.ix_(sub, sub)]
        rs = W.sum(axis=1, keepdims=True)
        rs[rs == 0] = 1.0
        return sub, W / rs

    # Per method: full-coverage methods keep one W (z's variance depends only
    # on W and on the permutation-invariant moments of x), so their permuted
    # z's are vectorised; partial-coverage methods rebuild W per permutation.
    perms = np.array([rng.permutation(n) for _ in range(n_perm)])
    z_obs, z_perm = [], np.zeros(n_perm)
    for x, a in zip(X, A):
        sub, W = _w(a)
        xc = x[sub] - x[sub].mean()
        den = float(xc @ xc)
        I_obs = float(xc @ (W @ xc)) / den
        z_obs.append(morans_analytic_z(W, xc, I_obs))
        if a.all():
            EI = -1.0 / (n - 1.0)
            sd = 1.0 / morans_analytic_z(W, xc, EI + 1.0)     # null SD of I on this W
            xcf = x - x.mean()
            P = xcf[perms]                                  # (B, n)
            I_p = np.einsum("ij,ij->i", P @ W.T, P) / den
            z_perm += (I_p - EI) / sd
        else:
            for b, pm in enumerate(perms):
                xp, ap = x[pm], a[pm]
                sub_p, Wp = _w(ap)
                xcp = xp[sub_p] - xp[sub_p].mean()
                dp = float(xcp @ xcp)
                Ip = float(xcp @ (Wp @ xcp)) / dp
                z_perm[b] += morans_analytic_z(Wp, xcp, Ip)
    k = len(used)
    Z = float(np.sum(z_obs)) / np.sqrt(k)
    Zp = z_perm / np.sqrt(k)
    p = (1 + int((Zp >= Z - 1e-12).sum())) / (1 + n_perm)
    out.update({"morans_z": round(Z, 4), "morans_p": round(float(p), 4)})
    if skipped:
        out["morans_note"] = "skipped:" + ",".join(skipped)
    return out


# name -> (run function, input kind). kind: 'continuous' = per-method float
# values; 'discrete' = per-method F1/F2/Amb labels; 'pair' = all methods'
# float values at once (one row per pair, method='COMBINED').
TESTS: Dict[str, Tuple[Callable[..., dict], str]] = {
    "morans": (run_morans, "continuous"),
    "fitch": (run_fitch, "discrete"),
    "nn": (run_nn, "discrete"),
    "dstat": (run_dstat, "discrete"),
    "clade_shift": (run_clade_shift, "continuous"),
    "clade_shift_combined": (run_clade_shift_combined, "pair"),
    "morans_combined": (run_morans_combined, "pair"),
}

# Per-test p-value columns that get a per-method BH adjustment (<col>_bh).
P_COLUMNS: Dict[str, Tuple[str, ...]] = {
    "morans": ("morans_p",),
    "fitch": ("parsimony_p",),
    "nn": ("nn_concordance_p", "nn3_concordance_p"),
    "dstat": (),
    "clade_shift": ("clade_p",),
    "clade_shift_combined": ("clade_p",),
    "morans_combined": ("morans_p",),
}


def run_test(name: str, tree, leaf_values, n_perm: int,
             rng: Optional[np.random.Generator] = None, **opts) -> dict:
    """Dispatch to a registered test by name."""
    if name not in TESTS:
        raise KeyError(f"unknown phylo test {name!r}; known: {sorted(TESTS)}")
    fn, _ = TESTS[name]
    if rng is None:
        rng = np.random.default_rng(opts.get("seed", 12345))
    return fn(tree, leaf_values, n_perm, rng, **opts)
