#!/usr/bin/env python3
"""Power simulation for the phylogenetic fold-preference tests (referee #5).

Question: on the REAL coarse cluster trees of the fold-switch pairs, how often
do Moran's I and the clade-shift scan detect a planted clade-structured
fold-preference signal, as a function of effect size, clade size, and how
noisy the per-cluster measurement is?

Trees
-----
For each pair, ``utils.ancestral_utils.build_coarse_tree`` (one leaf per MSA
cluster) followed by ``_name_internal_nodes`` -- exactly the tree that
``scripts/phylo_placement.py`` tests on. Default: all pairs with
>= ``--min-clusters`` (20) leaves, sub-sampled to ``--n-trees`` (20) trees
spread evenly over the size bins (cluster counts are read from
``docs/phylo_placement.csv`` when present so that only the selected trees are
built). ``--pairs A,B,..`` overrides the selection.

Simulated values (per replicate, per leaf c)
--------------------------------------------
true_c   = signal_c + e_c                       e_c ~ N(0, 1) iid ("biological"
                                                 scatter; the noise SD unit)
observed = true_c + m_c                          m_c ~ N(0, s_m^2)
           with s_m^2 = V_true (1 - r) / r so that the reliability
           r = Var(true) / Var(observed).  V_true is the EXPECTED variance of
           the true value (see below), fixed per scenario and tree.

Scenarios
---------
(a) clade : signal_c = delta * 1[c in clade]; the clade is the internal node
            (a split of the coarse tree, as enumerated by the clade scan)
            whose leaf fraction is closest to f.
            f in {0.05, 0.10, 0.25, 0.50}, delta in {0.5, 1, 1.5, 2, 3}
            (delta in SD units of e).  V_true = 1 + delta^2 p (1 - p),
            p = realised clade fraction.
(b) bm    : true_c = sqrt(h) * B_c + sqrt(1 - h) * e_c, where B is Brownian
            motion on the coarse tree, rescaled so the mean (over leaves)
            marginal variance is 1.  h in {0.1, 0.25, 0.5} is the BM variance
            fraction (a phylogenetic-heritability-like signal).  V_true = 1.
            Branch lengths follow the same convention as the Moran weights
            (``_patristic_matrix``): missing/zero length counts as 1.0.
(c) null  : true_c = e_c (r irrelevant: noise on pure noise is still noise).
Reliability axis r in {1, 0.4, 0.1} is crossed with every (a)/(b) scenario.
Motivation: the observed cross-method per-cluster fold-preference
correlations have median ~0.04 and max ~0.39; if two methods measured the
same true value with independent errors, corr = r, so realistic per-method
reliabilities are r ~ 0.04-0.4.

Tests (same implementations as the paper)
-----------------------------------------
* Moran's I: ``utils.phylo_tests.morans_I_test`` -- row-standardised inverse
  patristic weights, one-sided permutation p.  For speed the weight matrix is
  built once per tree and the I / p computation is replicated in
  ``morans_fast`` (identical numbers to ``morans_I_test`` given the same rng;
  checked in tests/test_power_negtraits.py). I_norm is not needed for power
  and is skipped.
* clade-shift scan: ``utils.phylo_tests._scan_with_null`` on the split
  matrix from ``clade_split_matrix`` -- i.e. ``run_clade_shift`` with no
  missing values (two-sided max-|z| Mann-Whitney over clades).
Permutations are reduced for speed (paper: Moran 1000, clade 10000);
defaults here --morans-perm 199, --clade-perm 999. With p = (1+k)/(1+B) the
p < 0.05 decision is exact for these B (both give attainable p-values with
B*0.05 integer-ish), so power is essentially unaffected; only the p-value
resolution changes.

Outputs
-------
--out-raw   per (pair, scenario, test): n_reps, n_reject, power, clade hit
            rate (scan's best split == planted clade, Jaccard >= 0.5)
--out       power per scenario x test x size bin (and bin 'ALL'), pooled
            over replicates of all trees in the bin
--fig       power curves (docs/figures/global/power_simulation.png)

Usage (repo root)
-----------------
python3 scripts/power_simulation_phylo.py --n-jobs 8
python3 scripts/power_simulation_phylo.py --dry-run           # plan + timing
python3 scripts/power_simulation_phylo.py --synthetic 6 --n-reps 20   # no data
python3 scripts/power_simulation_phylo.py --pairs 1dzlA_5keqF,2n54B_2hdmA
python3 scripts/power_simulation_phylo.py --plot-only          # redraw from --out
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-cache-" + os.environ.get("USER", "x"))
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)

import argparse
import multiprocessing as mp
import signal
import sys
import time
import zlib
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

from utils import phylo_tests as pt  # noqa: E402

CLADE_FRACS = (0.05, 0.10, 0.25, 0.50)
DELTAS = (0.5, 1.0, 1.5, 2.0, 3.0)
BM_H = (0.10, 0.25, 0.50)
RELIABILITIES = (1.0, 0.4, 0.1)
SIZE_BINS = ((20, 49, "20-49"), (50, 99, "50-99"), (100, 10 ** 9, "100+"))
TESTS = ("morans", "clade_shift")
ALPHA = 0.05


# ---------------------------------------------------------------------------
# Tree preparation (done once per tree)
# ---------------------------------------------------------------------------

def _bl(node) -> float:
    """Branch length with the ``_patristic_matrix`` convention."""
    return node.dist if (node.dist and node.dist > 0) else 1.0


def prepare_tree(tree, pair_id: str = "", min_clade: int = 2) -> dict:
    """Compact, picklable per-tree data for the simulation.

    names : leaf labels (all leaves of the coarse tree)
    W     : row-standardised inverse-patristic Moran weights (as morans_I_test)
    M     : (S, n) bool clade split matrix (as run_clade_shift)
    A     : (n, K) BM loading matrix, A[l, k] = sqrt(bl_k) if non-root node k
            is on the root->leaf-l path; B = A @ z / sqrt(mean depth)
    """
    names = [l.name for l in tree.iter_leaves()]
    n = len(names)
    eps = 1e-9
    W = 1.0 / (pt._patristic_matrix(tree, names) + eps)
    np.fill_diagonal(W, 0.0)
    rs = W.sum(axis=1, keepdims=True)
    rs[rs == 0] = 1.0
    W = W / rs
    M, _ = pt.clade_split_matrix(tree, names, min_clade)

    idx = {nm: j for j, nm in enumerate(names)}
    nodes = [nd for nd in tree.traverse("preorder") if not nd.is_root()]
    col = {id(nd): k for k, nd in enumerate(nodes)}
    A = np.zeros((n, len(nodes)))
    for leaf in tree.iter_leaves():
        j = idx[leaf.name]
        nd = leaf
        while nd is not None and not nd.is_root():
            A[j, col[id(nd)]] = np.sqrt(_bl(nd))
            nd = nd.up
    depth = (A ** 2).sum(axis=1)
    A = A / np.sqrt(depth.mean())               # mean marginal BM variance = 1
    return {"pair_id": pair_id, "names": names, "n": n, "W": W, "M": M,
            "A": A, "min_clade": min_clade}


def pick_clade(td: dict, frac: float) -> np.ndarray:
    """Split (clade) whose leaf fraction is closest to ``frac``; ties -> the
    first in clade_split_matrix (postorder) order. Returns a bool mask."""
    sizes = td["M"].sum(axis=1)
    s = int(np.argmin(np.abs(sizes - frac * td["n"])))
    return td["M"][s].copy()


# ---------------------------------------------------------------------------
# Scenarios and value simulation
# ---------------------------------------------------------------------------

def build_scenarios(fracs=CLADE_FRACS, deltas=DELTAS, hs=BM_H,
                    rels=RELIABILITIES) -> List[dict]:
    sc = [{"scenario": "null", "kind": "null", "f": np.nan, "delta": 0.0,
           "h": 0.0, "r": 1.0}]
    for r in rels:
        for f in fracs:
            for d in deltas:
                sc.append({"scenario": f"clade_f{f:.2f}_d{d:g}_r{r:g}",
                           "kind": "clade", "f": f, "delta": d, "h": np.nan,
                           "r": r})
        for h in hs:
            sc.append({"scenario": f"bm_h{h:.2f}_r{r:g}", "kind": "bm",
                       "f": np.nan, "delta": np.nan, "h": h, "r": r})
    return sc


def simulate(td: dict, sc: dict, rng: np.random.Generator,
             clade: Optional[np.ndarray] = None) -> np.ndarray:
    """One replicate of observed per-leaf values (see module docstring)."""
    n = td["n"]
    e = rng.standard_normal(n)
    if sc["kind"] == "clade":
        p = clade.mean()
        true = sc["delta"] * clade + e
        v_true = 1.0 + sc["delta"] ** 2 * p * (1 - p)
    elif sc["kind"] == "bm":
        h = sc["h"]
        B = td["A"] @ rng.standard_normal(td["A"].shape[1])
        true = np.sqrt(h) * B + np.sqrt(1 - h) * e
        v_true = 1.0
    else:
        return e
    r = sc["r"]
    if r >= 1.0:
        return true
    return true + rng.standard_normal(n) * np.sqrt(v_true * (1 - r) / r)


# ---------------------------------------------------------------------------
# Tests on precomputed tree data
# ---------------------------------------------------------------------------

def morans_fast(W: np.ndarray, x: np.ndarray, n_perm: int,
                rng: np.random.Generator, eps: float = 1e-9):
    """(I, p) exactly as ``pt.morans_I_test`` (same rng consumption), with the
    weight matrix precomputed. Returns (nan, nan) for zero variance."""
    xc = x - x.mean()
    den = float(xc @ xc)
    if den < eps:
        return np.nan, np.nan
    I_obs = float((xc @ (W @ xc)) / den)
    P = np.empty((n_perm, len(xc)))
    for k in range(n_perm):
        P[k] = rng.permutation(xc)
    I_perm = np.einsum("ij,ij->i", P @ W, P) / den
    return I_obs, (1 + int((I_perm >= I_obs).sum())) / (n_perm + 1)


def clade_fast(M: np.ndarray, x: np.ndarray, n_perm: int,
               rng: np.random.Generator, min_clade: int = 2):
    """(T, p, best split index) as ``pt.run_clade_shift`` (no missing values)."""
    n = len(x)
    r, tt = pt.average_ranks(x)
    if M.shape[0] == 0 or tt >= n ** 3 - n - 1e-9:
        return np.nan, np.nan, -1
    zc, _, T, p = pt._scan_with_null(M, r[None, :], np.ones((1, n), bool),
                                     np.array([tt]), min_clade, n_perm, rng)
    if not np.isfinite(T):
        return np.nan, np.nan, -1
    return T, p, int(np.nanargmax(np.abs(zc)))


def _jaccard(a: np.ndarray, b: np.ndarray) -> float:
    u = (a | b).sum()
    return float((a & b).sum() / u) if u else 0.0


def run_cell(td: dict, sc: dict, n_reps: int, seed: int,
             morans_perm: int, clade_perm: int,
             tests: Sequence[str] = TESTS) -> List[dict]:
    """All replicates of one (tree, scenario); one row per test."""
    ss = np.random.SeedSequence([seed, zlib.crc32(td["pair_id"].encode()),
                                 zlib.crc32(sc["scenario"].encode())])
    rng_val, rng_m, rng_c = [np.random.default_rng(s) for s in ss.spawn(3)]
    clade = pick_clade(td, sc["f"]) if sc["kind"] == "clade" else None
    rej = {t: 0 for t in tests}
    ok = {t: 0 for t in tests}
    hits = 0
    pvals = {t: [] for t in tests}
    for _ in range(n_reps):
        x = simulate(td, sc, rng_val, clade)
        if "morans" in tests:
            _, p = morans_fast(td["W"], x, morans_perm, rng_m)
            if np.isfinite(p):
                ok["morans"] += 1
                rej["morans"] += p < ALPHA
                pvals["morans"].append(p)
        if "clade_shift" in tests:
            _, p, s = clade_fast(td["M"], x, clade_perm, rng_c, td["min_clade"])
            if np.isfinite(p):
                ok["clade_shift"] += 1
                rej["clade_shift"] += p < ALPHA
                pvals["clade_shift"].append(p)
                if clade is not None and s >= 0 and p < ALPHA:
                    hits += _jaccard(td["M"][s], clade) >= 0.5
    rows = []
    for t in tests:
        rows.append({
            "pair_id": td["pair_id"], "n_leaves": td["n"],
            "size_bin": size_bin(td["n"]), **sc,
            "clade_frac_actual": float(clade.mean()) if clade is not None else np.nan,
            "test": t, "n_reps": ok[t], "n_reject": int(rej[t]),
            "power": rej[t] / ok[t] if ok[t] else np.nan,
            "median_p": float(np.median(pvals[t])) if pvals[t] else np.nan,
            "clade_hit_rate": (hits / ok[t] if (t == "clade_shift" and clade is not None
                                                and ok[t]) else np.nan),
        })
    return rows


def size_bin(n: int) -> str:
    for lo, hi, lab in SIZE_BINS:
        if lo <= n <= hi:
            return lab
    return f"<{SIZE_BINS[0][0]}"


# ---------------------------------------------------------------------------
# Tree loading (real data on the cluster; synthetic for local tests)
# ---------------------------------------------------------------------------

class _Node:
    """ete3-like stand-in used for --synthetic trees."""

    def __init__(self, name="", dist=1.0):
        self.name, self.dist, self.children, self.up = name, dist, [], None

    def add_child(self, c):
        c.up = self
        self.children.append(c)
        return c

    def is_leaf(self):
        return not self.children

    def is_root(self):
        return self.up is None

    def traverse(self, strategy="levelorder"):
        if strategy == "postorder":
            for c in self.children:
                yield from c.traverse("postorder")
            yield self
        else:
            yield self
            for c in self.children:
                yield from c.traverse("preorder")

    def iter_leaves(self):
        for nd in self.traverse("preorder"):
            if nd.is_leaf():
                yield nd


def random_tree(n_leaves: int, seed: int = 0) -> _Node:
    """Random rooted binary tree (random joining), exponential branch lengths."""
    rs = np.random.RandomState(seed)
    pool = [_Node(f"S{i:02d}", float(rs.exponential(0.3)) + 0.01)
            for i in range(n_leaves)]
    k = 0
    while len(pool) > 1:
        i, j = sorted(rs.choice(len(pool), 2, replace=False))
        b, a = pool.pop(j), pool.pop(i)
        p = _Node(f"N{k}", float(rs.exponential(0.3)) + 0.01)
        k += 1
        p.add_child(a)
        p.add_child(b)
        pool.append(p)
    pool[0].dist = 0.0
    return pool[0]


def _alarm(signum, frame):
    raise TimeoutError("coarse-tree build timed out")


def load_real_tree(args_tuple):
    """Worker: build + prepare one pair's coarse tree. Returns dict or str."""
    pair_id, min_clade, timeout = args_tuple
    try:
        from utils.ancestral_utils import build_coarse_tree, _name_internal_nodes
        signal.signal(signal.SIGALRM, _alarm)
        signal.alarm(timeout)
        try:
            tree, _, _ = build_coarse_tree(pair_id)
            _name_internal_nodes(tree)
        finally:
            signal.alarm(0)
        return prepare_tree(tree, pair_id, min_clade)
    except Exception as e:  # noqa: BLE001
        return f"{pair_id}: {type(e).__name__}: {e}"


def discover_pairs() -> List[str]:
    from config import DATA_DIR, PAIR_DIR_RE
    if not os.path.isdir(DATA_DIR):
        return []
    return sorted(d for d in os.listdir(DATA_DIR)
                  if PAIR_DIR_RE.match(d) and os.path.isdir(os.path.join(DATA_DIR, d)))


def select_pairs(pairs: List[str], n_trees: int, min_clusters: int,
                 sizes_csv: str, seed: int) -> List[str]:
    """Pre-select ~n_trees pairs spread evenly over SIZE_BINS using cluster
    counts from ``sizes_csv`` (docs/phylo_placement.csv, column n_clusters).
    Without the CSV, a seeded random subset (size filter applied after build)."""
    rng = np.random.default_rng(seed)
    if not os.path.exists(sizes_csv):
        print(f"[select] {sizes_csv} missing; random subset, size filter after build")
        if n_trees and len(pairs) > n_trees:
            # oversample to survive the size filter
            pairs = sorted(rng.choice(pairs, min(len(pairs), 2 * n_trees), replace=False))
        return list(pairs)
    sz = pd.read_csv(sizes_csv).groupby("pair_id")["n_clusters"].max()
    sz = sz[sz.index.isin(pairs) & (sz >= min_clusters)]
    return stratify(sz, n_trees, rng)


def stratify(sz: pd.Series, n_trees: int, rng: np.random.Generator) -> List[str]:
    """Pick n_trees ids from {id: n_leaves}, spread round-robin over SIZE_BINS."""
    if not n_trees or len(sz) <= n_trees:
        return sorted(sz.index)
    bins = {lab: sorted(sz[(sz >= lo) & (sz <= hi)].index) for lo, hi, lab in SIZE_BINS}
    bins = {k: v for k, v in bins.items() if v}
    quota = {k: 0 for k in bins}
    left = n_trees
    while left > 0:                       # round-robin, capped by bin size
        moved = False
        for k in bins:
            if left and quota[k] < len(bins[k]):
                quota[k] += 1
                left -= 1
                moved = True
        if not moved:
            break
    chosen = []
    for k, v in bins.items():
        chosen += [str(c) for c in rng.choice(v, quota[k], replace=False)]
    return sorted(chosen)


# ---------------------------------------------------------------------------
# Aggregation and figure
# ---------------------------------------------------------------------------

def aggregate(raw: pd.DataFrame) -> pd.DataFrame:
    keys = ["scenario", "kind", "f", "delta", "h", "r", "test"]
    parts = []
    for b, sub in list(raw.groupby("size_bin")) + [("ALL", raw)]:
        g = sub.groupby(keys, dropna=False).agg(
            n_trees=("pair_id", "nunique"), n_reps=("n_reps", "sum"),
            n_reject=("n_reject", "sum"),
            mean_clade_frac=("clade_frac_actual", "mean"),
            clade_hit_rate=("clade_hit_rate", "mean")).reset_index()
        g["power"] = g["n_reject"] / g["n_reps"]
        g["power_se"] = np.sqrt(g["power"] * (1 - g["power"]) / g["n_reps"])
        g.insert(0, "size_bin", b)
        parts.append(g)
    out = pd.concat(parts, ignore_index=True)
    order = {lab: i for i, (_, _, lab) in enumerate(SIZE_BINS)}
    order["ALL"] = 99
    out["_o"] = out["size_bin"].map(order).fillna(98)
    return out.sort_values(["_o", "kind", "r", "f", "delta", "h", "test"]
                           ).drop(columns="_o").reset_index(drop=True)


def plot_power(agg: pd.DataFrame, fig_path: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    col = {"morans": "#2a78d6", "clade_shift": "#eb6834"}
    lab = {"morans": "Moran's I", "clade_shift": "clade-shift scan"}
    ls = {1.0: "-", 0.4: "--", 0.1: ":"}
    mk = {1.0: "o", 0.4: "s", 0.1: "^"}
    a = agg[agg["size_bin"] == "ALL"]
    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    fracs = sorted(a.loc[a.kind == "clade", "f"].dropna().unique())
    for ax, f in zip(axes.flat[:4], fracs):
        sub = a[(a.kind == "clade") & np.isclose(a.f, f)]
        for t in TESTS:
            for r in sorted(sub.r.unique(), reverse=True):
                s = sub[(sub.test == t) & (sub.r == r)].sort_values("delta")
                ax.plot(s.delta, s.power, color=col[t], ls=ls.get(r, "-"),
                        marker=mk.get(r, "o"), ms=5, lw=2)
        fa = sub["mean_clade_frac"].mean()
        ax.set_title(f"single-clade shift, f = {f:.2f} (realised {fa:.2f})", fontsize=10)
        ax.set_xlabel("effect δ (noise SD units)")
    ax = axes.flat[4]
    sub = a[a.kind == "bm"]
    for t in TESTS:
        for r in sorted(sub.r.unique(), reverse=True):
            s = sub[(sub.test == t) & (sub.r == r)].sort_values("h")
            ax.plot(s.h, s.power, color=col[t], ls=ls.get(r, "-"),
                    marker=mk.get(r, "o"), ms=5, lw=2)
    ax.set_title("diffuse Brownian-motion trait", fontsize=10)
    ax.set_xlabel("BM variance fraction h")
    # null calibration by size bin
    ax = axes.flat[5]
    nl = agg[agg.kind == "null"]
    bins = [b for b in [s[2] for s in SIZE_BINS] + ["ALL"] if b in set(nl.size_bin)]
    xs = np.arange(len(bins))
    for k, t in enumerate(TESTS):
        v = [nl[(nl.size_bin == b) & (nl.test == t)].power.mean() for b in bins]
        e = [nl[(nl.size_bin == b) & (nl.test == t)].power_se.mean() for b in bins]
        ax.bar(xs + (k - 0.5) * 0.36, v, 0.34, yerr=e, color=col[t],
               edgecolor="white", linewidth=2, capsize=3)
    ax.axhline(ALPHA, color="#52514e", lw=1, ls="--")
    ax.text(-0.45, ALPHA + 0.004, "α = 0.05", ha="left", fontsize=8,
            color="#52514e")
    ax.set_xticks(xs)
    ax.set_xticklabels(bins)
    ax.set_title("null: false-positive rate by tree size", fontsize=10)
    ax.set_xlabel("clusters per tree")
    for k, ax in enumerate(axes.flat):
        ax.set_ylim(0, 1.02 if k < 5 else 0.2)
        ax.grid(axis="y", color="#e4e3df", lw=0.8)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[0, 0].set_ylabel("power (p < 0.05)")
    axes[1, 0].set_ylabel("power (p < 0.05)")
    axes[1, 2].set_ylabel("rejection rate under the null")
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=col[t], lw=2, label=lab[t]) for t in TESTS]
    handles += [Line2D([], [], color="#52514e", ls=ls[r], marker=mk[r], lw=1.5,
                       label=f"reliability r = {r:g}") for r in (1.0, 0.4, 0.1)]
    n_tr = int(a["n_trees"].max()) if len(a) else 0
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False,
               bbox_to_anchor=(0.5, -0.01))
    fig.suptitle(f"Power of phylogenetic fold-preference tests on real coarse "
                 f"cluster trees ({n_tr} trees)", fontsize=12)
    fig.tight_layout(rect=(0, 0.05, 1, 0.97))
    os.makedirs(os.path.dirname(os.path.abspath(fig_path)), exist_ok=True)
    fig.savefig(fig_path, dpi=200)
    plt.close(fig)
    print(f"[fig] {fig_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

_TREES: Dict[str, dict] = {}


def _init_worker(trees):
    global _TREES
    _TREES = trees


def _cell_worker(job):
    pair_id, sc, n_reps, seed, mp_, cp_ = job
    return run_cell(_TREES[pair_id], sc, n_reps, seed, mp_, cp_)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--pairs", default="ALL",
                    help="ALL (selection by size) or comma-separated pair ids")
    ap.add_argument("--n-trees", type=int, default=20,
                    help="trees to sample across size bins when --pairs ALL (0 = all)")
    ap.add_argument("--min-clusters", type=int, default=20)
    ap.add_argument("--max-pairs", type=int, default=0, help="cap #pairs (quick tests)")
    ap.add_argument("--sizes-csv", default=os.path.join(REPO, "docs", "phylo_placement.csv"))
    ap.add_argument("--synthetic", type=int, default=0,
                    help="use N synthetic random trees (sizes 20..100) instead of data")
    ap.add_argument("--n-reps", type=int, default=100)
    ap.add_argument("--morans-perm", type=int, default=199)
    ap.add_argument("--clade-perm", type=int, default=999)
    ap.add_argument("--min-clade", type=int, default=2)
    ap.add_argument("--fracs", default=",".join(map(str, CLADE_FRACS)))
    ap.add_argument("--deltas", default=",".join(map(str, DELTAS)))
    ap.add_argument("--hs", default=",".join(map(str, BM_H)))
    ap.add_argument("--rels", default=",".join(map(str, RELIABILITIES)))
    ap.add_argument("--seed", type=int, default=20261021)
    ap.add_argument("--n-jobs", type=int, default=1)
    ap.add_argument("--tree-timeout", type=int, default=300)
    ap.add_argument("--out", default=os.path.join(REPO, "docs", "power_simulation.csv"))
    ap.add_argument("--out-raw", default=os.path.join(REPO, "docs", "power_simulation_per_tree.csv"))
    ap.add_argument("--fig", default=os.path.join(REPO, "docs", "figures", "global",
                                                   "power_simulation.png"))
    ap.add_argument("--dry-run", action="store_true",
                    help="build trees, time one cell per tree, print plan; no outputs")
    ap.add_argument("--plot-only", action="store_true", help="redraw figure from --out")
    args = ap.parse_args(argv)

    if args.plot_only:
        plot_power(pd.read_csv(args.out), args.fig)
        return

    fl = lambda s: tuple(float(v) for v in s.split(",") if v)  # noqa: E731
    scenarios = build_scenarios(fl(args.fracs), fl(args.deltas), fl(args.hs), fl(args.rels))
    n_jobs = max(1, args.n_jobs)

    t0 = time.time()
    trees: Dict[str, dict] = {}
    if args.synthetic:
        sizes = np.linspace(20, 100, args.synthetic).round().astype(int)
        for i, n in enumerate(sizes):
            pid = f"synthetic_{i:02d}_n{n}"
            trees[pid] = prepare_tree(random_tree(int(n), seed=args.seed + i), pid,
                                      args.min_clade)
    else:
        if args.pairs.upper() == "ALL":
            pairs = select_pairs(discover_pairs(), args.n_trees, args.min_clusters,
                                 args.sizes_csv, args.seed)
        else:
            pairs = [p.strip() for p in args.pairs.split(",") if p.strip()]
        if args.max_pairs:
            pairs = pairs[: args.max_pairs]
        print(f"[trees] building {len(pairs)} coarse trees with {n_jobs} jobs...", flush=True)
        jobs = [(p, args.min_clade, args.tree_timeout) for p in pairs]
        if n_jobs > 1:
            with mp.Pool(n_jobs) as pool:
                res = pool.map(load_real_tree, jobs, chunksize=1)
        else:
            res = [load_real_tree(j) for j in jobs]
        for r in res:
            if isinstance(r, str):
                print(f"  [skip] {r}")
            elif r["n"] < args.min_clusters:
                print(f"  [skip] {r['pair_id']}: {r['n']} leaves < {args.min_clusters}")
            else:
                trees[r["pair_id"]] = r
        if args.pairs.upper() == "ALL" and args.n_trees and len(trees) > args.n_trees:
            keep = stratify(pd.Series({k: v["n"] for k, v in trees.items()}),
                            args.n_trees, np.random.default_rng(args.seed))
            trees = {k: trees[k] for k in keep}
    if not trees:
        sys.exit("[error] no usable trees")
    print(f"[trees] {len(trees)} trees ready in {time.time() - t0:.1f}s; sizes: "
          + ", ".join(f"{k}={v['n']}" for k, v in sorted(trees.items())), flush=True)

    n_cells = len(trees) * len(scenarios)
    print(f"[plan] {len(scenarios)} scenarios x {len(trees)} trees = {n_cells} cells, "
          f"{args.n_reps} reps each; perms Moran {args.morans_perm}, clade {args.clade_perm}")

    if args.dry_run:
        tt = []
        sc = [s for s in scenarios if s["kind"] == "clade"][0]
        for pid, td in sorted(trees.items()):
            t1 = time.time()
            run_cell(td, sc, 5, args.seed, args.morans_perm, args.clade_perm)
            tt.append((time.time() - t1) / 5)
            print(f"  {pid}: n={td['n']} splits={td['M'].shape[0]} "
                  f"{tt[-1] * 1000:.0f} ms/replicate (both tests)")
        est = np.mean(tt) * args.n_reps * n_cells / n_jobs
        print(f"[plan] estimated wall time ~{est / 60:.1f} min on {n_jobs} jobs")
        return

    jobs = [(pid, sc, args.n_reps, args.seed, args.morans_perm, args.clade_perm)
            for pid in sorted(trees) for sc in scenarios]
    rows: List[dict] = []
    t1 = time.time()
    if n_jobs > 1:
        with mp.Pool(n_jobs, initializer=_init_worker, initargs=(trees,)) as pool:
            for k, out in enumerate(pool.imap_unordered(_cell_worker, jobs, chunksize=2), 1):
                rows += out
                if k % max(1, len(jobs) // 20) == 0:
                    print(f"  {k}/{len(jobs)} cells, {time.time() - t1:.0f}s", flush=True)
    else:
        _init_worker(trees)
        for k, j in enumerate(jobs, 1):
            rows += _cell_worker(j)
            if k % max(1, len(jobs) // 20) == 0:
                print(f"  {k}/{len(jobs)} cells, {time.time() - t1:.0f}s", flush=True)
    raw = pd.DataFrame(rows).sort_values(["pair_id", "scenario", "test"])
    for path in (args.out, args.out_raw):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    raw.to_csv(args.out_raw, index=False)
    agg = aggregate(raw)
    agg.to_csv(args.out, index=False)
    print(f"[out] {args.out_raw}\n[out] {args.out}  ({time.time() - t0:.0f}s total)")

    a = agg[agg.size_bin == "ALL"]
    show = a[(a.kind == "null") | ((a.r == 1.0) & (a.kind == "bm"))
             | ((a.kind == "clade") & a.delta.isin([1.0, 2.0]))]
    with pd.option_context("display.width", 160, "display.max_rows", 200):
        print(show.pivot_table(index=["kind", "f", "delta", "h", "r"], columns="test",
                               values="power", dropna=False).dropna(how="all").round(3))
    plot_power(agg, args.fig)


if __name__ == "__main__":
    main()
