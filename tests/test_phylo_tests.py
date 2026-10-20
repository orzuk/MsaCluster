#!/usr/bin/env python3
"""Unit tests for utils/phylo_tests.py (tree-test registry + clade-shift scan).

Run either way (from the repo root):
    python3 -m pytest tests/test_phylo_tests.py -q
    python3 tests/test_phylo_tests.py

Trees are built with a tiny ete3-compatible stand-in (``MiniNode``) so the
tests run without ete3; when ete3 is installed, ``test_ete3_tree_agrees``
also checks that a real ete3 tree gives the same numbers.
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils import phylo_tests as pt  # noqa: E402


class SkipTest(Exception):
    pass


def _skip(msg):
    try:
        import pytest
        pytest.skip(msg)
    except ImportError:
        raise SkipTest(msg)


# ---------------------------------------------------------------------------
# Minimal ete3-like tree
# ---------------------------------------------------------------------------

class MiniNode:
    """Subset of the ete3.TreeNode API used by utils/phylo_tests.py."""

    def __init__(self, name="", dist=1.0):
        self.name = name
        self.dist = dist
        self.children = []
        self.up = None

    def add_child(self, child):
        child.up = self
        self.children.append(child)
        return child

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
        for n in self.traverse("preorder"):
            if n.is_leaf():
                yield n

    def write_newick(self):
        if self.is_leaf():
            return f"{self.name}:{self.dist}"
        inner = ",".join(c.write_newick() for c in self.children)
        return f"({inner}){self.name}:{self.dist}"


def random_tree(n_leaves, seed=0):
    """Random rooted binary tree (random joining) with exponential branch
    lengths; leaves L000.., internal nodes N0.."""
    rs = np.random.RandomState(seed)
    pool = [MiniNode(f"L{i:03d}", float(rs.exponential(0.3)) + 0.01)
            for i in range(n_leaves)]
    k = 0
    while len(pool) > 1:
        i, j = sorted(rs.choice(len(pool), 2, replace=False))
        b, a = pool.pop(j), pool.pop(i)
        p = MiniNode(f"N{k}", float(rs.exponential(0.3)) + 0.01)
        k += 1
        p.add_child(a)
        p.add_child(b)
        pool.append(p)
    root = pool[0]
    root.dist = 0.0
    return root


def _clade_of_size(tree, size, lo=4, hi=7):
    """Leaf names of an internal node with exactly ``size`` leaves (else the
    closest size in [lo, hi]); None if none."""
    best = None
    for node in tree.traverse("postorder"):
        if node.is_leaf() or node.is_root():
            continue
        leaves = [l.name for l in node.iter_leaves()]
        if lo <= len(leaves) <= hi:
            if best is None or abs(len(leaves) - size) < abs(len(best) - size):
                best = leaves
    return best


def _planted(seed, n=60, size=5, shift=1.5):
    tree = random_tree(n, seed=seed)
    clade = _clade_of_size(tree, size)
    rs = np.random.RandomState(1000 + seed)
    vals = {l.name: float(rs.normal()) for l in tree.iter_leaves()}
    for nm in clade:
        vals[nm] += shift
    return tree, vals, clade


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_vectorized_z_matches_scipy():
    """Matrix-product rank sums + tie-corrected z == scipy mannwhitneyu."""
    try:
        from scipy.stats import mannwhitneyu, norm
    except ImportError:
        _skip("scipy not installed")
    for seed in range(5):
        tree = random_tree(30, seed=seed)
        names = [l.name for l in tree.iter_leaves()]
        rs = np.random.RandomState(seed)
        x = np.round(rs.normal(size=len(names)), 1)      # rounding -> ties
        M, _ = pt.clade_split_matrix(tree, names, min_clade=2)
        r, tt = pt.average_ranks(x)
        n1 = M.sum(axis=1)
        z = pt.mann_whitney_z(M.astype(float) @ r, n1, len(names), tt)
        for s in range(M.shape[0]):
            res = mannwhitneyu(x[M[s]], x[~M[s]], alternative="two-sided",
                               method="asymptotic", use_continuity=False)
            U1 = r[M[s]].sum() - n1[s] * (n1[s] + 1) / 2
            assert abs(U1 - res.statistic) < 1e-9
            assert abs(2 * norm.sf(abs(z[s])) - res.pvalue) < 1e-9


def test_split_matrix_dedup():
    """Root's two children give one split (smaller side kept); sizes OK."""
    tree = random_tree(12, seed=3)
    names = [l.name for l in tree.iter_leaves()]
    M, nodes = pt.clade_split_matrix(tree, names, min_clade=2)
    keys = {tuple(row) if not row[0] else tuple(~row) for row in M}
    assert len(keys) == M.shape[0]
    assert (M.sum(axis=1) >= 2).all() and ((~M).sum(axis=1) >= 2).all()


# Pre-refactor outputs (scripts/phylo_placement.py @ eb3f845) on the inputs
# built in _legacy_inputs(); the registry must reproduce them exactly.
_LEGACY_REF = {
    1: {"morans": {"morans_I": 0.0059, "morans_p": 0.2236, "morans_z": 0.811,
                   "morans_lambda_max": 0.488, "morans_max_perm": 0.4337,
                   "morans_I_norm": 0.1035, "morans_n": 24, "morans_note": ""},
        "parsimony": (7, 0.7662), "nn1": (0.5714, 0.408), "D": 1.3508},
    2: {"morans": {"morans_I": -0.0478, "morans_p": 0.509, "morans_z": -0.039,
                   "morans_lambda_max": 0.4356, "morans_max_perm": 0.3899,
                   "morans_I_norm": -0.0153, "morans_n": 24, "morans_note": ""},
        "parsimony": (7, 0.6368), "nn1": (0.5, 0.4428), "D": 1.0773},
    3: {"morans": {"morans_I": -0.0787, "morans_p": 0.7126, "morans_z": -0.607,
                   "morans_lambda_max": 0.5033, "morans_max_perm": 0.3912,
                   "morans_I_norm": -0.1823, "morans_n": 24, "morans_note": ""},
        "parsimony": (6, 0.6866), "nn1": (0.35, 0.9353), "D": 1.0448},
}


def _legacy_inputs(s):
    t = random_tree(25, seed=s)
    rs = np.random.RandomState(100 + s)
    leaves = [l.name for l in t.iter_leaves()]
    vals = {l: float(rs.normal()) for l in leaves}
    vals[leaves[3]] = None
    states = {l: ("F1" if (v is not None and v > 0.3) else
                  ("F2" if (v is not None and v < -0.3) else "Amb"))
              for l, v in vals.items()}
    return t, vals, states


def test_legacy_outputs_unchanged():
    for s, ref in _LEGACY_REF.items():
        t, vals, states = _legacy_inputs(s)
        m = pt.run_test("morans", t, vals, 500,
                        rng=np.random.default_rng(12345 + 7))
        assert m == ref["morans"], (s, m)
        # direct call with a seed == registry call with default_rng(seed)
        assert pt.morans_I_test(t, vals, n_perm=500, seed=12345 + 7) == m
        f = pt.run_test("fitch", t, states, 200, seed=12345)
        assert (f["parsimony_score"], f["parsimony_p"]) == ref["parsimony"]
        nn = pt.run_test("nn", t, states, 200, seed=12346)
        assert (nn["nn_concordance"], nn["nn_concordance_p"]) == ref["nn1"]
        d = pt.run_test("dstat", t, states, 50, seed=12348)
        assert d["D_statistic"] == ref["D"]


def _planted_power(shift, seeds):
    p_clade, p_mor, n_found = [], [], 0
    for seed in seeds:
        tree, vals, clade = _planted(seed, shift=shift)
        r = pt.run_test("clade_shift", tree, vals, 2000,
                        rng=np.random.default_rng(seed))
        m = pt.run_test("morans", tree, vals, 1000,
                        rng=np.random.default_rng(seed))
        p_clade.append(r["clade_p"])
        p_mor.append(m["morans_p"])
        best = set(r["clade_best_members"].split(";"))
        n_found += len(best & set(clade)) >= 0.6 * len(clade)
    return np.asarray(p_clade), np.asarray(p_mor), n_found / len(seeds)


def test_planted_clade_detected():
    """One 5-leaf clade among 60 leaves shifted, N(0,1) noise.

    Measured (40 seeds): shift 1.5 SD -> clade p<.05 in 50%, p<.01 in 20%
    (Moran: 30% / 10%); shift 2.5 SD -> clade p<.01 in 78% (Moran 38%).
    A rank test on 5 vs 55 leaves caps |z| at ~3.7, so 1.5 SD is near the
    detection limit; the 2.5 SD case checks the "p<.01 in most seeds" goal.
    """
    seeds = range(20)
    pc, pm, found = _planted_power(1.5, seeds)
    assert np.mean(pc < 0.05) >= 0.4, pc
    assert np.mean(pc < 0.05) > np.mean(pm < 0.05), (pc, pm)
    assert np.median(pc) < np.median(pm), (np.median(pc), np.median(pm))
    pc, pm, found = _planted_power(2.5, seeds)
    assert np.mean(pc < 0.01) >= 0.7, pc
    assert found >= 0.8, found
    assert np.mean(pc < 0.01) > np.mean(pm < 0.01), (pc, pm)


def test_null_calibrated():
    """No shift: P(p < 0.05) over 50 seeds within [0.01, 0.12]."""
    ps = []
    for seed in range(50):
        tree = random_tree(40, seed=500 + seed)
        rs = np.random.RandomState(seed)
        vals = {l.name: float(rs.normal()) for l in tree.iter_leaves()}
        r = pt.run_test("clade_shift", tree, vals, 999,
                        rng=np.random.default_rng(seed))
        ps.append(r["clade_p"])
    frac = float(np.mean(np.asarray(ps) < 0.05))
    assert 0.01 <= frac <= 0.12, (frac, sorted(ps)[:6])


def test_combined_null_calibrated():
    """No shift, 4 correlated methods (one with a missing leaf): the joint
    leaf permutation keeps the COMBINED p calibrated (measured over 200
    seeds: P(p<.05) = 0.035)."""
    ps = []
    for seed in range(50):
        tree = random_tree(40, seed=900 + seed)
        names = [l.name for l in tree.iter_leaves()]
        rs = np.random.RandomState(seed)
        common = rs.normal(size=len(names))
        per = {m: {nm: float(0.7 * common[i] + 0.7 * rs.normal())
                   for i, nm in enumerate(names)} for m in "ABCD"}
        per["B"][names[2]] = None
        r = pt.run_test("clade_shift_combined", tree, per, 999,
                        rng=np.random.default_rng(seed))
        ps.append(r["clade_p"])
    frac = float(np.mean(np.asarray(ps) < 0.05))
    assert 0.0 <= frac <= 0.12, frac


def test_combined_and_missing():
    """COMBINED: 3 methods share a planted shift -> stronger than any single
    method; missing values / zero-variance methods handled."""
    tree, _, clade = _planted(4)
    names = [l.name for l in tree.iter_leaves()]
    rs = np.random.RandomState(7)
    per = {}
    for m in ("A", "B", "C"):
        v = {nm: float(rs.normal()) + (1.2 if nm in clade else 0.0)
             for nm in names}
        per[m] = v
    per["B"][names[0]] = None                 # missing value
    per["B"][names[1]] = float("nan")
    per["Z"] = {nm: 0.0 for nm in names}      # zero variance -> skipped
    comb = pt.run_test("clade_shift_combined", tree, per, 2000,
                       rng=np.random.default_rng(0))
    singles = [pt.run_test("clade_shift", tree, per[m], 2000,
                           rng=np.random.default_rng(0)) for m in "ABC"]
    assert comb["clade_n_methods"] == 3 and "Z(zero_variance)" in comb["clade_note"]
    assert comb["clade_T"] > max(s["clade_T"] for s in singles)
    assert comb["clade_p"] < 0.01
    assert comb["clade_best_direction"] == 1
    assert singles[1]["clade_n_leaves"] == len(names) - 2
    # zero-variance single method
    z = pt.run_test("clade_shift", tree, per["Z"], 100,
                    rng=np.random.default_rng(0))
    assert z["clade_note"] == "zero_variance" and np.isnan(z["clade_p"])


def test_members_mapped_to_tags():
    tree, vals, clade = _planted(0)
    l2t = {nm: f"ShallowMsa_{int(nm[1:]):03d}" for nm in vals}
    r = pt.run_test("clade_shift", tree, vals, 200,
                    rng=np.random.default_rng(0), label_to_tag=l2t)
    assert all(m.startswith("ShallowMsa_") for m in
               r["clade_best_members"].split(";"))


def test_ete3_tree_agrees():
    try:
        from ete3 import Tree
    except ImportError:
        _skip("ete3 not installed")
    tree, vals, _ = _planted(2)
    et = Tree(tree.write_newick() + ";", format=1)
    a = pt.run_test("clade_shift", tree, vals, 500, rng=np.random.default_rng(1))
    b = pt.run_test("clade_shift", et, vals, 500, rng=np.random.default_rng(1))
    assert a["clade_T"] == b["clade_T"]
    assert a["clade_best_members"] == b["clade_best_members"]
    ma = pt.run_test("morans", tree, vals, 200, rng=np.random.default_rng(1))
    mb = pt.run_test("morans", et, vals, 200, rng=np.random.default_rng(1))
    assert ma["morans_I"] == mb["morans_I"]


if __name__ == "__main__":
    import time
    n_fail = 0
    for name, fn in list(globals().items()):
        if not (name.startswith("test_") and callable(fn)):
            continue
        t0 = time.time()
        try:
            fn()
            print(f"PASS  {name} ({time.time() - t0:.1f}s)")
        except SkipTest as e:
            print(f"SKIP  {name}: {e}")
        except AssertionError as e:
            n_fail += 1
            print(f"FAIL  {name}: {e}")
    sys.exit(1 if n_fail else 0)
