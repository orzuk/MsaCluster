#!/usr/bin/env python3
"""Unit tests for utils/cluster_reps.py (coarse-tree representative leaves).

Run either way (from the repo root):
    python3 -m pytest tests/test_cluster_reps.py -q
    python3 tests/test_cluster_reps.py

Uses the ete3-free ``MiniNode`` stand-in from tests/test_phylo_tests.py.
"""

from __future__ import annotations

import itertools
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from utils.cluster_reps import (choose_cluster_reps, cluster_pieces,  # noqa: E402
                                piece_medoid)
from test_phylo_tests import MiniNode, random_tree  # noqa: E402


def T(name="", dist=1.0, *children):
    n = MiniNode(name, dist)
    for c in children:
        n.add_child(c)
    return n


def _groups(tree, leaf_to_tag):
    g = {}
    for lf in tree.iter_leaves():
        if lf.name in leaf_to_tag:
            g.setdefault(leaf_to_tag[lf.name], []).append(lf.name)
    return g


def _toy():
    """Cluster X: small piece {x1} listed first (sister of y1), large piece
    {x2..x5} second, containing an unassigned leaf u that must be ignored.
    Cluster Z is monophyletic. x3 is the medoid of the large X piece."""
    tree = T("root", 0.0,
             T("", 1.0, T("x1", 1.0), T("y1", 1.0)),
             T("", 1.0,
               T("", 0.1, T("x2", 0.5), T("x3", 0.1)),
               T("", 0.1, T("", 0.2, T("x4", 1.0), T("u", 0.3)), T("x5", 1.0))),
             T("", 1.0, T("y2", 1.0), T("y3", 1.0)),
             T("", 2.0, T("z1", 0.2), T("", 0.1, T("z2", 0.2), T("z3", 0.9))))
    l2t = {**{f"x{i}": "X" for i in range(1, 6)},
           **{f"y{i}": "Y" for i in range(1, 4)},
           **{f"z{i}": "Z" for i in range(1, 4)}}
    return tree, l2t


def test_first_vs_largest_piece_medoid():
    tree, l2t = _toy()
    groups = _groups(tree, l2t)
    assert groups["X"][0] == "x1"
    reps_f, info_f = choose_cluster_reps(tree, groups, "first")
    reps_m, info_m = choose_cluster_reps(tree, groups, "largest_piece_medoid")
    assert reps_f["X"] == "x1"
    assert info_f["X"]["rep_piece_size"] == 1
    assert abs(info_f["X"]["rep_piece_frac"] - 0.2) < 1e-12
    assert reps_m["X"] in {"x2", "x3", "x4", "x5"}
    assert reps_m["X"] == "x3"           # medoid of the large piece
    assert info_m["X"]["n_pieces"] == 2  # u ignored -> {x1}, {x2..x5}
    assert info_m["X"]["rep_piece_size"] == 4
    assert abs(info_m["X"]["rep_piece_frac"] - 0.8) < 1e-12
    # Y: pieces {y1}, {y2,y3}
    assert info_m["Y"]["n_pieces"] == 2 and reps_m["Y"] in {"y2", "y3"}
    # Z monophyletic: both strategies keep a Z leaf, whole cluster in piece
    for reps, info in ((reps_f, info_f), (reps_m, info_m)):
        assert reps["Z"] in {"z1", "z2", "z3"}
        assert info["Z"]["n_pieces"] == 1 and info["Z"]["rep_piece_frac"] == 1.0
    assert reps_m["Z"] == "z2"   # sums: z1 1.7, z2 1.6, z3 2.3
    _check_medoid_bruteforce(tree, l2t)


def _patristic(a, b):
    def path(n):
        out = []
        while n is not None:
            out.append(n)
            n = n.up
        return out
    pa, pb = path(a), path(b)
    sb = set(map(id, pb))
    lca = next(n for n in pa if id(n) in sb)
    d = 0.0
    for p in (pa, pb):
        for n in p:
            if n is lca:
                break
            d += n.dist
    return d


def _check_medoid_bruteforce(tree, l2t):
    pieces, _ = cluster_pieces(tree, l2t)
    for tag, plist in pieces.items():
        for p in plist:
            mem = [lf for lf in p.iter_leaves() if l2t.get(lf.name) == tag]
            sums = [sum(_patristic(a, b) for b in mem) for a in mem]
            best = min(sums)
            got = piece_medoid(p, {m.name for m in mem})
            got_sum = sums[[m.name for m in mem].index(got)]
            assert abs(got_sum - best) < 1e-9, (tag, got, got_sum, best)
            # every member of the piece is labelled tag; parent is not pure
            assert all(l2t.get(lf.name) in (tag, None) for lf in p.iter_leaves())


def test_random_trees_pieces_and_medoids():
    for seed in range(20):
        tree = random_tree(80, seed=seed)
        rs = np.random.RandomState(seed)
        names = [lf.name for lf in tree.iter_leaves()]
        l2t = {}
        for nm in names:
            r = rs.rand()
            if r < 0.15:
                continue                 # unassigned leaf
            l2t[nm] = f"C{rs.randint(4)}"
        groups = _groups(tree, l2t)
        pieces, n_assigned = cluster_pieces(tree, l2t)
        # pieces partition each cluster's members
        for tag, mem in groups.items():
            covered = list(itertools.chain.from_iterable(
                [lf.name for lf in p.iter_leaves() if lf.name in l2t]
                for p in pieces[tag]))
            assert sorted(covered) == sorted(mem)
        _check_medoid_bruteforce(tree, l2t)
        reps, info = choose_cluster_reps(tree, groups)
        for tag in groups:
            sizes = [n_assigned[p] for p in pieces[tag]]
            assert info[tag]["rep_piece_size"] == max(sizes)
            assert reps[tag] in groups[tag]
        reps_f, _ = choose_cluster_reps(tree, groups, "first")
        assert all(reps_f[t] == groups[t][0] for t in groups)


if __name__ == "__main__":
    n_fail = 0
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"PASS {name}")
            except Exception as e:  # noqa: BLE001
                n_fail += 1
                print(f"FAIL {name}: {type(e).__name__}: {e}")
    sys.exit(1 if n_fail else 0)
