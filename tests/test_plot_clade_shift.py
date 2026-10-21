#!/usr/bin/env python3
"""Smoke test for scripts/plot_xcl1_clade_shift.py on synthetic inputs.

Run (from the repo root):
    python3 -m pytest tests/test_plot_clade_shift.py -q
    python3 tests/test_plot_clade_shift.py [OUT_DIR]   # keeps the PNG/PDF

Builds a random 100-leaf coarse tree (ete3-like stand-in, as in
tests/test_phylo_tests.py) with a planted 4-leaf clade, a synthetic survey CSV
(7 methods, the clade shifted toward F1), a synthetic clade-shift CSV, and a
toy a3m directory for the % identity column. ete3 is stubbed if missing
(utils.phytree_utils imports it at module level; only its matplotlib
drawing helper is used here).
"""

from __future__ import annotations

import os
import sys
import tempfile
import types

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    import ete3  # noqa: F401
except ImportError:
    _stub = types.ModuleType("ete3")
    for _n in ("Tree", "TreeStyle", "NodeStyle", "faces", "RectFace", "TextFace"):
        setattr(_stub, _n, object)
    sys.modules["ete3"] = _stub

from scripts import plot_xcl1_clade_shift as px  # noqa: E402
from tests.test_phylo_tests import MiniNode, random_tree  # noqa: E402

HL = ["ShallowMsa_005", "ShallowMsa_037", "ShallowMsa_038", "ShallowMsa_075"]


class Node(MiniNode):
    def get_tree_root(self):
        n = self
        while n.up is not None:
            n = n.up
        return n

    def get_leaves(self):
        return list(self.iter_leaves())


def _as_node(m):
    n = Node(m.name, m.dist)
    for c in m.children:
        n.add_child(_as_node(c))
    return n


def synthetic_tree(n=100, seed=3):
    """Random tree, then graft a 4-leaf clade carrying the HL labels."""
    t = _as_node(random_tree(n - 4, seed=seed))
    other = [f"S{i:02d}" for i in range(n) if f"ShallowMsa_{i:03d}" not in HL]
    for lf, name in zip(list(t.iter_leaves()), other):
        lf.name = name
    a, b = Node("", 0.05), Node("", 0.05)
    a.add_child(Node("S05", 0.1)); a.add_child(Node("S37", 0.1))
    b.add_child(Node("S38", 0.1)); b.add_child(Node("S75", 0.1))
    clade = Node("", 0.6)
    clade.add_child(a); clade.add_child(b)
    root = Node("", 0.0)
    old = t
    old.dist = 0.4
    root.add_child(old); root.add_child(clade)
    tag_to_label = {f"ShallowMsa_{i:03d}": f"S{i:02d}" for i in range(n)}
    return root, tag_to_label


def synthetic_inputs(tmp, n=100, seed=0):
    rs = np.random.RandomState(seed)
    scale = {"AF2": 0.13, "AF3": 0.15, "ESM": 0.035, "Boltz2": 0.13,
             "BioEmu": 0.07, "DDG": 1.5, "S4PRED": 0.09}
    rows = []
    for i in range(n):
        tag = f"ShallowMsa_{i:03d}"
        for m, sd in scale.items():
            v = rs.normal(0, sd)
            if tag in HL and m != "ESM":
                v = abs(v) + 2.5 * sd
            if m == "Boltz2" and i == 17:
                v = np.nan                      # a missing cell
            rows.append({"pair_id": "2n54B_2hdmA", "cluster": tag,
                         "method": m, "TMdiff_centered": v})
    survey = os.path.join(tmp, "survey.csv")
    pd.DataFrame(rows).to_csv(survey, index=False)
    members = ";".join(HL)
    crow = []
    for m in list(scale) + ["COMBINED"]:
        same = m not in ("ESM", "DDG")
        crow.append({"pair_id": "2n54B_2hdmA", "method": m,
                     "clade_p": 1e-4 if same else 0.21,
                     "clade_best_members": members if same
                     else "ShallowMsa_010;ShallowMsa_011",
                     "clade_best_direction": 1 if same else -1})
    clade = os.path.join(tmp, "clade.csv")
    pd.DataFrame(crow).to_csv(clade, index=False)
    return survey, clade


def toy_pair_dir(tmp):
    """DeepMsa query + 2 cluster a3m files for the PID path."""
    pd_ = os.path.join(tmp, "pair")
    os.makedirs(os.path.join(pd_, "output_get_msa"))
    os.makedirs(os.path.join(pd_, "output_msa_cluster"))
    with open(os.path.join(pd_, "output_get_msa", "DeepMsa.a3m"), "w") as f:
        f.write(">q\nACDEFGHIKL\n>x\nACDEFGHIKL\n")
    with open(os.path.join(pd_, "output_msa_cluster", "ShallowMsa_005.a3m"), "w") as f:
        f.write(">q\nACDEFGHIKL\n>a\nACDEFGHIKV\n>b\nACDEFGHI--\n")   # 90%, 100%
    with open(os.path.join(pd_, "output_msa_cluster", "ShallowMsa_001.a3m"), "w") as f:
        f.write(">q\nACDEFGHIKL\n>a\nACWWWWWWWW\n")                     # 2/10
    return pd_


def test_helpers():
    assert px.parse_members("ShallowMsa_005;ShallowMsa_037") == HL[:2]
    assert px.parse_members("S05, 37") == HL[:2]
    assert px.short_label("ShallowMsa_075") == "S75"
    assert px.resolve_highlight(None, {"AF2": {"clade_best_members": ";".join(HL)}}) == HL
    assert px._runs([0, 1, 1, 0, 1]) == [(1, 2), (4, 4)]
    assert abs(px.robust_range(np.r_[np.zeros(50), 1.0, np.nan], 100) - 1) < 1e-12
    with tempfile.TemporaryDirectory() as tmp:
        pid = px.cluster_pid_from_a3m(toy_pair_dir(tmp))
        assert abs(pid["ShallowMsa_005"] - 95.0) < 1e-9, pid
        assert abs(pid["ShallowMsa_001"] - 20.0) < 1e-9, pid


def run_render(out_dir):
    survey, clade = synthetic_inputs(out_dir)
    tree, t2l = synthetic_tree()
    rows = px.load_clade_rows(clade, "2n54B_2hdmA")
    hl = px.resolve_highlight(None, rows)
    vals = px.load_values(survey, "2n54B_2hdmA", px.FIG_METHODS)
    rs = np.random.RandomState(1)
    pid = {t: (62 + rs.normal(0, 4)) if t in HL else 30 + rs.normal(0, 6)
           for t in t2l}
    out = os.path.join(out_dir, "synthetic_xcl1_clade_shift")
    data = px.make_figure(tree, t2l, vals, hl, rows, pid, "2n54B_2hdmA", out,
                          mark_label="XCL",
                          fold_names=("F1 2n54B", "F2 2hdmA"))
    return out, data


def test_render():
    with tempfile.TemporaryDirectory() as tmp:
        out, data = run_render(tmp)
        for ext in (".png", ".pdf", "_data.csv"):
            assert os.path.getsize(out + ext) > 0
        assert int(data["in_clade"].sum()) == 4 and len(data) == 100
        assert sorted(data.index[data["in_clade"]]) == HL


if __name__ == "__main__":
    test_helpers()
    if len(sys.argv) > 1:
        os.makedirs(sys.argv[1], exist_ok=True)
        print(run_render(sys.argv[1])[0] + ".png")
    else:
        test_render()
    print("ok")
