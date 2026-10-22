#!/usr/bin/env python3
"""Unit tests for scripts/trait_null_phylo.py (fold preference vs random
sequence traits) and the 4-state ancestral gain/loss test
(utils/ancestral_tests.py, scripts/run_ancestral_analysis.py).
Synthetic trees / toy a3m only; no ete3, no data.

Run from the repo root:
    python3 -m pytest tests/test_trait_null_ancestral.py -q
    python3 tests/test_trait_null_ancestral.py
"""

from __future__ import annotations

import ast
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

from utils import phylo_tests as pt  # noqa: E402
from utils import ancestral_tests as at  # noqa: E402
import power_simulation_phylo as ps  # noqa: E402
import negative_traits_phylo as nt  # noqa: E402
import trait_null_phylo as tn  # noqa: E402

_Node = ps._Node


def _bm_traits(td, K, h, rng):
    """K traits with BM variance fraction h (phylogenetically structured)."""
    B = td["A"] @ rng.standard_normal((td["A"].shape[1], K))      # (n, K)
    return np.sqrt(h) * B + np.sqrt(1 - h) * rng.standard_normal((td["n"], K))


# ---------------------------------------------------------------------------
# Task A: trait null
# ---------------------------------------------------------------------------

def test_fast_stats_match_registry():
    tree = ps.random_tree(40, seed=3)
    names = [l.name for l in tree.iter_leaves()]
    pre = tn.tree_precompute(tree, names)
    rng = np.random.default_rng(1)
    for _ in range(3):
        x = rng.normal(size=40)
        ref = pt.morans_I_test(tree, dict(zip(names, x)), n_perm=4000,
                               rng=np.random.default_rng(2))
        I, I_norm, z, T = tn.trait_stats(pre, x)
        assert round(I, 4) == ref["morans_I"]
        assert round(I_norm, 4) == ref["morans_I_norm"]
        assert abs(z - ref["morans_z"]) < 0.1, (z, ref["morans_z"])   # MC error only
        cs = pt.run_clade_shift(tree, dict(zip(names, x)), 9, np.random.default_rng(0))
        assert round(T, 4) == cs["clade_T"]


def test_planted_clade_beats_random_traits_and_null_uniform():
    tree = ps.random_tree(60, seed=1)
    td = ps.prepare_tree(tree, "t")
    names = td["names"]
    rng = np.random.default_rng(7)
    K = 99
    X = _bm_traits(td, K, 0.3, rng)                   # structured random traits
    rand = {nm: X[j] for j, nm in enumerate(names)}
    clade = ps.pick_clade(td, 0.25)
    fold = {"planted": dict(zip(names, 3.0 * clade + 0.3 * rng.normal(size=td["n"])))}
    # fold vectors drawn from the SAME distribution as the random traits
    Y = _bm_traits(td, 300, 0.3, rng)
    for k in range(Y.shape[1]):
        fold[f"null{k}"] = dict(zip(names, Y[:, k]))
    rows = {r["method"]: r for r in tn.trait_null_rows(tree, rand, fold)}
    pl = rows["planted"]
    assert pl["K"] == K and pl["n_clusters"] == 60
    assert pl["p_trait_inorm"] <= 0.05 and pl["p_trait_T"] <= 0.05, pl
    for s in ("inorm", "z", "T"):
        p = np.array([rows[f"null{k}"][f"p_trait_{s}"] for k in range(Y.shape[1])])
        assert 0.40 <= p.mean() <= 0.60, (s, p.mean())
        assert (p <= 0.05).mean() <= 0.12, (s, (p <= 0.05).mean())
        assert (p >= 0.5).mean() >= 0.35


def test_missing_values_and_small_trees():
    tree = ps.random_tree(15, seed=2)
    names = [l.name for l in tree.iter_leaves()]
    rng = np.random.default_rng(0)
    rand = {nm: rng.normal(size=20) for nm in names}
    rand[names[0]] = np.full(20, np.nan)              # cluster without sequences
    vals = dict(zip(names, rng.normal(size=15)))
    vals[names[1]] = None
    rows = tn.trait_null_rows(tree, rand, {"AF2": vals, "ESM": {names[2]: 1.0}},
                              min_clusters=10)
    r = {x["method"]: x for x in rows}
    assert r["AF2"]["n_clusters"] == 13 and np.isfinite(r["AF2"]["p_trait_inorm"])
    assert r["ESM"]["note"] == "n<10" and np.isnan(r["ESM"]["p_trait_inorm"])
    assert 1 / 21 <= r["AF2"]["p_trait_inorm"] <= 1


def test_family_table_bonferroni_bh():
    rows = []
    for pair, ps_ in (("A", [0.01, 0.2, 0.5]), ("B", [0.3, 0.6, np.nan]),
                      ("C", [0.002, 0.9, 0.9])):
        for m, p in zip(("AF2", "AF3", "ESM"), ps_):
            rows.append({"pair_id": pair, "labels": "raw", "method": m, "n_clusters": 20,
                         "I_norm": 0.1, "z": 1.0, "T": 2.0, "p_trait_inorm": p,
                         "p_trait_z": p, "p_trait_T": p})
    fam = tn.family_table(pd.DataFrame(rows)).set_index("pair_id")
    assert abs(fam.loc["A", "p_fam_inorm"] - 0.03) < 1e-12
    assert abs(fam.loc["B", "p_fam_inorm"] - 0.6) < 1e-12        # 2 methods with p
    assert abs(fam.loc["C", "p_fam_inorm"] - 0.006) < 1e-12
    assert fam.loc["A", "best_method_inorm"] == "AF2"
    q = at.bh_adjust([0.03, 0.6, 0.006])
    assert np.allclose(fam.loc[["A", "B", "C"], "q_inorm"].values, q)


def test_random_proj_vectorised_matches_loop():
    seqs = ["ACDEFGHIK", "ACDEFGHIK", "AC-EFaaGHIR", "KK--FGHIX", "----"]
    tabs = [nt.random_scores(k, 777) for k in (1, 2, 3)]
    R = nt.random_score_matrix(3, 777)
    v = nt.random_proj_values(seqs[1:], R)
    for k, tab in enumerate(tabs):
        loop = np.nanmean([nt._mean_score("".join(c for c in s if c.isupper()), tab)
                           for s in seqs[1:4]])
        assert abs(v[k] - loop) < 1e-12
    t = nt.cluster_traits(seqs[:4], None, tabs)
    assert abs(t["random_proj_2"] - v[1]) < 1e-12
    assert np.isnan(nt.random_proj_values(["--", "ab"], R)).all()


# ---------------------------------------------------------------------------
# Task B: Fitch events, discretisation, tests
# ---------------------------------------------------------------------------

def _tiny(states):
    """((a,b),(c,(d,e))) with leaf names a..e."""
    root = _Node("R", 0.0)
    X = root.add_child(_Node("X"))
    Y = root.add_child(_Node("Y"))
    X.add_child(_Node("a"))
    X.add_child(_Node("b"))
    Y.add_child(_Node("c"))
    Z = Y.add_child(_Node("Z"))
    Z.add_child(_Node("d"))
    Z.add_child(_Node("e"))
    return root, dict(zip("abcde", states))


def test_fitch_event_typing_tiny_tree():
    # loss: root resolves to "both" (tie-break), F1 / F2 lineages lose it
    tree, ls = _tiny(["F1", "F1", "both", "F2", "F2"])
    ti = at.index_tree(tree)
    st = dict(zip(ti["names"], at.fitch_reconstruct(ti, ls)))
    assert st["R"] == "both" and st["X"] == "F1" and st["Y"] == "both" and st["Z"] == "F2"
    ev = at.fitch_events(ti, [st[n] for n in ti["names"]])
    assert sorted((e["parent"], e["child"], e["event"]) for e in ev) == [
        ("R", "X", "loss"), ("Y", "Z", "loss")]
    # gain: a single derived "both" clade
    tree, ls = _tiny(["F1", "F1", "F1", "both", "both"])
    ti = at.index_tree(tree)
    ev = at.fitch_events(ti, at.fitch_reconstruct(ti, ls))
    assert [(e["parent"], e["child"], e["event"], e["n_leaves_below"]) for e in ev] == [
        ("Y", "Z", "gain", 2)]
    # swap + other
    tree, ls = _tiny(["F1", "F1", "none", "F2", "F2"])
    ti = at.index_tree(tree)
    summ, ev, states = at.ancestral_core(tree, ls, n_perm=50)
    assert summ["n_transitions"] == 2 and summ["n_swap"] == 1 and summ["n_other"] == 1
    assert summ["testable"] and summ["n_F1"] == 2 and summ["n_none"] == 1
    assert at.event_type("none", "both") == "gain" and at.event_type("both", "none") == "loss"
    assert at.event_type("F2", "F1") == "swap" and at.event_type("F1", "F1") is None
    # untestable: one informative state only
    tree, ls = _tiny(["F1", "F1", "none", "none", "F1"])
    summ, _, _ = at.ancestral_core(tree, ls, n_perm=50)
    assert not summ["testable"] and np.isnan(summ["p_value"])


def _load_reference_fitch():
    """Exec ancestral_utils' _fitch_bottom_up / _fitch_top_down without
    importing the module (it needs ete3)."""
    src = open(os.path.join(ROOT, "utils", "ancestral_utils.py"), encoding="utf-8").read()
    mod = ast.parse(src)
    fns = [n for n in mod.body if isinstance(n, ast.FunctionDef)
           and n.name in ("_fitch_bottom_up", "_fitch_top_down")]
    ns = {"_FITCH_PREF": list(at.FITCH_PREF)}
    exec(compile(ast.Module(body=fns, type_ignores=[]), "ref", "exec"), ns)
    return ns["_fitch_bottom_up"], ns["_fitch_top_down"]


def test_fast_fitch_matches_reference_implementation():
    bottom_up, top_down = _load_reference_fitch()
    rng = np.random.default_rng(0)
    for seed in range(25):
        tree = ps.random_tree(int(rng.integers(5, 40)), seed=seed)
        leaves = list(tree.iter_leaves())
        ls = {l.name: str(rng.choice(["F1", "F2", "both", "none"])) for l in leaves}
        for l in leaves:
            l._state = ls[l.name]
        bottom_up(tree)
        top_down(tree)
        ref = [nd._state for nd in tree.traverse("preorder")]
        ti = at.index_tree(tree)
        assert at.fitch_reconstruct(ti, ls) == ref
        ref_n = sum(1 for nd in tree.traverse("preorder")
                    if not nd.is_root() and nd._state != nd.up._state)
        assert at.count_changes(ti, ref) == ref_n


def test_discretise_preserves_frequencies_exactly():
    rng = np.random.default_rng(3)
    for trial in range(50):
        n = int(rng.integers(5, 80))
        labels = list(rng.choice(["F1", "F2", "both", "none"], size=n,
                                 p=rng.dirichlet(np.ones(4))))
        v = rng.normal(size=n)
        v[rng.random(n) < 0.1] = np.nan                  # missing values
        v[: n // 4] = np.round(v[: n // 4], 0)           # ties
        out = at.discretise_like(v, labels, np.random.default_rng(trial))
        assert sorted(out) == sorted(labels)
    # rank order: most frequent state gets the largest values
    labels = ["none"] * 5 + ["F1"] * 3 + ["both"] * 2
    out = at.discretise_like(np.arange(10.0), labels)
    assert out == ["both"] * 2 + ["F1"] * 3 + ["none"] * 5
    # frequency tie broken by FITCH_PREF order (both before F1)
    out = at.discretise_like([3.0, 2.0, 1.0, 0.0], ["F1", "F1", "both", "both"])
    assert out == ["both", "both", "F1", "F1"]


def test_permutation_and_trait_control():
    tree = ps.random_tree(40, seed=5)
    td = ps.prepare_tree(tree, "t")
    clade = ps.pick_clade(td, 0.3)
    ls = {nm: ("both" if c else "F1") for nm, c in zip(td["names"], clade)}
    rng = np.random.default_rng(0)
    iid = {nm: rng.normal(size=30) for nm in td["names"]}
    summ, ev, _ = at.ancestral_core(tree, ls, iid, n_perm=200, seed=1)
    assert summ["n_transitions"] == 1 and summ["n_gain"] + summ["n_loss"] == 1
    assert summ["p_value"] <= 0.01 and summ["p_trait_anc"] <= 1 / 31 + 1e-12
    assert summ["K_trait"] == 30 and summ["n_trait_leaves"] == 40
    # fold labels that ARE a discretised random trait -> p_trait_anc not small
    ti = at.index_tree(tree)
    names = ti["leaf_names"]
    X = _bm_traits(td, 31, 0.5, rng)
    lab_map = dict(zip(td["names"], ["F1"] * 25 + ["F2"] * 15))
    labs = at.discretise_like(X[:, 0], [lab_map[nm] for nm in td["names"]],
                              np.random.default_rng(1))
    obs_ls = dict(zip(td["names"], labs))
    obs = at.leaf_transitions(ti, names, [obs_ls[nm] for nm in names])
    traits = X[:, 1:].T                                       # (30, n) in td order
    order = [td["names"].index(nm) for nm in names]
    p, tr = at.trait_transition_null(ti, obs_ls, traits[:, order], obs, seed=2)
    assert len(tr) == 30 and p > 0.05, (p, obs, tr)


if __name__ == "__main__":
    fails = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"PASS {name}")
            except Exception as e:  # noqa: BLE001
                fails += 1
                import traceback
                traceback.print_exc()
                print(f"FAIL {name}: {type(e).__name__}: {e}")
    sys.exit(1 if fails else 0)
