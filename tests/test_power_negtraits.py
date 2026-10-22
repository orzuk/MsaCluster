#!/usr/bin/env python3
"""Unit tests for scripts/power_simulation_phylo.py and
scripts/negative_traits_phylo.py (synthetic trees / toy a3m; no ete3, no data).

Run from the repo root:
    python3 -m pytest tests/test_power_negtraits.py -q
    python3 tests/test_power_negtraits.py
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

from utils import phylo_tests as pt  # noqa: E402
import power_simulation_phylo as ps  # noqa: E402
import negative_traits_phylo as nt  # noqa: E402


def _td(n, seed=0):
    return ps.prepare_tree(ps.random_tree(n, seed=seed), f"t{n}_{seed}")


# ---------------------------------------------------------------------------
# Power simulation core
# ---------------------------------------------------------------------------

def test_fast_tests_match_registry():
    """morans_fast / clade_fast give the registry's numbers (same rng)."""
    for seed in range(3):
        tree = ps.random_tree(30, seed=seed)
        td = ps.prepare_tree(tree, "x")
        x = np.random.default_rng(seed).normal(size=td["n"])
        vals = dict(zip(td["names"], x))
        ref = pt.morans_I_test(tree, vals, n_perm=199, rng=np.random.default_rng(5))
        I, p = ps.morans_fast(td["W"], x, 199, np.random.default_rng(5))
        assert abs(round(I, 4) - ref["morans_I"]) < 1e-12
        assert abs(round(p, 4) - ref["morans_p"]) < 1e-12
        ref = pt.run_clade_shift(tree, vals, 499, np.random.default_rng(9))
        T, p, s = ps.clade_fast(td["M"], x, 499, np.random.default_rng(9))
        assert abs(round(T, 4) - ref["clade_T"]) < 1e-12
        assert abs(round(p, 6) - ref["clade_p"]) < 1e-12
        assert int(td["M"][s].sum()) == ref["clade_best_size"]


def test_planted_large_shift_detected():
    td = _td(60, seed=1)
    sc = {"scenario": "clade_big", "kind": "clade", "f": 0.25, "delta": 3.0,
          "h": np.nan, "r": 1.0}
    rows = ps.run_cell(td, sc, n_reps=20, seed=1, morans_perm=199, clade_perm=499)
    pw = {r["test"]: r["power"] for r in rows}
    assert pw["morans"] >= 0.9 and pw["clade_shift"] >= 0.9, pw
    hit = [r["clade_hit_rate"] for r in rows if r["test"] == "clade_shift"][0]
    assert hit >= 0.8, hit
    # heavy measurement noise must cost power
    sc_noisy = dict(sc, scenario="clade_big_noisy", r=0.1)
    rows_n = ps.run_cell(td, sc_noisy, n_reps=20, seed=1, morans_perm=199, clade_perm=499)
    assert all(rn["power"] < r["power"] for rn, r in zip(rows_n, rows))


def test_null_pvalues_roughly_uniform():
    try:
        from scipy.stats import kstest
    except ImportError:
        kstest = None
    td = _td(40, seed=2)
    rng = np.random.default_rng(3)
    null = {"scenario": "null", "kind": "null", "f": np.nan, "delta": 0.0,
            "h": 0.0, "r": 1.0}
    pm, pc = [], []
    for _ in range(300):
        x = ps.simulate(td, null, rng)
        pm.append(ps.morans_fast(td["W"], x, 99, rng)[1])
        pc.append(ps.clade_fast(td["M"], x, 99, rng)[1])
    for p in (np.array(pm), np.array(pc)):
        assert 0.02 <= (p < 0.05).mean() <= 0.09, (p < 0.05).mean()
        assert 0.40 <= np.median(p) <= 0.60, np.median(p)
        if kstest is not None:
            assert kstest(p, "uniform").pvalue > 0.001


def test_bm_scaling_and_reliability():
    td = _td(50, seed=4)
    assert abs((td["A"] ** 2).sum(axis=1).mean() - 1.0) < 1e-9
    rng = np.random.default_rng(0)
    clade = ps.pick_clade(td, 0.25)
    assert abs(clade.mean() - 0.25) <= 0.15
    sc = {"scenario": "c", "kind": "clade", "f": 0.25, "delta": 2.0, "h": np.nan, "r": 0.4}
    p = clade.mean()
    v_true = 1 + 4 * p * (1 - p)
    obs = np.array([ps.simulate(td, sc, rng, clade) for _ in range(4000)])
    # per-leaf variance across replicates: in-clade/out-clade leaves both
    # have var = 1/r of the e-variance (+ measurement); mean total var check:
    resid_var = obs.var(axis=0).mean()            # = 1 + s_m^2
    expect = 1 + v_true * (1 - 0.4) / 0.4
    assert abs(resid_var / expect - 1) < 0.05, (resid_var, expect)


def test_scenarios_and_aggregate():
    sc = ps.build_scenarios()
    assert len(sc) == 1 + 3 * (4 * 5 + 3)
    assert len({s["scenario"] for s in sc}) == len(sc)
    td = _td(25, seed=5)
    rows = []
    for s in sc[:3]:
        rows += ps.run_cell(td, s, 3, 0, 19, 19)
    agg = ps.aggregate(pd.DataFrame(rows))
    assert set(agg.size_bin) == {"20-49", "ALL"}
    assert (agg.power.between(0, 1)).all()


# ---------------------------------------------------------------------------
# Negative traits
# ---------------------------------------------------------------------------

TOY_A3M = """>query
ACDEFGHIK
>m1
ACDEFGHIK
>m2
AC-EFaaGHIR
>m3
KK--FGHIK
"""


def test_toy_a3m_traits(tmp_path=None):
    import tempfile
    d = tempfile.mkdtemp()
    f = os.path.join(d, "ShallowMsa_000.a3m")
    with open(f, "w") as fh:
        fh.write(TOY_A3M)
    seqs = nt.load_a3m_seqs(f)
    assert seqs[2] == "AC-EFaaGHIR"
    tab = {a: 1.0 for a in nt.AA20}
    t = nt.cluster_traits(seqs, None, [tab])
    # identity to row 0 over query match columns (lowercase in member skipped)
    pids = [1.0, 7 / 8, 5 / 7]
    assert abs(t["pid_to_query"] - np.mean(pids)) < 1e-12
    kd = nt.KYTE_DOOLITTLE
    hyd = [np.mean([kd[c] for c in "ACDEFGHIK"]),
           np.mean([kd[c] for c in "ACEFGHIR"]),
           np.mean([kd[c] for c in "KKFGHIK"])]
    assert abs(t["hydrophobicity"] - np.mean(hyd)) < 1e-12
    assert abs(t["gap_frac"] - np.mean([0, 1 / 9, 2 / 9])) < 1e-12
    assert abs(t["log_size"] - np.log(3)) < 1e-12
    # charges: m1 D,E,K -> -1; m2 E,R -> 0; m3 K,K,K -> +3
    assert abs(t["net_charge"] - np.mean([-1, 0, 3])) < 1e-12
    assert abs(t["random_proj_1"] - 1.0) < 1e-12
    # explicit query overrides row 0
    t2 = nt.cluster_traits(seqs, "ACDEFGHIR", [tab])
    assert abs(t2["pid_to_query"] - np.mean([8 / 9, 1.0, 4 / 7])) < 1e-12
    assert nt.cluster_traits(seqs[:1], None, [tab]) is None


def test_random_scores_fixed():
    a, b = nt.random_scores(1, 777), nt.random_scores(1, 777)
    assert a == b and a != nt.random_scores(2, 777)


def test_trait_morans_and_summary():
    tree = ps.random_tree(30, seed=6)
    td = ps.prepare_tree(tree, "p")
    rng = np.random.default_rng(0)
    clade = ps.pick_clade(td, 0.4)
    lv = {nm: {"hydrophobicity": float(rng.normal() + 4 * c),
               "random_proj_1": float(rng.normal())}
          for nm, c in zip(td["names"], clade)}
    rows = nt.morans_rows("p", tree, lv, ["hydrophobicity", "random_proj_1"], 1000, 12345)
    r = {x["trait"]: x for x in rows}
    assert r["hydrophobicity"]["morans_p"] < 0.01
    assert r["hydrophobicity"]["morans_I"] > r["random_proj_1"]["morans_I"]
    tr = pd.DataFrame(rows)
    fp = pd.DataFrame({"pair_id": ["p", "p"], "method": ["AF2", "ESM"],
                       "morans_I": [0.01, 0.2], "morans_I_norm": [0.1, 0.5],
                       "morans_p": [0.5, 0.01], "label_set": "corrected"})
    summ, fam = nt.summarise(tr, fp)
    assert set(summ.group) == {"hydrophobicity", "random_proj_1",
                               "FP_corrected:AF2", "FP_corrected:ESM"}
    assert fam.loc[0, "best_fp_I_corrected"] == 0.2
    assert fam.loc[0, "best_fp_method_corrected"] == "ESM"


if __name__ == "__main__":
    fails = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"PASS {name}")
            except Exception as e:  # noqa: BLE001
                fails += 1
                print(f"FAIL {name}: {type(e).__name__}: {e}")
    sys.exit(1 if fails else 0)
