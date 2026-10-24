#!/usr/bin/env python3
"""Unit tests for the calibration and matched-dispersion analyses:
scripts/cluster_composition.py (fold-free per-cluster features),
scripts/concordance_fp_calibration.py (concordance false-positive floor),
scripts/tm_subsample_matched.py (cluster-count-matched TM dispersion) and
scripts/ancestor_divergence_baseline.py (XCL1 ancestor baselines).
Synthetic inputs only; no Pipeline data, no cluster, no GPU.

Run from the repo root:
    python3 -m pytest tests/test_fp_calibration.py -q
    python3 tests/test_fp_calibration.py
"""

from __future__ import annotations

import os
import sys
import tempfile

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

import cluster_composition as cc            # noqa: E402
import concordance_fp_calibration as fp     # noqa: E402
import tm_subsample_matched as ts           # noqa: E402
import ancestor_divergence_baseline as ab   # noqa: E402


# ---------------------------------------------------------------------------
# cluster_composition: fold-free per-cluster features
# ---------------------------------------------------------------------------

def test_read_a3m_skips_hash_and_joins_wrapped():
    """DeepMsa.a3m carries a '#<len>\t<n>' line that is not a sequence, and
    rows may wrap over several lines."""
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "x.a3m")
        open(p, "w").write(">q\n#334\t1\nACDE\nFGHI\n>m1\nAC-E\n")
        seqs = cc.read_a3m(p)
    assert seqs == ["ACDEFGHI", "AC-E"], seqs


def test_features_frequencies_and_gap_fraction():
    """Frequencies are over match-state residues and sum to 1; '-' counts as a
    match-column gap; lower case is an insertion and is ignored entirely."""
    d = cc.features(["AAAA", "AA--"])
    freq = np.array([d[f"f_{a}"] for a in cc.AA20])
    assert abs(freq.sum() - 1.0) < 1e-12
    assert abs(d["f_A"] - 1.0) < 1e-12                    # only alanine present
    # member 1 has no gaps, member 2 is half gaps -> mean 0.25
    assert abs(d["gap_frac"] - 0.25) < 1e-12
    assert abs(d["log_size"] - np.log(2)) < 1e-12
    assert d["n_members"] == 2

    ins = cc.features(["AAAA", "AAAAwww"])                 # insertions ignored
    assert abs(ins["f_A"] - 1.0) < 1e-12
    assert abs(ins["gap_frac"]) < 1e-12


def test_features_charge_and_hydrophobicity():
    d = cc.features(["KKDD"])                              # +2 -2 = 0
    assert abs(d["net_charge"]) < 1e-12
    d2 = cc.features(["KKKK"])
    assert abs(d2["net_charge"] - 4.0) < 1e-12
    # all-isoleucine is the Kyte-Doolittle maximum, 4.5
    assert abs(cc.features(["IIII"])["hydrophobicity"] - 4.5) < 1e-12


def test_features_returns_none_without_residues():
    assert cc.features(["----"]) is None


def test_write_table_roundtrip():
    """A two-cluster toy pipeline tree produces one row per cluster, with row 0
    of each a3m dropped as the query."""
    import csv
    with tempfile.TemporaryDirectory() as d:
        cdir = os.path.join(d, "1abcA_2defB", "output_msa_cluster")
        os.makedirs(cdir)
        open(os.path.join(cdir, "ShallowMsa_000.a3m"), "w").write(
            ">q\nAAAA\n>m\nAAAA\n")
        open(os.path.join(cdir, "ShallowMsa_001.a3m"), "w").write(
            ">q\nKKKK\n>m\nKKKK\n>m2\nKKKK\n")
        open(os.path.join(cdir, "DeepMsa.a3m"), "w").write(">q\nAAAA\n")
        out = os.path.join(d, "comp.csv")
        n_ok, n_skip = cc.write_table(out, d)
        rows = list(csv.DictReader(open(out)))
    assert (n_ok, n_skip) == (2, 0), (n_ok, n_skip)        # DeepMsa not matched
    assert [r["cluster"] for r in rows] == ["ShallowMsa_000", "ShallowMsa_001"]
    assert [int(r["n_members"]) for r in rows] == [1, 2]
    assert abs(float(rows[1]["f_K"]) - 1.0) < 1e-12


# ---------------------------------------------------------------------------
# concordance_fp_calibration: the correction must match the real pipeline
# ---------------------------------------------------------------------------

def _ref_regress(y, x1, x2):
    """seq_divergence_correction.py's branch logic, transcribed."""
    ok = np.isfinite(y) & np.isfinite(x1)
    if ok.sum() < 3:
        return y
    yv, x1v, x2v = y[ok], x1[ok], x2[ok]
    has_F2 = not np.isnan(x2v).all()
    if has_F2 and np.std(x2v) >= 1e-9 and np.std(x1v) >= 1e-9:
        X = np.column_stack([np.ones(len(yv)), x1v, x2v])
        a, b1, b2 = np.linalg.lstsq(X, yv, rcond=None)[0]
        return y - (a + b1 * x1 + b2 * x2)
    if np.std(x1v) < 1e-9:
        return y
    b1, a = np.polyfit(x1v, yv, 1)
    return y - (b1 * x1 + a)


def test_regress_out_matches_reference_two_predictor():
    rng = np.random.default_rng(0)
    for _ in range(20):
        n = rng.integers(3, 30)
        x1, x2 = rng.random(n), rng.random(n)
        y = 2.0 * x1 - 1.5 * x2 + 0.3 + rng.standard_normal(n) * 0.1
        got, want = fp.regress_out(y, x1, x2), _ref_regress(y, x1, x2)
        assert np.allclose(got, want, atol=1e-10)


def test_regress_out_removes_the_divergence_gradient():
    """A response that is exactly linear in the two identities leaves no
    residual, which is the point of the correction."""
    rng = np.random.default_rng(1)
    x1, x2 = rng.random(40), rng.random(40)
    y = 3.1 * x1 - 0.7 * x2 + 2.0
    assert np.allclose(fp.regress_out(y, x1, x2), 0.0, atol=1e-8)


def test_regress_out_falls_back_when_second_reference_missing():
    """All-NaN x2 must use the single-predictor branch, not propagate NaN."""
    rng = np.random.default_rng(2)
    x1 = rng.random(20)
    x2 = np.full(20, np.nan)
    y = 1.3 * x1 + rng.standard_normal(20) * 0.05
    r = fp.regress_out(y, x1, x2)
    assert np.isfinite(r).all()
    assert np.allclose(r, _ref_regress(y, x1, x2), atol=1e-10)


def test_regress_out_partial_nan_second_reference_is_finite():
    """A partially missing x2 makes np.std NaN, so the single-predictor branch
    must take over. Using np.nanstd here would feed NaN into lstsq."""
    rng = np.random.default_rng(3)
    x1 = rng.random(15)
    x2 = rng.random(15)
    x2[::4] = np.nan
    y = 0.9 * x1 + rng.standard_normal(15) * 0.05
    r = fp.regress_out(y, x1, x2)
    assert np.isfinite(r).all(), r


def test_regress_out_noop_when_underdetermined_or_degenerate():
    y = np.array([1.0, 2.0])
    assert np.allclose(fp.regress_out(y, np.array([0.1, 0.2]),
                                      np.array([np.nan, np.nan])), y)
    y3 = np.array([1.0, 2.0, 3.0])
    const = np.full(3, 0.5)
    assert np.allclose(fp.regress_out(y3, const, np.full(3, np.nan)), y3)


# ---------------------------------------------------------------------------
# The pseudo-methods must be fold-free and land on chance
# ---------------------------------------------------------------------------

def test_pseudo_methods_are_centred_on_chance():
    """Seven random projections of one composition matrix share sequence
    structure but no fold, so their mean pairwise sign concordance must sit at
    0.5 rather than above it."""
    rng = np.random.default_rng(7)
    vals = []
    for _ in range(300):
        Z = rng.standard_normal((40, 24))
        Z = (Z - Z.mean(axis=0)) / Z.std(axis=0)
        cols = []
        for _m in range(7):
            v = Z @ rng.standard_normal(24)
            cols.append(v - np.median(v))
        signs = np.sign(np.column_stack(cols))
        a, _n = fp._mean_concordance_from_residuals(signs)
        vals.append(a)
    assert abs(float(np.mean(vals)) - 0.5) < 0.01, np.mean(vals)


def test_mean_concordance_matches_hand_count():
    """Two methods agreeing on 3 of 4 clusters give 0.75; a zero sign marks a
    cluster the method does not cover and is excluded pairwise."""
    signs = np.array([[1, 1], [1, 1], [1, -1], [-1, -1]], dtype=float)
    a, n = fp._mean_concordance_from_residuals(signs)
    assert abs(a - 0.75) < 1e-12 and n == 4
    signs[0, 1] = 0.0                                   # method 2 misses row 0
    a2, n2 = fp._mean_concordance_from_residuals(signs)
    assert abs(a2 - 2.0 / 3.0) < 1e-12 and n2 == 3


def test_bh_adjust_is_monotone_and_matches_formula():
    p = np.array([0.001, 0.008, 0.039, 0.041, 0.9])
    adj = fp.bh_adjust(p)
    assert np.all(np.diff(adj) >= -1e-12)               # monotone in p
    assert abs(adj[0] - 0.005) < 1e-12                  # 0.001 * 5 / 1
    assert adj.max() <= 1.0
    assert np.isnan(fp.bh_adjust(np.array([np.nan]))[0])


# ---------------------------------------------------------------------------
# tm_subsample_matched: matched cluster budgets
# ---------------------------------------------------------------------------

def test_target_stats_none_when_too_few_clusters():
    mat = np.full((3, 3), 0.8)
    np.fill_diagonal(mat, np.nan)
    assert ts.target_stats(mat, 5, 10, np.random.default_rng(0)) is None
    assert ts.target_stats(mat, 3, 10, np.random.default_rng(0)) is not None


def test_target_stats_exact_on_a_known_matrix():
    """With k equal to n there is only one subset, so the statistics are the
    plain min / mean / p10 over the off-diagonal entries."""
    mat = np.array([[np.nan, 0.2, 0.4],
                    [0.2, np.nan, 0.6],
                    [0.4, 0.6, np.nan]])
    st = ts.target_stats(mat, 3, 5, np.random.default_rng(0))
    vals = np.array([0.2, 0.4, 0.6])
    assert abs(st["min_tm"] - 0.2) < 1e-12
    assert abs(st["mean_tm"] - vals.mean()) < 1e-12
    assert abs(st["p10_tm"] - np.quantile(vals, 0.1)) < 1e-12


def test_subsampling_removes_the_extreme_value_advantage():
    """Two targets with the same TM distribution but different cluster counts
    differ in their raw minimum and agree once the budget is matched."""
    rng = np.random.default_rng(5)
    def mat_of(n):
        m = rng.uniform(0.5, 0.9, size=(n, n))
        m = (m + m.T) / 2
        np.fill_diagonal(m, np.nan)
        return m
    big, small = mat_of(60), mat_of(8)
    iu_b, iu_s = np.triu_indices(60, 1), np.triu_indices(8, 1)
    assert np.nanmin(big[iu_b]) < np.nanmin(small[iu_s])     # pure arithmetic
    sb = ts.target_stats(big, 8, 400, np.random.default_rng(1))
    ss = ts.target_stats(small, 8, 400, np.random.default_rng(1))
    assert abs(sb["min_tm"] - ss["min_tm"]) < 0.05, (sb, ss)


# ---------------------------------------------------------------------------
# ancestor_divergence_baseline
# ---------------------------------------------------------------------------

def test_rho_matches_scipy_including_ties():
    try:
        from scipy.stats import spearmanr
    except Exception:                                        # pragma: no cover
        return
    rng = np.random.default_rng(11)
    for _ in range(30):
        x = rng.integers(0, 4, 8).astype(float)              # ties on purpose
        y = rng.integers(0, 4, 8).astype(float)
        if np.std(x) < 1e-12 or np.std(y) < 1e-12:
            continue
        assert abs(ab.rho(x, y) - spearmanr(x, y).statistic) < 1e-10


def test_exact_p_is_over_all_120_orderings():
    pct = np.array([0.0, 0.0, 8.0, 90.7, 50.0])
    r, p, n = ab.exact_p(np.arange(5.0), pct)
    assert n == 120
    assert abs(p * 120 - round(p * 120)) < 1e-9              # a multiple of 1/120
    assert abs(r - 0.872081599272381) < 1e-9                 # the lineage ceiling
    assert abs(p - 0.05) < 1e-12


def test_lineage_ceiling_bounds_any_monotone_score():
    """No score increasing along the lineage can beat the lineage order, which
    is why reaching rho = 0.87 is not evidence of reading fold."""
    pct = np.array([0.0, 0.0, 8.0, 90.7, 50.0])
    ceiling = ab.rho(np.arange(5.0), pct)
    rng = np.random.default_rng(13)
    for _ in range(200):
        v = np.sort(rng.random(5))                           # strictly increasing
        assert ab.rho(v, pct) <= ceiling + 1e-12


def test_gapaware_pid_ignores_shared_gaps():
    assert abs(ab.gapaware_pid("AC-E", "AC-E") - 1.0) < 1e-12
    assert abs(ab.gapaware_pid("ACDE", "ACDQ") - 0.75) < 1e-12
    assert abs(ab.gapaware_pid("A--E", "AC-E") - 2.0 / 3.0) < 1e-12


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
