#!/usr/bin/env python3
"""Unit tests for the fine-cluster contact provenance guard and top-L scoring
(utils/contact_utils.py, Analysis/cmap_analysis.py). Synthetic inputs only.

Run from the repo root:
    python3 -m pytest tests/test_contacts_fine.py -q
    python3 tests/test_contacts_fine.py
"""

from __future__ import annotations

import os
import shutil
import sys
import tempfile
import time
import types

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from utils import contact_utils as cu  # noqa: E402


def _import_cmap_analysis():
    """Import Analysis.cmap_analysis; stub utils.protein_utils if torch etc.
    are not installed (only load_cmap_and_idx is needed from it)."""
    try:
        import utils.protein_utils  # noqa: F401
    except Exception:
        stub = types.ModuleType("utils.protein_utils")

        def load_cmap_and_idx(path):
            z = np.load(path)
            return z["cmap"], z["idx"], (z["keep_cols"] if "keep_cols" in z.files else None)

        stub.load_cmap_and_idx = load_cmap_and_idx
        sys.modules["utils.protein_utils"] = stub
    from Analysis import cmap_analysis
    return cmap_analysis


def _write_a3m(path, recs):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        for h, s in recs:
            fh.write(f">{h}\n{s}\n")


# ---------------------------------------------------------------------------
# top-L binarisation
# ---------------------------------------------------------------------------
def test_top_l_symmetric_k_entries_and_sep():
    rng = np.random.default_rng(0)
    n, k, sep = 40, 25, 6
    P = rng.random((n, n))            # deliberately asymmetric
    B = cu.top_l_binarize(P, k, sep_min=sep)
    assert B.dtype == np.uint8
    assert (B == B.T).all()
    assert int(np.triu(B, 1).sum()) == k
    assert np.diag(B).sum() == 0
    i, j = np.nonzero(np.triu(B, 1))
    assert (np.abs(i - j) >= sep).all()
    # chosen pairs are the k best eligible ones of the symmetrised map
    S = 0.5 * (P + P.T)
    iu, ju = np.triu_indices(n, k=sep)
    thr = np.sort(S[iu, ju])[::-1][k - 1]
    assert S[i, j].min() >= thr - 1e-12


def test_top_l_short_range_never_selected_even_if_high():
    n = 20
    P = np.zeros((n, n))
    for d in range(1, 6):              # huge short-range scores
        for a in range(n - d):
            P[a, a + d] = 10.0
    P[0, 15] = 1.0
    P[2, 12] = 0.5
    B = cu.top_l_binarize(P, 2, sep_min=6)
    assert {tuple(x) for x in np.argwhere(np.triu(B, 1))} == {(0, 15), (2, 12)}


def test_top_l_uses_idx_for_separation_and_caps_k():
    P = np.ones((4, 4))
    idx = np.array([0, 1, 10, 11])    # seed columns: only cross pairs are >= 6 apart
    B = cu.top_l_binarize(P, 100, sep_min=6, idx=idx)
    assert int(np.triu(B, 1).sum()) == 4      # (0,2),(0,3),(1,2),(1,3)
    assert B[0, 1] == 0 and B[2, 3] == 0
    P2 = np.full((5, 5), np.nan)
    assert cu.top_l_binarize(P2, 3, sep_min=1).sum() == 0


# ---------------------------------------------------------------------------
# manifest hash logic
# ---------------------------------------------------------------------------
def test_manifest_hash_match_and_mismatch():
    tmp = tempfile.mkdtemp()
    try:
        pair = os.path.join(tmp, "1abcA_2defB")
        a3m = os.path.join(pair, "output_msa_cluster", "ShallowMsa_000.a3m")
        _write_a3m(a3m, [("S1 desc", "ACDEF"), ("u1", "AC-EF"), ("u2", "ACDEG")])
        h1, n = cu.a3m_cluster_hash(a3m)
        assert n == 3
        # row order / header annotations do not change the hash
        _write_a3m(a3m, [("u2 other", "ACDEG"), ("S1", "ACDEF"), ("u1", "AC-EF")])
        assert cu.a3m_cluster_hash(a3m)[0] == h1

        mdir = cu.method_dir(pair, "ccmpred")
        out = "ShallowMsa_000.ccmpred.npz"
        # legacy output, no manifest -> stale
        os.makedirs(mdir, exist_ok=True)
        np.savez(os.path.join(mdir, out), cmap=np.zeros((5, 5)), idx=np.arange(5))
        assert cu.entry_status(pair, "ccmpred", out)[:2] == (False, "no_manifest_entry")

        cu.record_run(pair, "ccmpred", a3m, "fine")
        # output written BEFORE the recorded run -> stale
        old = time.time() - 3600
        os.utime(os.path.join(mdir, out), (old, old))
        assert cu.entry_status(pair, "ccmpred", out)[1] == "output_predates_run"
        np.savez(os.path.join(mdir, out), cmap=np.zeros((5, 5)), idx=np.arange(5))
        ok, reason, entry = cu.entry_status(pair, "ccmpred", out)
        assert ok and reason == "ok" and entry["cluster_hash"] == h1
        assert entry["input_a3m"] == "output_msa_cluster/ShallowMsa_000.a3m"
        assert cu.a3m_is_fresh(pair, "ccmpred", a3m)

        # re-clustering: same label, different members -> stale
        _write_a3m(a3m, [("S1", "ACDEF"), ("v9", "ACDEY")])
        assert cu.entry_status(pair, "ccmpred", out)[:2] == (False, "hash_mismatch")
        assert not cu.a3m_is_fresh(pair, "ccmpred", a3m)
        os.remove(a3m)
        assert cu.entry_status(pair, "ccmpred", out)[:2] == (False, "input_missing")
    finally:
        shutil.rmtree(tmp)


# ---------------------------------------------------------------------------
# scorer skips stale npz; top-L recorded
# ---------------------------------------------------------------------------
def _make_pair(tmp):
    pair_id = "1abcA_2defB"
    pair = os.path.join(tmp, pair_id)
    n = 30
    seed = "A" * n
    _write_a3m(os.path.join(pair, "_seed_both.a3m"), [("1abcA", seed), ("2defB", seed)])
    clus = os.path.join(pair, "output_msa_cluster")
    for k, tail in enumerate(["C", "D", "E"]):
        _write_a3m(os.path.join(clus, f"ShallowMsa_{k:03d}.a3m"),
                   [("1abcA", seed), (f"m{k}", seed[:-1] + tail)])
    _write_a3m(os.path.join(pair, "output_msa_cluster_coarse", "CoarseMsa_000.a3m"),
               [("1abcA", seed)])
    return pair_id, pair, n


def test_scorer_skips_stale_and_uses_top_l():
    ca = _import_cmap_analysis()
    tmp = tempfile.mkdtemp()
    old_dd, old_truth = ca.DATA_DIR, ca._truth_contact_map
    try:
        pair_id, pair, n = _make_pair(tmp)
        truth = np.zeros((n, n), dtype=bool)
        for a in range(0, n - 8, 2):
            truth[a, a + 8] = truth[a + 8, a] = True
        ca.DATA_DIR = tmp
        ca._truth_contact_map = lambda *a, **k: truth.copy()

        mdir = cu.method_dir(pair, "ccmpred")
        os.makedirs(mdir, exist_ok=True)
        rng = np.random.default_rng(1)

        def put(name):
            P = rng.random((n, n)) * 0.2       # coupling-scale scores (<0.4)
            np.savez(os.path.join(mdir, name), cmap=P.astype(np.float32),
                     idx=np.arange(n, dtype=np.int32))

        clus = os.path.join(pair, "output_msa_cluster")
        # 000: fresh; 001: recorded then re-clustered; 002: legacy (no entry)
        for k in (0, 1):
            cu.record_run(pair, "ccmpred", os.path.join(clus, f"ShallowMsa_{k:03d}.a3m"), "fine")
        for k in (0, 1, 2):
            put(f"ShallowMsa_{k:03d}.ccmpred.npz")
        _write_a3m(os.path.join(clus, "ShallowMsa_001.a3m"), [("1abcA", "A" * n), ("zz", "C" * n)])
        cu.record_run(pair, "ccmpred",
                      os.path.join(pair, "output_msa_cluster_coarse", "CoarseMsa_000.a3m"), "coarse")
        put("CoarseMsa_000.ccmpred.npz")

        df = ca.compute_cmap_metrics_for_pair(pair_id, cmap_method="ccmpred")
        assert list(df["cluster"]) == ["ShallowMsa_000"], df["cluster"].tolist()
        row = df.iloc[0]
        assert row["top_l_factor"] == 1.0 and row["top_l_k"] == n
        assert row["t1_contacts_pred"] == n          # exactly top-L pairs predicted
        assert row["cluster_hash"] == cu.a3m_cluster_hash(
            os.path.join(clus, "ShallowMsa_000.a3m"))[0]
        assert os.path.isfile(os.path.join(pair, "Analysis", "df_cmap_ccmpred.csv"))

        # coarse only when asked
        dfc = ca.compute_cmap_metrics_for_pair(pair_id, cmap_method="ccmpred", resolution="coarse")
        assert sorted(dfc["cluster"]) == ["CoarseMsa_000"]
        # threshold mode (old behaviour) on coupling-scale scores -> 0 contacts
        dft = ca.compute_cmap_metrics_for_pair(pair_id, cmap_method="ccmpred", top_l_factor=None)
        assert int(dft.iloc[0]["t1_contacts_pred"]) == 0
        # unsafe opt-out scores everything at fine resolution
        dfa = ca.compute_cmap_metrics_for_pair(pair_id, cmap_method="ccmpred",
                                               require_manifest=False)
        assert len(dfa) == 3

        # survey-side filter: fresh rows kept, legacy CSV ignored, re-cluster drops
        kept = cu.filter_fresh_cmap_rows(df, pair, "ccmpred", verbose=False)
        assert len(kept) == 1
        legacy = df.drop(columns=["cluster_hash"])
        assert cu.filter_fresh_cmap_rows(legacy, pair, "ccmpred", verbose=False).empty
        _write_a3m(os.path.join(clus, "ShallowMsa_000.a3m"), [("1abcA", "A" * n)])
        assert cu.filter_fresh_cmap_rows(df, pair, "ccmpred", verbose=False).empty
    finally:
        ca.DATA_DIR, ca._truth_contact_map = old_dd, old_truth
        shutil.rmtree(tmp)


def test_resolve_top_l_factor():
    assert cu.resolve_top_l_factor("ccmpred") == 1.0
    assert cu.resolve_top_l_factor("msa_transformer") is None
    assert cu.resolve_top_l_factor("msa_transformer", "0.5") == 0.5
    assert cu.resolve_top_l_factor("ccmpred", "none") is None
    assert cu.resolve_top_l_factor("ccmpred", None) is None


if __name__ == "__main__":
    n_ok = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
            n_ok += 1
    print(f"{n_ok} tests passed")
