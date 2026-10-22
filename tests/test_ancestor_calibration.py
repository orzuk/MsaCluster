#!/usr/bin/env python3
"""Unit tests for scripts/ancestor_calibration.py.

Covers the parts that can be wrong silently and that no cluster run would
catch:

  1. the Delta sign conventions (TM / DDG / SS) -- positive must mean
     "prefers F1 = the alternative 2n54B fold" for every method;
  2. ranks + Spearman rho (incl. ties) against scipy when available, and
     against hand-computed values when not;
  3. the EXACT permutation p: enumerates 5! = 120 orderings, so the
     attainable values are multiples of 1/120 and a perfectly ordered
     method hits the 2/120 floor imposed by the Anc.0 / Anc.2 tie;
  4. projecting an ancestor onto the DeepMsa match-state column frame;
  5. end-to-end staging + collection on a synthetic pseudo-pair (no
     cluster data, no GPU, no external binaries);
  6. --dry_run prints the planned cluster commands and touches nothing.

Run:  python3 tests/test_ancestor_calibration.py
"""

from __future__ import annotations

import itertools
import math
import os
import random
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

import ancestor_calibration as ac  # noqa: E402

FAILS: list[str] = []


def _stub_scoring_deps() -> list[str]:
    """Make ``utils.score_against_folds`` importable without the heavy
    scientific stack (torch, tmtools, mdtraj, ...).

    Imports it, and on each ModuleNotFoundError installs a permissive stub
    for the missing module and retries. TM numbers from the stubbed
    ``tmtools`` are meaningless -- only the plumbing (one row per cluster,
    the right columns, delta = tm_f1 - tm_f2) is being checked here. On a
    machine with the real stack nothing is stubbed.

    Returns the list of stubbed module names.
    """
    import types
    from unittest.mock import MagicMock

    class _TMRes:
        tm_norm_chain1 = 0.5
        tm_norm_chain2 = 0.5
        t = np.zeros(3)
        u = np.eye(3)
        rmsd = 0.0

    stubbed: list[str] = []
    for _ in range(20):
        try:
            from utils.score_against_folds import safe_tm_vs_truth  # noqa: F401
            return stubbed
        except ModuleNotFoundError as e:
            name = e.name
            if not name or name.startswith("utils"):
                raise
            if name == "tmtools":
                mod = types.ModuleType("tmtools")
                mod.tm_align = lambda *a, **k: _TMRes()
                sys.modules["tmtools"] = mod
                io_mod = types.ModuleType("tmtools.io")
                io_mod.get_structure = lambda *a, **k: None
                io_mod.get_residue_data = lambda *a, **k: (np.zeros((1, 3)), "A")
                sys.modules["tmtools.io"] = io_mod
            else:
                sys.modules[name] = MagicMock()
            stubbed.append(name)
    raise RuntimeError(f"could not import score_against_folds; stubbed {stubbed}")


def check(cond, msg):
    if cond:
        print(f"  ok   {msg}")
    else:
        print(f"  FAIL {msg}")
        FAILS.append(msg)


# ----------------------------------------------------------------------
# 1. Delta conventions
# ----------------------------------------------------------------------

def test_deltas():
    print("\n[1] delta conventions (positive = prefers F1 = alternative fold)")
    check(ac.delta_from_tm(0.80, 0.30) > 0,
          "TM: higher TM to F1 -> positive delta")
    check(ac.delta_from_tm(0.30, 0.80) < 0,
          "TM: higher TM to F2 -> negative delta")
    check(abs(ac.delta_from_tm(0.72, 0.41) - 0.31) < 1e-12,
          "TM: delta == tm_f1 - tm_f2 exactly")
    check(math.isnan(ac.delta_from_tm(float("nan"), 0.4)),
          "TM: NaN propagates")
    check(math.isnan(ac.delta_from_tm(None, 0.4)),
          "TM: None -> NaN")

    # DDG: bias = sum DDG(F2) - sum DDG(F1). A sequence destabilised MORE on
    # F2 (larger positive DDG on F2) prefers F1.
    check(ac.delta_from_ddg(1.0, 5.0) > 0,
          "DDG: more destabilised on F2 -> prefers F1 -> positive")
    check(ac.delta_from_ddg(5.0, 1.0) < 0,
          "DDG: more destabilised on F1 -> prefers F2 -> negative")
    check(abs(ac.delta_from_ddg(-2.5, 1.5) - 4.0) < 1e-12,
          "DDG: delta == sum_F2 - sum_F1 exactly (kcal/mol)")
    # must agree with run_DDG._score_sequence's own 'bias'
    c1 = np.array([0.5, -1.0, 2.0], dtype=np.float32)
    c2 = np.array([1.5, 0.0, 2.5], dtype=np.float32)
    bias_run_ddg = float(c2.sum() - c1.sum())
    check(abs(ac.delta_from_ddg(float(c1.sum()), float(c2.sum()))
              - bias_run_ddg) < 1e-6,
          "DDG: matches run_DDG._score_sequence's bias formula")

    check(ac.delta_from_ss(0.9, 0.4) > 0,
          "SS: closer to F1's (%H,%E) -> positive")
    check(abs(ac.delta_from_ss(0.9, 0.4) - 0.5) < 1e-12,
          "SS: delta == sim_f1 - sim_f2 exactly")
    check(math.isnan(ac.delta_from_ss(float("nan"), 0.4)),
          "SS: NaN propagates")


# ----------------------------------------------------------------------
# 2. ranks + Spearman
# ----------------------------------------------------------------------

def test_spearman():
    print("\n[2] rankdata / spearman_rho")
    check(list(ac.rankdata([10, 20, 30])) == [1.0, 2.0, 3.0],
          "rankdata: strictly increasing")
    check(list(ac.rankdata([5, 5, 9])) == [1.5, 1.5, 3.0],
          "rankdata: ties averaged")
    check(list(ac.rankdata([3, 1, 1, 1])) == [4.0, 2.0, 2.0, 2.0],
          "rankdata: 3-way tie averaged")

    check(abs(ac.spearman_rho([1, 2, 3, 4, 5], [1, 2, 3, 4, 5]) - 1.0) < 1e-12,
          "rho(x, x) == 1")
    check(abs(ac.spearman_rho([1, 2, 3, 4, 5], [5, 4, 3, 2, 1]) + 1.0) < 1e-12,
          "rho(x, reversed x) == -1")

    # the real ground-truth vector has a tie, so a perfect method cannot hit 1
    y = [g["pct_alt"] for g in ac.GROUND_TRUTH]
    check(y == [0.0, 0.0, 8.0, 90.7, 50.0],
          "GROUND_TRUTH %alt in evolutionary order == [0, 0, 8, 90.7, 50]")
    perfect = [0.0, 0.0, 0.2, 0.9, 0.5]
    r = ac.spearman_rho(perfect, y)
    check(abs(r - 1.0) < 1e-12,
          "rho == 1 when a method reproduces the order including the tie")

    try:
        from scipy.stats import spearmanr
    except Exception:
        print("  skip scipy cross-check (scipy not installed)")
        return
    rng = random.Random(7)
    worst = 0.0
    for _ in range(300):
        n = rng.randint(4, 12)
        a = [rng.choice([1, 2, 2, 3, 4, 5, 5, 6, 7]) for _ in range(n)]
        b = [rng.gauss(0, 1) for _ in range(n)]
        worst = max(worst, abs(ac.spearman_rho(a, b) - spearmanr(a, b).statistic))
    check(worst < 1e-10,
          f"spearman_rho matches scipy on 300 random vectors (max |diff|={worst:.2e})")


# ----------------------------------------------------------------------
# 3. exact permutation p
# ----------------------------------------------------------------------

def test_exact_p():
    print("\n[3] exact permutation p over 5! = 120 orderings")
    y = [g["pct_alt"] for g in ac.GROUND_TRUTH]

    p, n = ac.exact_perm_p([0.0, 0.1, 0.2, 0.9, 0.5], y, "greater")
    check(n == 120, "enumerates exactly 120 permutations")
    check(abs(p - 2.0 / 120.0) < 1e-12,
          f"perfectly ordered method -> p = 2/120 = {2/120:.4f} "
          f"(the Anc.0/Anc.2 tie makes 2 permutations equivalent); got {p:.4f}")

    p_rev, _ = ac.exact_perm_p([0.9, 0.5, 0.2, 0.0, 0.1], y, "greater")
    check(abs(p_rev - 1.0) < 1e-12,
          "exactly anti-ordered method -> one-sided p = 1.0")
    p_rev_l, _ = ac.exact_perm_p([0.9, 0.5, 0.2, 0.0, 0.1], y, "less")
    check(p_rev_l <= 2.0 / 120.0 + 1e-12,
          "exactly anti-ordered method -> 'less' p hits the floor")

    # p must be a multiple of 1/120 and in (0, 1]
    rng = random.Random(3)
    ok = True
    for _ in range(80):
        x = [rng.gauss(0, 1) for _ in range(5)]
        pv, _ = ac.exact_perm_p(x, y, "greater")
        k = pv * 120.0
        ok &= (abs(k - round(k)) < 1e-9) and (0 < pv <= 1.0)
    check(ok, "p is always a multiple of 1/120 and lies in (0, 1]")

    # uniformity: under the null, P(p <= a) <= a
    cnt = 0
    N = 1200
    for _ in range(N):
        x = [rng.gauss(0, 1) for _ in range(5)]
        pv, _ = ac.exact_perm_p(x, y, "greater")
        cnt += int(pv <= 0.10)
    frac = cnt / N
    check(frac <= 0.16,
          f"null calibration: P(p<=0.10) = {frac:.3f} is not inflated")

    # two-sided is never smaller than the better one-sided tail
    for _ in range(40):
        x = [rng.gauss(0, 1) for _ in range(5)]
        pg, _ = ac.exact_perm_p(x, y, "greater")
        pl, _ = ac.exact_perm_p(x, y, "less")
        pt, _ = ac.exact_perm_p(x, y, "two-sided")
        if pt < min(pg, pl) - 1e-12:
            check(False, "two-sided p >= min(one-sided tails)")
            return
    check(True, "two-sided p >= min(one-sided tails)")

    check(math.isnan(ac.exact_perm_p([1, 2, float("nan"), 4, 5], y)[0]),
          "NaN in the method vector -> NaN p (no silent partial test)")


def test_zscore():
    print("\n[4] zscore")
    z = ac.zscore([1.0, 2.0, 3.0])
    check(abs(float(np.mean(z))) < 1e-12 and abs(float(np.std(z)) - 1) < 1e-12,
          "zscore: mean 0, sd 1")
    check(np.isnan(ac.zscore([2.0, 2.0, 2.0])).all(),
          "zscore: no spread -> all NaN (never divides by zero)")


# ----------------------------------------------------------------------
# 5. projection into the DeepMsa column frame
# ----------------------------------------------------------------------

def test_projection():
    print("\n[5] project_to_query_frame")
    q = "ACDEFGHIKL"
    row, dropped = ac.project_to_query_frame("ACDEFGHIKL", q)
    check(row == "ACDEFGHIKL" and dropped == 0,
          "identical sequence -> identity projection")

    row, dropped = ac.project_to_query_frame("ACDEFGHIKLMNPQRSTVWY",
                                             "ACDEFGHIKLMNPQRSTVWY")
    check(len(row) == 20 and dropped == 0, "20-mer identity projection")

    # deletion relative to the query: the row keeps the query's width
    q2 = "ACDEFGHIKLMNPQRSTVWY"
    row, dropped = ac.project_to_query_frame("ACDEFGHIKLMNPQRSTVWY"[:10] +
                                             "MNPQRSTVWY", q2)
    check(len(row) == len(q2), "projected row always has the query's width")

    # a gapped query row: output width == query row width, residues land in
    # the query's occupied columns only
    qrow = "AC-DEF--GH"
    row, _ = ac.project_to_query_frame("ACDEFGH", qrow)
    check(len(row) == len(qrow), "gapped query row: width preserved")
    check(all(row[j] == "-" for j, c in enumerate(qrow) if c == "-"),
          "gapped query row: query-gap columns stay gaps")
    check(ac.ungap(row) == "ACDEFGH",
          "gapped query row: all target residues placed")


# ----------------------------------------------------------------------
# 6. end-to-end on a synthetic pseudo-pair
# ----------------------------------------------------------------------

CA_TMPL = ("ATOM  {idx:5d}  CA  ALA {ch}{res:4d}    "
           "{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n")


def _write_toy_pdb(path: Path, chain: str, n: int, helix: bool) -> None:
    lines = []
    for i in range(n):
        if helix:
            t = i * 100.0 * math.pi / 180.0
            x, y, z = 2.3 * math.cos(t), 2.3 * math.sin(t), 1.5 * i
        else:
            x, y, z = 3.3 * i, 0.0, 0.0
        lines.append(CA_TMPL.format(idx=i + 1, ch=chain, res=i + 1,
                                    x=x, y=y, z=z))
    lines.append("TER\nEND\n")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(lines))


def test_staging_and_collection():
    print("\n[6] staging + collection on a synthetic pseudo-pair")
    anc = ac.read_fasta(REPO / "data" / "ancestors" / "dishman_ancestors.fasta")
    names = [n for n, _ in anc]
    check(names == ac.ANCESTOR_NAMES,
          f"data/ancestors FASTA holds {ac.ANCESTOR_NAMES} in order")
    xcl1 = dict(anc)["XCL1"]

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        ref = td / "ref" / ac.DEFAULT_PAIR
        # synthetic truth structures + chain fastas
        _write_toy_pdb(ref / "chain_pdb_files" / "2n54B.pdb", "B", 60, False)
        _write_toy_pdb(ref / "chain_pdb_files" / "2hdmA.pdb", "A", 60, True)
        _write_toy_pdb(ref / "2n54.pdb", "B", 60, False)
        _write_toy_pdb(ref / "2hdm.pdb", "A", 60, True)
        (ref / "fasta_chain_files").mkdir(parents=True, exist_ok=True)
        (ref / "fasta_chain_files" / "2n54B.fasta").write_text(
            f">2n54B\n{xcl1}\n")
        (ref / "fasta_chain_files" / "2hdmA.fasta").write_text(
            f">2hdmA\n{xcl1}\n")
        # synthetic deep MSA: query = XCL1, plus 40 mutated family members
        rng = random.Random(11)
        rows = [(">XCL1_query", xcl1)]
        for i in range(40):
            s = list(xcl1)
            for _ in range(rng.randint(2, 12)):
                j = rng.randrange(len(s))
                s[j] = rng.choice("ACDEFGHIKLMNPQRSTVWY")
            rows.append((f">hit_{i:03d}", "".join(s)))
        deep = ref / "output_get_msa" / "DeepMsa.a3m"
        deep.parent.mkdir(parents=True, exist_ok=True)
        deep.write_text("".join(f"{h}\n{s}\n" for h, s in rows))

        single = td / "work" / "single" / "FoldPairs" / ac.DEFAULT_PAIR
        withmsa = td / "work" / "withmsa" / "FoldPairs" / ac.DEFAULT_PAIR
        ac.stage_pseudo_pair(ref, single, ac.DEFAULT_PAIR)
        ac.stage_pseudo_pair(ref, withmsa, ac.DEFAULT_PAIR)
        check((single / "chain_pdb_files" / "2n54B.pdb").is_file()
              and (single / "output_get_msa" / "DeepMsa.a3m").is_file()
              and (single / "2hdm.pdb").is_file(),
              "stage_pseudo_pair links truth PDBs, chain FASTAs and DeepMsa")
        check(sorted(p.name for p in ref.rglob("*")
                     if p.name.startswith("ShallowMsa")) == [],
              "production pair dir is never written to")

        manifest = ac.build_cluster_a3ms(
            single / "output_get_msa" / "DeepMsa.a3m", anc,
            n_controls=8, control_seed=0, msa_top_n=5,
            single_dir=single / "output_msa_cluster",
            withmsa_dir=withmsa / "output_msa_cluster")
        check(len(manifest) == 13, "5 ancestors + 8 controls -> 13 clusters")
        check([m["label"] for m in manifest[:5]] == ac.ANCESTOR_NAMES,
              "ancestors occupy ShallowMsa_000..004 in evolutionary order")

        singles = sorted((single / "output_msa_cluster").glob("*.a3m"))
        check(len(singles) == 13, "13 single-sequence cluster a3ms written")
        n_rows_single = {len(ac.read_fasta(p)) for p in singles}
        check(n_rows_single == {1},
              "every 'single' cluster a3m has exactly one sequence "
              "(so DDG/S4PRED score the ancestor itself, not its neighbours)")
        w = sorted((withmsa / "output_msa_cluster").glob("*.a3m"))
        check(all(len(ac.read_fasta(p)) == 6 for p in w),
              "every 'withmsa' cluster a3m has the query + 5 neighbours")
        widths = {len(s) for p in w for _, s in ac.read_fasta(p)}
        check(len(widths) == 1,
              "withmsa a3m preserves the DeepMsa column width "
              "(required by _write_pair_a3m_for_chain)")
        # The projected row lives in the DeepMsa match-state frame, so
        # residues that align to no query column are dropped (an a3m would
        # push them into lowercase insert columns, which every per-cluster
        # reader in this repo strips). Check it is a subsequence, that the
        # loss is small, and that the manifest records it.
        def _is_subseq(small, big):
            it = iter(big)
            return all(c in it for c in small)

        for m, p in zip(manifest[:5], singles[:5]):
            got = ac.ungap(ac.read_fasta(p)[0][1])
            full = dict(anc)[m["label"]]
            check(_is_subseq(got, full),
                  f"{m['label']}: projected row is a subsequence of the ancestor")
            check(len(full) - len(got) == m["n_dropped_insertions"],
                  f"{m['label']}: manifest n_dropped_insertions "
                  f"({m['n_dropped_insertions']}) matches the residues lost")
            check(len(got) >= 0.9 * len(full),
                  f"{m['label']}: >=90% of residues survive the projection")

        # ---- collection: fake ESMFold2 outputs with a planted trend ----
        for i, m in enumerate(manifest):
            # ancestors 0..4: progressively more "extended" (F1-like)
            frac = (i / 4.0) if m["kind"] == "ancestor" else rng.random()
            _write_toy_pdb(single / "output_esmfold2" / f"{m['cluster']}.pdb",
                           "A", 60, helix=(frac < 0.5))
        stubbed = _stub_scoring_deps()
        if stubbed:
            print(f"  note stubbed for the TM path: {', '.join(stubbed)}")
        rows_e = ac.collect_esmfold2(single, manifest)
        check(len(rows_e) == 13, "collect_esmfold2 returns one row per cluster")
        check(all(set(("tm_f1", "tm_f2", "delta")) <= set(r) for r in rows_e),
              "ESMFold2 rows carry tm_f1 / tm_f2 / delta")
        check(all(r["pred_pdb"].endswith(".pdb") for r in rows_e),
              "ESMFold2 rows record the prediction they were scored from")

        # ---- collection: synthetic df_ddg.csv ----
        (single / "Analysis").mkdir(parents=True, exist_ok=True)
        with open(single / "Analysis" / "df_ddg.csv", "w") as fh:
            fh.write("pair_id,cluster,n_seqs,bias_mean,bias_std,"
                     "frac_F1_preferring,sum_F1_mean,sum_F2_mean\n")
            for i, m in enumerate(manifest):
                s1, s2 = 10.0, 10.0 + (i if m["kind"] == "ancestor" else 0.0)
                fh.write(f"{ac.DEFAULT_PAIR},{m['cluster_num']},1,"
                         f"{s2 - s1},0.0,0.0,{s1},{s2}\n")
        rows_d = ac.collect_ddg(single, manifest)
        check(len(rows_d) == 13, "collect_ddg returns one row per cluster")
        d_by = {r["label"]: r["delta"] for r in rows_d}
        check(d_by["Anc.0"] == 0.0 and abs(d_by["XCL1"] - 4.0) < 1e-9,
              "collect_ddg reproduces sum_F2 - sum_F1 (kcal/mol)")

        # ---- stats on a planted, perfectly ordered method ----
        planted = []
        for i, m in enumerate(manifest):
            if m["kind"] == "ancestor":
                # tied at the two monomorphic ancestors, exactly like the
                # ground truth, so a perfect method can reach rho == 1
                val = {"Anc.0": 0.0, "Anc.2": 0.0, "Anc.3": 0.2,
                       "Anc.4": 0.9, "XCL1": 0.5}[m["label"]]
            else:
                val = rng.gauss(0.3, 0.3)
            planted.append({"label": m["label"], "kind": m["kind"],
                            "cluster": m["cluster"], "method": "PLANTED",
                            "mode": "single_sequence", "delta": val})
        st = ac.method_stats(planted, n_control_draws=500, seed=0)
        check(len(st) == 1 and abs(st[0]["rho"] - 1.0) < 1e-12,
              "method_stats: planted perfect ordering -> rho == 1")
        check(abs(st[0]["p_exact_greater"] - 2.0 / 120.0) < 1e-12,
              "method_stats: planted perfect ordering -> exact p = 2/120")
        check(st[0]["n_controls"] == 8 and 0.0 <= st[0]["p_control"] <= 1.0,
              "method_stats: control null uses all 8 control sequences")


# ----------------------------------------------------------------------
# 7. --dry_run
# ----------------------------------------------------------------------

def test_dry_run():
    print("\n[7] --dry_run")
    env = dict(os.environ)
    env["MPLBACKEND"] = "Agg"
    with tempfile.TemporaryDirectory() as td:
        r = subprocess.run(
            [sys.executable, str(REPO / "scripts" / "ancestor_calibration.py"),
             "--dry_run", "--methods", "all",
             "--work_root", str(Path(td) / "wr"),
             "--s4pred_dir", "/fake/s4pred",
             "--thermompnn_dir", "/fake/thermompnn"],
            capture_output=True, text=True, env=env)
        out = r.stdout + r.stderr
        check(r.returncode == 0, f"--dry_run exits 0 (got {r.returncode})")
        for needle, what in [
            ("run_ESMFold2.py", "plans ESMFold2"),
            ("run_DDG.py", "plans DDG"),
            ("run_s4pred_per_cluster.py", "plans S4PRED"),
            ("run_Boltz2.py", "plans Boltz-2"),
            ("run_bioemu.py", "plans BioEmu"),
            ("RunAF2_Colabfold.sh", "plans AF2"),
            ("MSACLUSTER_DATA_DIR", "passes MSACLUSTER_DATA_DIR to the "
                                    "pair-id-driven runners"),
            ("query_type=first", "builds AF2 a3m with query_type=first"),
        ]:
            check(needle in out, what)
        check(not (Path(td) / "wr").exists(),
              "--dry_run stages nothing on disk")


# ----------------------------------------------------------------------
# 8. full main() on synthetic outputs: CSV + stats + figure
# ----------------------------------------------------------------------

def test_end_to_end_collect(keep_dir: Path | None = None):
    print("\n[8] main() --collect_only: CSV, stats CSV and figure")
    import json
    td = Path(tempfile.mkdtemp()) if keep_dir is None else keep_dir
    td.mkdir(parents=True, exist_ok=True)
    work = td / "work"
    single = work / "single" / "FoldPairs" / ac.DEFAULT_PAIR
    (single / "Analysis").mkdir(parents=True, exist_ok=True)

    # manifest: 5 ancestors + 20 controls
    manifest = []
    for i, name in enumerate(ac.ANCESTOR_NAMES):
        manifest.append({"cluster": f"ShallowMsa_{i:03d}",
                         "cluster_num": f"{i:03d}", "label": name,
                         "kind": "ancestor", "seq": "A" * 66, "seq_len": 66,
                         "n_dropped_insertions": 0, "n_neighbours": 10,
                         "source": "dishman_fig_S1C"})
    rng = random.Random(5)
    for j in range(20):
        i = 5 + j
        manifest.append({"cluster": f"ShallowMsa_{i:03d}",
                         "cluster_num": f"{i:03d}", "label": f"ctrl_{j:02d}",
                         "kind": "control", "seq": "A" * 66, "seq_len": 66,
                         "n_dropped_insertions": 0, "n_neighbours": 10,
                         "source": "DeepMsa"})
    (work / "ancestor_manifest.json").write_text(json.dumps(manifest))

    # a DDG table that reproduces the measured ordering for the ancestors
    # and is pure noise for the controls
    planted = {"Anc.0": 0.0, "Anc.2": 0.0, "Anc.3": 1.0,
               "Anc.4": 6.0, "XCL1": 3.0}
    with open(single / "Analysis" / "df_ddg.csv", "w") as fh:
        fh.write("pair_id,cluster,n_seqs,bias_mean,bias_std,"
                 "frac_F1_preferring,sum_F1_mean,sum_F2_mean\n")
        for m in manifest:
            b = planted.get(m["label"], rng.gauss(0.0, 2.0))
            fh.write(f"{ac.DEFAULT_PAIR},{m['cluster_num']},1,{b},0,0,"
                     f"0.0,{b}\n")

    out_csv = td / "ancestor_calibration.csv"
    stats_csv = td / "ancestor_calibration_stats.csv"
    fig = td / "ancestor_calibration.png"
    rc = ac.main(["--collect_only", "--methods", "ddg",
                  "--work_root", str(work), "--n_controls", "20",
                  "--out", str(out_csv), "--stats_out", str(stats_csv),
                  "--fig", str(fig)])
    check(rc == 0, f"main(--collect_only) exits 0 (got {rc})")
    check(out_csv.is_file(), "writes the per-(sequence, method, mode) CSV")
    check(stats_csv.is_file(), "writes the per-method stats CSV")
    check(fig.is_file() and fig.stat().st_size > 10000,
          "writes a non-trivial calibration figure")
    check(fig.with_suffix(".pdf").is_file(), "writes the PDF companion")

    import csv as _csv
    with open(out_csv) as fh:
        rows = list(_csv.DictReader(fh))
    check(len(rows) == 25, f"25 rows (5 ancestors + 20 controls); got {len(rows)}")
    anc_rows = {r["label"]: r for r in rows if r["kind"] == "ancestor"}
    check(set(anc_rows) == set(ac.ANCESTOR_NAMES), "all 5 ancestors present")
    check(anc_rows["Anc.4"]["pct_alt_measured"] == "90.7"
          and anc_rows["Anc.3"]["pct_alt_measured_sd"] == "3.7",
          "the measured NMR values are carried into the CSV")
    check(abs(float(anc_rows["Anc.4"]["delta"]) - 6.0) < 1e-9,
          "delta column holds the DDG bias in kcal/mol")
    zs = [float(anc_rows[n]["delta_z"]) for n in ac.ANCESTOR_NAMES]
    check(abs(sum(zs) / 5) < 1e-9, "delta_z is z-scored over the 5 ancestors")

    with open(stats_csv) as fh:
        st = list(_csv.DictReader(fh))
    check(len(st) == 1 and st[0]["method"] == "DDG", "one stats row for DDG")
    check(abs(float(st[0]["rho"]) - 1.0) < 1e-9,
          "planted monotone DDG -> rho == 1")
    check(abs(float(st[0]["p_exact_greater"]) - 2.0 / 120.0) < 1e-9,
          "planted monotone DDG -> exact one-sided p = 2/120")
    check(int(st[0]["n_perms"]) == 120, "stats CSV records the 120 permutations")
    check(int(st[0]["n_controls"]) == 20, "all 20 controls enter the null")
    check(float(st[0]["p_control"]) < 0.25,
          f"random family sequences rarely reproduce the trend "
          f"(p_control={st[0]['p_control']})")
    if keep_dir is None:
        import shutil
        shutil.rmtree(td, ignore_errors=True)
    else:
        print(f"  kept outputs under {td}")


def main() -> int:
    print("=" * 70)
    print("tests/test_ancestor_calibration.py")
    print("=" * 70)
    test_deltas()
    test_spearman()
    test_exact_p()
    test_zscore()
    test_projection()
    test_staging_and_collection()
    test_dry_run()
    keep = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    test_end_to_end_collect(keep)
    print("\n" + "=" * 70)
    if FAILS:
        print(f"{len(FAILS)} FAILURE(S):")
        for f in FAILS:
            print(f"  - {f}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
