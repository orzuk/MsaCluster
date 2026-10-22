#!/usr/bin/env python3
"""Calibrate our per-sequence fold-preference predictors against the
EXPERIMENTALLY measured evolutionary trajectory of XCL1 (lymphotactin).

The experiment
--------------
Dishman et al. (Science 2021, 371:86-90) reconstructed ancestral XC
chemokines and measured, by NMR, what fraction of each ancestor populates
the canonical chemokine fold vs the alternative all-beta dimer fold. That
gives a ground-truth ordering along one real evolutionary lineage:

    Anc.0 ~ Anc.2  <  Anc.3  <  XCL1  <  Anc.4      (preference for ALT fold)

This script asks whether the per-cluster fold-preference methods used
throughout this repo reproduce that ordering when each ancestral sequence
is handed to them as a one-sequence "cluster". It is a calibration, not a
new method: if the predictors cannot recover a measured trajectory on the
one family where the ground truth exists, no amount of better MSA
clustering can support clade-level inference. If they can, the methods
have genuine per-sequence fold sensitivity.

Reference folds: the fold-switch pair ``2n54B_2hdmA``.
    F1 = 2N54 chain B  -- alternative (all-beta dimer) XCL1 fold
    F2 = 2HDM chain A  -- canonical chemokine fold
Repo-wide convention: ``Delta = score(F1) - score(F2)``, so
**positive Delta = prefers the ALTERNATIVE fold**.

How it reuses the pipeline
--------------------------
No method code is re-implemented. The script stages two FoldPairs-shaped
*pseudo-pair* trees under ``--work_root`` (the production pair dir is
never written to) and points the existing runners at them:

    <work_root>/single/FoldPairs/2n54B_2hdmA/     one sequence per "cluster"
    <work_root>/withmsa/FoldPairs/2n54B_2hdmA/    sequence + top-N neighbours

Truth PDBs, chain FASTAs and ``output_get_msa/DeepMsa.a3m`` are symlinked
in from the production pair dir. Each ancestor becomes
``output_msa_cluster/ShallowMsa_NNN.a3m``, written in the DeepMsa column
frame (the ancestor is projected onto the DeepMsa query by the project's
canonical pairwise aligner, ``utils.msa_utils.get_align_indexes``).

    method     runner                                reached via
    ---------  ------------------------------------  --------------------------
    ESMFold2   run_ESMFold2.py                       --input_dir / --output_dir
    DDG        run_DDG.py (ThermoMPNN)               MSACLUSTER_DATA_DIR
    S4PRED     scripts/run_s4pred_per_cluster.py     --fpairs_root
    Boltz-2    run_Boltz2.py                         MSACLUSTER_DATA_DIR
    BioEmu     run_bioemu.py                         MSACLUSTER_DATA_DIR
    AF2        scripts/shell/RunAF2_Colabfold.sh     a3m built by
               (a3m from run_foldswitch_pipeline      _write_pair_a3m_for_chain
                ._write_pair_a3m_for_chain)           with --query_type first

``MSACLUSTER_DATA_DIR`` is the one-line override added to ``config.py``;
unset, everything resolves to ``Pipeline/FoldPairs`` exactly as before.

Negative control
----------------
The same computation is run for ``--n_controls`` random XCL-family
sequences drawn from the production deep MSA. They should NOT show the
monotone trend, which rules out the trend being an artefact of sequence
length or of identity to the 2N54/2HDM reference chains.

Outputs
-------
    docs/ancestor_calibration.csv              one row per (sequence, method, mode)
    docs/ancestor_calibration_stats.csv        per-method rho + exact p + control null
    docs/figures/global/ancestor_calibration.png

Examples
--------
    # what would run, no cluster needed
    python3 scripts/ancestor_calibration.py --dry_run

    # cheap CPU/GPU methods only
    python3 scripts/ancestor_calibration.py --methods esmfold2,ddg,s4pred \
        --s4pred_dir $S4PRED_DIR --thermompnn_dir $THERMOMPNN_DIR

    # re-score / re-plot from outputs already on disk
    python3 scripts/ancestor_calibration.py --methods all --collect_only
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import os
import random
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


# ======================================================================
# Ground truth -- Dishman et al. 2021, Science 371:86-90 (fig. S1C +
# main-text NMR populations). See data/ancestors/README.md.
# ======================================================================

#: Evolutionary order, root -> modern tip. Anc.4 is the direct ancestor of
#: XCL1, so the lineage is Anc.0 -> Anc.2 -> Anc.3 -> Anc.4 -> XCL1.
#: ``pct_alt`` is the measured % of the ALTERNATIVE (all-beta, 2N54) fold
#: at 37 C / 150 mM NaCl; ``pct_alt_sd`` is the reported s.d. where one was
#: given. Anc.0 and Anc.2 are monomorphic chemokine-fold, i.e. 0 % alt.
GROUND_TRUTH: List[dict] = [
    {"name": "Anc.0", "evo_order": 0, "pct_chemokine": 100.0, "pct_chemokine_sd": float("nan"),
     "pct_alt": 0.0, "pct_alt_sd": float("nan"),
     "note": "monomorphic chemokine fold; NMR structure PDB 7JH1"},
    {"name": "Anc.2", "evo_order": 1, "pct_chemokine": 100.0, "pct_chemokine_sd": float("nan"),
     "pct_alt": 0.0, "pct_alt_sd": float("nan"),
     "note": "single fold at all tested conditions"},
    {"name": "Anc.3", "evo_order": 2, "pct_chemokine": 92.0, "pct_chemokine_sd": 3.7,
     "pct_alt": 8.0, "pct_alt_sd": 3.7,
     "note": "metamorphic, 92 +/- 3.7 % chemokine at 37 C / 150 mM NaCl"},
    {"name": "Anc.4", "evo_order": 3, "pct_chemokine": 9.3, "pct_chemokine_sd": 2.0,
     "pct_alt": 90.7, "pct_alt_sd": 2.0,
     "note": "metamorphic, 9.3 +/- 2.0 % chemokine at 37 C / 150 mM NaCl"},
    {"name": "XCL1", "evo_order": 4, "pct_chemokine": 50.0, "pct_chemokine_sd": float("nan"),
     "pct_alt": 50.0, "pct_alt_sd": float("nan"),
     "note": "modern human tip, near 50:50"},
]

ANCESTOR_NAMES = [g["name"] for g in GROUND_TRUTH]

DEFAULT_PAIR = "2n54B_2hdmA"
F1_TAG, F2_TAG = "2n54B", "2hdmA"          # F1 = alternative fold, F2 = chemokine fold

#: Methods whose "cluster" MSA must be exactly one sequence (per-sequence
#: scorers: averaging over neighbours would wash out the ancestor signal).
SINGLE_TREE_METHODS = ("esmfold2", "ddg", "s4pred", "boltz2", "af2_single")
#: Methods that want the sequence plus its nearest family neighbours.
WITHMSA_TREE_METHODS = ("af2", "bioemu", "bioemu2")
ALL_METHODS = SINGLE_TREE_METHODS + WITHMSA_TREE_METHODS


# ======================================================================
# Delta conventions -- the single place each method's number is turned
# into a comparable "preference for F1 (= the alternative fold)".
# ======================================================================

def delta_from_tm(tm_f1: float, tm_f2: float) -> float:
    """TM-based delta: TM(F1) - TM(F2). Positive = prefers F1 (alt fold).

    Used by every structure-predicting method (ESMFold2, AF2, Boltz-2,
    BioEmu), matching ``utils.score_against_folds`` ("TMdiff").
    """
    if tm_f1 is None or tm_f2 is None:
        return float("nan")
    if not (np.isfinite(tm_f1) and np.isfinite(tm_f2)):
        return float("nan")
    return float(tm_f1) - float(tm_f2)


def delta_from_ddg(sum_ddg_f1: float, sum_ddg_f2: float) -> float:
    """ThermoMPNN conformation-bias delta, in kcal/mol.

    ``bias = sum DDG(F2) - sum DDG(F1)``: a sequence threaded onto both
    truth backbones prefers F1 when it is destabilised MORE on F2. This is
    exactly ``run_DDG._score_sequence``'s ``bias`` and keeps the sign
    convention "positive = prefers F1".
    """
    if not (np.isfinite(sum_ddg_f1) and np.isfinite(sum_ddg_f2)):
        return float("nan")
    return float(sum_ddg_f2) - float(sum_ddg_f1)


def delta_from_ss(sim_f1: float, sim_f2: float) -> float:
    """S4PRED secondary-structure delta: sim(F1) - sim(F2).

    ``sim`` is ``scripts.aggregate_s4pred_into_survey.he_similarity``,
    1 - normalised Euclidean distance in the (%H, %E) plane between the
    predicted Q3 composition and each fold's DSSP composition.
    """
    if not (np.isfinite(sim_f1) and np.isfinite(sim_f2)):
        return float("nan")
    return float(sim_f1) - float(sim_f2)


# ======================================================================
# Spearman + EXACT permutation p (5 points -> 5! = 120 orderings)
# ======================================================================

def rankdata(x: Sequence[float]) -> np.ndarray:
    """Average ranks (1-based), ties averaged. Same as scipy's default."""
    a = np.asarray(x, dtype=float)
    order = np.argsort(a, kind="mergesort")
    ranks = np.empty(len(a), dtype=float)
    ranks[order] = np.arange(1, len(a) + 1, dtype=float)
    # average over tied groups
    sorted_a = a[order]
    i = 0
    while i < len(a):
        j = i
        while j + 1 < len(a) and sorted_a[j + 1] == sorted_a[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = ranks[order[i:j + 1]].mean()
        i = j + 1
    return ranks


def spearman_rho(x: Sequence[float], y: Sequence[float]) -> float:
    """Spearman rank correlation (Pearson on average ranks; ties handled)."""
    rx, ry = rankdata(x), rankdata(y)
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    denom = math.sqrt(float((rx ** 2).sum()) * float((ry ** 2).sum()))
    if denom == 0.0:
        return float("nan")
    return float((rx * ry).sum() / denom)


def exact_perm_p(x: Sequence[float], y: Sequence[float],
                 alternative: str = "greater") -> Tuple[float, int]:
    """EXACT permutation p for Spearman rho, enumerating all n! orderings.

    With n = 5 that is 120 permutations, so the p-value is exact, not
    sampled. Returns ``(p, n_perms)``.

    ``alternative``:
      "greater"   -- H1: rho > 0 (the pre-registered direction here: a
                     method that works should rank the ancestors in the
                     measured order)
      "less"      -- H1: rho < 0
      "two-sided" -- H1: |rho| >= |rho_obs|

    Note the attainable floor: the ground-truth vector has a tie
    (Anc.0 = Anc.2 = 0 % alternative fold), so a perfectly ordered method
    is matched by the 2 permutations that swap those two positions and the
    smallest possible one-sided p is 2/120 = 0.0167.
    """
    x = list(x)
    y = list(y)
    if len(x) != len(y):
        raise ValueError("x and y must be the same length")
    if any(not np.isfinite(v) for v in x) or any(not np.isfinite(v) for v in y):
        return float("nan"), 0
    obs = spearman_rho(x, y)
    if not np.isfinite(obs):
        return float("nan"), 0
    perms = list(itertools.permutations(range(len(x))))
    hits = 0
    for p in perms:
        r = spearman_rho([x[i] for i in p], y)
        if alternative == "greater":
            ok = r >= obs - 1e-12
        elif alternative == "less":
            ok = r <= obs + 1e-12
        elif alternative == "two-sided":
            ok = abs(r) >= abs(obs) - 1e-12
        else:
            raise ValueError(f"bad alternative: {alternative!r}")
        hits += int(ok)
    return hits / len(perms), len(perms)


def zscore(v: Sequence[float]) -> np.ndarray:
    """Z-score, NaN-safe; returns all-NaN if the vector has no spread."""
    a = np.asarray(v, dtype=float)
    m = np.nanmean(a)
    s = np.nanstd(a)
    if not np.isfinite(s) or s == 0:
        return np.full_like(a, np.nan)
    return (a - m) / s


# ======================================================================
# FASTA / A3M helpers
# ======================================================================

STD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def read_fasta(path) -> List[Tuple[str, str]]:
    out: List[Tuple[str, str]] = []
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(chunks)))
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        out.append((name, "".join(chunks)))
    return out


def read_a3m_match_state(path) -> List[Tuple[str, str]]:
    """[(header, match-state row)] -- insertions (lowercase/'.') dropped.

    Thin wrapper over ``utils.cluster_representative._read_a3m_match_state``
    so this script and the per-cluster runners read a3m identically.
    """
    from utils.cluster_representative import _read_a3m_match_state
    return _read_a3m_match_state(Path(path))


def project_to_query_frame(target: str, query_row: str) -> Tuple[str, int]:
    """Place ``target`` into the column frame of a DeepMsa match-state row.

    ``query_row`` is the a3m's row-1 match-state string (length = number of
    match columns, L). Returns ``(row, n_dropped)`` where ``row`` has length
    L and ``n_dropped`` counts target residues that aligned to no query
    position (insertions relative to the query -- they cannot be
    represented in the match-state frame and are discarded, exactly as an
    a3m would push them into lowercase insert columns).
    """
    q_cols = [j for j, c in enumerate(query_row) if c.isalpha()]
    q_seq = "".join(query_row[j] for j in q_cols)
    target = "".join(c for c in target.upper() if c.isalpha())
    row = ["-"] * len(query_row)
    if target == q_seq:
        for k, j in enumerate(q_cols):
            row[j] = target[k]
        return "".join(row), 0
    from utils.msa_utils import get_align_indexes
    idx_t, idx_q = get_align_indexes(target, q_seq)
    for t, q in zip(idx_t, idx_q):
        row[q_cols[q]] = target[t]
    return "".join(row), len(target) - len(idx_t)


def ungap(s: str) -> str:
    return "".join(c for c in s if c.isalpha()).upper()


def _row_distance(row_a: str, row_b: str, keep: Sequence[bool]) -> int:
    """Mismatch count between two match-state rows over ``keep`` columns."""
    return sum(1 for a, b, k in zip(row_a, row_b, keep)
               if k and a.upper() != b.upper())


# ======================================================================
# Pseudo-pair staging
# ======================================================================

STATIC_DIRS = ("chain_pdb_files", "fasta_chain_files")


def _link_or_copy(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        return
    try:
        dst.symlink_to(src.resolve())
    except OSError:
        import shutil
        shutil.copy2(src, dst)


def stage_pseudo_pair(ref_pair_dir: Path, dest_pair_dir: Path,
                      pair_id: str) -> Path:
    """Create a FoldPairs-shaped dir that borrows the production pair's
    truth structures / chain FASTAs / deep MSA by symlink.

    Nothing is ever written into ``ref_pair_dir``.
    """
    dest_pair_dir.mkdir(parents=True, exist_ok=True)
    for sub in ("output_msa_cluster", "output_get_msa", "Analysis",
                "output_ddg", "output_esmfold2", "tmp_msa_files"):
        (dest_pair_dir / sub).mkdir(parents=True, exist_ok=True)

    for sub in STATIC_DIRS:
        src_dir = ref_pair_dir / sub
        if not src_dir.is_dir():
            continue
        for f in sorted(src_dir.iterdir()):
            if f.is_file():
                _link_or_copy(f, dest_pair_dir / sub / f.name)

    # Root-level truth PDBs: scripts/aggregate_s4pred_into_survey.find_pdb_path
    # looks for <pair_dir>/<pdbid>*.pdb, and score_against_folds falls back to
    # <pair_dir>/<pdb>.pdb / <pdb>_cif.pdb.
    for f in sorted(ref_pair_dir.glob("*.pdb")):
        _link_or_copy(f, dest_pair_dir / f.name)
    for f in sorted(ref_pair_dir.glob("*.cif")):
        _link_or_copy(f, dest_pair_dir / f.name)

    deep = ref_pair_dir / "output_get_msa" / "DeepMsa.a3m"
    if deep.is_file():
        _link_or_copy(deep, dest_pair_dir / "output_get_msa" / "DeepMsa.a3m")
    return dest_pair_dir


def pick_controls(deep_rows: List[Tuple[str, str]], query_row: str,
                  n: int, seed: int, min_cov: float = 0.7
                  ) -> List[Tuple[str, str]]:
    """Random family sequences from the deep MSA, as the negative control.

    Filters to rows with >= ``min_cov`` coverage of the query columns and
    only standard residues (so they are comparable in length/alphabet to
    the ancestors but otherwise arbitrary), then samples ``n`` without
    replacement using ``seed``.
    """
    L_q = sum(1 for c in query_row if c.isalpha())
    pool: List[Tuple[str, str]] = []
    seen = set()
    for hdr, row in deep_rows[1:]:
        if len(row) != len(query_row):
            continue
        s = ungap(row)
        if not s or len(s) < min_cov * L_q:
            continue
        if any(c not in STD_AA for c in s):
            continue
        if s in seen:
            continue
        seen.add(s)
        pool.append((hdr, row))
    rng = random.Random(seed)
    if len(pool) <= n:
        return pool
    return [pool[i] for i in sorted(rng.sample(range(len(pool)), n))]


def check_deep_query(deep_a3m: Path, pair_dir: Path) -> None:
    """Warn loudly if row 1 of the deep MSA is not one of the pair's chains.

    Everything downstream assumes the a3m's match-state frame is defined by
    row 1 (the colabfold query). If it is not the 2n54B / 2hdmA chain, the
    ancestors would be projected onto the wrong frame and every number would
    be quietly wrong, so this is checked rather than assumed.
    """
    rows = read_a3m_match_state(deep_a3m)
    if not rows:
        return
    q_hdr, q_row = rows[0]
    q = ungap(q_row)
    for tag in (F1_TAG, F2_TAG):
        fa = pair_dir / "fasta_chain_files" / f"{tag}.fasta"
        if not fa.is_file():
            continue
        chain = ungap(read_fasta(fa)[0][1]) if read_fasta(fa) else ""
        if chain and chain == q:
            print(f"[anc] DeepMsa query row 1 = {q_hdr} = chain {tag} "
                  f"({len(q)} aa), {len(q_row)} match columns -- as expected")
            return
    print(f"[anc] WARNING: DeepMsa row 1 ({q_hdr}, {len(q)} aa) matches "
          f"neither {F1_TAG} nor {F2_TAG}'s chain FASTA. The ancestors are "
          f"projected onto row 1's match-state frame, so check that frame "
          f"before trusting the numbers.", file=sys.stderr)


def build_cluster_a3ms(deep_a3m: Path, ancestors: List[Tuple[str, str]],
                       n_controls: int, control_seed: int,
                       msa_top_n: int,
                       single_dir: Path, withmsa_dir: Path) -> List[dict]:
    """Write one ``ShallowMsa_NNN.a3m`` per sequence into both trees.

    ``single_dir`` gets a one-row a3m (the sequence alone, in the DeepMsa
    column frame). ``withmsa_dir`` gets the same row 1 plus its
    ``msa_top_n`` nearest deep-MSA neighbours, so the a3m keeps the DeepMsa
    column structure that ``_write_pair_a3m_for_chain`` expects.

    Returns the manifest: one dict per cluster.
    """
    deep_rows = read_a3m_match_state(deep_a3m)
    if not deep_rows:
        raise SystemExit(f"[anc] empty or unreadable DeepMsa: {deep_a3m}")
    query_hdr, query_row = deep_rows[0]
    L = len(query_row)

    targets: List[dict] = []
    for name, seq in ancestors:
        row, dropped = project_to_query_frame(seq, query_row)
        targets.append({"label": name, "kind": "ancestor", "seq": ungap(seq),
                        "row": row, "n_dropped_insertions": dropped,
                        "source": "dishman_fig_S1C"})
    for hdr, row in pick_controls(deep_rows, query_row, n_controls,
                                  control_seed):
        targets.append({"label": hdr, "kind": "control", "seq": ungap(row),
                        "row": row, "n_dropped_insertions": 0,
                        "source": "DeepMsa"})

    single_dir.mkdir(parents=True, exist_ok=True)
    withmsa_dir.mkdir(parents=True, exist_ok=True)
    manifest: List[dict] = []
    for i, t in enumerate(targets):
        tag = f"ShallowMsa_{i:03d}"
        keep = [c.isalpha() for c in t["row"]]
        with open(single_dir / f"{tag}.a3m", "w") as fh:
            fh.write(f">{t['label']}\n{t['row']}\n")
        # nearest neighbours in the deep MSA, measured on the target's own
        # occupied columns (same metric as _write_pair_a3m_for_chain's
        # msa_top_n subsample)
        neigh: List[Tuple[int, int]] = []
        for j, (hdr, row) in enumerate(deep_rows):
            if len(row) != L or ungap(row) == t["seq"]:
                continue
            neigh.append((_row_distance(t["row"], row, keep), j))
        neigh.sort()
        with open(withmsa_dir / f"{tag}.a3m", "w") as fh:
            fh.write(f">{t['label']}\n{t['row']}\n")
            for _, j in neigh[:max(0, msa_top_n)]:
                fh.write(f">{deep_rows[j][0]}\n{deep_rows[j][1]}\n")
        manifest.append({"cluster": tag, "cluster_num": f"{i:03d}",
                         "label": t["label"], "kind": t["kind"],
                         "seq": t["seq"], "seq_len": len(t["seq"]),
                         "n_dropped_insertions": t["n_dropped_insertions"],
                         "n_neighbours": min(len(neigh), max(0, msa_top_n)),
                         "source": t["source"]})
    return manifest


# ======================================================================
# Running the pipeline methods
# ======================================================================

class Runner:
    """Shell runner that can print instead of execute (``--dry_run``)."""

    def __init__(self, dry_run: bool, cwd: Path):
        self.dry_run = dry_run
        self.cwd = cwd
        self.planned: List[str] = []

    def run(self, cmd: str, env: Optional[Dict[str, str]] = None,
            label: str = "") -> int:
        full_env = dict(os.environ)
        if env:
            full_env.update(env)
        prefix = " ".join(f"{k}={shlex.quote(v)}" for k, v in (env or {}).items())
        shown = f"{prefix} {cmd}".strip()
        self.planned.append(shown)
        print(f"[anc:{label or 'cmd'}] {shown}", flush=True)
        if self.dry_run:
            return 0
        return subprocess.call(cmd, shell=True, cwd=str(self.cwd), env=full_env)


def run_esmfold2(rn: Runner, pair_dir: Path, args) -> None:
    inner = (f"export XDG_CACHE_HOME={shlex.quote(args.esmfold2_hf_cache)} && "
             f"export HF_HOME={shlex.quote(args.esmfold2_hf_cache + '/huggingface')} && "
             f"export HF_HUB_OFFLINE=1 && "
             f"source {shlex.quote(args.esmfold2_venv)}/bin/activate && "
             f"python3 run_ESMFold2.py "
             f"--input_dir {shlex.quote(str(pair_dir / 'output_msa_cluster'))} "
             f"--output_dir {shlex.quote(str(pair_dir / 'output_esmfold2'))} "
             f"--representative_method first --device {shlex.quote(args.device)} "
             f"--skip_existing")
    rn.run(f"bash -c {shlex.quote(inner)}", label="esmfold2")


def run_ddg(rn: Runner, data_root: Path, pair_id: str, args) -> None:
    rn.run(f"python3 ./run_DDG.py -input {shlex.quote(pair_id)} "
           f"--thermompnn_dir {shlex.quote(args.thermompnn_dir)} "
           f"--device {shlex.quote(args.device)} --cluster_sample_n 1",
           env={"MSACLUSTER_DATA_DIR": str(data_root)}, label="ddg")


def run_s4pred(rn: Runner, data_root: Path, pair_id: str, work_root: Path,
               args) -> None:
    if not args.s4pred_dir:
        print("[anc:s4pred] skipped: --s4pred_dir / $S4PRED_DIR not set",
              file=sys.stderr)
        return
    rn.run(f"python3 scripts/run_s4pred_per_cluster.py "
           f"--pairs {shlex.quote(pair_id)} "
           f"--fpairs_root {shlex.quote(str(data_root))} "
           f"--s4pred_dir {shlex.quote(args.s4pred_dir)} "
           f"--work_root {shlex.quote(str(work_root / '_s4pred_work'))} "
           f"--n_reps 1 --device {'gpu' if args.device == 'cuda' else 'cpu'} "
           f"--out {shlex.quote(str(work_root / 's4pred_ancestors.csv'))}",
           label="s4pred")


def run_boltz2(rn: Runner, data_root: Path, pair_id: str, args) -> None:
    rn.run(f"python3 ./run_Boltz2.py -input {shlex.quote(pair_id)} --mode apo",
           env={"MSACLUSTER_DATA_DIR": str(data_root)}, label="boltz2")


def run_bioemu(rn: Runner, data_root: Path, pair_id: str, args) -> None:
    rn.run(f"python3 -u ./run_bioemu.py -input {shlex.quote(pair_id)} "
           f"--query_type first --af_msa_top_n {int(args.msa_top_n)} "
           f"--num_samples {int(args.bioemu_num_samples)}",
           env={"MSACLUSTER_DATA_DIR": str(data_root)}, label="bioemu")


def run_bioemu2(rn: Runner, data_root: Path, pair_id: str, args) -> None:
    """BioEmu2 on the same per-ancestor inputs BioEmu1 gets (see run_bioemu.py
    --engine): a separate model, so a separate column in the calibration."""
    rn.run(f"python3 -u ./run_bioemu.py -input {shlex.quote(pair_id)} "
           f"--engine bioemu2 "
           f"--query_type first --af_msa_top_n {int(args.msa_top_n)} "
           f"--num_samples {int(args.bioemu_num_samples)}",
           env={"MSACLUSTER_DATA_DIR": str(data_root)}, label="bioemu2")


def run_af2(rn: Runner, pair_dir: Path, manifest: List[dict], args,
            single_seq: bool) -> None:
    """AF2 per ancestral sequence.

    The per-cluster a3m is built by the pipeline's own
    ``_write_pair_a3m_for_chain`` with ``query_type="first"`` (row 1 of our
    cluster a3m IS the ancestor), then handed to the same
    ``scripts/shell/RunAF2_Colabfold.sh`` the production runs use.
    """
    tag_suffix = "_single" if single_seq else ""
    tmp_dir = pair_dir / f"tmp_msa_files{tag_suffix}"
    deep_a3m = pair_dir / "output_get_msa" / "DeepMsa.a3m"

    if not rn.dry_run:
        tmp_dir.mkdir(parents=True, exist_ok=True)
        from run_foldswitch_pipeline import _write_pair_a3m_for_chain
    for m in manifest:
        cl_a3m = pair_dir / "output_msa_cluster" / f"{m['cluster']}.a3m"
        out_a3m = tmp_dir / f"tmp_{m['cluster']}__{F1_TAG}.a3m"
        if rn.dry_run:
            print(f"[anc:af2] would build {out_a3m} "
                  f"(query_type=first, msa_top_n="
                  f"{0 if single_seq else args.msa_top_n})")
        else:
            ok = _write_pair_a3m_for_chain(
                str(cl_a3m), str(deep_a3m), F1_TAG, str(out_a3m),
                query_type="first",
                msa_top_n=0 if single_seq else int(args.msa_top_n))
            if not ok:
                print(f"[anc:af2] a3m build FAILED for {m['cluster']}",
                      file=sys.stderr)
                continue
            if single_seq:
                # keep row 1 only -> genuine single-sequence AF2
                recs = read_fasta(out_a3m)[:1]
                with open(out_a3m, "w") as fh:
                    for nm, sq in recs:
                        fh.write(f">{nm}\n{sq}\n")
        out_dir = pair_dir / f"output_AF{tag_suffix}" / m["cluster"] / F1_TAG
        rn.run(f"bash ./scripts/shell/RunAF2_Colabfold.sh "
               f"{shlex.quote(str(out_a3m))} {shlex.quote(str(out_dir))} "
               f"--num-models 5 --num-recycle 1 --model-type alphafold2_ptm",
               label=f"af2{tag_suffix}")


# ======================================================================
# Collecting results
# ======================================================================

def truth_paths(pair_dir: Path) -> Tuple[str, str]:
    from utils.score_against_folds import truth_chain_pdb
    f1 = truth_chain_pdb(pair_dir, F1_TAG[:-1], F1_TAG[-1])
    f2 = truth_chain_pdb(pair_dir, F2_TAG[:-1], F2_TAG[-1])
    return f1, f2


def collect_tm_method(pair_dir: Path, manifest: List[dict], method: str,
                      mode: str, pdb_for: Dict[str, Optional[Path]]
                      ) -> List[dict]:
    from utils.score_against_folds import safe_tm_vs_truth
    t1, t2 = truth_paths(pair_dir)
    rows = []
    for m in manifest:
        pred = pdb_for.get(m["cluster"])
        tm1 = safe_tm_vs_truth(pred, t1, F1_TAG[-1]) if pred else float("nan")
        tm2 = safe_tm_vs_truth(pred, t2, F2_TAG[-1]) if pred else float("nan")
        rows.append({"label": m["label"], "kind": m["kind"],
                     "cluster": m["cluster"], "method": method, "mode": mode,
                     "tm_f1": tm1, "tm_f2": tm2,
                     "delta": delta_from_tm(tm1, tm2),
                     "pred_pdb": str(pred) if pred else ""})
    return rows


def collect_esmfold2(pair_dir: Path, manifest: List[dict]) -> List[dict]:
    preds = pair_dir / "output_esmfold2"
    # ESMFold2 emits mmCIF when the chain id is multi-character; the shared
    # converter in scripts/score_esmfold2.py turns those into .pdb in place
    # (cached, idempotent) so TM-align sees what it expects.
    try:
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "_score_esmfold2", REPO_ROOT / "scripts" / "score_esmfold2.py")
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        mod.convert_cif_preds_to_pdb(preds, label="ancestors")
    except Exception as e:
        print(f"[anc:esmfold2] cif->pdb conversion unavailable: {e}",
              file=sys.stderr)
    pdbs: Dict[str, Optional[Path]] = {}
    for m in manifest:
        p = preds / f"{m['cluster']}.pdb"
        if not p.is_file():
            c = preds / f"{m['cluster']}.cif"
            p = c if c.is_file() else None
        pdbs[m["cluster"]] = p
    return collect_tm_method(pair_dir, manifest, "ESM", "single_sequence",
                             pdbs)


def collect_af2(pair_dir: Path, manifest: List[dict], single_seq: bool
                ) -> List[dict]:
    from run_foldswitch_pipeline import _find_best_af2_pdb
    suffix = "_single" if single_seq else ""
    pdbs: Dict[str, Optional[Path]] = {}
    for m in manifest:
        d = pair_dir / f"output_AF{suffix}" / m["cluster"] / F1_TAG
        best = _find_best_af2_pdb(str(d)) if d.is_dir() else None
        pdbs[m["cluster"]] = Path(best) if best else None
    mode = "single_sequence" if single_seq else "deepmsa_topN"
    return collect_tm_method(pair_dir, manifest, "AF2", mode, pdbs)


def collect_ddg(pair_dir: Path, manifest: List[dict]) -> List[dict]:
    import pandas as pd
    csv_path = pair_dir / "Analysis" / "df_ddg.csv"
    if not csv_path.is_file():
        return []
    df = pd.read_csv(csv_path, dtype={"cluster": str})
    by = {str(r["cluster"]).zfill(3): r for _, r in df.iterrows()}
    rows = []
    for m in manifest:
        r = by.get(m["cluster_num"])
        if r is None:
            continue
        s1 = float(r.get("sum_F1_mean", float("nan")))
        s2 = float(r.get("sum_F2_mean", float("nan")))
        delta = delta_from_ddg(s1, s2)
        # run_DDG already reports exactly this as bias_mean; cross-check.
        bias = float(r.get("bias_mean", float("nan")))
        if np.isfinite(bias) and np.isfinite(delta) and abs(bias - delta) > 1e-3:
            print(f"[anc:ddg] WARN {m['cluster']}: bias_mean={bias:.4f} != "
                  f"sum_F2-sum_F1={delta:.4f}", file=sys.stderr)
        rows.append({"label": m["label"], "kind": m["kind"],
                     "cluster": m["cluster"], "method": "DDG",
                     "mode": "single_sequence",
                     "ddg_sum_f1": s1, "ddg_sum_f2": s2,
                     "ddg_bias_kcal_mol": bias,
                     "delta": delta if np.isfinite(delta) else bias})
    return rows


def collect_s4pred(pair_dir: Path, manifest: List[dict]) -> List[dict]:
    import importlib.util
    import pandas as pd
    csv_path = pair_dir / "Analysis" / "df_s4pred.csv"
    if not csv_path.is_file():
        return []
    spec = importlib.util.spec_from_file_location(
        "_agg_s4pred", REPO_ROOT / "scripts" / "aggregate_s4pred_into_survey.py")
    agg = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(agg)

    t1, t2 = truth_paths(pair_dir)
    h1, e1, _, n1 = agg.true_ss_from_dssp(t1, F1_TAG[-1])
    h2, e2, _, n2 = agg.true_ss_from_dssp(t2, F2_TAG[-1])
    if not n1 or not n2:
        print("[anc:s4pred] DSSP failed on a truth fold; no rows",
              file=sys.stderr)
        return []
    df = pd.read_csv(csv_path)
    by = {str(r["cluster"]): r for _, r in df.iterrows()}
    rows = []
    for m in manifest:
        r = by.get(m["cluster"])
        if r is None:
            continue
        h, e = float(r.get("frac_H", np.nan)), float(r.get("frac_E", np.nan))
        s1 = agg.he_similarity(h, e, h1, e1)
        s2 = agg.he_similarity(h, e, h2, e2)
        rows.append({"label": m["label"], "kind": m["kind"],
                     "cluster": m["cluster"], "method": "S4PRED",
                     "mode": "single_sequence",
                     "ss_frac_H": h, "ss_frac_E": e,
                     "ss_sim_f1": s1, "ss_sim_f2": s2,
                     "delta": delta_from_ss(s1, s2)})
    return rows


def _collect_simple_csv(pair_dir: Path, manifest: List[dict], fname: str,
                        method: str, mode: str, c1: str, c2: str) -> List[dict]:
    import pandas as pd
    csv_path = pair_dir / "Analysis" / fname
    if not csv_path.is_file():
        return []
    df = pd.read_csv(csv_path, dtype={"cluster": str})
    by = {str(r["cluster"]): r for _, r in df.iterrows()}
    rows = []
    for m in manifest:
        # `or` on a pandas Series raises; pick explicitly.
        r = by.get(m["cluster"])
        if r is None:
            r = by.get(m["cluster_num"])
        if r is None:
            continue
        tm1, tm2 = float(r.get(c1, np.nan)), float(r.get(c2, np.nan))
        rows.append({"label": m["label"], "kind": m["kind"],
                     "cluster": m["cluster"], "method": method, "mode": mode,
                     "tm_f1": tm1, "tm_f2": tm2,
                     "delta": delta_from_tm(tm1, tm2)})
    return rows


def collect_boltz2(pair_dir, manifest):
    return _collect_simple_csv(pair_dir, manifest, "df_boltz2.csv", "Boltz2",
                               "single_sequence", "rank1_tm_fold1",
                               "rank1_tm_fold2")


def collect_bioemu(pair_dir, manifest):
    return _collect_simple_csv(pair_dir, manifest, "df_bioemu.csv", "BioEmu",
                               "deepmsa_topN", "mean_tm_fold1",
                               "mean_tm_fold2")


def collect_bioemu2(pair_dir, manifest):
    return _collect_simple_csv(pair_dir, manifest, "df_bioemu2.csv", "BioEmu2",
                               "deepmsa_topN", "mean_tm_fold1",
                               "mean_tm_fold2")


# ======================================================================
# Statistics + figure
# ======================================================================

def method_stats(rows: List[dict], n_control_draws: int = 2000,
                 seed: int = 0) -> List[dict]:
    """Per (method, mode): Spearman rho vs the measured % alternative fold,
    exact permutation p, and the control null.

    The control null answers "could any 5 family sequences have done this?":
    we draw 5 of the ``--n_controls`` random XCL-family sequences, pair them
    with the same experimental vector in a random order, and record rho.
    ``p_control`` is the fraction of draws reaching the ancestors' rho.
    """
    truth = {g["name"]: g["pct_alt"] for g in GROUND_TRUTH}
    y = [truth[n] for n in ANCESTOR_NAMES]
    rng = random.Random(seed)
    out = []
    keys = []
    for r in rows:
        k = (r["method"], r["mode"])
        if k not in keys:
            keys.append(k)
    for method, mode in keys:
        sub = [r for r in rows if r["method"] == method and r["mode"] == mode]
        by_label = {r["label"]: r["delta"] for r in sub}
        x = [by_label.get(n, float("nan")) for n in ANCESTOR_NAMES]
        ctrl = [r["delta"] for r in sub if r["kind"] == "control"
                and np.isfinite(r["delta"])]
        n_ok = int(sum(1 for v in x if np.isfinite(v)))
        if n_ok < len(ANCESTOR_NAMES):
            out.append({"method": method, "mode": mode, "n_ancestors": n_ok,
                        "rho": float("nan"), "p_exact_greater": float("nan"),
                        "p_exact_two_sided": float("nan"), "n_perms": 0,
                        "n_controls": len(ctrl),
                        "control_delta_mean": float(np.mean(ctrl)) if ctrl else float("nan"),
                        "control_delta_sd": float(np.std(ctrl)) if ctrl else float("nan"),
                        "p_control": float("nan"),
                        "deltas": ";".join(f"{v:.4f}" for v in x)})
            continue
        rho = spearman_rho(x, y)
        p_g, n_perm = exact_perm_p(x, y, "greater")
        p_t, _ = exact_perm_p(x, y, "two-sided")
        p_ctrl = float("nan")
        ctrl_rhos: List[float] = []
        if len(ctrl) >= 5:
            hits = 0
            for _ in range(n_control_draws):
                s = rng.sample(ctrl, 5)
                rng.shuffle(s)
                r_c = spearman_rho(s, y)
                if np.isfinite(r_c):
                    ctrl_rhos.append(r_c)
                if np.isfinite(r_c) and r_c >= rho - 1e-12:
                    hits += 1
            p_ctrl = hits / n_control_draws
        out.append({"control_rhos": ctrl_rhos,
                    "method": method, "mode": mode, "n_ancestors": n_ok,
                    "rho": rho, "p_exact_greater": p_g,
                    "p_exact_two_sided": p_t, "n_perms": n_perm,
                    "n_controls": len(ctrl),
                    "control_delta_mean": float(np.mean(ctrl)) if ctrl else float("nan"),
                    "control_delta_sd": float(np.std(ctrl)) if ctrl else float("nan"),
                    "p_control": p_ctrl,
                    "deltas": ";".join(f"{v:.4f}" for v in x)})
    return out


def _method_label(method: str, mode: str) -> str:
    """Display label, using the repo-wide method names so this figure reads
    the same as every other cross-method figure."""
    try:
        from utils.method_labels import disp
        name = disp(method)
    except Exception:
        name = method
    return name if mode == "single_sequence" else f"{name} ({mode})"


def make_figure(rows: List[dict], stats: List[dict], out_png: Path) -> None:
    """Two panels.

    Left: the trajectory. x = ancestor in evolutionary order, y = each
    method's Delta z-scored over the five ancestors (so methods on
    TM-score, kcal/mol and SS-similarity scales are comparable), with the
    measured % alternative fold on the right-hand axis.

    Right: the negative control. For each method, the ancestors' Spearman
    rho against the distribution of rho obtained by drawing five random
    XCL-family sequences from the deep MSA in a random order. A real
    signal sits far to the right of its own control cloud.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    truth = {g["name"]: g for g in GROUND_TRUTH}
    xs = np.arange(len(ANCESTOR_NAMES))

    keys: List[Tuple[str, str]] = []
    for r in rows:
        k = (r["method"], r["mode"])
        if k not in keys and r["kind"] == "ancestor":
            keys.append(k)
    try:
        from utils.method_labels import order_methods
        rank = {m: i for i, m in enumerate(
            order_methods({k[0] for k in keys}))}
        keys.sort(key=lambda k: (rank.get(k[0], 99), k[1]))
    except Exception:
        pass

    has_ctrl = any(s.get("control_rhos") for s in stats)
    fig, axes = plt.subplots(
        1, 2 if has_ctrl else 1,
        figsize=((11.4, 5.2) if has_ctrl else (8.0, 5.2)),
        gridspec_kw={"width_ratios": [2.1, 1.0]} if has_ctrl else None)
    ax = axes[0] if has_ctrl else axes
    axc = axes[1] if has_ctrl else None

    cmap = plt.get_cmap("tab10")
    colour = {k: cmap(i % 10) for i, k in enumerate(keys)}
    plotted = 0
    for k in keys:
        method, mode = k
        by_label = {r["label"]: r["delta"] for r in rows
                    if r["method"] == method and r["mode"] == mode}
        y = np.array([by_label.get(n, np.nan) for n in ANCESTOR_NAMES],
                     dtype=float)
        if not np.isfinite(y).any():
            continue
        z = zscore(y)
        st = next((s for s in stats
                   if s["method"] == method and s["mode"] == mode), None)
        lbl = _method_label(method, mode)
        if st and np.isfinite(st["rho"]):
            lbl += (f"   $\\rho$={st['rho']:.2f}, "
                    f"p={st['p_exact_greater']:.3f}")
        ax.plot(xs, z, "o-", color=colour[k], lw=1.8, ms=6, label=lbl,
                zorder=3)
        plotted += 1

    ax.axhline(0.0, color="0.78", lw=0.8, zorder=1)
    ax.set_xticks(xs)
    ax.set_xticklabels(ANCESTOR_NAMES)
    ax.set_xlabel("reconstructed ancestor (evolutionary order; "
                  "XCL1 = modern tip)")
    ax.set_ylabel("predicted preference for the alternative fold\n"
                  r"$\Delta$ = score(2n54B) $-$ score(2hdmA),"
                  "\nz-scored over the five ancestors")

    ax2 = ax.twinx()
    pct = np.array([truth[n]["pct_alt"] for n in ANCESTOR_NAMES], dtype=float)
    err = np.array([truth[n]["pct_alt_sd"] for n in ANCESTOR_NAMES],
                   dtype=float)
    err = np.where(np.isfinite(err), err, 0.0)
    ax2.errorbar(xs, pct, yerr=err, color="black", lw=2.4, ls="--",
                 marker="s", ms=7, capsize=4, zorder=2,
                 label="measured % alternative fold (NMR)")
    ax2.set_ylabel("measured % alternative fold (NMR, Dishman 2021)")
    ax2.set_ylim(-8, 108)

    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=7.5, loc="upper left",
              framealpha=0.92)
    ax.set_title("Measured XCL1 ancestral trajectory\nvs our fold-preference "
                 "predictors", fontsize=11)

    if axc is not None:
        shown = [k for k in keys
                 if any(s["method"] == k[0] and s["mode"] == k[1]
                        and s.get("control_rhos") for s in stats)]
        for i, k in enumerate(shown):
            st = next(s for s in stats
                      if s["method"] == k[0] and s["mode"] == k[1])
            cr = np.asarray(st["control_rhos"], dtype=float)
            jitter = (np.random.default_rng(0).normal(0, 0.055, size=len(cr)))
            axc.scatter(cr, np.full(len(cr), i) + jitter, s=3, alpha=0.12,
                        color="0.45", linewidths=0, zorder=2)
            axc.plot([np.percentile(cr, 2.5), np.percentile(cr, 97.5)],
                     [i, i], color="0.3", lw=1.6, zorder=3)
            if np.isfinite(st["rho"]):
                axc.plot([st["rho"]], [i], marker="D", ms=8,
                         color=colour[k], mec="black", mew=0.8, zorder=4)
        axc.set_yticks(range(len(shown)))
        axc.set_yticklabels([_method_label(*k) for k in shown], fontsize=8)
        axc.set_ylim(-0.7, len(shown) - 0.3)
        axc.set_xlim(-1.08, 1.08)
        axc.axvline(0.0, color="0.78", lw=0.8)
        axc.set_xlabel(r"Spearman $\rho$ vs measured % alternative fold")
        axc.set_title("Negative control:\n5 random XCL-family sequences",
                      fontsize=10)
        from matplotlib.lines import Line2D
        axc.legend(handles=[
            Line2D([], [], marker="D", ls="", color="0.2", ms=7,
                   label="ancestors (observed)"),
            Line2D([], [], marker="o", ls="", color="0.45", ms=4,
                   label="random 5-sequence draws"),
            Line2D([], [], color="0.3", lw=1.6, label="control 95 % range"),
        ], fontsize=7, loc="lower left", framealpha=0.92)

    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=200)
    fig.savefig(out_png.with_suffix(".pdf"))
    plt.close(fig)
    print(f"[anc] wrote {out_png} ({plotted} method curves)")


# ======================================================================
# CLI
# ======================================================================

CSV_FIELDS = ["label", "kind", "cluster", "method", "mode", "evo_order",
              "pct_alt_measured", "pct_alt_measured_sd", "seq_len",
              "n_dropped_insertions",
              "tm_f1", "tm_f2", "delta", "delta_z",
              "ddg_sum_f1", "ddg_sum_f2", "ddg_bias_kcal_mol",
              "ss_frac_H", "ss_frac_E", "ss_sim_f1", "ss_sim_f2",
              "pred_pdb"]

STATS_FIELDS = ["method", "mode", "n_ancestors", "rho", "p_exact_greater",
                "p_exact_two_sided", "n_perms", "n_controls",
                "control_delta_mean", "control_delta_sd", "p_control",
                "deltas"]


def parse_methods(s: str) -> List[str]:
    if s.strip().lower() in ("all", ""):
        return list(ALL_METHODS)
    out = []
    for tok in s.replace(",", " ").split():
        tok = tok.strip().lower()
        if tok not in ALL_METHODS:
            raise SystemExit(f"[anc] unknown method {tok!r}; "
                             f"choose from {','.join(ALL_METHODS)} or 'all'")
        out.append(tok)
    return out


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ancestors_fasta",
                    default=str(REPO_ROOT / "data" / "ancestors" /
                                "dishman_ancestors.fasta"))
    ap.add_argument("--pair_id", default=DEFAULT_PAIR,
                    help="Reference fold-switch pair (F1_F2).")
    ap.add_argument("--ref_pair_dir", default="",
                    help="Production pair dir to borrow truth PDBs / chain "
                         "FASTAs / DeepMsa.a3m from (read-only). "
                         "Default: $MSACLUSTER_DATA_DIR or "
                         "Pipeline/FoldPairs, / <pair_id>.")
    ap.add_argument("--work_root", default="Pipeline/AncestorCalibration",
                    help="Where the two pseudo-pair trees are staged.")
    ap.add_argument("--methods", default="esmfold2,ddg,s4pred",
                    help=f"Comma list from {','.join(ALL_METHODS)}, or 'all'.")
    ap.add_argument("--n_controls", type=int, default=20,
                    help="Random XCL-family sequences from the deep MSA used "
                         "as the negative control.")
    ap.add_argument("--control_seed", type=int, default=0)
    ap.add_argument("--msa_top_n", type=int, default=10,
                    help="Neighbours per sequence in the 'withmsa' tree "
                         "(AF2 / BioEmu). Matches --af_msa_top_n.")
    ap.add_argument("--device", default="cuda", choices=["cuda", "cpu"])
    ap.add_argument("--thermompnn_dir",
                    default=os.environ.get(
                        "THERMOMPNN_DIR",
                        "/sci/labs/orzuk/orzuk/github/ThermoMPNN"))
    ap.add_argument("--s4pred_dir", default=os.environ.get("S4PRED_DIR", ""))
    ap.add_argument("--esmfold2_venv",
                    default="/sci/labs/orzuk/orzuk/venvs/esmfold2-venv")
    ap.add_argument("--esmfold2_hf_cache",
                    default="/sci/labs/orzuk/orzuk/tmp/xdg-cache")
    ap.add_argument("--bioemu_num_samples", type=int, default=100)
    ap.add_argument("--out", default="docs/ancestor_calibration.csv")
    ap.add_argument("--stats_out", default="docs/ancestor_calibration_stats.csv")
    ap.add_argument("--fig", default="docs/figures/global/ancestor_calibration.png")
    ap.add_argument("--dry_run", action="store_true",
                    help="Stage nothing, run nothing: print the planned "
                         "commands and exit.")
    ap.add_argument("--stage_only", action="store_true",
                    help="Build the pseudo-pair trees and the per-ancestor "
                         "a3ms, then stop.")
    ap.add_argument("--collect_only", action="store_true",
                    help="Skip every runner; just re-score / re-plot from "
                         "outputs already on disk.")
    ap.add_argument("--no_figure", action="store_true")
    args = ap.parse_args(argv)

    methods = parse_methods(args.methods)
    work_root = Path(args.work_root)
    if not work_root.is_absolute():
        work_root = REPO_ROOT / work_root
    single_root = work_root / "single" / "FoldPairs"
    withmsa_root = work_root / "withmsa" / "FoldPairs"
    single_dir = single_root / args.pair_id
    withmsa_dir = withmsa_root / args.pair_id

    if args.ref_pair_dir:
        ref_pair_dir = Path(args.ref_pair_dir)
    else:
        base = os.environ.get("MSACLUSTER_DATA_DIR") or str(
            REPO_ROOT / "Pipeline" / "FoldPairs")
        ref_pair_dir = Path(base) / args.pair_id

    ancestors = read_fasta(args.ancestors_fasta)
    missing = [n for n in ANCESTOR_NAMES
               if n not in {a for a, _ in ancestors}]
    if missing:
        raise SystemExit(f"[anc] {args.ancestors_fasta} is missing "
                         f"{missing}; expected {ANCESTOR_NAMES}")
    # keep GROUND_TRUTH's evolutionary order
    by_name = dict(ancestors)
    ancestors = [(n, by_name[n]) for n in ANCESTOR_NAMES]

    rn = Runner(args.dry_run, REPO_ROOT)
    print(f"[anc] pair={args.pair_id}  F1={F1_TAG} (alternative fold)  "
          f"F2={F2_TAG} (chemokine fold)")
    print(f"[anc] reference pair dir (read-only): {ref_pair_dir}")
    print(f"[anc] work root: {work_root}")
    print(f"[anc] methods: {', '.join(methods)}")

    manifest_path = work_root / "ancestor_manifest.json"
    manifest: List[dict]

    if args.dry_run:
        manifest = [{"cluster": f"ShallowMsa_{i:03d}",
                     "cluster_num": f"{i:03d}", "label": n,
                     "kind": "ancestor", "seq": s, "seq_len": len(s),
                     "n_dropped_insertions": 0, "n_neighbours": args.msa_top_n,
                     "source": "dishman_fig_S1C"}
                    for i, (n, s) in enumerate(ancestors)]
        print(f"[anc] (dry run) + {args.n_controls} control clusters "
              f"ShallowMsa_{len(ancestors):03d}.."
              f"ShallowMsa_{len(ancestors) + args.n_controls - 1:03d}")
    elif args.collect_only:
        if not manifest_path.is_file():
            raise SystemExit(f"[anc] --collect_only needs {manifest_path}; "
                             f"run the staging step first")
        manifest = json.loads(manifest_path.read_text())
    else:
        if not ref_pair_dir.is_dir():
            raise SystemExit(f"[anc] reference pair dir not found: "
                             f"{ref_pair_dir} (run this on the cluster, or "
                             f"pass --ref_pair_dir)")
        for dest in (single_dir, withmsa_dir):
            stage_pseudo_pair(ref_pair_dir, dest, args.pair_id)
        deep = single_dir / "output_get_msa" / "DeepMsa.a3m"
        if not deep.is_file():
            raise SystemExit(f"[anc] DeepMsa.a3m not found under "
                             f"{ref_pair_dir / 'output_get_msa'}")
        check_deep_query(deep, single_dir)
        manifest = build_cluster_a3ms(
            deep, ancestors, args.n_controls, args.control_seed,
            args.msa_top_n,
            single_dir / "output_msa_cluster",
            withmsa_dir / "output_msa_cluster")
        manifest_path.write_text(json.dumps(manifest, indent=2))
        n_anc = sum(1 for m in manifest if m["kind"] == "ancestor")
        print(f"[anc] staged {len(manifest)} clusters "
              f"({n_anc} ancestors + {len(manifest) - n_anc} controls) "
              f"-> {manifest_path}")
        for m in manifest:
            if m["kind"] == "ancestor" and m["n_dropped_insertions"]:
                print(f"[anc] WARN {m['label']}: "
                      f"{m['n_dropped_insertions']} residue(s) aligned to no "
                      f"DeepMsa match column and were dropped")
        if args.stage_only:
            return 0

    if not args.collect_only:
        if "esmfold2" in methods:
            run_esmfold2(rn, single_dir, args)
        if "ddg" in methods:
            run_ddg(rn, single_root, args.pair_id, args)
        if "s4pred" in methods:
            run_s4pred(rn, single_root, args.pair_id, work_root, args)
        if "boltz2" in methods:
            run_boltz2(rn, single_root, args.pair_id, args)
        if "bioemu" in methods:
            run_bioemu(rn, withmsa_root, args.pair_id, args)
        if "bioemu2" in methods:
            run_bioemu2(rn, withmsa_root, args.pair_id, args)
        if "af2_single" in methods:
            run_af2(rn, single_dir, manifest, args, single_seq=True)
        if "af2" in methods:
            run_af2(rn, withmsa_dir, manifest, args, single_seq=False)

    if args.dry_run:
        print(f"\n[anc] dry run: {len(rn.planned)} command(s) planned; "
              f"nothing staged or executed.")
        return 0

    rows: List[dict] = []
    if "esmfold2" in methods:
        rows += collect_esmfold2(single_dir, manifest)
    if "ddg" in methods:
        rows += collect_ddg(single_dir, manifest)
    if "s4pred" in methods:
        rows += collect_s4pred(single_dir, manifest)
    if "boltz2" in methods:
        rows += collect_boltz2(single_dir, manifest)
    if "bioemu" in methods:
        rows += collect_bioemu(withmsa_dir, manifest)
    if "bioemu2" in methods:
        rows += collect_bioemu2(withmsa_dir, manifest)
    if "af2_single" in methods:
        rows += collect_af2(single_dir, manifest, single_seq=True)
    if "af2" in methods:
        rows += collect_af2(withmsa_dir, manifest, single_seq=False)

    if not rows:
        print("[anc] no results collected (did the runners produce output?)",
              file=sys.stderr)
        return 2

    truth = {g["name"]: g for g in GROUND_TRUTH}
    len_by_label = {m["label"]: m["seq_len"] for m in manifest}
    drop_by_label = {m["label"]: m["n_dropped_insertions"] for m in manifest}
    # per-method z-score, over the 5 ancestors only (so controls are scored
    # on the ancestors' scale and stay comparable)
    for method, mode in {(r["method"], r["mode"]) for r in rows}:
        sub = [r for r in rows if r["method"] == method and r["mode"] == mode]
        anc = [r["delta"] for r in sub if r["kind"] == "ancestor"
               and np.isfinite(r["delta"])]
        mu = float(np.mean(anc)) if anc else float("nan")
        sd = float(np.std(anc)) if anc else float("nan")
        for r in sub:
            r["delta_z"] = ((r["delta"] - mu) / sd
                            if (sd and np.isfinite(sd) and sd > 0
                                and np.isfinite(r["delta"])) else float("nan"))
    for r in rows:
        g = truth.get(r["label"])
        r["evo_order"] = g["evo_order"] if g else ""
        r["pct_alt_measured"] = g["pct_alt"] if g else ""
        r["pct_alt_measured_sd"] = (g["pct_alt_sd"] if g and
                                    np.isfinite(g["pct_alt_sd"]) else "")
        r["seq_len"] = len_by_label.get(r["label"], "")
        # Residues of the ancestor that aligned to no DeepMsa match column.
        # Every per-cluster reader in the pipeline works in the match-state
        # frame, so these residues are invisible to all methods -- recorded
        # here so the effect is auditable rather than silent.
        r["n_dropped_insertions"] = drop_by_label.get(r["label"], "")

    rows.sort(key=lambda r: (r["method"], r["mode"], r["kind"] != "ancestor",
                             r.get("evo_order", 99) if r.get("evo_order") != "" else 99,
                             r["label"]))
    out_csv = Path(args.out)
    if not out_csv.is_absolute():
        out_csv = REPO_ROOT / out_csv
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_FIELDS, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"[anc] wrote {out_csv} ({len(rows)} rows)")

    stats = method_stats(rows, seed=args.control_seed)
    stats_csv = Path(args.stats_out)
    if not stats_csv.is_absolute():
        stats_csv = REPO_ROOT / stats_csv
    with open(stats_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=STATS_FIELDS,
                           extrasaction="ignore")   # drops control_rhos
        w.writeheader()
        w.writerows(stats)
    print(f"[anc] wrote {stats_csv}")

    print("\n  method            mode              rho    p_exact(>)  "
          "p_two_sided  p_control  n_ctrl")
    for s in stats:
        print(f"  {s['method']:<16} {s['mode']:<16} "
              f"{s['rho']:>6.3f}  {s['p_exact_greater']:>9.4f}  "
              f"{s['p_exact_two_sided']:>10.4f}  "
              f"{s['p_control']:>9.4f}  {s['n_controls']:>5d}")
    print(f"\n  exact p enumerates all {math.factorial(len(ANCESTOR_NAMES))} "
          f"orderings of {len(ANCESTOR_NAMES)} points; the attainable floor is "
          f"2/120 = 0.0167 because Anc.0 and Anc.2 tie at 0 % alternative fold.")

    if not args.no_figure:
        fig_path = Path(args.fig)
        if not fig_path.is_absolute():
            fig_path = REPO_ROOT / fig_path
        make_figure(rows, stats, fig_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
