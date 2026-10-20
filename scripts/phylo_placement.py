#!/usr/bin/env python3
"""Test whether per-cluster fold preference (F1c/F2c/Amc) is phylogenetically
structured for each fold-switching pair, on the COARSE cluster-level tree.

For every pair p and every method m in {AF2, ESM, MSAT, DDG}:

  1. Build a COARSE cluster-level tree with one leaf per ShallowMsa_NNN
     cluster (plus DeepMsa, when present), via
     ``utils.ancestral_utils.build_coarse_tree``. That helper
       a. loads the per-sequence Newick at
          Pipeline/FoldPairs/<p>/output_phytree/DeepMsa_tree.nwk,
       b. assigns each leaf to its cluster via the per-cluster a3m files
          at Pipeline/FoldPairs/<p>/output_msa_cluster/ShallowMsa_*.a3m,
       c. picks one representative leaf per cluster and prunes the tree,
       d. renames each surviving leaf 'S00', 'S01', ..., 'D' (DeepMsa).
  2. Read each (pair, method, cluster)'s centered preference from
     docs/fold_diversity_survey.csv (column ``pref_centered``, with a
     fallback to recomputing it on the fly if the CSV lacks that column).
  3. Compute three phylogenetic-structure statistics on the F1c vs F2c
     binary partition of the COARSE tree (Amc leaves treated as missing):
       a. Fitch parsimony score (count of state changes along edges)
       b. k=1 / k=3 nearest-neighbour concordance via cophenetic distance
       c. Phylogenetic D statistic (Fritz & Purvis 2010), if computable
     Each statistic gets an empirical permutation-null p-value. The null
     shuffles the F1c/F2c labels among the cluster-level leaves; the
     coarse-tree shuffle is more conservative than a sequence-level shuffle
     because the number of "free" labels equals the number of clusters
     (typically 8-15), not the number of sequences (hundreds).
  4. Output one row per (pair, method) to docs/phylo_placement.csv.

Per-pair figures show the COARSE tree (10-20 nodes) with each leaf
coloured by F1c / F2c / Amc.

WHY COARSE?
-----------
The fold-preference label is defined PER CLUSTER, so the natural unit of
analysis is the cluster, not the sequence. Running the same statistics on
a fine per-sequence tree
   * makes within-cluster leaves trivially share a label, inflating
     NN-concordance toward 1.0;
   * leaves Fitch parsimony invariant (a constant subtree contributes 0),
     but makes the permutation null too easy to reject because most
     "shuffles" of sequence-level labels still preserve the cluster
     structure that the survey already imposed.
The coarse tree fixes both issues and produces a much cleaner figure.

Usage
-----
    python3 scripts/phylo_placement.py
    python3 scripts/phylo_placement.py --pairs 2namA_1uxmK,2n54B_2hdmA
    python3 scripts/phylo_placement.py --pairs ALL --make-figures all
    # add the clade-shift scan (per method + COMBINED row per pair):
    python3 scripts/phylo_placement.py --pairs ALL \
        --tests morans,clade_shift,clade_shift_combined --clade-n-perm 10000

Tests live in utils/phylo_tests.py (registry ``TESTS``); ``--tests`` picks
them. The default reproduces the pre-registry output exactly (Moran's I, plus
the discrete Fitch/NN/D tests with --discrete-perms).
"""

from __future__ import annotations

import os
# Headless rendering env vars MUST be set BEFORE importing matplotlib / ete3.
# Without these, on a headless compute node ete3.TreeStyle hangs waiting on Qt
# and matplotlib fails to find a writable config dir.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-cache-orzuk")
os.environ.setdefault("PYTHONUNBUFFERED", "1")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)

import argparse
import random
import re
import signal
import sys
import time
from collections import Counter
from contextlib import contextmanager
from typing import Dict, List, Optional, Sequence, Tuple


@contextmanager
def _per_pair_timeout(seconds: int):
    """SIGALRM-based per-pair timeout. Skips a single hung pair without
    killing the whole run. Linux-only; uses SIGALRM."""
    def _handler(signum, frame):
        raise TimeoutError(f"per-pair timeout after {seconds}s")
    old = signal.signal(signal.SIGALRM, _handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)

import numpy as np
import pandas as pd

# Allow running from repo root and from anywhere
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

print(f"[phylo] importing modules...", flush=True)
_t0 = time.time()
from config import DATA_DIR, PAIR_DIR_RE, TABLES_RES  # noqa: E402
from utils.ancestral_utils import build_coarse_tree, _name_internal_nodes  # noqa: E402
print(f"[phylo] imports done in {time.time()-_t0:.1f}s", flush=True)


# Methods we score. Methods missing from the survey CSV are silently
# skipped (n_pairs=0). Add new methods here when fold_diversity_survey.py
# starts emitting them.
SCORED_METHODS: Tuple[str, ...] = (
    "AF2", "AF3", "Boltz2",          # MSA-conditioned structure predictors
    "BioEmu",                         # generative equilibrium-ensemble (9th method)
    "ESM",                            # single-sequence ESMFold (v1) or ESMFold2 (v2)
    "MSAT", "CCMpred",                # contact-map predictors
    "DDG",                            # Rosetta ΔΔG per cluster
    "S4PRED",                         # secondary-structure-based fold preference
)


# ---------------------------------------------------------------------------
# Cluster-tag normalization (matches scripts/fold_diversity_survey.py)
# ---------------------------------------------------------------------------

_MSA_TAG_RE = re.compile(r"(?:Shallow)?Msa_(\d+)", re.IGNORECASE)


def normalize_cluster(s: str) -> str:
    """Normalize a cluster tag to either 'ShallowMsa_NNN' or 'DeepMsa'."""
    s = str(s).strip()
    s = re.sub(r"^(msa_t_)?(clusters_|cmap_)?", "", s)
    if s.startswith("DeepMsa") or s.lower() in (
        "deep", "deepmsa", "msa_deep", "deep_msa", "query"
    ):
        return "DeepMsa"
    m = _MSA_TAG_RE.fullmatch(s)
    if m:
        return f"ShallowMsa_{int(m.group(1)):03d}"
    m2 = re.search(r"ShallowMsa_(\d+)", s)
    if m2:
        return f"ShallowMsa_{int(m2.group(1)):03d}"
    if s.isdigit():
        return f"ShallowMsa_{int(s):03d}"
    return s


# ---------------------------------------------------------------------------
# Fold-preference ingestion (from fold_diversity_survey.csv)
# ---------------------------------------------------------------------------

def bh_adjust(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (NaN-safe, monotone)."""
    p = np.asarray(pvals, dtype=float)
    mask = np.isfinite(p)
    out = np.full_like(p, np.nan, dtype=float)
    if not mask.any():
        return out
    pv = p[mask]
    n = len(pv)
    order = np.argsort(pv)
    ranks = np.arange(1, n + 1)
    adj_sorted = pv[order] * n / ranks
    adj_sorted = np.minimum.accumulate(adj_sorted[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.clip(adj_sorted, 0.0, 1.0)
    out[mask] = adj
    return out


def load_corrected_pref(survey_csv: str) -> Dict[Tuple[str, str, str], str]:
    """Build (pair, method, cluster) -> residual-corrected preference label.

    Reads the 'pref_corrected' column produced by
    scripts/seq_divergence_correction.py from
    docs/fold_diversity_survey_corrected.csv.
    Used when --labels corrected is selected; corresponds to F1c/F2c/Amc
    labels computed from the regression residuals (PID-to-query stripped).
    """
    df = pd.read_csv(survey_csv)
    if "pref_corrected" not in df.columns:
        raise SystemExit(
            f"Expected column 'pref_corrected' in {survey_csv}. "
            "Run scripts/seq_divergence_correction.py first."
        )
    cluster_col = "cluster_norm" if "cluster_norm" in df.columns else "cluster"
    df[cluster_col] = df[cluster_col].astype(str).apply(normalize_cluster)

    out: Dict[Tuple[str, str, str], str] = {}
    for _, r in df.iterrows():
        label = r.get("pref_corrected")
        if not isinstance(label, str) or label not in ("F1", "F2", "Amb"):
            label = "Amb"
        key = (str(r["pair_id"]), str(r["method"]), str(r[cluster_col]))
        out[key] = label
    return out


def load_fold_value(survey_csv: str, value_col: str
                    ) -> Dict[Tuple[str, str, str], Optional[float]]:
    """Build (pair, method, cluster) -> CONTINUOUS signed fold preference.

    Used by the Moran's-I phylogenetic-signal test (no thresholding). The
    value is the method's native signed preference; Moran's I is invariant
    to per-method affine rescaling, so no normalization is needed.

      raw labels       -> value_col = 'TMdiff_centered' (median-centred)
      corrected labels -> value_col = 'TMdiff_residual' (seq-divergence
                          regression residual; produced by
                          scripts/seq_divergence_correction.py)
    """
    df = pd.read_csv(survey_csv)
    if value_col not in df.columns:
        raise SystemExit(
            f"Expected continuous column '{value_col}' in {survey_csv}."
        )
    cluster_col = "cluster_norm" if "cluster_norm" in df.columns else "cluster"
    df[cluster_col] = df[cluster_col].astype(str).apply(normalize_cluster)
    out: Dict[Tuple[str, str, str], Optional[float]] = {}
    for _, r in df.iterrows():
        v = pd.to_numeric(r.get(value_col), errors="coerce")
        key = (str(r["pair_id"]), str(r["method"]), str(r[cluster_col]))
        out[key] = float(v) if pd.notna(v) else None
    return out


def load_fold_pref(survey_csv: str, delta_tm: float = 0.05,
                   delta_ddg: float = 2.0
                   ) -> Dict[Tuple[str, str, str], str]:
    """Build (pair, method, cluster) -> centered preference (F1/F2/Amb).

    If the survey CSV already has the centered column ``pref_centered``,
    use it directly. Otherwise we recompute centering on the fly so older
    survey CSVs still work.
    """
    df = pd.read_csv(survey_csv)
    df["cluster"] = df["cluster"].apply(normalize_cluster)

    if "pref_centered" not in df.columns:
        df["pref_centered"] = "Amb"
        for (pair_id, method), grp in df.groupby(["pair_id", "method"]):
            shallow = grp[grp["cluster"] != "DeepMsa"]
            tdvals = pd.to_numeric(shallow.get("TMdiff_max"), errors="coerce")
            tdvals = tdvals.dropna()
            if len(tdvals) >= 2:
                med = float(np.median(tdvals))
            else:
                med = 0.0
            delta = delta_ddg if method == "DDG" else delta_tm
            for idx, r in grp.iterrows():
                tdm = pd.to_numeric(r.get("TMdiff_max"), errors="coerce")
                if pd.isna(tdm):
                    continue
                centered = float(tdm) - (med if r["cluster"] != "DeepMsa" else 0.0)
                if centered > delta:
                    df.at[idx, "pref_centered"] = "F1"
                elif centered < -delta:
                    df.at[idx, "pref_centered"] = "F2"
                else:
                    df.at[idx, "pref_centered"] = "Amb"

    out: Dict[Tuple[str, str, str], str] = {}
    for _, r in df.iterrows():
        key = (str(r["pair_id"]), str(r["method"]), str(r["cluster"]))
        pc = r.get("pref_centered")
        if not isinstance(pc, str) or pc not in ("F1", "F2", "Amb"):
            pc = "Amb"
        out[key] = pc
    return out


# ---------------------------------------------------------------------------
# Phylogenetic-structure statistics: see utils/phylo_tests.py (test registry).
# Re-exported here so `from phylo_placement import morans_I_test` etc. keep
# working.
# ---------------------------------------------------------------------------

from utils.phylo_tests import (  # noqa: E402,F401
    TESTS, P_COLUMNS, run_test,
    _fitch_parsimony_binary, _cophenetic_neighbours, nn_concordance,
    d_statistic, permutation_p, _patristic_matrix, morans_I_test,
)


# ---------------------------------------------------------------------------
# Which tests to run (--tests, with the legacy flags as aliases)
# ---------------------------------------------------------------------------

DISCRETE_TESTS: Tuple[str, ...] = ("fitch", "nn", "dstat")


def resolve_tests(tests_arg: Optional[str], discrete_perms: bool = False,
                  do_d_stat: bool = True) -> List[str]:
    """Resolve the list of registered tests to run.

    Default (``tests_arg`` None/''): the pre-registry behaviour -- 'morans',
    plus 'fitch,nn,dstat' when ``--discrete-perms``; ``--no-d-statistic``
    drops 'dstat'. An explicit comma list is used as given ('discrete' is
    shorthand for fitch,nn,dstat), still honouring the two legacy flags.
    """
    if not tests_arg:
        tests = ["morans"] + (list(DISCRETE_TESTS) if discrete_perms else [])
    else:
        tests = []
        for t in tests_arg.split(","):
            t = t.strip()
            if not t:
                continue
            for tt in (DISCRETE_TESTS if t == "discrete" else (t,)):
                if tt not in TESTS:
                    raise SystemExit(f"[error] unknown --tests entry {tt!r}; "
                                     f"known: {','.join(TESTS)}")
                if tt not in tests:
                    tests.append(tt)
        if discrete_perms:
            tests += [t for t in DISCRETE_TESTS if t not in tests]
    if not do_d_stat:
        tests = [t for t in tests if t != "dstat"]
    return tests


# ---------------------------------------------------------------------------
# Per-pair, per-method analysis on the COARSE tree
# ---------------------------------------------------------------------------

def analyze_pair_method(pair_id: str, method: str,
                        tree, tag_to_label: Dict[str, str],
                        cluster_pref: Dict[str, str],
                        n_perm: int, seed: int,
                        do_d_stat: bool = True,
                        cluster_value: Optional[Dict[str, Optional[float]]] = None,
                        do_discrete_perms: bool = False,
                        tests: Optional[Sequence[str]] = None,
                        clade_n_perm: int = 10000,
                        min_clade: int = 2,
                        ) -> Optional[dict]:
    """Compute phylo-structure stats for one (pair, method) on the COARSE tree.

    ``tag_to_label`` maps each cluster tag (ShallowMsa_NNN / DeepMsa) to
    the short leaf label used in the coarse tree (e.g. 'S03', 'D').
    ``cluster_pref`` is {cluster_tag: 'F1'|'F2'|'Amb'}.
    ``tests`` lists registered tests (utils/phylo_tests.TESTS) to run; None
    derives it from the legacy ``do_discrete_perms`` / ``do_d_stat`` flags.
    """
    # Map coarse-tree leaf name -> centered preference
    leaf_states: Dict[str, str] = {}
    for tag, label in tag_to_label.items():
        leaf_states[label] = cluster_pref.get(tag, "Amb")
    for leaf in tree.iter_leaves():
        leaf_states.setdefault(leaf.name, "Amb")

    counts = Counter(leaf_states.values())
    n_f1 = counts.get("F1", 0)
    n_f2 = counts.get("F2", 0)
    n_amb = counts.get("Amb", 0)
    n_total = sum(counts.values())

    notes: List[str] = []
    if n_f1 == 0 or n_f2 == 0:
        notes.append("no_F1F2_split")

    out = {
        "pair_id": pair_id,
        "method": method,
        "n_leaves": n_total,                 # = n_clusters in coarse tree
        "n_F1c_leaves": n_f1,
        "n_F2c_leaves": n_f2,
        "n_Amc_leaves": n_amb,
        "parsimony_score": np.nan,
        "parsimony_p": np.nan,
        "nn_concordance": np.nan,
        "nn_concordance_p": np.nan,
        "nn3_concordance": np.nan,
        "nn3_concordance_p": np.nan,
        "D_statistic": np.nan,
        "morans_I": np.nan,
        "morans_p": np.nan,
        "morans_z": np.nan,
        "morans_lambda_max": np.nan,
        "morans_max_perm": np.nan,
        "morans_I_norm": np.nan,
        "morans_n": 0,
        "morans_note": "",
        "notes": "",
    }

    if tests is None:
        tests = resolve_tests(None, do_discrete_perms, do_d_stat)

    # ---- CONTINUOUS tests (no thresholding): Moran's I (primary), clade scan ----
    # Use the raw/corrected continuous preference per cluster; every cluster
    # participates (no F1/F2/Amb gate), so coverage is method-independent.
    if cluster_value is not None:
        leaf_values: Dict[str, Optional[float]] = {}
        for tag, label in tag_to_label.items():
            leaf_values[label] = cluster_value.get(tag)
        if "morans" in tests:
            mt = run_test("morans", tree, leaf_values, n_perm=min(n_perm, 1000),
                          rng=np.random.default_rng(seed + 7))
            for col in ("morans_I", "morans_p", "morans_z", "morans_lambda_max",
                        "morans_max_perm", "morans_I_norm", "morans_n",
                        "morans_note"):
                out[col] = mt.get(col, np.nan)
        if "clade_shift" in tests:
            label_to_tag = {lab: tag for tag, lab in tag_to_label.items()}
            out.update(run_test("clade_shift", tree, leaf_values,
                                n_perm=clade_n_perm,
                                rng=np.random.default_rng(seed + 11),
                                min_clade=min_clade,
                                label_to_tag=label_to_tag))

    # ---- SECONDARY (discrete; thresholded labels) ----
    # OFF by default: Moran's I is the primary phylo test.
    # The discrete Fitch/NN permutation tests are slow (1000 tree traversals
    # per stat per pair-method) and superseded. Code retained; enable with
    # --discrete-perms (or --tests ...,fitch,nn,dstat) for the legacy
    # parsimony/NN/D columns.
    if n_f1 > 0 and n_f2 > 0:
        if "fitch" in tests:
            out.update(run_test("fitch", tree, leaf_states, n_perm=n_perm,
                                seed=seed))
        if "nn" in tests:
            out.update(run_test("nn", tree, leaf_states,
                                n_perm=min(n_perm, 1000), seed=seed + 1))
        if "dstat" in tests:
            try:
                out.update(run_test("dstat", tree, leaf_states, n_perm=200,
                                    seed=seed + 3))
            except Exception as e:
                notes.append(f"D_failed:{type(e).__name__}")
    out["notes"] = ";".join(notes)
    out["_leaf_states"] = leaf_states  # stash for figure pass; popped before CSV
    return out


def analyze_pair_combined(pair_id: str, tree, tag_to_label: Dict[str, str],
                          cluster_value_per_method: Dict[str, Dict[str, Optional[float]]],
                          n_perm: int, seed: int, min_clade: int = 2) -> dict:
    """COMBINED multi-method clade-shift scan: one row per pair
    (method='COMBINED'); see utils.phylo_tests.run_clade_shift_combined."""
    values = {m: {lab: cv.get(tag) for tag, lab in tag_to_label.items()}
              for m, cv in cluster_value_per_method.items()}
    label_to_tag = {lab: tag for tag, lab in tag_to_label.items()}
    row = {"pair_id": pair_id, "method": "COMBINED",
           "n_leaves": len(tag_to_label)}
    row.update(run_test("clade_shift_combined", tree, values, n_perm=n_perm,
                        rng=np.random.default_rng(seed + 13),
                        min_clade=min_clade, label_to_tag=label_to_tag))
    return row


# ---------------------------------------------------------------------------
# Coarse-tree figure (per pair, per method)
# ---------------------------------------------------------------------------

_STATE_COLORS = {"F1": "#d62728", "F2": "#1f77b4", "Amb": "#bbbbbb"}
_STATE_DISPLAY = {"F1": "F1c", "F2": "F2c", "Amb": "Amc"}


def render_phylo_placement(pair_id: str, method: str,
                           tree, leaf_states: Dict[str, str],
                           output_png: str,
                           parsimony_score: Optional[int] = None,
                           parsimony_p: Optional[float] = None,
                           nn_concordance_val: Optional[float] = None,
                           ) -> None:
    """Save a PNG of the COARSE tree, with leaves coloured by F1c/F2c/Amc
    and labelled by short cluster tag (S03, D, ...).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    leaves = list(tree.iter_leaves())
    n_leaves = len(leaves)
    if n_leaves == 0:
        return

    node_x: Dict[int, float] = {}
    node_y: Dict[int, float] = {}

    def _depth(n, d):
        node_x[id(n)] = d
        for c in n.children:
            cd = c.dist if (c.dist and c.dist > 0) else 0.1
            _depth(c, d + cd)
    _depth(tree, 0.0)

    leaf_idx = 0

    def _yset(n):
        nonlocal leaf_idx
        if n.is_leaf():
            node_y[id(n)] = leaf_idx
            leaf_idx += 1
        else:
            ys = []
            for c in n.children:
                _yset(c)
                ys.append(node_y[id(c)])
            node_y[id(n)] = float(np.mean(ys))
    _yset(tree)

    max_x = max(node_x.values()) if node_x else 1.0
    if max_x <= 0:
        max_x = 1.0

    fig_h = max(3.0, 0.45 * n_leaves + 1.0)
    fig_w = 7.0
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    for node in tree.traverse("preorder"):
        if node.is_root():
            continue
        x0, y0 = node_x[id(node.up)], node_y[id(node.up)]
        x1, y1 = node_x[id(node)], node_y[id(node)]
        ax.plot([x0, x0], [y0, y1], color="#444", linewidth=1.1, zorder=1)
        ax.plot([x0, x1], [y1, y1], color="#444", linewidth=1.1, zorder=1)

    for leaf in leaves:
        s = leaf_states.get(leaf.name, "Amb")
        col = _STATE_COLORS.get(s, _STATE_COLORS["Amb"])
        x, y = node_x[id(leaf)], node_y[id(leaf)]
        ax.scatter(x, y, c=col, s=140, zorder=3,
                   edgecolors="black", linewidths=0.5)
        ax.text(x + max_x * 0.03, y,
                f" {leaf.name} ({_STATE_DISPLAY.get(s, s)})",
                va="center", ha="left", fontsize=9)

    # Internal-node markers (small grey squares)
    for node in tree.traverse("preorder"):
        if node.is_leaf():
            continue
        ax.scatter(node_x[id(node)], node_y[id(node)],
                   c="#cccccc", s=28, zorder=2, marker="s",
                   edgecolors="black", linewidths=0.4)

    title = f"{pair_id} | {method}  (coarse cluster-level tree)"
    if parsimony_score is not None and not (isinstance(parsimony_score, float)
                                            and np.isnan(parsimony_score)):
        title += f"\nparsimony={int(parsimony_score)}"
    if parsimony_p is not None and not (isinstance(parsimony_p, float)
                                        and np.isnan(parsimony_p)):
        title += f"  p={parsimony_p:.3f}"
    if nn_concordance_val is not None and not (isinstance(nn_concordance_val, float)
                                               and np.isnan(nn_concordance_val)):
        title += f"  NN={nn_concordance_val:.2f}"
    ax.set_title(title, fontsize=10)

    legend_handles = [
        mpatches.Patch(color=_STATE_COLORS["F1"], label="F1c"),
        mpatches.Patch(color=_STATE_COLORS["F2"], label="F2c"),
        mpatches.Patch(color=_STATE_COLORS["Amb"], label="Amc"),
    ]
    ax.legend(handles=legend_handles, loc="upper right", fontsize=8,
              framealpha=0.85)
    ax.set_yticks([])
    ax.set_xlabel("Branch length (cumulative)", fontsize=8)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)

    plt.tight_layout()
    os.makedirs(os.path.dirname(output_png), exist_ok=True)
    plt.savefig(output_png, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# XCL1 (2n54B_2hdmA) Dishman cross-check at the cluster level
# ---------------------------------------------------------------------------

def xcl1_cross_check(pair_id: str,
                     tag_to_label: Dict[str, str],
                     cluster_pref_per_method: Dict[str, Dict[str, str]],
                     out_dir: str) -> Optional[dict]:
    """For 2n54B_2hdmA: print cluster-level F1c/F2c/Amc by method.

    The Dishman 2021 hypothesis places fold-switching gain on the XCL1
    lineage. At the coarse-tree level the per-cluster centred preferences
    are the right object to inspect against that hypothesis.
    """
    if pair_id != "2n54B_2hdmA":
        return None

    tags = sorted(tag_to_label.keys())
    print(f"\n[XCL1 cross-check] {pair_id} (coarse tree, {len(tags)} clusters)")
    header = f"  {'cluster':<18s} {'leaf':<6s}"
    for m in SCORED_METHODS:
        header += f" {m:>8s}"
    print(header)
    rows = []
    for tag in tags:
        label = tag_to_label[tag]
        prefs = []
        for m in SCORED_METHODS:
            p = cluster_pref_per_method.get(m, {}).get(tag, "Amb")
            prefs.append(_STATE_DISPLAY.get(p, p))
        line = f"  {tag:<18s} {label:<6s}"
        for p in prefs:
            line += f" {p:>8s}"
        print(line)
        rows.append({"cluster": tag, "leaf": label,
                     **{f"pref_{m}": prefs[i]
                        for i, m in enumerate(SCORED_METHODS)}})

    out_csv = os.path.join(out_dir, "phylo_placement_xcl1_clusters.csv")
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    print(f"  Wrote {out_csv}")
    return {"pair_id": pair_id, "n_clusters": len(tags)}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def discover_pairs() -> List[str]:
    if not os.path.isdir(DATA_DIR):
        return []
    return sorted(d for d in os.listdir(DATA_DIR)
                  if PAIR_DIR_RE.match(d)
                  and os.path.isdir(os.path.join(DATA_DIR, d)))


def _submit_self_sbatch(args, script_path: str) -> None:
    """Re-invoke this script under sbatch with --run_job_mode local.

    Mirrors the run_foldswitch_pipeline.py pattern: --run_job_mode sbatch
    submits a single sbatch job that runs the actual computation inline.
    """
    import shlex
    import subprocess
    jobs_dir = "Pipeline/FoldPairs/jobs"
    os.makedirs(jobs_dir, exist_ok=True)
    log = os.path.join(jobs_dir, f"phylo_placement_{args.labels}.out")
    # Build the inner command — same args minus --run_job_mode
    inner_args = []
    for k, v in vars(args).items():
        if k == "run_job_mode":
            continue
        if v is None or v is False:
            continue
        flag = "--" + k.replace("_", "-")
        if v is True:
            inner_args.append(flag)
        else:
            inner_args.extend([flag, str(v)])
    inner = " ".join([shlex.quote(sys.executable), script_path,
                      "--run-job-mode", "local"] + inner_args)
    sbatch_cmd = [
        "sbatch", "--parsable",
        f"--job-name=phylo_placement_{args.labels}",
        f"--cpus-per-task={args.sbatch_cpus}",
        f"--mem={args.sbatch_mem}",
        f"--time={args.sbatch_time}",
        f"--output={log}",
        "--wrap", shlex.quote(inner),
    ]
    if args.sbatch_partition:
        sbatch_cmd.insert(1, f"--partition={args.sbatch_partition}")
    cmd_str = " ".join(sbatch_cmd)
    print(f"[sbatch] {cmd_str}")
    jid = subprocess.check_output(cmd_str, shell=True, text=True).strip()
    print(f"[sbatch] submitted phylo_placement ({args.labels}) job {jid}")
    print(f"[sbatch] log: {log}")


def main():
    p = argparse.ArgumentParser(
        description="Phylogenetic-structure test for cluster fold-preference "
                    "on the COARSE cluster-level tree")
    p.add_argument("--pairs", default="ALL",
                   help="Comma-separated pair IDs, or ALL")
    p.add_argument("--labels", choices=["centered", "corrected"],
                   default="centered",
                   help="centered = use pref_centered from survey CSV (default); "
                        "corrected = use pref_corrected from "
                        "fold_diversity_survey_corrected.csv (residual-corrected, "
                        "PID-to-query stripped via seq_divergence_correction.py)")
    p.add_argument("--survey-csv", default=None,
                   help="Path to survey CSV (default depends on --labels)")
    p.add_argument("--output", default=None,
                   help="Output CSV path (default depends on --labels)")
    p.add_argument("--make-figures", default="0",
                   help="'all' = render every (pair, method) figure; "
                        "an integer N = top-N by parsimony_p per method; "
                        "0 = no figures (default).")
    p.add_argument("--figure-dir", default=None,
                   help="Where to write figures "
                        "(default: <docs>/figures/phylo_placement)")
    p.add_argument("--n-perm", type=int, default=1000)
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--delta-tm", type=float, default=0.05)
    p.add_argument("--delta-ddg", type=float, default=2.0,
                   help="kcal/mol threshold for discrete DDG labels (default 2.0, "
                        "matching the survey and the correction)")
    p.add_argument("--no-d-statistic", action="store_true",
                   help="Skip D statistic (slowest stat)")
    p.add_argument("--discrete-perms", action="store_true",
                   help="Also run the legacy discrete Fitch/NN permutation "
                        "tests (parsimony_p, nn_concordance_p, nn3, D). OFF by "
                        "default - Moran's I is the primary test and these are "
                        "slow. Discrete F1/F2/Amb label counts are always kept "
                        "(for figures); only their permutation p-values are gated.")
    p.add_argument("--tests", default=None,
                   help="Comma list of registered tests (utils/phylo_tests.py): "
                        "morans, fitch, nn, dstat ('discrete' = all three), "
                        "clade_shift, clade_shift_combined. Default: morans "
                        "(+ fitch,nn,dstat with --discrete-perms) -- i.e. the "
                        "pre-registry behaviour. --discrete-perms / "
                        "--no-d-statistic still apply as aliases.")
    p.add_argument("--clade-n-perm", type=int, default=10000,
                   help="Permutations for clade_shift / clade_shift_combined "
                        "(independent of --n-perm; default 10000)")
    p.add_argument("--min-clade", type=int, default=2,
                   help="Minimum clade AND rest size (in leaves) for a "
                        "clade-shift split (default 2)")
    p.add_argument("--clade-methods", default="ALL",
                   help="Comma list of methods entering clade_shift_combined "
                        "(default ALL scored methods with values)")
    # sbatch self-submission (matches run_foldswitch_pipeline.py pattern)
    p.add_argument("--run-job-mode", choices=["local", "sbatch"], default="local",
                   help="local = run inline (default); sbatch = self-submit "
                        "one sbatch job that runs the full analysis inline")
    p.add_argument("--sbatch-partition", default=None)
    p.add_argument("--sbatch-cpus", default=2)
    p.add_argument("--sbatch-mem", default="8G")
    p.add_argument("--sbatch-time", default="04:00:00")
    args = p.parse_args()

    # --- sbatch self-submission ---
    if args.run_job_mode == "sbatch":
        _submit_self_sbatch(args, os.path.abspath(__file__))
        return

    # --- mode-aware paths ---
    if args.labels == "corrected":
        default_survey = os.path.join(TABLES_RES, "fold_diversity_survey_corrected.csv")
        default_out = os.path.join(TABLES_RES, "phylo_placement_corrected.csv")
    else:
        default_survey = os.path.join(TABLES_RES, "fold_diversity_survey.csv")
        default_out = os.path.join(TABLES_RES, "phylo_placement.csv")
    survey_csv = args.survey_csv or default_survey
    out_csv = args.output or default_out
    figure_dir = args.figure_dir or os.path.join(
        TABLES_RES, "figures", "phylo_placement"
    )

    if not os.path.isfile(survey_csv):
        print(f"[error] survey CSV not found: {survey_csv}")
        sys.exit(2)

    if args.pairs.upper() == "ALL":
        pair_ids = discover_pairs()
    else:
        pair_ids = [s.strip() for s in args.pairs.split(",") if s.strip()]
    if not pair_ids:
        print("[error] no pairs to process")
        sys.exit(2)

    tests = resolve_tests(args.tests, args.discrete_perms,
                          not args.no_d_statistic)
    if args.clade_methods.upper() == "ALL":
        combined_methods = list(SCORED_METHODS)
    else:
        combined_methods = [m.strip() for m in args.clade_methods.split(",")
                            if m.strip()]
    print(f"Tests: {','.join(tests)}"
          + (f" (clade n_perm={args.clade_n_perm}, min_clade={args.min_clade})"
             if any(t.startswith("clade") for t in tests) else ""))

    if args.make_figures.lower() == "all":
        fig_mode = "all"
        fig_topn = 0
    else:
        try:
            fig_topn = int(args.make_figures)
        except ValueError:
            print(f"[error] --make-figures must be 'all' or an integer, "
                  f"got {args.make_figures!r}")
            sys.exit(2)
        fig_mode = "topn" if fig_topn > 0 else "none"

    print(f"Loading fold preferences ({args.labels}): {survey_csv}")
    if args.labels == "corrected":
        pref_lookup = load_corrected_pref(survey_csv)
        value_col = "TMdiff_residual"
    else:
        pref_lookup = load_fold_pref(survey_csv,
                                     delta_tm=args.delta_tm,
                                     delta_ddg=args.delta_ddg)
        value_col = "TMdiff_centered"
    # Continuous per-cluster preference for the Moran's-I test (primary stat).
    value_lookup = load_fold_value(survey_csv, value_col)
    pairs_in_survey = set(k[0] for k in pref_lookup)
    print(f"  {len(pairs_in_survey)} pairs in survey")

    cached_states: Dict[Tuple[str, str], Tuple[object, Dict[str, str], dict]] = {}

    rows: List[dict] = []
    t_total = time.time()
    per_pair_timeout_s = int(os.environ.get("PHYLO_PAIR_TIMEOUT", "120"))
    for i, pair_id in enumerate(sorted(pair_ids), 1):
        t_pair = time.time()
        print(f"\n[{i}/{len(pair_ids)}] {pair_id}: building coarse tree "
              f"(timeout {per_pair_timeout_s}s)...", flush=True)
        try:
            with _per_pair_timeout(per_pair_timeout_s):
                tree, tag_to_label, _ = build_coarse_tree(pair_id)
                _name_internal_nodes(tree)
        except TimeoutError as e:
            print(f"  TIMEOUT after {per_pair_timeout_s}s on coarse-tree build; "
                  f"skipping {pair_id}", flush=True)
            continue
        except FileNotFoundError as e:
            print(f"  skipped: {e}", flush=True)
            continue
        except Exception as e:
            print(f"  coarse-tree build failed: {e}", flush=True)
            continue

        n_clusters = len(tag_to_label)
        if n_clusters < 2:
            print(f"  only {n_clusters} cluster(s) -- skip ({time.time()-t_pair:.1f}s)", flush=True)
            continue
        print(f"  coarse tree built with {n_clusters} cluster-leaves "
              f"({time.time()-t_pair:.1f}s)", flush=True)

        cluster_pref_per_method: Dict[str, Dict[str, str]] = {}
        cluster_value_per_method: Dict[str, Dict[str, Optional[float]]] = {}
        for method in SCORED_METHODS:
            t_method = time.time()
            cp: Dict[str, str] = {
                tag: pref_lookup.get((pair_id, method, tag), "Amb")
                for tag in tag_to_label
            }
            cluster_pref_per_method[method] = cp
            cv: Dict[str, Optional[float]] = {
                tag: value_lookup.get((pair_id, method, tag))
                for tag in tag_to_label
            }
            if method in combined_methods:
                cluster_value_per_method[method] = cv

            try:
                with _per_pair_timeout(per_pair_timeout_s):
                    row = analyze_pair_method(
                        pair_id, method, tree, tag_to_label, cp,
                        n_perm=args.n_perm, seed=args.seed,
                        do_d_stat=(not args.no_d_statistic),
                        cluster_value=cv,
                        do_discrete_perms=args.discrete_perms,
                        tests=tests, clade_n_perm=args.clade_n_perm,
                        min_clade=args.min_clade,
                    )
            except TimeoutError:
                print(f"    {method}: TIMEOUT after {per_pair_timeout_s}s; skipped", flush=True)
                continue
            except Exception as e:
                print(f"    {method}: failed: {e}", flush=True)
                continue
            print(f"    {method}: {time.time()-t_method:.1f}s", flush=True)
            if row is None:
                continue
            row["n_clusters"] = n_clusters
            leaf_states = row.pop("_leaf_states")
            cached_states[(pair_id, method)] = (tree, leaf_states, row)
            rows.append(row)
        if "clade_shift_combined" in tests:
            try:
                t_method = time.time()
                crow = analyze_pair_combined(
                    pair_id, tree, tag_to_label, cluster_value_per_method,
                    n_perm=args.clade_n_perm, seed=args.seed,
                    min_clade=args.min_clade)
                crow["n_clusters"] = n_clusters
                rows.append(crow)
                print(f"    COMBINED: {time.time()-t_method:.1f}s "
                      f"(T={crow.get('clade_T')}, p={crow.get('clade_p')})",
                      flush=True)
            except Exception as e:
                print(f"    COMBINED: failed: {e}", flush=True)
        print(f"  pair total: {time.time()-t_pair:.1f}s  (running {time.time()-t_total:.0f}s)", flush=True)

        if pair_id == "2n54B_2hdmA":
            try:
                xcl1_cross_check(pair_id, tag_to_label,
                                 cluster_pref_per_method, TABLES_RES)
            except Exception as e:
                print(f"  [warn] XCL1 cross-check failed: {e}")

    df = pd.DataFrame(rows)
    if df.empty:
        print("[warn] no rows produced")
        return

    # BH FDR across the screen for each p-value column (always-on)
    for pcol in ("morans_p", "parsimony_p", "nn_concordance_p", "nn3_concordance_p",
                 "clade_p"):
        if pcol in df.columns:
            # BH within each method, across pairs (the per-method screen).
            df[pcol + "_bh"] = np.nan
            for _m in df["method"].unique():
                mask = df["method"] == _m
                df.loc[mask, pcol + "_bh"] = bh_adjust(
                    df.loc[mask, pcol].to_numpy(dtype=float))

    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"\nWrote {out_csv} ({len(df)} rows; labels={args.labels})")

    # Per-method BH-significant counts on the PRIMARY stat (Moran's I).
    print("\n=== Moran's I phylogenetic signal: BH-significant pairs per method ===")
    for method in SCORED_METHODS:
        sub = df[(df.method == method) & df.morans_p.notna()].copy()
        if sub.empty:
            continue
        nsig = int((sub["morans_p_bh"] < 0.05).sum()) if "morans_p_bh" in sub else 0
        print(f"  {method:9} tested={len(sub):>3}  "
              f"raw_p<.05={int((sub.morans_p<0.05).sum()):>3}  BH<.05={nsig:>3}  "
              f"median_I={sub.morans_I.median():.3f}")

    if "clade_p" in df.columns:
        print("\n=== Clade-shift scan (max |z| Mann-Whitney over clades): "
              "significant pairs per method ===")
        for method in list(SCORED_METHODS) + ["COMBINED"]:
            sub = df[(df.method == method) & df.clade_p.notna()]
            if sub.empty:
                continue
            print(f"  {method:9} tested={len(sub):>3}  "
                  f"raw_p<.05={int((sub.clade_p<0.05).sum()):>3}  "
                  f"BH<.05={int((sub.clade_p_bh<0.05).sum()):>3}  "
                  f"median_T={sub.clade_T.median():.2f}")

    for method in SCORED_METHODS:
        sub = df[(df.method == method) & df.morans_p.notna()].copy()
        if sub.empty:
            continue
        # p floors at 1/(n_perm+1); break ties by SIGNED I (clustering first)
        sub = sub.sort_values(["morans_p", "morans_I"], ascending=[True, False])
        print(f"\nTop 10 phylogenetically-structured pairs ({method}, by p then signed I):")
        cols = ["pair_id", "morans_n", "morans_I", "morans_p", "morans_p_bh",
                "n_F1c_leaves", "n_F2c_leaves", "nn_concordance_p"]
        cols = [c for c in cols if c in sub.columns]
        print(sub[cols].head(10).to_string(index=False))

    if fig_mode != "none":
        os.makedirs(figure_dir, exist_ok=True)
        if fig_mode == "all":
            print(f"\nRendering ALL coarse-tree figures -> {figure_dir}")
            targets = [(k, v) for k, v in cached_states.items()]
        else:
            print(f"\nRendering top {fig_topn} coarse-tree figures per method "
                  f"-> {figure_dir}")
            targets = []
            for method in SCORED_METHODS:
                sub = df[(df.method == method) & df.parsimony_p.notna()].copy()
                if sub.empty:
                    continue
                sub = sub.sort_values("parsimony_p")
                for _, r in sub.head(fig_topn).iterrows():
                    key = (r.pair_id, method)
                    if key in cached_states:
                        targets.append((key, cached_states[key]))
        for (pair_id, method), (tree, leaf_states, row) in targets:
            out_png = os.path.join(
                figure_dir,
                f"phylo_placement_{pair_id}_{method}.png"
            )
            try:
                render_phylo_placement(
                    pair_id, method, tree, leaf_states, out_png,
                    parsimony_score=row.get("parsimony_score"),
                    parsimony_p=row.get("parsimony_p"),
                    nn_concordance_val=row.get("nn_concordance"),
                )
            except Exception as e:
                print(f"  [warn] figure failed for {pair_id}/{method}: {e}")

    print("\nDone.")


if __name__ == "__main__":
    main()
