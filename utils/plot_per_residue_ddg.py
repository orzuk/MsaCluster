#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_per_residue_ddg.py — per-residue ThermoMPNN ΔΔG figure for per-pair pages.

Replaces the old, unreliable raw-Rosetta `<pair>_ddg_aligned.png` (E1-E2 on
unrelaxed whole-assembly structures, dominated by clashes) with the ThermoMPNN
conformation-biasing signal computed by `run_DDG.py`.

For each aligned (seed match-state) position we compute, averaged over the
pair's cluster homolog sequences:

    delta(pos) = mean_s [ contrib_F2(pos, s) - contrib_F1(pos, s) ]

This is the per-position contribution to the conformation-bias score
`bias = sum_F2 - sum_F1`. Positive => the position destabilizes fold 2 more
than fold 1 => favors fold 1 (red); negative => favors fold 2 (blue).

The public entry point never raises: it returns the output PNG path on success
or `None` when any required input is missing/unreadable.

Per-position contributions are sourced in this priority order:
  1. `Analysis/df_ddg.csv` columns, IF a future version of run_DDG ever stores
     per-position contrib vectors there (it currently does not — see below).
  2. The per-cluster JSON files `output_ddg/cluster_*.json`, which DO store
     per-sequence `contribs_F1` / `contribs_F2` (written by run_DDG._score_sequence).
     This is the current, canonical source.
  3. Recompute from the cached `output_ddg/ddg_matrix_<tag>.npz` matrices plus
     the cluster a3m homolog sequences, mirroring run_DDG._score_sequence.
If none is available, return None gracefully.
"""
from __future__ import annotations

import os
import re
import glob
import json
from typing import List, Optional, Tuple

import numpy as np

import matplotlib
matplotlib.use("Agg")  # headless / no display
import matplotlib.pyplot as plt  # noqa: E402


# Amino-acid ordering must match run_DDG / ProteinMPNN.
AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")
AA_TO_IDX = {aa: i for i, aa in enumerate(AA_LIST)}

# Candidate names (lowercased) for per-position contrib columns in df_ddg.csv,
# should a future run_DDG ever serialize them per cluster.
_CONTRIB_F1_KEYS = ("contribs_f1", "contrib_f1", "per_pos_contribs_f1")
_CONTRIB_F2_KEYS = ("contribs_f2", "contrib_f2", "per_pos_contribs_f2")


# ----------------------------------------------------------------------
# Path helpers (kept local so this module has no hard dependency on config
# at import time; config is imported lazily inside _data_dir()).
# ----------------------------------------------------------------------

def _data_dir() -> Optional[str]:
    """Return the FoldPairs data directory, or None if config is unavailable."""
    try:
        from config import DATA_DIR
        return DATA_DIR
    except Exception:
        return None


def _pair_dir(pair_id: str) -> Optional[str]:
    dd = _data_dir()
    if dd is None:
        return None
    return os.path.join(dd, pair_id)


# ----------------------------------------------------------------------
# Contrib parsing helpers
# ----------------------------------------------------------------------

def _coerce_vector(val) -> Optional[np.ndarray]:
    """Coerce a cell/value into a 1-D float array, or None on failure.

    Accepts python lists, numpy arrays, or string reprs like "[0.1, -0.2, ...]"
    / "0.1 -0.2 ..." (as may appear after a CSV round-trip).
    """
    if val is None:
        return None
    try:
        if isinstance(val, (list, tuple, np.ndarray)):
            arr = np.asarray(val, dtype=np.float64)
            return arr if arr.ndim == 1 and arr.size else None
        if isinstance(val, str):
            s = val.strip()
            if not s or s.lower() in ("nan", "none", "[]"):
                return None
            # Try JSON first, then loose whitespace/comma split.
            try:
                parsed = json.loads(s)
                arr = np.asarray(parsed, dtype=np.float64)
                return arr if arr.ndim == 1 and arr.size else None
            except Exception:
                toks = re.split(r"[,\s]+", s.strip("[]() "))
                vals = [float(t) for t in toks if t not in ("", "nan", "None")]
                return np.asarray(vals, dtype=np.float64) if vals else None
    except Exception:
        return None
    return None


def _accumulate(delta_sum: Optional[np.ndarray], delta_cnt: Optional[np.ndarray],
                c1: np.ndarray, c2: np.ndarray
                ) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """Add one (contrib_F1, contrib_F2) sequence pair into running sums.

    Returns (delta_sum, delta_cnt). NaN positions are skipped per-position so a
    single bad value never poisons the whole column.
    """
    n = min(len(c1), len(c2))
    if n == 0:
        return delta_sum, delta_cnt
    c1 = np.asarray(c1[:n], dtype=np.float64)
    c2 = np.asarray(c2[:n], dtype=np.float64)
    d = c2 - c1
    valid = np.isfinite(d)
    if delta_sum is None:
        delta_sum = np.zeros(n, dtype=np.float64)
        delta_cnt = np.zeros(n, dtype=np.float64)
    elif len(delta_sum) < n:
        delta_sum = np.concatenate([delta_sum, np.zeros(n - len(delta_sum))])
        delta_cnt = np.concatenate([delta_cnt, np.zeros(n - len(delta_cnt))])
    delta_sum[:n][valid] += d[valid]
    delta_cnt[:n][valid] += 1.0
    return delta_sum, delta_cnt


# ----------------------------------------------------------------------
# Source 1: df_ddg.csv columns (only if a future version stores them)
# ----------------------------------------------------------------------

def _delta_from_df(pair_id: str) -> Optional[np.ndarray]:
    pd_dir = _pair_dir(pair_id)
    if pd_dir is None:
        return None
    df_path = os.path.join(pd_dir, "Analysis", "df_ddg.csv")
    if not os.path.isfile(df_path):
        return None
    try:
        import pandas as pd
        df = pd.read_csv(df_path)
    except Exception:
        return None
    if df is None or df.empty:
        return None

    lc = {str(c).lower(): c for c in df.columns}
    f1col = next((lc[k] for k in _CONTRIB_F1_KEYS if k in lc), None)
    f2col = next((lc[k] for k in _CONTRIB_F2_KEYS if k in lc), None)
    if f1col is None or f2col is None:
        return None  # current schema: not present -> fall through to JSON source

    delta_sum = delta_cnt = None
    try:
        for _, row in df.iterrows():
            c1 = _coerce_vector(row[f1col])
            c2 = _coerce_vector(row[f2col])
            if c1 is None or c2 is None:
                continue
            delta_sum, delta_cnt = _accumulate(delta_sum, delta_cnt, c1, c2)
    except Exception:
        return None
    return _finalize_mean(delta_sum, delta_cnt)


# ----------------------------------------------------------------------
# Source 2: per-cluster JSON files (canonical current source)
# ----------------------------------------------------------------------

def _delta_from_cluster_jsons(pair_id: str) -> Optional[np.ndarray]:
    pd_dir = _pair_dir(pair_id)
    if pd_dir is None:
        return None
    jsons = sorted(glob.glob(os.path.join(pd_dir, "output_ddg", "cluster_*.json")))
    if not jsons:
        return None

    delta_sum = delta_cnt = None
    for jp in jsons:
        try:
            with open(jp, "r") as fh:
                obj = json.load(fh)
        except Exception:
            continue
        seqs = obj.get("sequences") if isinstance(obj, dict) else None
        if not seqs:
            continue
        for s in seqs:
            try:
                c1 = _coerce_vector(s.get("contribs_F1"))
                c2 = _coerce_vector(s.get("contribs_F2"))
            except Exception:
                continue
            if c1 is None or c2 is None:
                continue
            delta_sum, delta_cnt = _accumulate(delta_sum, delta_cnt, c1, c2)
    return _finalize_mean(delta_sum, delta_cnt)


# ----------------------------------------------------------------------
# Source 3: recompute from cached ddg matrices + cluster a3m sequences
# ----------------------------------------------------------------------

def _strip_insertions(a3m_seq: str) -> str:
    """Match-state projection of an a3m sequence (drop lowercase insertions)."""
    return "".join(c for c in a3m_seq if c.isupper() or c == "-")


def _load_matrix(pd_dir: str, tag: str) -> Optional[Tuple[np.ndarray, str]]:
    cache = os.path.join(pd_dir, "output_ddg", f"ddg_matrix_{tag}.npz")
    if not os.path.isfile(cache):
        return None
    try:
        z = np.load(cache, allow_pickle=False)
        ddg = np.asarray(z["ddg"], dtype=np.float64)
        seq = str(z["seq"])
        if ddg.ndim != 2 or ddg.shape[1] != len(AA_LIST):
            return None
        return ddg, seq
    except Exception:
        return None


def _delta_from_recompute(pair_id: str) -> Optional[np.ndarray]:
    """Recompute the per-position mean delta from cached matrices + clusters.

    Mirrors run_DDG._score_sequence: for each homolog match-state sequence and
    each fold, contrib[k] = ddg[pdb_idx[k], AA_TO_IDX[aa_k]]. We then average
    (contrib_F2 - contrib_F1) over all homolog sequences and clusters.
    """
    pd_dir = _pair_dir(pair_id)
    if pd_dir is None:
        return None

    # Cached ddg matrices (need tags).
    try:
        from utils.utils import pair_str_to_tuple
        tagA, tagB = pair_str_to_tuple(pair_id)
    except Exception:
        # Fall back to discovering tags from the matrix filenames.
        mats = sorted(glob.glob(os.path.join(pd_dir, "output_ddg", "ddg_matrix_*.npz")))
        tags = [os.path.basename(m)[len("ddg_matrix_"):-len(".npz")] for m in mats]
        if len(tags) < 2:
            return None
        tagA, tagB = tags[0], tags[1]

    m1 = _load_matrix(pd_dir, tagA)
    m2 = _load_matrix(pd_dir, tagB)
    if m1 is None or m2 is None:
        return None
    ddg1, seq1 = m1
    ddg2, seq2 = m2

    # Seed match-state -> (msa column, pdb residue index) mapping per fold.
    deep_a3m = os.path.join(pd_dir, "output_get_msa", "DeepMsa.a3m")
    if not os.path.isfile(deep_a3m):
        return None
    try:
        from Analysis.cmap_analysis import seed_residue_indices
        msa_cols_F1, seed_match_F1 = seed_residue_indices(
            deep_a3m, tagA, fallback_index=0, return_seq=True)
        msa_cols_F2, seed_match_F2 = seed_residue_indices(
            deep_a3m, tagB, fallback_index=1, return_seq=True)
    except Exception:
        return None

    pdb_idx_F1 = _seed_to_pdb_map(seed_match_F1, seq1)
    pdb_idx_F2 = _seed_to_pdb_map(seed_match_F2, seq2)
    if pdb_idx_F1 is None or pdb_idx_F2 is None:
        return None
    msa_cols_F1 = np.asarray(msa_cols_F1, dtype=np.int64)
    msa_cols_F2 = np.asarray(msa_cols_F2, dtype=np.int64)
    if len(msa_cols_F1) != len(pdb_idx_F1) or len(msa_cols_F2) != len(pdb_idx_F2):
        return None

    # Cluster a3m sequences.
    cluster_a3ms = sorted(glob.glob(
        os.path.join(pd_dir, "output_msa_cluster", "ShallowMsa_*.a3m")))
    if not cluster_a3ms:
        return None

    try:
        from utils.protein_utils import read_msa, greedy_select
    except Exception:
        read_msa = None
        greedy_select = None
    if read_msa is None:
        return None

    def _contribs(match_seq: str, ddg: np.ndarray,
                  msa_cols: np.ndarray, pdb_idx: np.ndarray) -> np.ndarray:
        L = len(msa_cols)
        c = np.full(L, np.nan, dtype=np.float64)
        for k in range(L):
            p = int(pdb_idx[k])
            if p < 0 or p >= ddg.shape[0]:
                continue
            col = int(msa_cols[k])
            if col < 0 or col >= len(match_seq):
                continue
            aa = match_seq[col]
            if aa == "-" or aa not in AA_TO_IDX:
                continue
            c[k] = ddg[p, AA_TO_IDX[aa]]
        return c

    delta_sum = delta_cnt = None
    for a3m in cluster_a3ms:
        try:
            entries = read_msa(a3m)
        except Exception:
            continue
        if not entries:
            continue
        try:
            picked = greedy_select(entries, min(10, len(entries))) if greedy_select else entries[:10]
        except Exception:
            picked = entries[:10]
        if not picked:
            continue
        for item in picked:
            try:
                a3m_seq = item[1]
            except Exception:
                continue
            match = _strip_insertions(a3m_seq)
            c1 = _contribs(match, ddg1, msa_cols_F1, pdb_idx_F1)
            c2 = _contribs(match, ddg2, msa_cols_F2, pdb_idx_F2)
            delta_sum, delta_cnt = _accumulate(delta_sum, delta_cnt, c1, c2)

    return _finalize_mean(delta_sum, delta_cnt)


def _seed_to_pdb_map(seed_match_seq: str, pdb_seq: str) -> Optional[np.ndarray]:
    """Map each seed match-state position to a PDB residue index (or -1).

    Mirrors run_DDG._seed_to_pdb_map: fast path for identical sequences,
    otherwise the project's canonical pairwise aligner.
    """
    try:
        if seed_match_seq == pdb_seq:
            return np.arange(len(seed_match_seq), dtype=np.int64)
        from utils.msa_utils import get_align_indexes
        idx_seed, idx_pdb = get_align_indexes(seed_match_seq, pdb_seq)
        seed_to_pdb = np.full(len(seed_match_seq), -1, dtype=np.int64)
        for s, p in zip(idx_seed, idx_pdb):
            seed_to_pdb[s] = p
        return seed_to_pdb
    except Exception:
        return None


# ----------------------------------------------------------------------
# Finalize / plot
# ----------------------------------------------------------------------

def _finalize_mean(delta_sum: Optional[np.ndarray],
                   delta_cnt: Optional[np.ndarray]) -> Optional[np.ndarray]:
    if delta_sum is None or delta_cnt is None:
        return None
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(delta_cnt > 0, delta_sum / np.maximum(delta_cnt, 1.0), np.nan)
    if not np.any(np.isfinite(mean)):
        return None
    return mean


def _compute_delta(pair_id: str) -> Optional[np.ndarray]:
    """Per-position mean of (contrib_F2 - contrib_F1) across cluster homologs.

    Tries df_ddg.csv columns, then per-cluster JSONs, then recompute.
    """
    for fn in (_delta_from_df, _delta_from_cluster_jsons, _delta_from_recompute):
        try:
            delta = fn(pair_id)
        except Exception:
            delta = None
        if delta is not None and np.any(np.isfinite(delta)):
            return delta
    return None


def _default_out_path(pair_id: str) -> str:
    base = _data_dir()
    # Anchor docs/ next to the repo root (parent of Pipeline), falling back to cwd.
    repo_root = None
    if base is not None:
        # DATA_DIR = <MAIN_DIR>/Pipeline/FoldPairs -> MAIN_DIR is two levels up.
        repo_root = os.path.dirname(os.path.dirname(base))
    if not repo_root or not os.path.isdir(repo_root):
        repo_root = os.getcwd()
    return os.path.join(repo_root, "docs", "HTML", "figs", pair_id,
                        f"{pair_id}_per_residue_ddg.png")


def _chain_pdb_path(pair_id: str, tag: str):
    """Pipeline/FoldPairs/<pair>/chain_pdb_files/<tag>.pdb (or None)."""
    base = _data_dir()
    if not base or not tag:
        return None
    p = os.path.join(base, pair_id, "chain_pdb_files", f"{tag}.pdb")
    return p if os.path.isfile(p) else None


def _dssp_ss_seq(pdb_path):
    """(ss, seq) for a single-chain PDB via mdtraj DSSP (Q3 H/E/C). None on fail."""
    try:
        import mdtraj as md
        t = md.load(pdb_path)
        ss = "".join(md.compute_dssp(t, simplified=True)[0])     # 'H'/'E'/'C'/'NA'
        seq = t.topology.to_fasta()[0]
        ss = ss.replace("NA", "C")
        if len(ss) and len(seq) and abs(len(ss) - len(seq)) <= 2:
            return ss, seq
    except Exception as e:
        print(f"[ddg-ss] DSSP failed for {pdb_path}: {e}")
    return None, None


def _map_ss_onto_frame(seq_ref, seq_other, ss_other):
    """Map ss_other (on seq_other) onto seq_ref's positions via a pairwise
    alignment (align_utils' Biopython pairwise2). Returns a list len(seq_ref)
    of SS chars ('-' where seq_other has a gap). Best-effort; None on failure."""
    try:
        from Bio import pairwise2
        aln = pairwise2.align.globalms(seq_ref, seq_other, 2, -1, -5, -0.5,
                                       one_alignment_only=True)[0]
        a_ref, a_oth = aln.seqA, aln.seqB
        out = []
        j = 0  # index into seq_other / ss_other
        for ca, co in zip(a_ref, a_oth):
            if ca != "-":  # a reference position -> emit aligned other SS
                out.append(ss_other[j] if (co != "-" and j < len(ss_other)) else "-")
            if co != "-":
                j += 1
        return out
    except Exception as e:
        print(f"[ddg-ss] sequence alignment failed: {e}")
        return None


def _draw_ss_track(axt, ss_chars, n, label):
    """Colour an SS track (helix/strand/coil) along residues 1..len on axt."""
    from matplotlib.patches import Rectangle
    SS_COL = {"H": "#e07b39", "G": "#e07b39", "I": "#e07b39",   # helix shades
              "E": "#f4c430", "B": "#f4c430",                    # strand
              "C": "#e8e8e8", "T": "#e8e8e8", "S": "#e8e8e8", "-": "#ffffff"}
    for i, c in enumerate(ss_chars):
        axt.add_patch(Rectangle((i + 0.5, 0.0), 1.0, 1.0,
                                color=SS_COL.get(c, "#e8e8e8"), linewidth=0))
    axt.set_xlim(0.5, n + 0.5)
    axt.set_ylim(0.0, 1.0)
    axt.set_xticks([]); axt.set_yticks([])
    axt.set_ylabel(label, rotation=0, ha="right", va="center", fontsize=7)


def per_residue_ddg_fig(pair_id: str, out_path: str | None = None) -> str | None:
    """Build the per-residue ThermoMPNN ΔΔG figure for a fold-switch pair.

    For each aligned position, plots the mean over the pair's cluster homolog
    sequences of (contrib_F2 - contrib_F1) -- the per-position contribution to
    the conformation-bias score. Positive (blue) favors fold 1; negative (red)
    favors fold 2.

    Parameters
    ----------
    pair_id : str
        Fold-switch pair id, e.g. "4dxrA_4dxtA".
    out_path : str | None
        Output PNG path. Defaults to
        docs/HTML/figs/<pair_id>/<pair_id>_per_residue_ddg.png.

    Returns
    -------
    str | None
        Path to the written PNG, or None if inputs are missing/unusable.
        Never raises.
    """
    try:
        if not pair_id or not isinstance(pair_id, str):
            return None

        delta = _compute_delta(pair_id)
        if delta is None:
            return None

        if out_path is None:
            out_path = _default_out_path(pair_id)
        try:
            os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
        except Exception:
            return None

        # Parse fold tags for the legend/title when available.
        tagA = tagB = None
        try:
            from utils.utils import pair_str_to_tuple
            tagA, tagB = pair_str_to_tuple(pair_id)
        except Exception:
            pass

        n = len(delta)
        positions = np.arange(1, n + 1)
        vals = np.nan_to_num(delta, nan=0.0)
        # Unified convention: fold1/pdb1 = blue, fold2/pdb2 = red (matches the
        # tree, heatmap, contact-map ΔΔG strip and 3D viewers).
        colors = np.where(vals >= 0, "#2c6fbb", "#c0392b")  # blue fold1 / red fold2

        fig, ax = plt.subplots(figsize=(max(8.0, min(22.0, 0.06 * n + 6.0)), 4.2))
        ax.bar(positions, vals, color=list(colors), width=1.0, linewidth=0)
        ax.axhline(0.0, color="black", linewidth=0.8)

        # Annotate the few largest-magnitude positions.
        finite = np.where(np.isfinite(delta))[0]
        if finite.size:
            order = finite[np.argsort(-np.abs(delta[finite]))]
            n_annot = min(5, order.size)
            ymax = np.nanmax(np.abs(delta[finite])) or 1.0
            for idx in order[:n_annot]:
                d = float(delta[idx])
                if not np.isfinite(d) or d == 0.0:
                    continue
                ax.text(idx + 1, d + (0.02 * ymax if d >= 0 else -0.02 * ymax),
                        str(idx + 1), ha="center",
                        va="bottom" if d >= 0 else "top", fontsize=8)

        f1_label = f"favors fold 1 ({tagA})" if tagA else "favors fold 1"
        f2_label = f"favors fold 2 ({tagB})" if tagB else "favors fold 2"
        handles = [
            plt.Rectangle((0, 0), 1, 1, color="#2c6fbb"),
            plt.Rectangle((0, 0), 1, 1, color="#c0392b"),
        ]
        ax.legend(handles, [f1_label, f2_label], loc="best", fontsize=8, frameon=False)

        ax.set_xlabel("residue position")
        ax.set_ylabel("ΔΔG contribution to fold preference (kcal/mol)")
        ax.set_title(f"Per-residue fold-preference contribution ({pair_id})")
        ax.set_xlim(0.5, n + 0.5)
        ax.margins(x=0)

        # --- secondary-structure tracks (Porter Fig 1a/b style) ---
        # fold1 SS on top, fold2 SS on the bottom, on the same residue axis as
        # the ΔΔG bars, so ΔΔG peaks line up with helix<->strand changes between
        # the two folds. SS from DSSP (mdtraj); the two folds' residues are
        # matched with align_utils' pairwise alignment (fold2 mapped onto fold1).
        try:
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            ss1, seq1 = _dssp_ss_seq(_chain_pdb_path(pair_id, tagA)) if tagA else (None, None)
            ss2, seq2 = _dssp_ss_seq(_chain_pdb_path(pair_id, tagB)) if tagB else (None, None)
            if ss1 and ss2 and seq1 and seq2:
                ss2_on1 = _map_ss_onto_frame(seq1, seq2, ss2) or list(ss2)
                div = make_axes_locatable(ax)
                ax_top = div.append_axes("top", size="8%", pad=0.06, sharex=ax)
                ax_bot = div.append_axes("bottom", size="8%", pad=0.45, sharex=ax)
                _draw_ss_track(ax_top, list(ss1), n, f"{tagA or 'fold1'} SS ")
                _draw_ss_track(ax_bot, ss2_on1, n, f"{tagB or 'fold2'} SS ")
                ax_top.set_title("secondary structure: "
                                 "■ helix  ■ strand  ■ coil",
                                 fontsize=7, pad=2)
                ax_bot.set_xlabel("residue position (fold 1 frame)")
            else:
                print(f"[ddg-ss] no SS tracks for {pair_id} "
                      f"(DSSP/seq unavailable for one fold)")
        except Exception as e:
            print(f"[ddg-ss] SS tracks skipped for {pair_id}: {e}")

        fig.tight_layout()
        try:
            fig.savefig(out_path, dpi=170)
        finally:
            plt.close(fig)
        return out_path
    except Exception:
        try:
            plt.close("all")
        except Exception:
            pass
        return None


if __name__ == "__main__":
    import sys
    pid = sys.argv[1] if len(sys.argv) > 1 else "xxxx_yyyy"
    print(per_residue_ddg_fig(pid))
