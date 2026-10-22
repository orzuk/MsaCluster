"""Provenance manifest + top-L binarisation for per-cluster contact predictions.

Background
----------
MSA-Transformer (MSAT) and CCMpred outputs live in
``Pipeline/FoldPairs/<pair>/output_cmaps/{msa_transformer,ccmpred}/`` and are
named after the cluster file they were computed from (``ShallowMsa_007`` etc.).
When the fine clustering is re-run, the labels are re-used for DIFFERENT
sequence sets, so a file name alone says nothing about which cluster an output
belongs to. Old outputs were therefore silently joined to new clusters.

The fix: every contact run records, per output file, a hash of the input A3M
it was computed from (``output_cmaps/<method>/_manifest.json``). Scoring and
the survey only use outputs whose recorded hash equals the hash of the CURRENT
cluster file, and whose file was (re)written after the run was recorded.
Outputs without a manifest entry (legacy) are treated as stale.

Manifest layout::

    {"version": 1, "method": "ccmpred",
     "entries": {"ShallowMsa_007.ccmpred.npz": {
         "input_a3m": "output_msa_cluster/ShallowMsa_007.a3m",   # relative to pair dir
         "cluster_hash": "<sha256>", "n_seqs": 57,
         "resolution": "fine", "run_started": 1760000000.0}}}
"""
from __future__ import annotations

import hashlib
import json
import os
import time
from typing import Iterable, Optional, Tuple

import numpy as np

MANIFEST_NAME = "_manifest.json"
MANIFEST_VERSION = 1
# Allowed clock skew (seconds) between the recorded run start and the output
# file's mtime (NFS clocks on compute nodes can differ slightly).
MTIME_TOLERANCE_S = 120.0

METHOD_SUBDIRS = {"msa_transformer": "msa_transformer", "ccmpred": "ccmpred"}

# Top-L factor used when scoring each method (None => probability threshold).
# CCMpred scores are APC-corrected couplings (~0.02-0.3), not probabilities,
# so a fixed 0.4 threshold yields no contacts; rank-based top-L is standard.
DEFAULT_TOP_L_FACTOR = {"ccmpred": 1.0, "msa_transformer": None}


# ---------------------------------------------------------------------------
# Hashing
# ---------------------------------------------------------------------------
def _iter_a3m_records(path: str) -> Iterable[Tuple[str, str]]:
    hdr, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if hdr is not None:
                    yield hdr, "".join(buf)
                hdr, buf = line[1:], []
            elif hdr is not None:
                buf.append(line)
    if hdr is not None:
        yield hdr, "".join(buf)


def a3m_cluster_hash(path: str) -> Tuple[str, int]:
    """Return (sha256, n_seqs) of the cluster's membership.

    The hash is over the SORTED set of "<first header token>\\t<sequence>"
    lines, so it identifies the sequence set independent of row order and of
    header annotations after the ID. Not cached: a re-clustering can rewrite a
    file within the filesystem's mtime granularity and keep its size.
    """
    items = sorted(f"{(h.split() or [''])[0]}\t{s}" for h, s in _iter_a3m_records(path))
    h = hashlib.sha256("\n".join(items).encode()).hexdigest()
    return h, len(items)


# ---------------------------------------------------------------------------
# Output naming
# ---------------------------------------------------------------------------
def expected_output_names(method: str, stem: str) -> list:
    """Output file name(s) a contact run writes for input A3M <stem>.a3m."""
    if method == "ccmpred":
        return [f"{stem}.ccmpred.npz"]
    if method == "msa_transformer":
        if stem == "DeepMsa":
            # task_cmap_msa_transformer runs run_MSATrans.py --keyword DeepMsa
            # on an alias copy named DeepMsa.a3m.
            return ["msa_t_DeepMsa_DeepMsa.npz"]
        return [f"msa_t_clusters_{stem}.npz"]
    raise ValueError(f"unknown contact method {method!r}")


def output_resolution(fname: str) -> str:
    """'fine' | 'coarse' | 'deep' | 'other' from an output file name."""
    if "DeepMsa" in fname or "MSA_deep" in fname:
        return "deep"
    if "ShallowMsa_" in fname:
        return "fine"
    if "CoarseMsa_" in fname:
        return "coarse"
    return "other"


# ---------------------------------------------------------------------------
# Manifest I/O
# ---------------------------------------------------------------------------
def method_dir(pair_dir: str, method: str) -> str:
    return os.path.join(pair_dir, "output_cmaps", METHOD_SUBDIRS[method])


def load_manifest(mdir: str) -> dict:
    path = os.path.join(mdir, MANIFEST_NAME)
    if not os.path.isfile(path):
        return {"version": MANIFEST_VERSION, "entries": {}}
    try:
        with open(path) as fh:
            man = json.load(fh)
    except (OSError, ValueError):
        return {"version": MANIFEST_VERSION, "entries": {}}
    man.setdefault("entries", {})
    return man


def save_manifest(mdir: str, man: dict) -> None:
    os.makedirs(mdir, exist_ok=True)
    path = os.path.join(mdir, MANIFEST_NAME)
    tmp = f"{path}.tmp.{os.getpid()}"
    with open(tmp, "w") as fh:
        json.dump(man, fh, indent=1, sort_keys=True)
    os.replace(tmp, path)


def record_run(pair_dir: str, method: str, a3m_path: str, resolution: str,
               man: Optional[dict] = None, run_started: Optional[float] = None,
               save: bool = True) -> dict:
    """Record (BEFORE running) that <a3m_path> is being processed by <method>.

    Writes one entry per expected output name. Returns the manifest.
    """
    mdir = method_dir(pair_dir, method)
    if man is None:
        man = load_manifest(mdir)
    man["version"] = MANIFEST_VERSION
    man["method"] = method
    h, n = a3m_cluster_hash(a3m_path)
    stem = os.path.splitext(os.path.basename(a3m_path))[0]
    rel = os.path.relpath(os.path.abspath(a3m_path), os.path.abspath(pair_dir))
    t0 = time.time() if run_started is None else float(run_started)
    for out_name in expected_output_names(method, stem):
        man["entries"][out_name] = {
            "input_a3m": rel.replace(os.sep, "/"),
            "cluster_hash": h,
            "n_seqs": n,
            "resolution": resolution,
            "run_started": t0,
        }
    if save:
        save_manifest(mdir, man)
    return man


def entry_status(pair_dir: str, method: str, out_name: str,
                 man: Optional[dict] = None) -> Tuple[bool, str, Optional[dict]]:
    """Is output <out_name> fresh w.r.t. the CURRENT cluster file?

    Returns (fresh, reason, entry). reason is 'ok' or one of
    'no_manifest_entry' | 'input_missing' | 'hash_mismatch' |
    'output_missing' | 'output_predates_run'.
    """
    mdir = method_dir(pair_dir, method)
    if man is None:
        man = load_manifest(mdir)
    entry = man.get("entries", {}).get(out_name)
    if not entry:
        return False, "no_manifest_entry", None
    in_path = os.path.join(pair_dir, entry.get("input_a3m", ""))
    if not os.path.isfile(in_path):
        return False, "input_missing", entry
    h, _ = a3m_cluster_hash(in_path)
    if h != entry.get("cluster_hash"):
        return False, "hash_mismatch", entry
    out_path = os.path.join(mdir, out_name)
    if not (os.path.isfile(out_path) and os.path.getsize(out_path) > 0):
        return False, "output_missing", entry
    if os.path.getmtime(out_path) < float(entry.get("run_started", 0)) - MTIME_TOLERANCE_S:
        return False, "output_predates_run", entry
    return True, "ok", entry


def a3m_is_fresh(pair_dir: str, method: str, a3m_path: str,
                 man: Optional[dict] = None) -> bool:
    """True iff every output expected for <a3m_path> is fresh and was computed
    from exactly this file (path + content)."""
    stem = os.path.splitext(os.path.basename(a3m_path))[0]
    rel = os.path.relpath(os.path.abspath(a3m_path),
                          os.path.abspath(pair_dir)).replace(os.sep, "/")
    if man is None:
        man = load_manifest(method_dir(pair_dir, method))
    for out_name in expected_output_names(method, stem):
        ok, _, entry = entry_status(pair_dir, method, out_name, man)
        if not ok or entry.get("input_a3m") != rel:
            return False
    return True


def filter_fresh_cmap_rows(df, pair_dir: str, method: str, verbose: bool = True):
    """Keep only df_cmap rows whose output is fresh AND whose recorded
    cluster_hash equals the current one. Legacy CSVs (no cluster_hash column)
    yield an empty frame."""
    if df is None or len(df) == 0:
        return df
    if "cluster_hash" not in df.columns or "file" not in df.columns:
        if verbose:
            print(f"[contacts] {os.path.basename(pair_dir)} {method}: legacy df_cmap "
                  f"(no cluster_hash) -> ignored; re-score with --force_cmap TRUE")
        return df.iloc[0:0]
    man = load_manifest(method_dir(pair_dir, method))
    keep = []
    for fn, h in zip(df["file"].astype(str), df["cluster_hash"].astype(str)):
        ok, _, entry = entry_status(pair_dir, method, fn, man)
        keep.append(bool(ok and entry.get("cluster_hash") == h))
    n_drop = len(keep) - sum(keep)
    if n_drop and verbose:
        print(f"[contacts] {os.path.basename(pair_dir)} {method}: dropped {n_drop} "
              f"stale row(s) (cluster changed since contacts were computed)")
    return df[np.asarray(keep, dtype=bool)]


# ---------------------------------------------------------------------------
# Top-L binarisation
# ---------------------------------------------------------------------------
def top_l_binarize(P: np.ndarray, k: int, sep_min: int = 6,
                   idx: Optional[np.ndarray] = None) -> np.ndarray:
    """Binary symmetric map with the k highest-scoring pairs (i<j) set.

    P is symmetrised (mean of P and P.T) first; only pairs with sequence
    separation >= sep_min are eligible. Separation is |i-j| on P's axis, or
    |idx[i]-idx[j]| when the seed-frame column index ``idx`` is given.
    NaNs are never selected. Ties are broken by (i, j) order (stable).
    """
    P = np.asarray(P, dtype=float)
    n = P.shape[0]
    out = np.zeros((n, n), dtype=np.uint8)
    if n == 0 or k <= 0:
        return out
    S = 0.5 * (P + P.T)
    iu, ju = np.triu_indices(n, k=1)
    pos = np.arange(n) if idx is None else np.asarray(idx)[:n]
    sep_ok = np.abs(pos[ju] - pos[iu]) >= max(1, int(sep_min))
    vals = S[iu, ju]
    ok = sep_ok & np.isfinite(vals)
    iu, ju, vals = iu[ok], ju[ok], vals[ok]
    if vals.size == 0:
        return out
    k = min(int(k), vals.size)
    order = np.argsort(-vals, kind="stable")[:k]
    out[iu[order], ju[order]] = 1
    out[ju[order], iu[order]] = 1
    return out


def resolve_top_l_factor(method: str, top_l_factor="auto"):
    if isinstance(top_l_factor, str) and top_l_factor.lower() == "auto":
        return DEFAULT_TOP_L_FACTOR.get(method)
    if top_l_factor is None:
        return None
    if isinstance(top_l_factor, str) and top_l_factor.lower() in ("none", ""):
        return None
    f = float(top_l_factor)
    return f if f > 0 else None

