"""Single source of truth for extracting a representative sequence
from a cluster MSA.

Used by AF2/AF3 (via run_foldswitch_pipeline.py), Boltz-2 (run_Boltz2.py),
ESMFold2 (run_ESMFold2.py), and any other method that takes ONE
sequence per cluster and feeds it to a structure/sequence predictor.

API:
    get_cluster_representative(a3m_path, method='medoid')
        -> (header, sequence)

Internally:
    representative_from_rows(rows, method='medoid') -> sequence

`rows` is a List[Tuple[header, aligned_sequence]] as returned by
utils.protein_utils.read_msa.

Supported methods:
    'medoid'    : row with minimum sum-of-Hamming-distances to others.
                  Default. Biologically meaningful: the most
                  "central" cluster member.
    'consensus' : per-column majority vote (match-state columns only).
                  Returns a chimera that may not exist in the cluster.
    'first'     : first row in the file. Useful when the upstream
                  clustering has already placed the chosen query at
                  row 0 (e.g., after --query_type medoid swap).
"""
from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np


def _ungap_upper(s: str) -> str:
    """Strip insertions/gaps; uppercase."""
    return "".join(ch for ch in (s or "") if ch.isalpha()).upper()


def representative_from_rows(rows: List[Tuple[str, str]],
                              method: str = "medoid",
                              max_rows_for_medoid: int = 200,
                              ) -> Optional[str]:
    """Compute a representative sequence from a list of (header, aligned_seq).

    Returns the ungapped, uppercased representative sequence, or None if
    rows is empty or all sequences have inconsistent lengths.

    For 'medoid', sub-samples to max_rows_for_medoid (default 200) for
    tractability on large clusters; the result is deterministic
    (chooses the row whose sum-of-Hamming to the others in the subsample
    is minimum, breaking ties by row index).
    """
    if not rows:
        return None
    aligned = [aln for _, aln in rows if aln]
    if not aligned:
        return None
    L = len(aligned[0])
    aligned = [a for a in aligned if len(a) == L]
    if not aligned:
        return None

    if method == "first":
        return _ungap_upper(aligned[0])

    if method == "consensus":
        cols: list[str] = []
        for j in range(L):
            col_counts: dict[str, int] = {}
            for a in aligned:
                c = a[j]
                if c == "-" or c == ".":
                    continue
                col_counts[c] = col_counts.get(c, 0) + 1
            if not col_counts:
                cols.append("-")
            else:
                cols.append(max(col_counts, key=col_counts.get))
        return "".join(c for c in cols if c.isalpha()).upper()

    if method == "medoid":
        # Sub-sample for tractability on huge clusters
        sub = aligned[:max_rows_for_medoid]
        arr = np.array([list(s) for s in sub])
        n = arr.shape[0]
        dist = np.zeros((n, n), dtype=np.int32)
        for i in range(n):
            dist[i] = (arr[i] != arr).sum(axis=1)
        idx_best = int(dist.sum(axis=1).argmin())
        return _ungap_upper(sub[idx_best])

    raise ValueError(f"Unknown representative method: {method!r}. "
                     f"Expected one of: 'medoid', 'consensus', 'first'.")


def get_cluster_representative(a3m_path,
                                method: str = "medoid",
                                ) -> Tuple[Optional[str], Optional[str]]:
    """Extract a single representative from a cluster .a3m file.

    Returns (header, sequence). The header is the FIRST record's id from
    the file (provided as a stable label even though the representative
    sequence may come from a different row, especially for 'consensus' or
    when 'medoid' picks an inner row). The sequence is ungapped and
    uppercased.

    Raises FileNotFoundError if the file doesn't exist. Returns
    (None, None) if the file is empty or unparseable.
    """
    from utils.protein_utils import read_msa  # local import to avoid cycle

    a3m_path = Path(a3m_path)
    if not a3m_path.is_file():
        raise FileNotFoundError(f"Cluster a3m not found: {a3m_path}")
    rows = read_msa(str(a3m_path))
    if not rows:
        return None, None
    header = rows[0][0].split()[0]
    seq = representative_from_rows(rows, method=method)
    return header, seq
