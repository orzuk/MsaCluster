"""Shared per-cluster scoring vs the two truth folds (F1, F2).

Single source of truth for "given a list of predicted PDBs for one pair,
TM-align each against the two truth structures and return a standardized
DataFrame." Used by every v2-era method scorer (AF2/AF3, ESMFold2,
Boltz-2-rescoring, future methods).

Canonical output schema (matches existing pilot CSVs):
    mode, version, cluster, chain, TM_F1, TM_F2, TMdiff,
    TMdiff_centered, pred_pdb, error

Also home of the canonical single-prediction helpers ``truth_chain_pdb``
(locate one fold's truth PDB) and ``tm_vs_truth`` / ``safe_tm_vs_truth``
(TM-score of a prediction against one truth chain, normalised by the truth
length), used by the per-method runners and the AF2/AF3 scorer.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd

from utils.align_utils import compute_tmscore_align


def truth_chain_pdb(pair_dir, pdbid: str, chain: str) -> str:
    """Truth PDB for one fold (pdbid, chain) of a pair.

    Tries, in order:  <pair_dir>/chain_pdb_files/<pdbid><chain>.pdb
                      <pair_dir>/<pdbid>.pdb
                      <pair_dir>/<pdbid>_cif.pdb
    Returns "" if none exists. The full deposited PDB can hold both folds'
    chains (same-PDB pairs), so callers must still pass ``chain`` to the TM
    call (``tm_vs_truth`` requires it).
    """
    pair_dir = Path(pair_dir)
    for c in (pair_dir / "chain_pdb_files" / f"{pdbid}{chain}.pdb",
              pair_dir / f"{pdbid}.pdb",
              pair_dir / f"{pdbid}_cif.pdb"):
        if c.is_file():
            return str(c)
    return ""


def tm_vs_truth(pred_pdb, truth_pdb, truth_chain: str,
                pred_chain: Optional[str] = None) -> float:
    """TM-score of a prediction against one truth chain, normalised by the
    truth chain's length.

    The truth goes first (TM-align chain 1) restricted to ``truth_chain``, so
    the ``by_1`` normalisation is by the truth length. ``pred_chain``
    restricts the prediction to one chain (holo outputs are complexes; the
    protein copy is chain A); None takes the prediction's first chain.
    Raises on failure; see ``safe_tm_vs_truth`` for a NaN-returning variant.
    """
    if not truth_chain:
        raise ValueError("tm_vs_truth: truth_chain is required")
    return float(compute_tmscore_align(str(truth_pdb), str(pred_pdb),
                                       chain1=truth_chain, chain2=pred_chain,
                                       norm="by_1"))


def safe_tm_vs_truth(pred_pdb, truth_pdb, truth_chain: str,
                     pred_chain: Optional[str] = None) -> float:
    """``tm_vs_truth`` that returns NaN when either file is missing or the
    alignment fails."""
    if not truth_pdb or not pred_pdb or not Path(pred_pdb).is_file():
        return float("nan")
    try:
        return tm_vs_truth(pred_pdb, truth_pdb, truth_chain, pred_chain)
    except Exception:
        return float("nan")


def find_truth_pdbs(pair_id: str, pair_dir: Path
                     ) -> tuple[tuple[Path, Path], tuple[str, str], tuple[str, str]]:
    """Locate the two truth PDB files for a fold pair.

    Returns ((truth_F1_path, truth_F2_path), (chain_F1, chain_F2),
             (pdb_id_F1, pdb_id_F2)).

    Tries, in order:  <pair_dir>/<pdb>.pdb
                     <pair_dir>/<pdb>_cif.pdb
                     <pair_dir>/chain_pdb_files/<pdb>.pdb

    Raises ValueError on bad pair_id, FileNotFoundError if either PDB
    is missing.
    """
    m = re.match(r"^([0-9a-zA-Z]{4})([A-Za-z])_([0-9a-zA-Z]{4})([A-Za-z])$",
                 pair_id)
    if not m:
        raise ValueError(f"Cannot parse pair_id {pair_id!r} "
                         f"(expected <pdb><chain>_<pdb><chain>)")
    pdb1, ch1, pdb2, ch2 = m.groups()
    found: list[Path] = []
    for pdb in (pdb1, pdb2):
        for candidate in (pair_dir / f"{pdb}.pdb",
                          pair_dir / f"{pdb}_cif.pdb",
                          pair_dir / "chain_pdb_files" / f"{pdb}.pdb"):
            if candidate.is_file():
                found.append(candidate)
                break
        else:
            raise FileNotFoundError(f"Truth PDB for {pdb} not found "
                                    f"under {pair_dir}")
    return (tuple(found), (ch1, ch2), (pdb1, pdb2))


def _tm_pair(pred_pdb: Path, truth_pdb: Path, truth_chain: str) -> float:
    """TM-score of a prediction vs one truth chain, normalised by the truth
    length (``tm_vs_truth``). Raises on failure; the caller records the
    error in the row."""
    return tm_vs_truth(pred_pdb, truth_pdb, truth_chain)


def score_predictions_vs_truth(
    pair_id: str,
    pred_paths: list[Path],
    cluster_ids: list[str],
    mode_label: str,
    version_label: str,
    pair_dir: Optional[Path] = None,
    chain_label: Optional[str] = None,
    add_centered: bool = True,
    repo_root: Optional[Path] = None,
) -> pd.DataFrame:
    """TM-align each prediction against the pair's two truth folds.

    Parameters:
      pair_id        : "5jytA_2qkeE" etc.
      pred_paths     : prediction PDBs (one per cluster) -- parallel to cluster_ids
      cluster_ids    : cluster tags (e.g., "ShallowMsa_000")
      mode_label     : value written to the 'mode' column (e.g., "AF2_newmethod_top10")
      version_label  : "AF2" / "AF3" / "ESMFold2" / "Boltz-2" / ...
      pair_dir       : optional; defaults to Pipeline/FoldPairs/<pair_id>
      chain_label    : value for the 'chain' column. Single-sequence methods
                       (ESMFold2) pass the F1 chain tag for compatibility;
                       multi-chain AF methods pass the chain they predicted
      add_centered   : whether to add TMdiff_centered column
      repo_root      : repo root for making pred_pdb paths relative (display only)

    Returns DataFrame with canonical schema.
    """
    if len(pred_paths) != len(cluster_ids):
        raise ValueError("pred_paths and cluster_ids must be parallel lists")

    if pair_dir is None:
        # config.DATA_DIR = Pipeline/FoldPairs unless MSACLUSTER_DATA_DIR overrides it
        from config import DATA_DIR
        pair_dir = Path(DATA_DIR) / pair_id

    (truth_F1, truth_F2), (ch1, ch2), (pdb1, pdb2) = find_truth_pdbs(pair_id, pair_dir)
    nominal_chain = chain_label if chain_label is not None else f"{pdb1}{ch1}"

    rows = []
    for cid, pred in zip(cluster_ids, pred_paths):
        pred = Path(pred)
        if not pred.is_file():
            rows.append({"mode": mode_label, "version": version_label,
                         "cluster": cid, "chain": nominal_chain,
                         "TM_F1": None, "TM_F2": None, "TMdiff": None,
                         "pred_pdb": str(pred), "error": "missing_pred_pdb"})
            continue
        try:
            tm_F1 = _tm_pair(pred, truth_F1, ch1)
            tm_F2 = _tm_pair(pred, truth_F2, ch2)
            pred_label = (str(pred.relative_to(repo_root))
                          if repo_root and pred.is_absolute()
                          else str(pred))
            rows.append({"mode": mode_label, "version": version_label,
                         "cluster": cid, "chain": nominal_chain,
                         "TM_F1": tm_F1, "TM_F2": tm_F2,
                         "TMdiff": tm_F1 - tm_F2,
                         "pred_pdb": pred_label, "error": ""})
        except Exception as e:
            rows.append({"mode": mode_label, "version": version_label,
                         "cluster": cid, "chain": nominal_chain,
                         "TM_F1": None, "TM_F2": None, "TMdiff": None,
                         "pred_pdb": str(pred), "error": str(e)[:120]})

    df = pd.DataFrame(rows)
    if add_centered and not df.empty and "TMdiff" in df.columns:
        # Center per (mode, chain) so AF2 / AF3 / ESMFold2 are each centered separately
        df["TMdiff_centered"] = df.groupby(["mode", "chain"])["TMdiff"].transform(
            lambda g: g - g.median()
        )
    return df
