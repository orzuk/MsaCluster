"""Shared per-cluster scoring vs the two truth folds (F1, F2).

Single source of truth for "given a list of predicted PDBs for one pair,
TM-align each against the two truth structures and return a standardized
DataFrame." Used by every v2-era method scorer (AF2/AF3, ESMFold2,
Boltz-2-rescoring, future methods).

Canonical output schema (matches existing pilot CSVs):
    mode, version, cluster, chain, TM_F1, TM_F2, TMdiff,
    TMdiff_centered, pred_pdb, error
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd

from utils.align_utils import tmalign_unified


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
    """Run TM-align and return max(by_1, by_2). Returns nan on failure."""
    res = tmalign_unified(str(pred_pdb), str(truth_pdb),
                         chain1=None, chain2=truth_chain)
    return max(res.get("tm_by_1") or 0.0, res.get("tm_by_2") or 0.0)


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
        # Assume CWD is repo root or that Pipeline/FoldPairs is on the path
        pair_dir = Path("Pipeline") / "FoldPairs" / pair_id

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
