#!/usr/bin/env python3
"""
Compute intra-cluster prediction pairwise TM-score for fold-switching pairs,
matching the per-protein df_diversity.csv schema used by run_controls_pipeline.py
(task_postprocess).

For each fold-switching pair, collect the per-cluster AF2 ShallowMsa predictions
(at output_AF/AF2/ShallowMsa_*.pdb, one PDB per cluster per chain), compute
pairwise TM-score between every prediction pair, and write summary statistics
matching the controls schema:
  pair_id, n_clusters, n_total_pdbs, mean_pairwise_TM, median_pairwise_TM,
  min_pairwise_TM, max_pairwise_TM, std_pairwise_TM

Notes on schema parity vs. controls/df_diversity.csv:
  - Controls has columns: protein_id, control_type, n_clusters, n_total_pdbs,
    mean_pairwise_TM, median_pairwise_TM, min_pairwise_TM, max_pairwise_TM,
    std_pairwise_TM, cluster_mean_TM, cluster_median_TM.
  - Here we use `pair_id` instead of `protein_id` and omit `control_type`.
  - We also omit `cluster_mean_TM` / `cluster_median_TM` — for fold-switching
    pairs we restrict input to ShallowMsa only (no DeepMsa, no ESM), so those
    would be identical to mean/median anyway.

Output: docs/foldpair_intra_diversity.csv (one row per pair).
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import DATA_DIR, PAIR_DIR_RE
from utils.align_utils import compute_tmscore_align


def collect_pair_pdbs(pair_dir: Path) -> list:
    """Return list of (label, pdb_path) for AF2 ShallowMsa cluster predictions.

    Filenames at top level of output_AF/AF2/ look like
    'ShallowMsa_NNN__<chain>.pdb' (one per cluster per chain).
    DeepMsa predictions are intentionally skipped (the diversity question is
    about per-cluster predictions); this matches `shallow_labels` filtering in
    run_controls_pipeline.task_postprocess.
    """
    af2_dir = pair_dir / "output_AF" / "AF2"
    if not af2_dir.is_dir():
        return []
    out = []
    for pdb in sorted(af2_dir.glob("ShallowMsa_*.pdb")):
        out.append((pdb.stem, str(pdb)))
    return out


def compute_pairwise_tm(pdbs: list) -> np.ndarray:
    """Pairwise TM-score matrix (symmetric, diagonal=1)."""
    n = len(pdbs)
    M = np.eye(n, dtype=np.float32)
    for i in range(n):
        for j in range(i + 1, n):
            try:
                tm = compute_tmscore_align(pdbs[i][1], pdbs[j][1],
                                           chain1=None, chain2=None)
                M[i, j] = M[j, i] = float(tm)
            except Exception as e:
                print(f"  WARN: TM-align {pdbs[i][0]} vs {pdbs[j][0]}: {e}",
                      flush=True)
                M[i, j] = M[j, i] = np.nan
    return M


def summarize(pair_id: str, pdbs: list, M: np.ndarray) -> dict:
    n = len(pdbs)
    base = {
        "pair_id": pair_id,
        "n_clusters": n,        # match controls: count of ShallowMsa PDBs (not unique tags)
        "n_total_pdbs": n,
        "mean_pairwise_TM": np.nan,
        "median_pairwise_TM": np.nan,
        "min_pairwise_TM": np.nan,
        "max_pairwise_TM": np.nan,
        "std_pairwise_TM": np.nan,
    }
    if n < 2:
        base["n_clusters"] = 0
        return base
    mask = np.triu(np.ones_like(M, dtype=bool), k=1)
    vals = M[mask]
    vals = vals[~np.isnan(vals)]
    if len(vals) == 0:
        return base
    base.update({
        "mean_pairwise_TM": float(np.mean(vals)),
        "median_pairwise_TM": float(np.median(vals)),
        "min_pairwise_TM": float(np.min(vals)),
        "max_pairwise_TM": float(np.max(vals)),
        "std_pairwise_TM": float(np.std(vals)),
    })
    return base


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pairs", default="ALL",
                   help="Comma-separated pair IDs or ALL")
    p.add_argument("--output", default="docs/foldpair_intra_diversity.csv")
    args = p.parse_args()

    if args.pairs == "ALL":
        pair_ids = [d.name for d in Path(DATA_DIR).iterdir()
                    if d.is_dir() and PAIR_DIR_RE.match(d.name)]
    else:
        pair_ids = [s.strip() for s in args.pairs.split(",") if s.strip()]

    rows = []
    for pid in sorted(pair_ids):
        pair_dir = Path(DATA_DIR) / pid
        pdbs = collect_pair_pdbs(pair_dir)
        if len(pdbs) < 2:
            print(f"[skip] {pid}: <2 cluster PDBs", flush=True)
            continue
        print(f"[run] {pid}: {len(pdbs)} cluster PDBs", flush=True)
        M = compute_pairwise_tm(pdbs)
        rows.append(summarize(pid, pdbs, M))

    df = pd.DataFrame(rows)
    out = args.output
    out_dir = os.path.dirname(out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"\nWrote {out} ({len(df)} pairs)")
    if not df.empty:
        print(df[["pair_id", "n_clusters", "mean_pairwise_TM",
                  "median_pairwise_TM"]].head(20).to_string(index=False))


if __name__ == "__main__":
    main()
