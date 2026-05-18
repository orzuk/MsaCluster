#!/usr/bin/env python3
"""run_CFrandom.py -- per-cluster CF-random ensemble prediction.

For each ShallowMsa_XXX cluster of a fold-switching pair, run ColabFold
many times with RANDOMLY SAMPLED shallow MSA depths (max_seq, max_extra_seq).
This is the Lee/Schafer/Porter "CF-random" idea applied at the cluster level:
instead of subsampling the whole deep MSA, we subsample within each clade-
representing cluster and ask whether different clades produce different
conformational ensembles.

For each cluster we get a distribution over conformations rather than a
single point estimate. Per-cluster summary statistics (mean/median TM to
F1/F2, ensemble bimodality, fraction favoring each fold) feed downstream
fold-diversity / concordance / phylogenetic-placement analyses.

Inputs
  Pipeline/<pair>/output_msa_cluster/ShallowMsa_*.a3m

Outputs
  Pipeline/<pair>/output_cfrandom/<cluster_tag>/run_<i>_max_<N>__<pdbid><chain>.pdb
      one PDB per (cluster, sample i, max_seq depth N, target fold reference chain)
  Pipeline/<pair>/output_cfrandom/<cluster_tag>/manifest.csv
      per-sample log: cluster, run_idx, max_seq, max_extra_seq, mean_plddt,
      tm_to_fold1, tm_to_fold2
  Pipeline/<pair>/Analysis/df_cfrandom.csv
      per-cluster ensemble summary (one row per cluster per fold reference)

CF-random subsampling schedule (defaults match Lee 2025 Figure 2b):
  (max_seq, max_extra_seq) drawn uniformly from a small set of shallow
  combinations; pass --depths "2:4,4:8,8:16,16:32" to customise. Default is
  five settings spanning very-shallow to moderate.

Venv: my-python-venv (reuses colabfold 1.2.0 already installed). No new
installs required.

Usage
  python3 run_CFrandom.py -input 6c6sD_2ougC --n_samples 100
  python3 run_CFrandom.py -input 6c6sD_2ougC --n_samples 50 --depths 2:4,4:8
"""

from __future__ import annotations

import argparse
import csv
import os
import random
import re
import shutil
import subprocess
import sys
import time
from glob import glob
from pathlib import Path
from typing import List, Tuple

import numpy as np

# Repo root importable
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _THIS_DIR)

from config import DATA_DIR  # noqa: E402
from utils.protein_utils import compute_tmscore_align  # noqa: E402


def _parse_depths(s: str) -> List[Tuple[int, int]]:
    """'2:4,4:8,8:16,16:32,32:64' -> [(2,4),(4,8),(8,16),(16,32),(32,64)]"""
    out = []
    for tok in s.split(","):
        tok = tok.strip()
        if not tok:
            continue
        a, b = tok.split(":")
        out.append((int(a), int(b)))
    return out


def _truth_chain_pdb(pair_dir: Path, pdbid: str, chain: str) -> str:
    """Best-effort path to truth PDB for one fold's reference chain."""
    cand = [
        pair_dir / f"{pdbid}.pdb",
        pair_dir / "chain_pdb_files" / f"{pdbid}{chain}.pdb",
        pair_dir / f"{pdbid}_cif.pdb",
    ]
    for c in cand:
        if c.is_file():
            return str(c)
    return ""


def _run_colabfold(a3m_in: Path, out_dir: Path, max_seq: int,
                   max_extra_seq: int, num_recycle: int = 3) -> Path:
    """Invoke ColabFold once with the given shallow subsampling.

    Returns the path to the best PDB written (or None on failure).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    # Reuse the same colabfold_batch entry point the rest of the pipeline uses
    cmd = (
        f"colabfold_batch "
        f"--max-seq {max_seq} --max-extra-seq {max_extra_seq} "
        f"--num-recycle {num_recycle} "
        f"--model-type alphafold2_ptm "
        f"--num-models 1 "
        f"{a3m_in} {out_dir}"
    )
    try:
        subprocess.run(cmd, shell=True, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(f"    [cfrandom] colabfold failed for {a3m_in.name}: "
              f"{e.stderr.decode()[:200]}", file=sys.stderr)
        return None
    # Pick the rank_001 PDB output
    pdbs = sorted(out_dir.glob("*_unrelaxed_rank_001_*.pdb"))
    if not pdbs:
        pdbs = sorted(out_dir.glob("*.pdb"))
    return pdbs[0] if pdbs else None


def _mean_plddt(pdb_path: Path) -> float:
    """Compute mean B-factor (= pLDDT for AF2 outputs) over CA atoms."""
    vals = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    vals.append(float(line[60:66]))
                except ValueError:
                    pass
    return float(np.mean(vals)) if vals else float("nan")


def _safe_tm(pred_pdb: Path, truth_pdb: str) -> float:
    if not truth_pdb or not pred_pdb or not Path(pred_pdb).is_file():
        return float("nan")
    try:
        tm = compute_tmscore_align(str(pred_pdb), str(truth_pdb))
        if isinstance(tm, (list, tuple)):
            tm = float(tm[0])
        return float(tm)
    except Exception:
        return float("nan")


def run_one_cluster(pair_id: str, cluster_a3m: Path, pdbid1: str, chain1: str,
                    pdbid2: str, chain2: str, n_samples: int,
                    depths: List[Tuple[int, int]], pair_dir: Path,
                    seed: int) -> List[dict]:
    """Run n_samples CF-random predictions on one cluster MSA; return rows."""
    rng = random.Random(seed)
    tag = cluster_a3m.stem
    out_root = pair_dir / "output_cfrandom" / tag
    out_root.mkdir(parents=True, exist_ok=True)

    truth1 = _truth_chain_pdb(pair_dir, pdbid1, chain1)
    truth2 = _truth_chain_pdb(pair_dir, pdbid2, chain2)

    rows = []
    print(f"  [cfrandom] {tag}: {n_samples} samples")
    for i in range(n_samples):
        ms, mes = rng.choice(depths)
        sub_out = out_root / f"run_{i:04d}_max_{ms}_{mes}"
        if (sub_out / "manifest.flag").exists():
            # Already done
            continue
        pdb = _run_colabfold(cluster_a3m, sub_out, ms, mes)
        if not pdb:
            continue
        plddt = _mean_plddt(pdb)
        tm1 = _safe_tm(pdb, truth1)
        tm2 = _safe_tm(pdb, truth2)
        rows.append({
            "pair_id": pair_id,
            "cluster": tag,
            "run_idx": i,
            "max_seq": ms,
            "max_extra_seq": mes,
            "mean_plddt": plddt,
            "tm_fold1": tm1,
            "tm_fold2": tm2,
            "TMdiff": (tm1 - tm2) if (tm1 == tm1 and tm2 == tm2) else float("nan"),
            "pdb_path": str(pdb),
        })
        # Mark sub-run done
        (sub_out / "manifest.flag").touch()
    # Per-cluster manifest
    if rows:
        with open(out_root / "manifest.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
    return rows


def summarize_per_cluster(rows: List[dict]) -> List[dict]:
    """Group manifest rows -> per-cluster ensemble statistics."""
    by_cluster: dict = {}
    for r in rows:
        by_cluster.setdefault(r["cluster"], []).append(r)
    out = []
    for cluster, items in by_cluster.items():
        tm1 = np.array([x["tm_fold1"] for x in items if x["tm_fold1"] == x["tm_fold1"]])
        tm2 = np.array([x["tm_fold2"] for x in items if x["tm_fold2"] == x["tm_fold2"]])
        diff = np.array([x["TMdiff"] for x in items if x["TMdiff"] == x["TMdiff"]])
        if len(diff) == 0:
            continue
        frac_f1 = float((diff > 0.05).sum()) / len(diff)
        frac_f2 = float((diff < -0.05).sum()) / len(diff)
        out.append({
            "pair_id": items[0]["pair_id"],
            "cluster": cluster,
            "n_samples": len(items),
            "mean_tm_fold1": float(np.mean(tm1)) if len(tm1) else float("nan"),
            "mean_tm_fold2": float(np.mean(tm2)) if len(tm2) else float("nan"),
            "median_TMdiff": float(np.median(diff)),
            "std_TMdiff": float(np.std(diff)),
            "frac_F1_pref": frac_f1,
            "frac_F2_pref": frac_f2,
            "frac_ambiguous": 1.0 - frac_f1 - frac_f2,
            "mean_plddt": float(np.mean([x["mean_plddt"] for x in items
                                         if x["mean_plddt"] == x["mean_plddt"]])),
        })
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-input", "--input", dest="input", required=True,
                    help="Pair id, e.g. 6c6sD_2ougC")
    ap.add_argument("--n_samples", type=int, default=100,
                    help="CF-random samples per cluster (default 100)")
    ap.add_argument("--depths", default="2:4,4:8,8:16,16:32,32:64",
                    help="Shallow MSA depth combinations as comma-separated "
                         "max_seq:max_extra_seq pairs.")
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--force", action="store_true",
                    help="Recompute even if outputs exist.")
    args = ap.parse_args()

    pair_id = args.input
    pair_dir = Path(DATA_DIR) / pair_id
    if not pair_dir.is_dir():
        print(f"[error] no pair dir: {pair_dir}", file=sys.stderr)
        sys.exit(2)
    # Parse pair_id -> pdb1/chain1, pdb2/chain2
    m = re.match(r"^([0-9a-zA-Z]{4})([A-Za-z0-9])_([0-9a-zA-Z]{4})([A-Za-z0-9])$",
                 pair_id)
    if not m:
        print(f"[error] cannot parse pair_id={pair_id}", file=sys.stderr)
        sys.exit(2)
    pdbid1, chain1, pdbid2, chain2 = m.groups()

    depths = _parse_depths(args.depths)
    cluster_a3ms = sorted((pair_dir / "output_msa_cluster").glob("ShallowMsa_*.a3m"))
    if not cluster_a3ms:
        print(f"[error] no cluster MSAs under {pair_dir}/output_msa_cluster",
              file=sys.stderr)
        sys.exit(2)
    print(f"[cfrandom] pair={pair_id}: {len(cluster_a3ms)} clusters x "
          f"{args.n_samples} samples = {len(cluster_a3ms)*args.n_samples} "
          f"colabfold runs (depths={depths})")

    t0 = time.time()
    all_rows: List[dict] = []
    for i, a3m in enumerate(cluster_a3ms):
        seed_i = args.seed + i * 1009
        rows = run_one_cluster(pair_id, a3m, pdbid1, chain1, pdbid2, chain2,
                               args.n_samples, depths, pair_dir, seed_i)
        all_rows.extend(rows)
        print(f"    [cfrandom] cluster {a3m.stem}: {len(rows)} samples "
              f"completed at t={time.time()-t0:.0f}s")

    # Global manifest + per-cluster summary
    analysis_dir = pair_dir / "Analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    if all_rows:
        with open(analysis_dir / "df_cfrandom_raw.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
            w.writeheader()
            w.writerows(all_rows)
    summary = summarize_per_cluster(all_rows)
    if summary:
        with open(analysis_dir / "df_cfrandom.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
            w.writeheader()
            w.writerows(summary)
    print(f"[cfrandom] done in {time.time()-t0:.0f}s; "
          f"summary -> {analysis_dir / 'df_cfrandom.csv'}")


if __name__ == "__main__":
    main()
