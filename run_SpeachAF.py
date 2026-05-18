#!/usr/bin/env python3
"""run_SpeachAF.py -- per-cluster SPEACH-AF prediction.

SPEACH-AF (Stein 2022) escapes AlphaFold2 mode-collapse by perturbing the
INPUT MSA: it masks/mutates patches of the alignment to remove
coevolutionary information from specific regions, forcing AF2 to predict
alternative conformations. We run SPEACH-AF per cluster: for each
ShallowMsa cluster MSA we generate K perturbed variants of the MSA (random
masking patches) and call ColabFold on each.

Inputs
  Pipeline/<pair>/output_msa_cluster/ShallowMsa_*.a3m

Outputs
  Pipeline/<pair>/output_speachaf/<cluster_tag>/perturb_<i>__<pdbid><chain>.pdb
  Pipeline/<pair>/output_speachaf/<cluster_tag>/manifest.csv
  Pipeline/<pair>/Analysis/df_speachaf.csv         (per-cluster summary)

Perturbation strategy:
  Random column-masking. For each perturbation, replace a random window
  of W columns (default W=8) of K randomly-chosen sequences (default
  fraction=0.5 of sequences in the cluster) with the gap character '-'.
  Window position is drawn uniformly per perturbation; multiple
  perturbations can mask different windows.

Venv: my-python-venv (reuses colabfold 1.2.0). No new installs.

Usage
  python3 run_SpeachAF.py -input 6c6sD_2ougC --n_perturbations 20
  python3 run_SpeachAF.py -input 6c6sD_2ougC --n_perturbations 10 \
      --mask_window 8 --mask_fraction 0.5
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
from pathlib import Path
from typing import List

import numpy as np

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _THIS_DIR)

from config import DATA_DIR  # noqa: E402
from utils.protein_utils import compute_tmscore_align  # noqa: E402


def _read_a3m(path: Path) -> List[tuple]:
    entries = []
    name = None
    seq_parts: List[str] = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    entries.append((name, "".join(seq_parts)))
                name = line[1:]
                seq_parts = []
            elif line:
                seq_parts.append(line)
    if name is not None:
        entries.append((name, "".join(seq_parts)))
    return entries


def _perturb_msa(entries: List[tuple], rng: random.Random,
                 window: int, mask_fraction: float) -> List[tuple]:
    """Mask a random window of length `window` in `mask_fraction` of
    NON-query sequences. Returns a new list of (name, seq) entries.
    The query (first sequence) is left intact."""
    if len(entries) < 2:
        return entries
    L = len(entries[0][1])  # query length (a3m: query has no inserts)
    # Allow window to be larger than L for very small clusters
    if L < window + 1:
        window = max(1, L // 4)
    start = rng.randint(0, L - window) if L > window else 0
    end = start + window
    out = [entries[0]]
    n_target = max(1, int((len(entries) - 1) * mask_fraction))
    targets = set(rng.sample(range(1, len(entries)), n_target))
    for i, (name, seq) in enumerate(entries[1:], start=1):
        if i in targets:
            # Replace columns [start:end] with '-' in non-insert chars
            new_chars = []
            col = 0
            for ch in seq:
                if ch.islower() or ch == '.':
                    new_chars.append(ch)  # insert; preserve
                else:
                    if start <= col < end:
                        new_chars.append('-')
                    else:
                        new_chars.append(ch)
                    col += 1
            out.append((name, "".join(new_chars)))
        else:
            out.append((name, seq))
    return out


def _write_a3m(entries: List[tuple], path: Path) -> None:
    with open(path, "w") as f:
        for name, seq in entries:
            f.write(f">{name}\n{seq}\n")


def _run_colabfold(a3m_in: Path, out_dir: Path,
                   num_recycle: int = 3) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = (
        f"colabfold_batch "
        f"--num-recycle {num_recycle} "
        f"--model-type alphafold2_ptm "
        f"--num-models 1 "
        f"{a3m_in} {out_dir}"
    )
    try:
        subprocess.run(cmd, shell=True, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(f"    [speachaf] colabfold failed: {e.stderr.decode()[:200]}",
              file=sys.stderr)
        return None
    pdbs = sorted(out_dir.glob("*_unrelaxed_rank_001_*.pdb"))
    if not pdbs:
        pdbs = sorted(out_dir.glob("*.pdb"))
    return pdbs[0] if pdbs else None


def _mean_plddt(pdb_path: Path) -> float:
    vals = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    vals.append(float(line[60:66]))
                except ValueError:
                    pass
    return float(np.mean(vals)) if vals else float("nan")


def _safe_tm(pred_pdb, truth_pdb):
    if not truth_pdb or not pred_pdb or not Path(pred_pdb).is_file():
        return float("nan")
    try:
        tm = compute_tmscore_align(str(pred_pdb), str(truth_pdb))
        if isinstance(tm, (list, tuple)):
            tm = float(tm[0])
        return float(tm)
    except Exception:
        return float("nan")


def _truth_chain_pdb(pair_dir: Path, pdbid: str, chain: str) -> str:
    for c in (pair_dir / f"{pdbid}.pdb",
              pair_dir / "chain_pdb_files" / f"{pdbid}{chain}.pdb",
              pair_dir / f"{pdbid}_cif.pdb"):
        if c.is_file():
            return str(c)
    return ""


def run_one_cluster(pair_id, cluster_a3m, pdbid1, chain1, pdbid2, chain2,
                    n_perturbations, mask_window, mask_fraction, pair_dir,
                    seed):
    rng = random.Random(seed)
    tag = cluster_a3m.stem
    out_root = pair_dir / "output_speachaf" / tag
    out_root.mkdir(parents=True, exist_ok=True)
    truth1 = _truth_chain_pdb(pair_dir, pdbid1, chain1)
    truth2 = _truth_chain_pdb(pair_dir, pdbid2, chain2)

    entries = _read_a3m(cluster_a3m)
    if len(entries) < 2:
        print(f"  [speachaf] {tag}: too few sequences ({len(entries)}), skip")
        return []

    rows = []
    print(f"  [speachaf] {tag}: {n_perturbations} perturbations "
          f"(window={mask_window}, frac={mask_fraction})")
    for i in range(n_perturbations):
        sub_out = out_root / f"perturb_{i:04d}"
        if (sub_out / "manifest.flag").exists():
            continue
        perturbed = _perturb_msa(entries, rng, mask_window, mask_fraction)
        perturbed_a3m = sub_out.with_suffix(".a3m")
        _write_a3m(perturbed, perturbed_a3m)
        pdb = _run_colabfold(perturbed_a3m, sub_out)
        if not pdb:
            continue
        tm1 = _safe_tm(pdb, truth1)
        tm2 = _safe_tm(pdb, truth2)
        rows.append({
            "pair_id": pair_id, "cluster": tag, "perturb_idx": i,
            "mask_window": mask_window, "mask_fraction": mask_fraction,
            "mean_plddt": _mean_plddt(pdb),
            "tm_fold1": tm1, "tm_fold2": tm2,
            "TMdiff": (tm1 - tm2) if (tm1 == tm1 and tm2 == tm2) else float("nan"),
            "pdb_path": str(pdb),
        })
        (sub_out / "manifest.flag").touch()
    if rows:
        with open(out_root / "manifest.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
    return rows


def summarize_per_cluster(rows):
    by = {}
    for r in rows:
        by.setdefault(r["cluster"], []).append(r)
    out = []
    for tag, items in by.items():
        diff = np.array([x["TMdiff"] for x in items if x["TMdiff"] == x["TMdiff"]])
        if len(diff) == 0:
            continue
        out.append({
            "pair_id": items[0]["pair_id"],
            "cluster": tag,
            "n_perturbations": len(items),
            "mean_tm_fold1": float(np.nanmean([x["tm_fold1"] for x in items])),
            "mean_tm_fold2": float(np.nanmean([x["tm_fold2"] for x in items])),
            "median_TMdiff": float(np.median(diff)),
            "std_TMdiff": float(np.std(diff)),
            "frac_F1_pref": float((diff > 0.05).sum()) / len(diff),
            "frac_F2_pref": float((diff < -0.05).sum()) / len(diff),
        })
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-input", "--input", dest="input", required=True)
    ap.add_argument("--n_perturbations", type=int, default=20)
    ap.add_argument("--mask_window", type=int, default=8)
    ap.add_argument("--mask_fraction", type=float, default=0.5)
    ap.add_argument("--seed", type=int, default=12345)
    args = ap.parse_args()

    pair_id = args.input
    pair_dir = Path(DATA_DIR) / pair_id
    m = re.match(r"^([0-9a-zA-Z]{4})([A-Za-z0-9])_([0-9a-zA-Z]{4})([A-Za-z0-9])$",
                 pair_id)
    if not m:
        print(f"[error] cannot parse pair_id={pair_id}", file=sys.stderr)
        sys.exit(2)
    pdbid1, chain1, pdbid2, chain2 = m.groups()
    cluster_a3ms = sorted((pair_dir / "output_msa_cluster").glob("ShallowMsa_*.a3m"))
    if not cluster_a3ms:
        print(f"[error] no clusters under {pair_dir}", file=sys.stderr)
        sys.exit(2)

    t0 = time.time()
    all_rows = []
    for i, a3m in enumerate(cluster_a3ms):
        rows = run_one_cluster(pair_id, a3m, pdbid1, chain1, pdbid2, chain2,
                               args.n_perturbations, args.mask_window,
                               args.mask_fraction, pair_dir,
                               seed=args.seed + i * 1009)
        all_rows.extend(rows)

    ana = pair_dir / "Analysis"
    ana.mkdir(parents=True, exist_ok=True)
    if all_rows:
        with open(ana / "df_speachaf_raw.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
            w.writeheader()
            w.writerows(all_rows)
    summary = summarize_per_cluster(all_rows)
    if summary:
        with open(ana / "df_speachaf.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
            w.writeheader()
            w.writerows(summary)
    print(f"[speachaf] done in {time.time()-t0:.0f}s; "
          f"summary -> {ana / 'df_speachaf.csv'}")


if __name__ == "__main__":
    main()
