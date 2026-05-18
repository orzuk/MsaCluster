#!/usr/bin/env python3
"""run_AlphaFlow.py -- per-cluster AlphaFlow ensemble prediction.

AlphaFlow (Jing, Berger, Jaakkola 2024) is a flow-matching diffusion
model trained on AF2 outputs that produces conformational ENSEMBLES per
input. We run it per ShallowMsa cluster: for each cluster we sample N
conformations and compute per-cluster ensemble statistics, then ask
whether different clusters produce different ensemble distributions.

Inputs
  Pipeline/<pair>/output_msa_cluster/ShallowMsa_*.a3m
  (We feed AlphaFlow the cluster consensus sequence; it works
  single-sequence and ensembles via the diffusion sampler.)

Outputs
  Pipeline/<pair>/output_alphaflow/<cluster_tag>/sample_<i>__<pdbid><chain>.pdb
  Pipeline/<pair>/output_alphaflow/<cluster_tag>/manifest.csv
  Pipeline/<pair>/Analysis/df_alphaflow.csv      (per-cluster summary)

Venv: torch2-venv (separate from my-python-venv). The actual AlphaFlow
inference call goes through scripts/shell/RunAlphaFlow.sh which sources
the torch2-venv; this Python wrapper stays in my-python-venv and shells
out to that wrapper.

Usage
  python3 run_AlphaFlow.py -input 6c6sD_2ougC --samples 8
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path
from typing import List

import numpy as np

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _THIS_DIR)

from config import DATA_DIR, MAIN_DIR  # noqa: E402
from utils.protein_utils import compute_tmscore_align  # noqa: E402


SH_WRAPPER = os.path.join(MAIN_DIR, "scripts", "shell", "RunAlphaFlow.sh")


def _consensus_from_a3m(a3m: Path) -> str:
    """Compute the column-wise plurality consensus of an A3M cluster MSA,
    treating only uppercase letters (skip A3M lowercase inserts and gaps).
    Returns one-letter sequence. The query (first seq) is preferred when
    column has no clear plurality."""
    # Read query length first
    query = None
    seqs: List[str] = []
    name: str = ""
    cur: List[str] = []
    with open(a3m) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name and cur:
                    seqs.append("".join(cur))
                name, cur = line[1:], []
            elif line:
                cur.append(line)
        if name and cur:
            seqs.append("".join(cur))
    if not seqs:
        return ""
    # Reduce each sequence to non-insert columns (uppercase + '-')
    reduced = []
    for s in seqs:
        reduced.append("".join(ch for ch in s if (ch.isupper() or ch == "-")))
    L = max(len(s) for s in reduced)
    out_chars = []
    for col in range(L):
        chars = [s[col] for s in reduced if col < len(s) and s[col] != "-"]
        if not chars:
            out_chars.append(reduced[0][col] if col < len(reduced[0]) else "X")
            continue
        most, _ = Counter(chars).most_common(1)[0]
        out_chars.append(most)
    return "".join(out_chars)


def _truth_chain_pdb(pair_dir: Path, pdbid: str, chain: str) -> str:
    for c in (pair_dir / f"{pdbid}.pdb",
              pair_dir / "chain_pdb_files" / f"{pdbid}{chain}.pdb",
              pair_dir / f"{pdbid}_cif.pdb"):
        if c.is_file():
            return str(c)
    return ""


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


def run_one_cluster(pair_id, cluster_a3m, pdbid1, chain1, pdbid2, chain2,
                    samples, mode, pair_dir):
    tag = cluster_a3m.stem
    out_root = pair_dir / "output_alphaflow" / tag
    out_root.mkdir(parents=True, exist_ok=True)
    truth1 = _truth_chain_pdb(pair_dir, pdbid1, chain1)
    truth2 = _truth_chain_pdb(pair_dir, pdbid2, chain2)

    # Build single-line input CSV for AlphaFlow
    consensus = _consensus_from_a3m(cluster_a3m)
    if not consensus:
        print(f"  [alphaflow] {tag}: empty consensus, skip")
        return []
    in_csv = out_root / "input.csv"
    with open(in_csv, "w") as f:
        f.write("name,sequence\n")
        f.write(f"{tag},{consensus}\n")

    # Invoke shell wrapper (handles torch2-venv activation)
    cmd = (f"bash {shlex.quote(SH_WRAPPER)} "
           f"{shlex.quote(str(in_csv))} {shlex.quote(str(out_root))} "
           f"{samples} {shlex.quote(mode)}")
    print(f"  [alphaflow] {tag}: {samples} samples ({cmd})")
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"    [alphaflow] wrapper failed for {tag}: {e}", file=sys.stderr)
        return []

    # Collect PDB outputs (AlphaFlow predict.py writes one PDB per sample)
    pdbs = sorted(out_root.glob("*.pdb"))
    rows = []
    for i, pdb in enumerate(pdbs):
        tm1 = _safe_tm(pdb, truth1)
        tm2 = _safe_tm(pdb, truth2)
        rows.append({
            "pair_id": pair_id, "cluster": tag, "sample_idx": i,
            "mode": mode,
            "tm_fold1": tm1, "tm_fold2": tm2,
            "TMdiff": (tm1 - tm2) if (tm1 == tm1 and tm2 == tm2) else float("nan"),
            "pdb_path": str(pdb),
        })
    if rows:
        with open(out_root / "manifest.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
    return rows


def summarize(rows):
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
            "n_samples": len(items),
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
    ap.add_argument("--samples", type=int, default=8)
    ap.add_argument("--mode", default="esmflow_md_distilled",
                    help="AlphaFlow inference mode; see RunAlphaFlow.sh")
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

    if not os.path.isfile(SH_WRAPPER):
        print(f"[error] missing shell wrapper {SH_WRAPPER}", file=sys.stderr)
        sys.exit(2)

    t0 = time.time()
    all_rows = []
    for a3m in cluster_a3ms:
        rows = run_one_cluster(pair_id, a3m, pdbid1, chain1, pdbid2, chain2,
                               args.samples, args.mode, pair_dir)
        all_rows.extend(rows)

    ana = pair_dir / "Analysis"
    ana.mkdir(parents=True, exist_ok=True)
    if all_rows:
        with open(ana / "df_alphaflow_raw.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
            w.writeheader()
            w.writerows(all_rows)
    summary = summarize(all_rows)
    if summary:
        with open(ana / "df_alphaflow.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
            w.writeheader()
            w.writerows(summary)
    print(f"[alphaflow] done in {time.time()-t0:.0f}s; "
          f"summary -> {ana / 'df_alphaflow.csv'}")


if __name__ == "__main__":
    main()
