#!/usr/bin/env python3
"""run_Boltz2.py -- per-cluster Boltz-2 fold prediction.

Boltz-2 (Wohlwend et al. 2024-2025, MIT) is an open-source AlphaFold3-class
all-atom predictor with binding-affinity head. We invoke it per
ShallowMsa cluster on the cluster consensus sequence. Two modes:

  --partner none      : predict the protein alone, per cluster. Cf. AF3.
  --partner <yaml>    : include a binding partner (another protein chain,
                         a ligand SMILES, an RNA sequence) -- specifically
                         useful for binding-/context-triggered fold-
                         switchers like RfaH + RNAP or KaiB + KaiC.

Inputs
  Pipeline/<pair>/output_msa_cluster/ShallowMsa_*.a3m

Outputs
  Pipeline/<pair>/output_boltz2/<cluster_tag>/
      Boltz-2 prediction dir per cluster (PDBs, plddt, optionally affinity)
  Pipeline/<pair>/Analysis/df_boltz2.csv
      per-cluster summary

Venv: torch2-venv. Inference goes through scripts/shell/RunBoltz2.sh.

Usage
  python3 run_Boltz2.py -input 6c6sD_2ougC
  python3 run_Boltz2.py -input 6c6sD_2ougC --partner_yaml data/rfah_rnap.yaml
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shlex
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
from utils.align_utils import compute_tmscore_align  # noqa: E402


SH_WRAPPER = os.path.join(MAIN_DIR, "scripts", "shell", "RunBoltz2.sh")


def _consensus_from_a3m(a3m: Path) -> str:
    seqs: List[str] = []
    name, cur = "", []
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
    reduced = ["".join(ch for ch in s if (ch.isupper() or ch == "-"))
               for s in seqs]
    L = max(len(s) for s in reduced)
    out = []
    for col in range(L):
        chars = [s[col] for s in reduced if col < len(s) and s[col] != "-"]
        if not chars:
            out.append(reduced[0][col] if col < len(reduced[0]) else "X")
            continue
        out.append(Counter(chars).most_common(1)[0][0])
    return "".join(out)


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


def _write_boltz_yaml(tag: str, sequence: str, partner_yaml: str,
                      out_yaml: Path) -> None:
    """Compose the Boltz YAML input for a single cluster.

    If partner_yaml is provided, prepend/merge its content (caller's
    responsibility to format the partner block); we just write the
    protein chain block.
    """
    if partner_yaml and os.path.isfile(partner_yaml):
        partner_block = open(partner_yaml).read().rstrip() + "\n"
    else:
        partner_block = ""
    with open(out_yaml, "w") as f:
        f.write("version: 1\n")
        f.write("sequences:\n")
        f.write(f"  - protein:\n")
        f.write(f"      id: A\n")
        f.write(f"      sequence: {sequence}\n")
        if partner_block:
            f.write(partner_block)


def run_one_cluster(pair_id, cluster_a3m, pdbid1, chain1, pdbid2, chain2,
                    pair_dir, partner_yaml):
    tag = cluster_a3m.stem
    out_root = pair_dir / "output_boltz2" / tag
    out_root.mkdir(parents=True, exist_ok=True)
    truth1 = _truth_chain_pdb(pair_dir, pdbid1, chain1)
    truth2 = _truth_chain_pdb(pair_dir, pdbid2, chain2)

    consensus = _consensus_from_a3m(cluster_a3m)
    if not consensus:
        print(f"  [boltz2] {tag}: empty consensus, skip")
        return []
    in_yaml = out_root / "input.yaml"
    _write_boltz_yaml(tag, consensus, partner_yaml, in_yaml)

    cmd = (f"bash {shlex.quote(SH_WRAPPER)} "
           f"{shlex.quote(str(in_yaml))} {shlex.quote(str(out_root))}")
    print(f"  [boltz2] {tag}: {cmd}")
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"    [boltz2] failed for {tag}: {e}", file=sys.stderr)
        return []

    # Boltz writes outputs into a subdir; pick the rank_001 PDB
    pdbs = sorted(out_root.rglob("*rank_001*.pdb"))
    if not pdbs:
        pdbs = sorted(out_root.rglob("*.pdb"))
    rows = []
    for i, pdb in enumerate(pdbs):
        tm1 = _safe_tm(pdb, truth1)
        tm2 = _safe_tm(pdb, truth2)
        rows.append({
            "pair_id": pair_id, "cluster": tag, "rank": i,
            "tm_fold1": tm1, "tm_fold2": tm2,
            "TMdiff": (tm1 - tm2) if (tm1 == tm1 and tm2 == tm2) else float("nan"),
            "pdb_path": str(pdb),
            "partner_yaml": partner_yaml or "",
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
        # Boltz-2 produces ranked models per cluster; we take rank-1 (or
        # full ranking if available) as the cluster prediction.
        top = items[0]
        out.append({
            "pair_id": top["pair_id"], "cluster": tag,
            "n_models": len(items),
            "rank1_tm_fold1": float(top["tm_fold1"]),
            "rank1_tm_fold2": float(top["tm_fold2"]),
            "rank1_TMdiff": float(top["TMdiff"]),
            "mean_tm_fold1": float(np.nanmean([x["tm_fold1"] for x in items])),
            "mean_tm_fold2": float(np.nanmean([x["tm_fold2"] for x in items])),
            "partner_yaml": top.get("partner_yaml", ""),
        })
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-input", "--input", dest="input", required=True)
    ap.add_argument("--partner_yaml", default="",
                    help="Optional Boltz-2 partner block (extra protein/ligand). "
                         "If provided, appended to per-cluster input YAML.")
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
                               pair_dir, args.partner_yaml)
        all_rows.extend(rows)

    ana = pair_dir / "Analysis"
    ana.mkdir(parents=True, exist_ok=True)
    if all_rows:
        with open(ana / "df_boltz2_raw.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
            w.writeheader()
            w.writerows(all_rows)
    summary = summarize(all_rows)
    if summary:
        with open(ana / "df_boltz2.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
            w.writeheader()
            w.writerows(summary)
    print(f"[boltz2] done in {time.time()-t0:.0f}s; "
          f"summary -> {ana / 'df_boltz2.csv'}")


if __name__ == "__main__":
    main()
