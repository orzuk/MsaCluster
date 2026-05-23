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
from typing import List, Optional

import numpy as np

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _THIS_DIR)

from config import DATA_DIR, MAIN_DIR  # noqa: E402
from utils.align_utils import compute_tmscore_align  # noqa: E402
from utils.trigger_utils import (  # noqa: E402
    get_trigger_class,
    ligand_codes_for_holo,
)


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
                      out_yaml: Path,
                      msa_path: Optional[Path] = None,
                      ligand_ccds: Optional[List[str]] = None) -> None:
    """Compose the Boltz YAML input for a single cluster.

    Parameters
    ----------
    msa_path
        Optional path to the cluster's ShallowMsa A3M file. When given,
        it is referenced from the YAML as the chain's MSA so Boltz does
        NOT fall back to fetching a generic colabfold MSA. This is the
        cluster-specific evolutionary signal that makes our per-cluster
        prediction meaningful; passing it preserves the same MSA the
        AF2/AF3 per-cluster pipeline already uses.
    ligand_ccds
        Optional list of PDB chemical-component IDs (e.g. ['ATP', 'MG'])
        to append as ``ligand:`` blocks for a holo prediction. Codes are
        taken from ``docs/triggers_from_pdb.csv`` (columns
        ``ligands_pdb1`` / ``ligands_pdb2``).
    partner_yaml
        Path to an existing YAML fragment to append (protein partner,
        SMILES ligand, RNA partner). Used for ad-hoc holo experiments
        outside the trigger CSV.
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
        if msa_path is not None:
            # Boltz reads custom MSA when this is set (overrides server fetch)
            f.write(f"      msa: {msa_path}\n")
        # Append holo ligand blocks (PDB CCD codes) from trigger CSV
        if ligand_ccds:
            next_id = ord("B")
            for code in ligand_ccds:
                f.write(f"  - ligand:\n")
                f.write(f"      id: {chr(next_id)}\n")
                f.write(f"      ccd: {code}\n")
                next_id += 1
        if partner_block:
            f.write(partner_block)


def run_one_cluster(pair_id, cluster_a3m, pdbid1, chain1, pdbid2, chain2,
                    pair_dir, partner_yaml,
                    ligand_ccds: Optional[List[str]] = None):
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
    # Pass the cluster's own MSA so Boltz does not re-fetch a generic
    # one via --use_msa_server. This is the whole point of per-cluster
    # prediction (same as AF2/AF3 per-cluster).
    _write_boltz_yaml(tag, consensus, partner_yaml, in_yaml,
                      msa_path=cluster_a3m, ligand_ccds=ligand_ccds)

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
    ap.add_argument(
        "--mode", choices=["apo", "holo", "auto"], default="apo",
        help="apo: predict each cluster with no ligand/partner (default). "
             "holo: predict each cluster WITH the trigger's ligand block "
             "(only implemented for trigger_class='ligand' so far; other "
             "classes raise an error). auto: behave as holo when the "
             "pair is ligand-triggered per docs/triggers_from_pdb.csv, "
             "else apo. The result CSV's `mode` column records what ran.")
    ap.add_argument("--partner_yaml", default="",
                    help="Optional Boltz-2 partner block (extra protein/"
                         "ligand/RNA). If provided, appended to every "
                         "cluster's input YAML; useful for ad-hoc holo "
                         "experiments that aren't captured by the trigger "
                         "CSV (e.g. specific RfaH-RNAP runs).")
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

    # Resolve mode + ligand list from the trigger CSV
    trigger_class = get_trigger_class(pair_id)
    ligand_ccds: List[str] = []
    effective_mode = args.mode
    if effective_mode == "auto":
        effective_mode = "holo" if trigger_class == "ligand" else "apo"
    if effective_mode == "holo":
        if trigger_class != "ligand":
            # Other holo subclasses (protein_binding / oligomerization /
            # mutation) are not yet implemented because they need partner
            # sequence extraction or homo-multimer scaffolding.
            print(f"[boltz2] holo mode requested for {pair_id} but "
                  f"trigger_class={trigger_class!r} — only 'ligand' is "
                  f"supported. Falling back to apo. Edit run_Boltz2.py "
                  f"to add protein_binding/oligomerization holo paths.",
                  file=sys.stderr)
            effective_mode = "apo"
        else:
            ligand_ccds = ligand_codes_for_holo(pair_id)
            if not ligand_ccds:
                print(f"[boltz2] {pair_id}: trigger_class=ligand but no "
                      f"ligand codes in CSV; falling back to apo.",
                      file=sys.stderr)
                effective_mode = "apo"

    if effective_mode == "holo":
        print(f"[boltz2] {pair_id}: HOLO mode, ligand CCDs = {ligand_ccds}")
    else:
        print(f"[boltz2] {pair_id}: APO mode "
              f"(trigger_class={trigger_class}, requested={args.mode})")

    t0 = time.time()
    all_rows = []
    for a3m in cluster_a3ms:
        rows = run_one_cluster(pair_id, a3m, pdbid1, chain1, pdbid2, chain2,
                               pair_dir, args.partner_yaml,
                               ligand_ccds=ligand_ccds if effective_mode == "holo" else None)
        # Record which mode produced these rows
        for r in rows:
            r["mode"] = effective_mode
            r["ligand_ccds"] = ",".join(ligand_ccds) if ligand_ccds else ""
        all_rows.extend(rows)

    ana = pair_dir / "Analysis"
    ana.mkdir(parents=True, exist_ok=True)
    # Distinguish apo / holo outputs so both can coexist for the same pair
    raw_fn = "df_boltz2_raw.csv" if effective_mode == "apo" else "df_boltz2_holo_raw.csv"
    sum_fn = "df_boltz2.csv" if effective_mode == "apo" else "df_boltz2_holo.csv"
    if all_rows:
        with open(ana / raw_fn, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
            w.writeheader()
            w.writerows(all_rows)
    summary = summarize(all_rows)
    if summary:
        # add mode/ligand cols to summary too
        for s in summary:
            s["mode"] = effective_mode
            s["ligand_ccds"] = ",".join(ligand_ccds) if ligand_ccds else ""
        with open(ana / sum_fn, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
            w.writeheader()
            w.writerows(summary)
    print(f"[boltz2] done in {time.time()-t0:.0f}s ({effective_mode}); "
          f"summary -> {ana / sum_fn}")


if __name__ == "__main__":
    main()
