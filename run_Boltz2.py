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
from utils.trigger_utils import condition_spec  # noqa: E402


SH_WRAPPER = os.path.join(MAIN_DIR, "scripts", "shell", "RunBoltz2.sh")


def _representative_from_a3m(a3m: Path, method: str = "medoid") -> str:
    """Extract one representative sequence from the cluster a3m.

    Delegates to utils.cluster_representative.get_cluster_representative so
    AF2/AF3/Boltz-2/ESMFold2 all use identical logic for picking the
    cluster's representative. Default is 'medoid' (the cluster's central
    member). Was 'consensus' historically; medoid is more biologically
    meaningful and matches the AF2/AF3 v2 methodology.
    """
    from utils.cluster_representative import get_cluster_representative
    _hdr, seq = get_cluster_representative(a3m, method=method)
    return seq or ""


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
                      ligand_ccds: Optional[List[str]] = None,
                      n_copies: int = 1) -> None:
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
    n_copies
        Number of identical protein chains to emit (homo-oligomer holo for
        the 'oligomerization' trigger class). ``1`` = monomer (default);
        ``n>1`` emits chains A, B, ... all sharing ``sequence`` and the
        same per-cluster ``msa_path``. The chain IDs consumed by copies
        are skipped when numbering subsequent ligand blocks.
    partner_yaml
        Path to an existing YAML fragment to append (protein partner,
        SMILES ligand, RNA partner). Used for the 'protein_binding' class
        (whose partner sequence is not in the trigger CSV) and other
        ad-hoc holo experiments.
    """
    if partner_yaml and os.path.isfile(partner_yaml):
        partner_block = open(partner_yaml).read().rstrip() + "\n"
    else:
        partner_block = ""
    n_copies = max(1, int(n_copies))
    with open(out_yaml, "w") as f:
        f.write("version: 1\n")
        f.write("sequences:\n")
        # Emit n_copies identical protein chains (A, B, ... for homo-oligomer).
        next_id = ord("A")
        for _ in range(n_copies):
            f.write(f"  - protein:\n")
            f.write(f"      id: {chr(next_id)}\n")
            f.write(f"      sequence: {sequence}\n")
            if msa_path is not None:
                # Boltz reads custom MSA when this is set (overrides server fetch)
                f.write(f"      msa: {msa_path}\n")
            next_id += 1
        # Append holo ligand blocks (PDB CCD codes) from trigger CSV, using
        # the next free chain IDs after the protein copies.
        if ligand_ccds:
            for code in ligand_ccds:
                f.write(f"  - ligand:\n")
                f.write(f"      id: {chr(next_id)}\n")
                f.write(f"      ccd: {code}\n")
                next_id += 1
        if partner_block:
            f.write(partner_block)


def run_one_cluster(pair_id, cluster_a3m, pdbid1, chain1, pdbid2, chain2,
                    pair_dir, partner_yaml,
                    ligand_ccds: Optional[List[str]] = None,
                    n_copies: int = 1,
                    force_rerun: bool = False):
    tag = cluster_a3m.stem
    out_root = pair_dir / "output_boltz2" / tag
    out_root.mkdir(parents=True, exist_ok=True)
    truth1 = _truth_chain_pdb(pair_dir, pdbid1, chain1)
    truth2 = _truth_chain_pdb(pair_dir, pdbid2, chain2)

    # Skip if we already have a valid rank_001 PDB for this cluster.
    # An empty manifest.csv or zero-PDB output is treated as a failed
    # previous attempt and re-run.
    existing_rank1 = sorted(out_root.rglob("*rank_001*.pdb"))
    if not existing_rank1:
        existing_rank1 = sorted(out_root.rglob("*.pdb"))
    if existing_rank1 and not force_rerun:
        # Check the PDB has nonzero content (filters out zero-byte files
        # from sbatch failures)
        if any(p.stat().st_size > 100 for p in existing_rank1):
            print(f"  [boltz2] {tag}: skipping ({len(existing_rank1)} PDB(s) on disk)")
            # Still emit rows from the on-disk PDBs so the summary CSV
            # picks them up.
            rows = []
            for i, pdb in enumerate(existing_rank1):
                tm1 = _safe_tm(pdb, truth1)
                tm2 = _safe_tm(pdb, truth2)
                rows.append({
                    "pair_id": pair_id, "cluster": tag, "rank": i,
                    "tm_fold1": tm1, "tm_fold2": tm2,
                    "TMdiff": (tm1 - tm2) if (tm1 == tm1 and tm2 == tm2) else float("nan"),
                    "pdb_path": str(pdb),
                    "partner_yaml": partner_yaml or "",
                })
            return rows

    # v2 methodology: use medoid (the cluster's central member) as the
    # query, matching what AF2/AF3 use. Was 'consensus' previously.
    seq = _representative_from_a3m(cluster_a3m, method="medoid")
    if not seq:
        print(f"  [boltz2] {tag}: empty representative, skip")
        return []
    # Strip any dash characters before passing to Boltz (defensive).
    # Dashes occur when an MSA column is all-gap; Boltz parses them as
    # invalid CCD components and aborts ("CCD component '-' error" seen in
    # 49/93 pair logs in the original consensus mode). Replace with 'X'.
    consensus = seq.replace("-", "X").replace(".", "X")
    in_yaml = out_root / "input.yaml"
    # Pass the cluster's own MSA so Boltz does not re-fetch a generic
    # one via --use_msa_server. This is the whole point of per-cluster
    # prediction (same as AF2/AF3 per-cluster).
    _write_boltz_yaml(tag, consensus, partner_yaml, in_yaml,
                      msa_path=cluster_a3m, ligand_ccds=ligand_ccds,
                      n_copies=n_copies)

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
             "holo: predict each cluster WITH the trigger encoded, per "
             "docs/triggers_from_pdb.csv -- ligand CCD block (ligand class), "
             "homo-oligomer copies (oligomerization class), or an explicit "
             "partner block (protein_binding class, via --partner_yaml). "
             "auto: holo when the pair is in a representable triggered class, "
             "else apo. The result CSV's `mode` column records what ran.")
    ap.add_argument("--partner_yaml", default="",
                    help="Optional Boltz-2 partner block (extra protein/"
                         "ligand/RNA). Required for protein_binding holo "
                         "(the partner sequence is not in the trigger CSV); "
                         "also useful for ad-hoc holo experiments (e.g. "
                         "specific RfaH-RNAP runs).")
    ap.add_argument("--copies", type=int, default=0,
                    help="Homo-oligomer chain copies for oligomerization "
                         "holo. 0 (default) = derive from the trigger CSV "
                         "chain counts, falling back to 2 (dimer) when the "
                         "CSV does not record the stoichiometry. Ignored in "
                         "apo mode.")
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

    # Resolve mode + condition (ligand / oligomer copies / partner) from the
    # trigger CSV via the shared condition_spec (utils.trigger_utils).
    spec = condition_spec(pair_id)
    trigger_class = spec["trigger_class"]
    ligand_ccds: List[str] = []
    n_copies = 1
    effective_mode = args.mode
    if effective_mode == "auto":
        effective_mode = "holo" if spec["representable"] else "apo"

    if effective_mode == "holo":
        if not spec["representable"]:
            print(f"[boltz2] holo requested for {pair_id} but "
                  f"trigger_class={trigger_class!r} is not representable as a "
                  f"holo condition ({spec['note']}). Falling back to apo.",
                  file=sys.stderr)
            effective_mode = "apo"
        elif trigger_class == "ligand":
            ligand_ccds = spec["ligand_ccds"]
            if not ligand_ccds:
                print(f"[boltz2] {pair_id}: ligand class but no CCD codes; "
                      f"falling back to apo.", file=sys.stderr)
                effective_mode = "apo"
        elif trigger_class == "oligomerization":
            # Homo-oligomer: repeat the chain. The CSV rarely records the true
            # stoichiometry (chain counts are often 1/1), so honor --copies,
            # else the CSV-derived count, else default to a dimer.
            n_copies = args.copies or spec["n_copies"] or 2
            if not args.copies and spec["n_copies"] is None:
                print(f"[boltz2] {pair_id}: oligomer stoichiometry unknown "
                      f"from CSV; defaulting to a dimer (--copies to override).",
                      file=sys.stderr)
        elif trigger_class == "protein_binding":
            # The partner SEQUENCE is not in the trigger CSV; it must be
            # supplied via --partner_yaml. Without it we cannot build the holo.
            if not (args.partner_yaml and os.path.isfile(args.partner_yaml)):
                print(f"[boltz2] {pair_id}: protein_binding holo needs a "
                      f"--partner_yaml (partner sequence not in CSV). "
                      f"Falling back to apo.", file=sys.stderr)
                effective_mode = "apo"

    if effective_mode == "holo":
        print(f"[boltz2] {pair_id}: HOLO mode (class={trigger_class}, "
              f"ligand CCDs={ligand_ccds}, copies={n_copies}, "
              f"partner={'yes' if args.partner_yaml else 'no'})")
    else:
        print(f"[boltz2] {pair_id}: APO mode "
              f"(trigger_class={trigger_class}, requested={args.mode})")

    t0 = time.time()
    all_rows = []
    for a3m in cluster_a3ms:
        rows = run_one_cluster(pair_id, a3m, pdbid1, chain1, pdbid2, chain2,
                               pair_dir, args.partner_yaml,
                               ligand_ccds=ligand_ccds if effective_mode == "holo" else None,
                               n_copies=n_copies if effective_mode == "holo" else 1)
        # Record which mode produced these rows
        for r in rows:
            r["mode"] = effective_mode
            r["ligand_ccds"] = ",".join(ligand_ccds) if ligand_ccds else ""
            r["n_copies"] = n_copies if effective_mode == "holo" else 1
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
            s["n_copies"] = n_copies if effective_mode == "holo" else 1
        with open(ana / sum_fn, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
            w.writeheader()
            w.writerows(summary)
    print(f"[boltz2] done in {time.time()-t0:.0f}s ({effective_mode}); "
          f"summary -> {ana / sum_fn}")


if __name__ == "__main__":
    main()
