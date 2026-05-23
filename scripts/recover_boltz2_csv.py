#!/usr/bin/env python3
"""Recover df_boltz2.csv (apo) and df_boltz2_holo.csv (holo) summary CSVs
from per-pair output_boltz2/ trees AND/OR from existing df_boltz2_raw.csv
files, without re-running Boltz inference.

Why this exists
---------------
The original run_Boltz2.py writes a raw per-PDB CSV (df_boltz2_*_raw.csv)
and then computes a per-cluster summary (df_boltz2[_holo].csv). The summary
step drops clusters with all-NaN TMdiff, which has been silently producing
empty summary CSVs when TM-align fails for a pair (e.g. due to chain
mismatch or path resolution). This script rebuilds the summary from disk:

  1. If df_boltz2[_holo]_raw.csv exists, parse it and re-emit the
     summary CSV. This is the fast path.
  2. Otherwise, walk output_boltz2/<cluster>/.../*.pdb and (re)compute
     per-PDB TM-scores via tmtools. Slower but fully recovers from
     disk-resident predictions.

Usage
-----
    python3 scripts/recover_boltz2_csv.py             # all pairs
    python3 scripts/recover_boltz2_csv.py --pairs 6c6sD_2ougC
    python3 scripts/recover_boltz2_csv.py --pairs ALL --rescan    # ignore raw, recompute
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from pathlib import Path

import numpy as np

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_THIS_DIR))

from config import DATA_DIR  # noqa: E402

PAIR_RE = re.compile(r"^([0-9a-zA-Z]{4})([A-Z0-9])_([0-9a-zA-Z]{4})([A-Z0-9])$")


def parse_pair(pair_id: str):
    m = PAIR_RE.match(pair_id)
    if not m:
        return None
    return m.groups()


def truth_pdb(pair_dir: Path, pdbid: str, chain: str) -> Path | None:
    for cand in (
        pair_dir / "chain_pdb_files" / f"{pdbid}{chain}.pdb",
        pair_dir / f"{pdbid}.pdb",
    ):
        if cand.exists() and cand.stat().st_size > 0:
            return cand
    return None


def safe_tm(pdb1: Path, pdb2: Path | None, chain2: str | None) -> float:
    if pdb2 is None or not pdb2.exists():
        return float("nan")
    try:
        from utils.align_utils import compute_tmscore_align
        tm = compute_tmscore_align(str(pdb1), str(pdb2), chain2=chain2)
        if isinstance(tm, (list, tuple)) and tm:
            tm = tm[0]
        return float(tm)
    except Exception:
        return float("nan")


def summarize_rows(rows: list[dict]) -> list[dict]:
    """Group per-PDB rows by cluster, keep rank-1 TM-scores per cluster."""
    by_cluster: dict[str, list[dict]] = {}
    for r in rows:
        by_cluster.setdefault(r["cluster"], []).append(r)
    out: list[dict] = []
    for tag, items in by_cluster.items():
        # rank-1 = first item (or smallest rank if present)
        items_sorted = sorted(items, key=lambda x: int(x.get("rank", 0)))
        top = items_sorted[0]
        try:
            tm1 = float(top["tm_fold1"])
            tm2 = float(top["tm_fold2"])
        except (KeyError, TypeError, ValueError):
            continue
        # Compute means across the cluster's PDBs
        tms1 = [float(x["tm_fold1"]) for x in items
                if x.get("tm_fold1") not in (None, "", "nan")]
        tms2 = [float(x["tm_fold2"]) for x in items
                if x.get("tm_fold2") not in (None, "", "nan")]
        mean1 = float(np.nanmean(tms1)) if tms1 else float("nan")
        mean2 = float(np.nanmean(tms2)) if tms2 else float("nan")
        tmdiff = (tm1 - tm2) if (tm1 == tm1 and tm2 == tm2) else float("nan")
        out.append({
            "pair_id": top.get("pair_id", ""),
            "cluster": tag,
            "n_models": len(items),
            "rank1_tm_fold1": tm1,
            "rank1_tm_fold2": tm2,
            "rank1_TMdiff": tmdiff,
            "mean_tm_fold1": mean1,
            "mean_tm_fold2": mean2,
            "partner_yaml": top.get("partner_yaml", ""),
            "mode": top.get("mode", ""),
            "ligand_ccds": top.get("ligand_ccds", ""),
        })
    return out


def write_summary(rows: list[dict], out_csv: Path) -> int:
    """Write summary rows to CSV. Returns number of rows written."""
    if not rows:
        return 0
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    return len(rows)


def recover_from_raw(raw_csv: Path) -> list[dict]:
    """Parse df_boltz2[_holo]_raw.csv and return rows."""
    rows: list[dict] = []
    with open(raw_csv) as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append(r)
    return rows


def rescan_from_disk(pair_id: str, mode_tag: str) -> list[dict]:
    """Scan output_boltz2/<cluster>/.../*.pdb and compute per-PDB rows."""
    parts = parse_pair(pair_id)
    if parts is None:
        return []
    pdbid1, ch1, pdbid2, ch2 = parts
    pair_dir = Path(DATA_DIR) / pair_id
    out_dir = pair_dir / "output_boltz2"
    if not out_dir.exists():
        return []
    truth1 = truth_pdb(pair_dir, pdbid1, ch1)
    truth2 = truth_pdb(pair_dir, pdbid2, ch2)
    rows: list[dict] = []
    for cluster_dir in sorted(out_dir.glob("ShallowMsa_*")):
        # New Boltz writes predictions/input/input_model_*.pdb
        pdbs = sorted(cluster_dir.rglob("*model*.pdb"))
        # Fallback to any PDB in the tree
        if not pdbs:
            pdbs = sorted(cluster_dir.rglob("*.pdb"))
        for i, pdb in enumerate(pdbs):
            tm1 = safe_tm(pdb, truth1, ch1)
            tm2 = safe_tm(pdb, truth2, ch2)
            rows.append({
                "pair_id": pair_id,
                "cluster": cluster_dir.name,
                "rank": i,
                "tm_fold1": tm1,
                "tm_fold2": tm2,
                "TMdiff": (tm1 - tm2) if (tm1 == tm1 and tm2 == tm2) else float("nan"),
                "pdb_path": str(pdb),
                "partner_yaml": "",
                "mode": mode_tag,
                "ligand_ccds": "",
            })
    return rows


def recover_pair(pair_id: str, rescan: bool = False) -> dict:
    """Recover both apo and holo CSVs for one pair.

    Returns dict with counts for logging.
    """
    pair_dir = Path(DATA_DIR) / pair_id
    ana = pair_dir / "Analysis"
    result = {"pair": pair_id, "apo_rows": 0, "holo_rows": 0,
              "apo_path": "", "holo_path": "", "skipped": ""}

    if not (pair_dir / "output_boltz2").exists():
        result["skipped"] = "no_output_boltz2_dir"
        return result

    for mode_tag, raw_fn, sum_fn in (
        ("apo", "df_boltz2_raw.csv", "df_boltz2.csv"),
        ("holo", "df_boltz2_holo_raw.csv", "df_boltz2_holo.csv"),
    ):
        raw_path = ana / raw_fn
        sum_path = ana / sum_fn

        # Decide which input rows to use
        rows: list[dict] = []
        if raw_path.exists() and not rescan:
            try:
                rows = recover_from_raw(raw_path)
            except Exception as e:
                print(f"  [{pair_id}/{mode_tag}] failed to parse {raw_path}: {e}",
                      file=sys.stderr)
                rows = []

        # If raw is empty/missing and we're in apo mode, fall back to disk scan.
        # For holo, only rescan if the user explicitly asks (the trigger CSV
        # ligand mapping isn't reproduced here).
        if not rows and (mode_tag == "apo" or rescan):
            rows = rescan_from_disk(pair_id, mode_tag)

        if not rows:
            continue

        summary = summarize_rows(rows)
        n_written = write_summary(summary, sum_path)
        result[f"{mode_tag}_rows"] = n_written
        result[f"{mode_tag}_path"] = str(sum_path) if n_written else ""
    return result


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs", default="ALL",
                    help="Comma- or space-separated pair IDs, or ALL. Default: ALL.")
    ap.add_argument("--rescan", action="store_true",
                    help="Ignore existing df_boltz2_raw.csv and recompute "
                         "from PDB files on disk (slower; needed when the raw "
                         "CSV is missing or stale).")
    args = ap.parse_args()

    if args.pairs.upper() == "ALL":
        from utils.utils import list_protein_pairs
        pairs = []
        for p in list_protein_pairs():
            if isinstance(p, tuple):
                pairs.append("_".join(p))
            else:
                pairs.append(p)
    else:
        pairs = [s.strip() for s in args.pairs.replace(",", " ").split() if s.strip()]

    n_apo_rec = n_holo_rec = 0
    skipped = []
    for pid in pairs:
        res = recover_pair(pid, rescan=args.rescan)
        if res["skipped"]:
            skipped.append(pid)
            continue
        if res["apo_rows"]:
            n_apo_rec += 1
            print(f"  {pid}: apo {res['apo_rows']} clusters -> {res['apo_path']}")
        if res["holo_rows"]:
            n_holo_rec += 1
            print(f"  {pid}: holo {res['holo_rows']} clusters -> {res['holo_path']}")

    print(f"\nSummary: recovered apo on {n_apo_rec} pairs, holo on {n_holo_rec} pairs.")
    if skipped:
        print(f"  Skipped (no output_boltz2 dir): {len(skipped)} pairs")


if __name__ == "__main__":
    main()
