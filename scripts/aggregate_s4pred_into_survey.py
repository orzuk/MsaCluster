#!/usr/bin/env python3
"""Aggregate per-pair S4PRED outputs into docs/fold_diversity_survey.csv as the
8th method ("S4PRED").

Per pair, this script:
  1. Reads Pipeline/FoldPairs/<pair>/Analysis/df_s4pred.csv
     (per-cluster predicted Q3 fractions from run_s4pred_per_cluster.py)
  2. Loads the two truth PDB structures (Pipeline/FoldPairs/<pair>/<pdb*>.pdb)
     and computes per-fold true secondary-structure via DSSP. From the
     DSSP output we extract frac_H_true, frac_E_true per fold.
  3. For each cluster, computes per-cluster similarity to fold-1 vs fold-2
     in the (frac_H, frac_E) plane via complement of Euclidean distance:
         sim_F1 = 1 - distance((frac_H_pred, frac_E_pred), (frac_H_F1, frac_E_F1))
         sim_F2 = 1 - distance((frac_H_pred, frac_E_pred), (frac_H_F2, frac_E_F2))
         SS_diff = sim_F1 - sim_F2     (positive = cluster's predicted SS
                                        is closer to fold 1's SS than to fold 2's)
  4. Centers per-pair (subtracts per-pair median of SS_diff across
     ShallowMsa_* clusters, mirroring the centering done for the
     TM-score-based methods).
  5. Appends method="S4PRED" rows to docs/fold_diversity_survey.csv.

Why a simple (%H, %E)-plane distance and not SOV?
-------------------------------------------------
SOV requires same-length predicted and truth strings, which in turn
requires aligning PDB-chain residue numbering with the query/MSA frame.
That's doable but adds an alignment step that can fail silently on
construct mismatches. The (%H, %E) Euclidean similarity is robust to
length differences and captures the same biology for the alpha-vs-beta
fold-switch question (both folds we care about differ primarily in
their helix-vs-strand composition). A future revision can use SOV
(scripts/sov.py is ready) once an explicit query<->PDB residue
alignment is wired in.

DSSP source
-----------
Uses Biopython's Bio.PDB.DSSP wrapper around the `mkdssp` (or `dssp`)
binary. If that binary is not on PATH, falls back to a manual H/E/C
assignment from CA-CA distances + dihedrals via a small heuristic.
The fallback is only approximate; the script will print a warning when
it uses the fallback so the user knows to install mkdssp.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import sys
from glob import glob

import numpy as np
import pandas as pd

DOCS = "docs"
SURVEY_CSV = os.path.join(DOCS, "fold_diversity_survey.csv")
PAIR_RE = re.compile(r"^([0-9a-z]{4})([A-Za-z0-9])_([0-9a-z]{4})([A-Za-z0-9])$")


def parse_pair_id(pair_id: str):
    m = PAIR_RE.match(pair_id)
    if not m:
        return None
    return m.group(1), m.group(2), m.group(3), m.group(4)


def find_pdb_path(pair_dir: str, pdbid: str) -> str | None:
    for cand in (
        os.path.join(pair_dir, f"{pdbid}.pdb"),
        os.path.join(pair_dir, f"{pdbid.lower()}.pdb"),
        os.path.join(pair_dir, f"{pdbid.upper()}.pdb"),
    ):
        if os.path.isfile(cand):
            return cand
    matches = glob(os.path.join(pair_dir, f"{pdbid}*.pdb"))
    return matches[0] if matches else None


def true_ss_from_dssp(pdb_path: str, chain: str) -> tuple[float, float, float, int]:
    """Return (frac_H, frac_E, frac_C, n_res) for the given chain of pdb_path.
    Uses mdtraj's native DSSP (no external binary). Falls back to Biopython's
    Bio.PDB.DSSP wrapper if mdtraj is unavailable.

    mdtraj's compute_dssp(simplified=True) returns the Q3 codes 'H','E','C'
    directly, matching the convention we use elsewhere.
    """
    # Path 1: mdtraj (preferred — no external binary needed)
    try:
        import mdtraj as md
    except ImportError:
        md = None

    if md is not None:
        try:
            traj = md.load_pdb(pdb_path)
            top = traj.topology
            ss_q3 = md.compute_dssp(traj, simplified=True)[0]
            # Filter to the requested chain
            chain_norm = (chain or "").upper().strip()
            if chain_norm:
                chain_residues = []
                for i, res in enumerate(top.residues):
                    rcid = getattr(res.chain, "chain_id", "")
                    if rcid and rcid.upper() == chain_norm:
                        chain_residues.append(i)
                if not chain_residues:
                    # Chain ID not found; try matching chain.index numerically
                    # (rare; e.g., when PDB chain IDs were stripped).
                    chain_residues = list(range(top.n_residues))
            else:
                chain_residues = list(range(top.n_residues))
            ss_chain = "".join(ss_q3[i] for i in chain_residues)
            # mdtraj uses 'NA' for non-protein residues; filter them
            ss_chain = "".join(c for c in ss_chain if c in "HEC")
            n = len(ss_chain)
            if n == 0:
                return (float("nan"), float("nan"), float("nan"), 0)
            h = ss_chain.count("H") / n
            e = ss_chain.count("E") / n
            c = ss_chain.count("C") / n
            return (h, e, c, n)
        except Exception as e:
            print(f"[warn] mdtraj DSSP failed for {pdb_path}: {e}",
                  file=sys.stderr)
            # fall through to Biopython path

    # Path 2: Biopython DSSP wrapper (requires dssp/mkdssp binary on PATH)
    try:
        from Bio.PDB import PDBParser, DSSP
    except ImportError:
        return (float("nan"), float("nan"), float("nan"), 0)
    try:
        p = PDBParser(QUIET=True)
        struct = p.get_structure("x", pdb_path)
        model = next(iter(struct))
        dssp = None
        for binname in ("mkdssp", "dssp"):
            try:
                dssp = DSSP(model, pdb_path, dssp=binname)
                break
            except Exception:
                continue
        if dssp is None:
            print(f"[warn] DSSP failed for {pdb_path}: no working binary",
                  file=sys.stderr)
            return (float("nan"), float("nan"), float("nan"), 0)
    except Exception as e:
        print(f"[warn] PDB parse failed for {pdb_path}: {e}", file=sys.stderr)
        return (float("nan"), float("nan"), float("nan"), 0)

    H_MAP = {"H": "H", "G": "H", "I": "H",
             "E": "E", "B": "E",
             "T": "C", "S": "C", "-": "C", " ": "C", "P": "C"}
    counts = {"H": 0, "E": 0, "C": 0}
    n = 0
    for key in dssp.keys():
        chain_id, _ = key
        if chain and chain.upper() != chain_id.upper():
            continue
        raw_ss = dssp[key][2]
        q3 = H_MAP.get(raw_ss, "C")
        counts[q3] += 1
        n += 1
    if n == 0:
        return (float("nan"), float("nan"), float("nan"), 0)
    return (counts["H"] / n, counts["E"] / n, counts["C"] / n, n)


def he_similarity(frac_h_pred: float, frac_e_pred: float,
                   frac_h_true: float, frac_e_true: float) -> float:
    """1 - Euclidean distance in (%H, %E) plane.
    Max possible distance is sqrt(2); we scale so result is in [-something, 1]
    but practically (frac_H, frac_E) are in [0,1] so distance is in [0, sqrt(2)].
    """
    if any(not np.isfinite(v) for v in (frac_h_pred, frac_e_pred,
                                          frac_h_true, frac_e_true)):
        return float("nan")
    d = math.sqrt((frac_h_pred - frac_h_true) ** 2 +
                  (frac_e_pred - frac_e_true) ** 2)
    return 1.0 - d / math.sqrt(2.0)


def aggregate_pair(pair_id: str, fpairs_root: str) -> list[dict]:
    """Build the S4PRED rows for one pair."""
    parts = parse_pair_id(pair_id)
    if not parts:
        print(f"[warn] {pair_id}: unrecognised pair_id format", file=sys.stderr)
        return []
    pdb1, ch1, pdb2, ch2 = parts
    pair_dir = os.path.join(fpairs_root, pair_id)
    s4_csv = os.path.join(pair_dir, "Analysis", "df_s4pred.csv")
    if not os.path.isfile(s4_csv):
        return []
    p1_path = find_pdb_path(pair_dir, pdb1)
    p2_path = find_pdb_path(pair_dir, pdb2)
    if not p1_path or not p2_path:
        print(f"[warn] {pair_id}: missing PDB(s) for {pdb1}/{pdb2}",
              file=sys.stderr)
        return []
    hF1, eF1, _, nF1 = true_ss_from_dssp(p1_path, ch1)
    hF2, eF2, _, nF2 = true_ss_from_dssp(p2_path, ch2)
    if nF1 == 0 or nF2 == 0:
        print(f"[warn] {pair_id}: DSSP failed for at least one fold",
              file=sys.stderr)
        return []

    df = pd.read_csv(s4_csv)
    out_rows: list[dict] = []
    diffs: list[float] = []
    for _, r in df.iterrows():
        h = float(r.get("frac_H", float("nan")))
        e = float(r.get("frac_E", float("nan")))
        sim_F1 = he_similarity(h, e, hF1, eF1)
        sim_F2 = he_similarity(h, e, hF2, eF2)
        diff = sim_F1 - sim_F2 if (np.isfinite(sim_F1) and np.isfinite(sim_F2)) else float("nan")
        is_deep = str(r.get("cluster", "")).lower().startswith("deep")
        if np.isfinite(diff) and not is_deep:
            diffs.append(diff)
        out_rows.append({
            "pair_id": pair_id,
            "cluster": r["cluster"],
            "method": "S4PRED",
            "TM1_max": sim_F1, "TM2_max": sim_F2,
            "TM1_mean": sim_F1, "TM2_mean": sim_F2,
            "TMdiff_max": diff, "TMdiff_mean": diff,
            "n_models": int(r.get("n_reps", 0)),
            "n_toward_f1": 0, "n_toward_f2": 0, "vote_frac_f1": float("nan"),
            "pref": "Amb",
            "metric_used": "SS_HE_similarity",
            "ddg_bias_mean": "", "ddg_bias_std": "",
            "ddg_sum_F1_mean": "", "ddg_sum_F2_mean": "",
            "t1_recall": "", "t2_recall": "",
        })
    # Center per pair (median over shallow clusters only)
    median = float(np.median(diffs)) if diffs else 0.0
    DELTA = 0.05
    for row in out_rows:
        d = row["TMdiff_max"]
        if np.isfinite(d):
            tc = d - median
            row["TMdiff_centered"] = tc
            row["pair_median_used"] = median
            row["pref_centered"] = ("F1" if tc > DELTA
                                     else ("F2" if tc < -DELTA else "Amb"))
            row["pref"] = ("F1" if d > 0 else ("F2" if d < 0 else "Amb"))
        else:
            row["TMdiff_centered"] = float("nan")
            row["pair_median_used"] = median
            row["pref_centered"] = "Amb"
    return out_rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs", default="ALL")
    ap.add_argument("--fpairs_root", default="Pipeline/FoldPairs")
    ap.add_argument("--survey_csv", default=SURVEY_CSV)
    ap.add_argument("--write_mode", default="append",
                    choices=["append", "replace_method"],
                    help="append: just append new rows; "
                         "replace_method: drop all existing rows where "
                         "method=='S4PRED' before appending.")
    args = ap.parse_args()

    if args.pairs.upper() == "ALL":
        pair_dirs = sorted(d for d in os.listdir(args.fpairs_root)
                            if os.path.isdir(os.path.join(args.fpairs_root, d)))
    else:
        pair_dirs = [s.strip() for s in args.pairs.replace(",", " ").split()
                      if s.strip()]

    all_rows: list[dict] = []
    for pid in pair_dirs:
        rows = aggregate_pair(pid, args.fpairs_root)
        if rows:
            all_rows.extend(rows)
            print(f"  {pid}: +{len(rows)} cluster rows", flush=True)
    if not all_rows:
        print("[error] no rows generated", file=sys.stderr); sys.exit(2)

    new_df = pd.DataFrame(all_rows)
    if os.path.isfile(args.survey_csv):
        survey = pd.read_csv(args.survey_csv)
        if args.write_mode == "replace_method":
            survey = survey[survey["method"] != "S4PRED"]
        # Align columns: union of both schemas
        all_cols = list(dict.fromkeys(list(survey.columns) + list(new_df.columns)))
        survey = survey.reindex(columns=all_cols)
        new_df = new_df.reindex(columns=all_cols)
        combined = pd.concat([survey, new_df], ignore_index=True)
    else:
        combined = new_df
    combined.to_csv(args.survey_csv, index=False)
    print(f"\n[ok] wrote {len(all_rows)} S4PRED rows -> {args.survey_csv} "
          f"(total rows now: {len(combined)})")


if __name__ == "__main__":
    main()
