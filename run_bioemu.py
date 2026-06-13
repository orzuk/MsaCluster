#!/usr/bin/env python3
"""run_bioemu.py -- per-cluster BioEmu equilibrium-ensemble fold preference.

BioEmu (Lewis et al., Science 2025; MIT, github.com/microsoft/bioemu) is a
generative diffusion model that samples the *approximate equilibrium ensemble*
of a single protein chain. Unlike AF2/AF3/Boltz-2 (which give a point estimate
-- the single most likely structure per cluster), BioEmu samples MANY
conformations, so for a fold-switching family we can read off the *population
balance* between the two folds, not just a hard call. That is the novelty:
a cluster AF calls cleanly fold-1 may still be, say, 70/30 fold-1/fold-2 in the
BioEmu ensemble -- a direct, physical measure of bistability.

This is the 9th per-cluster method, wired exactly like the others (Boltz-2 in
particular -- see run_Boltz2.py): per ShallowMsa cluster we run BioEmu on the
SAME per-cluster MSA AF2/AF3 use, sample N structures, TM-score each sample
against both truth folds (via utils.align_utils -- never hand-rolled), and
aggregate to a signed per-cluster ensemble fold preference that drops straight
into docs/fold_diversity_survey.csv as method "BioEmu".

Inputs
  Pipeline/<pair>/output_msa_cluster/ShallowMsa_*.a3m   (same as AF/Boltz-2)
Outputs
  Pipeline/<pair>/output_bioemu/<cluster_tag>/          (samples, per-frame PDBs)
  Pipeline/<pair>/Analysis/df_bioemu_raw.csv            (one row per sample)
  Pipeline/<pair>/Analysis/df_bioemu.csv                (one row per cluster)

Venv: a DEDICATED bioemu-venv (see INSTALL_bioemu.md). Inference is shelled out
through scripts/shell/RunBioEmu.sh so the heavy GPU deps stay isolated.

Usage
  python3 run_bioemu.py -input 2n54B_2hdmA
  python3 run_bioemu.py -input 2n54B_2hdmA --num_samples 100
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
from pathlib import Path

import numpy as np

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _THIS_DIR)

from config import DATA_DIR, MAIN_DIR  # noqa: E402
from utils.align_utils import compute_tmscore_align  # noqa: E402
from utils.msa_utils import write_query_anchored_a3m  # noqa: E402

SH_WRAPPER = os.path.join(MAIN_DIR, "scripts", "shell", "RunBioEmu.sh")
DEFAULT_NUM_SAMPLES = 50


def _truth_chain_pdb(pair_dir: Path, pdbid: str, chain: str) -> str:
    """Locate a truth-fold backbone PDB (same resolution order as run_Boltz2)."""
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


def _split_ensemble_to_pdbs(out_root: Path) -> list:
    """Turn BioEmu's (topology.pdb + samples.xtc) into per-frame PDBs.

    BioEmu writes the ensemble as a trajectory; we split it into individual
    PDBs (samples/frame_XXXX.pdb) so each can be TM-scored. If per-frame PDBs
    already exist we reuse them. Returns the sorted list of PDB paths (possibly
    empty). Never raises.
    """
    samples_dir = out_root / "samples"
    existing = sorted(samples_dir.glob("frame_*.pdb")) if samples_dir.is_dir() else []
    if existing:
        return existing
    top = out_root / "topology.pdb"
    xtc = out_root / "samples.xtc"
    # Some BioEmu versions already dump per-model PDBs; pick those up too.
    if not (top.is_file() and xtc.is_file()):
        direct = sorted(out_root.rglob("*.pdb"))
        direct = [p for p in direct if p.name != "topology.pdb"]
        return direct
    try:
        import mdtraj as md
        traj = md.load(str(xtc), top=str(top))
        samples_dir.mkdir(parents=True, exist_ok=True)
        paths = []
        for i in range(traj.n_frames):
            p = samples_dir / f"frame_{i:04d}.pdb"
            traj[i].save_pdb(str(p))
            paths.append(p)
        return paths
    except Exception as e:
        print(f"    [bioemu] could not split ensemble in {out_root}: {e}",
              file=sys.stderr)
        return []


def run_one_cluster(pair_id, cluster_a3m, pdbid1, chain1, pdbid2, chain2,
                    pair_dir, num_samples, force_rerun=False):
    tag = cluster_a3m.stem
    out_root = pair_dir / "output_bioemu" / tag
    out_root.mkdir(parents=True, exist_ok=True)
    truth1 = _truth_chain_pdb(pair_dir, pdbid1, chain1)
    truth2 = _truth_chain_pdb(pair_dir, pdbid2, chain2)

    frames = _split_ensemble_to_pdbs(out_root)
    need_sample = force_rerun or not any(p.stat().st_size > 100 for p in frames)
    if need_sample:
        # BioEmu needs the query (a3m row 0) ungapped; cluster a3m keeps it
        # aligned/gapped. Project onto query columns -> ungapped-query a3m.
        bioemu_a3m = out_root / "bioemu_query.a3m"
        try:
            write_query_anchored_a3m(cluster_a3m, bioemu_a3m)
        except Exception as e:
            print(f"    [bioemu] could not build query-anchored a3m for {tag}: {e}",
                  file=sys.stderr, flush=True)
            return []
        cmd = (f"bash {shlex.quote(SH_WRAPPER)} "
               f"{shlex.quote(str(bioemu_a3m))} {shlex.quote(str(out_root))} "
               f"{int(num_samples)}")
        print(f"  [bioemu] {tag}: {cmd}", flush=True)
        try:
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"    [bioemu] failed for {tag}: {e}", file=sys.stderr, flush=True)
            return []
        frames = _split_ensemble_to_pdbs(out_root)
    else:
        print(f"  [bioemu] {tag}: skipping ({len(frames)} samples on disk)", flush=True)

    rows = []
    for i, pdb in enumerate(frames):
        tm1 = _safe_tm(pdb, truth1)
        tm2 = _safe_tm(pdb, truth2)
        rows.append({
            "pair_id": pair_id, "cluster": tag, "sample": i,
            "tm_fold1": tm1, "tm_fold2": tm2,
            "TMdiff": (tm1 - tm2) if (tm1 == tm1 and tm2 == tm2) else float("nan"),
            "pdb_path": str(pdb),
        })
    return rows


def summarize(rows):
    """One row per cluster: ensemble means + the population balance frac_f1."""
    by = {}
    for r in rows:
        by.setdefault(r["cluster"], []).append(r)
    out = []
    for tag, items in by.items():
        tm1 = np.array([x["tm_fold1"] for x in items], dtype=float)
        tm2 = np.array([x["tm_fold2"] for x in items], dtype=float)
        diff = tm1 - tm2
        valid = np.isfinite(diff)
        if valid.sum() == 0:
            continue
        n = int(valid.sum())
        # Population balance: fraction of the ensemble closer to fold 1. This is
        # the BioEmu-specific signal (bistability), not available from a single
        # AF/Boltz prediction.
        frac_f1 = float(np.mean(diff[valid] > 0))
        out.append({
            "pair_id": items[0]["pair_id"], "cluster": tag,
            "n_samples": n,
            "mean_tm_fold1": float(np.nanmean(tm1)),
            "mean_tm_fold2": float(np.nanmean(tm2)),
            "mean_TMdiff": float(np.nanmean(diff[valid])),
            "best_tm_fold1": float(np.nanmax(tm1)) if np.isfinite(tm1).any() else float("nan"),
            "best_tm_fold2": float(np.nanmax(tm2)) if np.isfinite(tm2).any() else float("nan"),
            "frac_toward_f1": round(frac_f1, 3),
        })
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-input", "--input", dest="input", required=True)
    ap.add_argument("--num_samples", type=int, default=DEFAULT_NUM_SAMPLES,
                    help=f"BioEmu samples per cluster (default {DEFAULT_NUM_SAMPLES}).")
    ap.add_argument("--max_clusters", type=int, default=0,
                    help="If >0, only process the first N ShallowMsa clusters "
                         "(pilot mode). 0 = all clusters.")
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

    n_total = len(cluster_a3ms)
    if args.max_clusters and args.max_clusters > 0 and n_total > args.max_clusters:
        cluster_a3ms = cluster_a3ms[: args.max_clusters]
        print(f"[bioemu] PILOT: capping to first {args.max_clusters} of "
              f"{n_total} clusters (--max_clusters)", flush=True)

    print(f"[bioemu] {pair_id}: {len(cluster_a3ms)} clusters, "
          f"{args.num_samples} samples/cluster", flush=True)
    t0 = time.time()
    all_rows = []
    for k, a3m in enumerate(cluster_a3ms, 1):
        print(f"[bioemu] cluster {k}/{len(cluster_a3ms)}: {a3m.stem} "
              f"({time.time()-t0:.0f}s elapsed)", flush=True)
        all_rows.extend(run_one_cluster(pair_id, a3m, pdbid1, chain1,
                                        pdbid2, chain2, pair_dir,
                                        args.num_samples))

    ana = pair_dir / "Analysis"
    ana.mkdir(parents=True, exist_ok=True)
    if all_rows:
        with open(ana / "df_bioemu_raw.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
            w.writeheader()
            w.writerows(all_rows)
    summary = summarize(all_rows)
    if summary:
        with open(ana / "df_bioemu.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
            w.writeheader()
            w.writerows(summary)
    print(f"[bioemu] done in {time.time()-t0:.0f}s; "
          f"summary -> {ana / 'df_bioemu.csv'} ({len(summary)} clusters)")


if __name__ == "__main__":
    main()
