"""Regenerate case-study cluster figures against current AHC clusters.

The paper's case-study subsections (SUN2 = 4dxrA_4dxtA and Phytochrome
= 4o01D_4o0pA) reference figures whose filenames bake in HDBSCAN
cluster IDs (cluster_47.png, cluster_18.png, cluster47vs18.png,
cluster45vs66.png). Those cluster IDs no longer exist under the current
AHC clustering. This script:

  1. Finds the best-F1-favoring and best-F2-favoring AHC cluster IDs
     for each case-study pair from df_af.csv.
  2. Reuses the existing PyMOL-based renderer
     utils.plot_3d._render_true_vs_best_models_generic (or its
     align_and_visualize_proteins helper) to produce the truth-vs-best
     overlays for each pair, both folds.
  3. Renames outputs to generic AHC-agnostic filenames the paper now
     references (best_F1.png, best_F2.png, bestF1_vs_bestF2.png).
  4. Writes a small JSON manifest with the AHC cluster IDs used so the
     paper can quote them.

Outputs (per pair, all under figs/case_study/<pair_id>/):
  best_F1.png         AHC best-F1 cluster AF2 prediction overlaid on fold1 truth
  best_F2.png         AHC best-F2 cluster AF2 prediction overlaid on fold2 truth
  bestF1_vs_bestF2.png  the two predictions overlaid on each other
  manifest.json       {best_f1_cluster, best_f2_cluster, best_f1_tm, best_f2_tm}

Usage:
    python3 scripts/regen_case_study_figs.py
    python3 scripts/regen_case_study_figs.py --pairs 4dxrA_4dxtA
"""

from __future__ import annotations

import os
# Headless rendering env vars set BEFORE matplotlib / pymol import.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-cache-orzuk")
os.environ.setdefault("PYTHONUNBUFFERED", "1")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)

import argparse
import json
import re
import sys
import time
from pathlib import Path

# Allow running from repo root
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

print("[case-study-figs] importing modules...", flush=True)
_t0 = time.time()
import pandas as pd
from utils.plot_3d import (
    _best_cluster_from_df,
    _pred_pdb_path,
    align_and_visualize_proteins,
    _add_png_legend,
)
from config import DATA_DIR
print(f"[case-study-figs] imports done in {time.time()-_t0:.1f}s", flush=True)

CASE_STUDIES = {
    "4dxrA_4dxtA": ("4dxrA", "4dxtA", "SUN2"),
    "4o01D_4o0pA": ("4o01D", "4o0pA", "Phytochrome"),
}


def regen_pair(pair_id: str, pdbid1_full: str, pdbid2_full: str, label: str,
               out_root: Path) -> dict:
    pdbid1, chain1 = pdbid1_full[:4], pdbid1_full[4]
    pdbid2, chain2 = pdbid2_full[:4], pdbid2_full[4]

    out_dir = out_root / pair_id
    out_dir.mkdir(parents=True, exist_ok=True)

    best_f1 = _best_cluster_from_df(pair_id, model_ver="2", fold_idx=0)
    best_f2 = _best_cluster_from_df(pair_id, model_ver="2", fold_idx=1)
    print(f"  [{pair_id}] best_F1 cluster={best_f1}  best_F2 cluster={best_f2}", flush=True)

    pair_dir = Path(DATA_DIR) / pair_id
    truth_f1 = pair_dir / f"{pdbid1_full}.pdb"
    truth_f2 = pair_dir / f"{pdbid2_full}.pdb"

    manifest: dict = {
        "pair_id": pair_id,
        "label": label,
        "fold1": pdbid1_full,
        "fold2": pdbid2_full,
    }

    if best_f1 is not None:
        pred_f1 = _pred_pdb_path(pair_id, "AF", "2", best_f1, pdbid1, chain1)
        if Path(pred_f1).is_file() and truth_f1.is_file():
            out_png = out_dir / "best_F1.png"
            try:
                align_and_visualize_proteins(
                    str(truth_f1), str(pred_f1), str(out_png),
                    open_environment=False,
                )
                _add_png_legend(
                    str(out_png), pdbid1_full, f"AF2(cluster_{best_f1})",
                    title=f"{label} best-F1 cluster vs {pdbid1_full} truth",
                )
                print(f"    wrote {out_png}", flush=True)
            except Exception as e:
                print(f"    [WARN] best_F1.png failed: {e}", flush=True)
        else:
            print(f"    [WARN] missing files for best_F1: truth={truth_f1} pred={pred_f1}", flush=True)
        manifest["best_f1_cluster"] = best_f1

    if best_f2 is not None:
        pred_f2 = _pred_pdb_path(pair_id, "AF", "2", best_f2, pdbid2, chain2)
        if Path(pred_f2).is_file() and truth_f2.is_file():
            out_png = out_dir / "best_F2.png"
            try:
                align_and_visualize_proteins(
                    str(truth_f2), str(pred_f2), str(out_png),
                    open_environment=False,
                )
                _add_png_legend(
                    str(out_png), pdbid2_full, f"AF2(cluster_{best_f2})",
                    title=f"{label} best-F2 cluster vs {pdbid2_full} truth",
                )
                print(f"    wrote {out_png}", flush=True)
            except Exception as e:
                print(f"    [WARN] best_F2.png failed: {e}", flush=True)
        else:
            print(f"    [WARN] missing files for best_F2: truth={truth_f2} pred={pred_f2}", flush=True)
        manifest["best_f2_cluster"] = best_f2

    # bestF1 vs bestF2: overlay the two predictions on each other (no truth)
    if best_f1 is not None and best_f2 is not None:
        pred_f1 = _pred_pdb_path(pair_id, "AF", "2", best_f1, pdbid1, chain1)
        pred_f2 = _pred_pdb_path(pair_id, "AF", "2", best_f2, pdbid2, chain2)
        if Path(pred_f1).is_file() and Path(pred_f2).is_file():
            out_png = out_dir / "bestF1_vs_bestF2.png"
            try:
                align_and_visualize_proteins(
                    pred_f1, pred_f2, str(out_png), open_environment=False,
                )
                _add_png_legend(
                    str(out_png), f"AF2(cluster_{best_f1})",
                    f"AF2(cluster_{best_f2})",
                    title=f"{label} best-F1 vs best-F2 cluster predictions",
                )
                print(f"    wrote {out_png}", flush=True)
            except Exception as e:
                print(f"    [WARN] bestF1_vs_bestF2.png failed: {e}", flush=True)

    with open(out_dir / "manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"    wrote {out_dir / 'manifest.json'}", flush=True)
    return manifest


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pairs", default="ALL",
                    help="comma-separated pair ids, or ALL (default: all case studies)")
    ap.add_argument("--out-dir", default="figs/case_study",
                    help="output root (default: figs/case_study)")
    args = ap.parse_args()

    if args.pairs.upper() == "ALL":
        targets = list(CASE_STUDIES.items())
    else:
        wanted = {p.strip() for p in args.pairs.split(",") if p.strip()}
        targets = [(p, v) for p, v in CASE_STUDIES.items() if p in wanted]

    out_root = Path(args.out_dir)
    for pair_id, (f1, f2, label) in targets:
        print(f"\n=== {pair_id} ({label}) ===", flush=True)
        regen_pair(pair_id, f1, f2, label, out_root)


if __name__ == "__main__":
    main()
