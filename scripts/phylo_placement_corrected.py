"""Phylogenetic placement on REGRESSION-RESIDUAL-CORRECTED F1c/F2c labels.

Companion to scripts/phylo_placement.py. The original script reads the
centered-classification labels from docs/fold_diversity_survey.csv (i.e.
labels driven by per-pair-median-centered TMdiff). This wrapper instead
reads the residual-corrected labels produced by
scripts/seq_divergence_correction.py:

  docs/fold_diversity_survey_corrected.csv  -- has column 'pref_corrected'
  with F1/F2/Amb labels computed from regression residuals.

Why: The phylogenetic-placement test on raw centered labels is NOT
automatically immune to a per-cluster sequence-divergence (PID-to-query)
gradient, because the sequence tree itself is built from sequence
distances that are correlated with PID-to-query. If the F1c/F2c labels
are driven by a monotone PID gradient, the labels can cluster
phylogenetically without reflecting independent biological signal. The
residual-corrected labels (with the per-pair linear PID dependence
removed per method) are the natural input to a deconfounded
phylogenetic test.

This script also applies Benjamini-Hochberg FDR adjustment across the
phylogenetic screen so that the screen-wide multiple-testing
question is answered explicitly.

Outputs:
  docs/phylo_placement_corrected.csv  -- one row per (pair, method) with:
    pair_id, method,
    n_leaves, n_F1c_leaves, n_F2c_leaves, n_Amc_leaves,
    parsimony_score, parsimony_p,
    nn_concordance, nn_concordance_p,
    nn3_concordance, nn3_concordance_p,
    parsimony_p_bh, nn_concordance_p_bh, nn3_concordance_p_bh,
    notes

Usage:
  python3 scripts/phylo_placement_corrected.py
  python3 scripts/phylo_placement_corrected.py --pairs 2namA_1uxmK,1zk9A_3jv6A
  python3 scripts/phylo_placement_corrected.py --n-perm 10000
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from collections import Counter
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# Make repo root importable when invoked from anywhere
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_THIS_DIR))

from config import DATA_DIR, PAIR_DIR_RE, TABLES_RES  # noqa: E402
from utils.ancestral_utils import build_coarse_tree, _name_internal_nodes  # noqa: E402
from scripts.phylo_placement import (  # noqa: E402
    analyze_pair_method,
    normalize_cluster,
    _per_pair_timeout,
    discover_pairs,
)


CORRECTED_SURVEY = os.path.join(TABLES_RES, "fold_diversity_survey_corrected.csv")
OUTPUT_CSV = os.path.join(TABLES_RES, "phylo_placement_corrected.csv")

METHOD_ORDER = ["AF2", "ESM", "MSAT", "DDG"]


def load_corrected_pref(survey_csv: str) -> Dict[Tuple[str, str, str], str]:
    """Build (pair, method, cluster) -> residual-corrected preference label.

    Uses the 'pref_corrected' column produced by
    scripts/seq_divergence_correction.py. Falls back to 'Amb' for
    rows where the residual is missing or the cell is empty.
    """
    df = pd.read_csv(survey_csv)
    if "pref_corrected" not in df.columns:
        raise SystemExit(
            f"Expected column 'pref_corrected' in {survey_csv}. "
            "Run scripts/seq_divergence_correction.py first."
        )
    # Use cluster_norm if present (already normalised), else fall back to cluster
    cluster_col = "cluster_norm" if "cluster_norm" in df.columns else "cluster"
    df[cluster_col] = df[cluster_col].astype(str).apply(normalize_cluster)

    out: Dict[Tuple[str, str, str], str] = {}
    for _, r in df.iterrows():
        label = r.get("pref_corrected")
        if not isinstance(label, str) or label not in ("F1", "F2", "Amb"):
            label = "Amb"
        key = (str(r["pair_id"]), str(r["method"]), str(r[cluster_col]))
        out[key] = label
    return out


def bh_adjust(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (NaN-safe, monotone)."""
    p = np.asarray(pvals, dtype=float)
    mask = np.isfinite(p)
    out = np.full_like(p, np.nan, dtype=float)
    if not mask.any():
        return out
    pv = p[mask]
    n = len(pv)
    order = np.argsort(pv)
    ranks = np.arange(1, n + 1)
    adj_sorted = pv[order] * n / ranks
    # Enforce monotonicity from the right
    adj_sorted = np.minimum.accumulate(adj_sorted[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.clip(adj_sorted, 0.0, 1.0)
    out[mask] = adj
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pairs", default="ALL",
                    help="Comma-separated pair IDs, or ALL.")
    ap.add_argument("--input", default=CORRECTED_SURVEY,
                    help=f"Corrected survey CSV (default: {CORRECTED_SURVEY})")
    ap.add_argument("--output", default=OUTPUT_CSV,
                    help=f"Output CSV (default: {OUTPUT_CSV})")
    ap.add_argument("--n-perm", type=int, default=1000,
                    help="Permutations for empirical null (default: 1000). "
                         "Increase to ~10000 for stable BH p-values.")
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--per-pair-timeout", type=int, default=120,
                    help="Seconds per pair before skipping (default: 120)")
    args = ap.parse_args()

    if not os.path.isfile(args.input):
        print(f"[error] Corrected survey CSV not found: {args.input}")
        print("        Run scripts/seq_divergence_correction.py first.")
        sys.exit(2)

    print(f"Loading residual-corrected labels: {args.input}")
    pref_lookup = load_corrected_pref(args.input)
    pairs_in_survey = sorted({k[0] for k in pref_lookup})
    print(f"  {len(pairs_in_survey)} pairs in corrected survey")

    if args.pairs.upper() == "ALL":
        pair_ids = discover_pairs()
    else:
        pair_ids = [s.strip() for s in args.pairs.split(",") if s.strip()]

    rows: List[dict] = []
    t_total = time.time()
    for i, pair_id in enumerate(sorted(pair_ids), 1):
        t_pair = time.time()
        print(f"\n[{i}/{len(pair_ids)}] {pair_id}: building coarse tree "
              f"(timeout {args.per_pair_timeout}s)...", flush=True)
        try:
            with _per_pair_timeout(args.per_pair_timeout):
                tree, tag_to_label, _ = build_coarse_tree(pair_id)
                _name_internal_nodes(tree)
        except TimeoutError:
            print(f"  TIMEOUT after {args.per_pair_timeout}s on coarse-tree "
                  f"build; skipping {pair_id}", flush=True)
            continue
        except FileNotFoundError as e:
            print(f"  skipped: {e}", flush=True)
            continue
        except Exception as e:
            print(f"  failed coarse-tree build: {e}", flush=True)
            continue
        print(f"  tree built in {time.time()-t_pair:.1f}s, "
              f"{len(tag_to_label)} leaves", flush=True)

        for method in METHOD_ORDER:
            cluster_pref: Dict[str, str] = {
                cluster: lab
                for (pid, m, cluster), lab in pref_lookup.items()
                if pid == pair_id and m == method
            }
            if not cluster_pref:
                continue
            try:
                with _per_pair_timeout(args.per_pair_timeout):
                    row = analyze_pair_method(
                        pair_id, method,
                        tree, tag_to_label,
                        cluster_pref,
                        n_perm=args.n_perm, seed=args.seed,
                        do_d_stat=False,  # D-stat is slow; skip
                    )
            except TimeoutError:
                print(f"    TIMEOUT on {pair_id} / {method}; skipping",
                      flush=True)
                continue
            except Exception as e:
                print(f"    failed {pair_id} / {method}: {e}", flush=True)
                continue
            if row is None:
                continue
            rows.append(row)

    if not rows:
        print("[error] no rows produced")
        sys.exit(1)

    df = pd.DataFrame(rows)
    # Apply BH FDR per p-value column over the full screen
    for pcol in ("parsimony_p", "nn_concordance_p", "nn3_concordance_p"):
        if pcol in df.columns:
            df[pcol + "_bh"] = bh_adjust(df[pcol].to_numpy(dtype=float))

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    df.to_csv(args.output, index=False)
    print(f"\nWrote {args.output}: {len(df)} (pair, method) rows "
          f"in {time.time()-t_total:.1f}s total")

    # Headline summary
    print("\n=== BH-significant (post-correction) at alpha = 0.05 ===")
    for pcol in ("parsimony_p_bh", "nn_concordance_p_bh", "nn3_concordance_p_bh"):
        if pcol not in df.columns:
            continue
        passing = df[df[pcol] <= 0.05].sort_values(pcol)
        print(f"\n-- {pcol} <= 0.05 ({len(passing)} rows) --")
        for _, r in passing.iterrows():
            print(f"  {r['pair_id']:<18s} {r['method']:<5s}  "
                  f"raw_p={r[pcol.replace('_bh', '')]:.4f}  "
                  f"BH_p={r[pcol]:.4f}")


if __name__ == "__main__":
    main()
