#!/usr/bin/env python3
"""Run ancestral fold-preference reconstruction across fold-switching pairs.

Usage:
  python scripts/run_ancestral_analysis.py --pairs ALL --delta 0.05
  python scripts/run_ancestral_analysis.py --pairs 2qqjA_4qdsA,1fzpD_2frhA
  python scripts/run_ancestral_analysis.py --pairs ALL --summary-only
"""

import os
# Headless rendering env vars MUST be set BEFORE importing matplotlib / ete3.
# Without these, on a headless compute node ete3.TreeStyle hangs waiting on Qt
# and matplotlib fails to find a writable config dir.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-cache-orzuk")
os.environ.setdefault("PYTHONUNBUFFERED", "1")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)

import argparse
import sys
import time

import pandas as pd

# Allow running from repo root
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

print(f"[ancestral] importing modules...", flush=True)
_t0 = time.time()
from config import DATA_DIR, PAIR_DIR_RE, TABLES_RES
from utils.ancestral_utils import analyze_pair
print(f"[ancestral] imports done in {time.time()-_t0:.1f}s", flush=True)


def discover_eligible_pairs() -> list[str]:
    """Find pairs that have both a tree file and df_af.csv."""
    eligible = []
    if not os.path.isdir(DATA_DIR):
        print(f"[error] DATA_DIR not found: {DATA_DIR}")
        return eligible

    for name in sorted(os.listdir(DATA_DIR)):
        if not PAIR_DIR_RE.match(name):
            continue
        pair_dir = os.path.join(DATA_DIR, name)
        tree_file = os.path.join(pair_dir, "output_phytree", "DeepMsa_tree.nwk")
        af_csv = os.path.join(pair_dir, "Analysis", "df_af.csv")
        cluster_dir = os.path.join(pair_dir, "output_msa_cluster")
        if (os.path.isfile(tree_file) and os.path.isfile(af_csv)
                and os.path.isdir(cluster_dir)):
            eligible.append(name)
    return eligible


def main():
    parser = argparse.ArgumentParser(
        description="Ancestral fold-preference reconstruction")
    parser.add_argument("--pairs", type=str, default="ALL",
                        help="Comma-separated pair IDs, or ALL")
    parser.add_argument("--delta", type=float, default=0.05,
                        help="Threshold for F1/F2/Amb assignment (default 0.05)")
    parser.add_argument("--n-perm", type=int, default=1000,
                        help="Number of permutations for significance test")
    parser.add_argument("--summary-only", action="store_true",
                        help="Only print/save summary, skip per-pair figures")
    parser.add_argument("--output-csv", type=str, default=None,
                        help="Path for summary CSV (default: docs/ancestral_summary.csv)")
    args = parser.parse_args()

    # Discover pairs
    if args.pairs.upper() == "ALL":
        pairs = discover_eligible_pairs()
        print(f"Found {len(pairs)} eligible pairs (have tree + df_af.csv)")
    else:
        pairs = [p.strip() for p in args.pairs.split(",") if p.strip()]

    if not pairs:
        print("[error] No pairs to process.")
        return

    # Process each pair
    results = []
    t0 = time.time()
    for i, pair_id in enumerate(pairs):
        t_pair = time.time()
        print(f"\n[{i+1}/{len(pairs)}] {pair_id} ... starting", flush=True)
        try:
            r = analyze_pair(pair_id, delta=args.delta, n_perm=args.n_perm)
            results.append(r)
            pref_str = f"F1={r['n_f1']} F2={r['n_f2']} Amb={r['n_amb']}"
            print(f"  {r['n_clusters']} clusters | {pref_str} | "
                  f"pattern={r['pattern']} | root={r['root_state']} | "
                  f"transitions={r['n_fold_transitions']}"
                  + (f" | p={r['p_value']:.3f}" if r['p_value'] is not None else "")
                  + f"  ({time.time()-t_pair:.1f}s)", flush=True)
        except Exception as e:
            print(f"  [FAILED] {e}  ({time.time()-t_pair:.1f}s)", flush=True)
            results.append({"pair_id": pair_id, "error": str(e)})

    elapsed = time.time() - t0
    print(f"\n{'='*60}")
    print(f"Processed {len(pairs)} pairs in {elapsed:.1f}s")

    # Build summary DataFrame
    df = pd.DataFrame(results)
    if df.empty:
        print("[warn] No results to summarize.")
        return

    # Summary statistics
    valid = df.dropna(subset=["pattern"]) if "pattern" in df.columns else df
    if not valid.empty and "pattern" in valid.columns:
        print(f"\nPattern distribution (n={len(valid)}):")
        for pat, count in valid["pattern"].value_counts().items():
            print(f"  {pat:25s} {count:3d} ({100*count/len(valid):.0f}%)")

        print(f"\nRoot state distribution:")
        for st, count in valid["root_state"].value_counts().items():
            print(f"  {st:5s} {count:3d}")

        has_pval = valid["p_value"].notna()
        if has_pval.any():
            sig = (valid.loc[has_pval, "p_value"] < 0.05).sum()
            tested = has_pval.sum()
            print(f"\nPermutation test: {sig}/{tested} pairs significant (p<0.05)")

    # Save CSV
    csv_path = args.output_csv or os.path.join(TABLES_RES, "ancestral_summary.csv")
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    df.to_csv(csv_path, index=False)
    print(f"\nSaved summary: {csv_path}")


if __name__ == "__main__":
    main()
