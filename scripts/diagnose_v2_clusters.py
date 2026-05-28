#!/usr/bin/env python3
"""Diagnose the v2 clustering output for all pairs.

Scans Pipeline/FoldPairs/*/output_msa_cluster/ for ShallowMsa_*.a3m
files, reports cluster-count distribution, identifies which pairs
pass the K>=10 viability threshold, and writes the viable list to
docs/v2_viable_pairs.txt for downstream stages.

Also accepts --resolution coarse to diagnose the coarse-resolution
clusters in output_msa_cluster_coarse/CoarseMsa_*.a3m.

Usage:
    python3 scripts/diagnose_v2_clusters.py
    python3 scripts/diagnose_v2_clusters.py --resolution coarse
    python3 scripts/diagnose_v2_clusters.py --k-min 5
"""
from __future__ import annotations

import argparse
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--resolution", default="fine",
                    choices=["fine", "coarse"],
                    help="Which resolution to scan (default fine)")
    ap.add_argument("--k-min", type=int, default=10,
                    help="Viability threshold (default 10 fine, 5 coarse)")
    ap.add_argument("--output", default=None,
                    help="Output viable-pairs file (default depends on resolution)")
    args = ap.parse_args()

    if args.resolution == "fine":
        cdir_name = "output_msa_cluster"
        glob_pat = "ShallowMsa_*.a3m"
        default_out = "docs/v2_viable_pairs.txt"
        default_kmin = 10
    else:
        cdir_name = "output_msa_cluster_coarse"
        glob_pat = "CoarseMsa_*.a3m"
        default_out = "docs/v2_viable_pairs_coarse.txt"
        default_kmin = 5

    k_min = args.k_min if args.k_min != 10 else default_kmin
    out_path = Path(args.output or default_out)

    counts = []
    for d in sorted(Path("Pipeline/FoldPairs").glob("*/")):
        pid = d.name
        if pid in ("_s4pred_work", "jobs"):
            continue
        cdir = d / cdir_name
        n = sum(1 for _ in cdir.glob(glob_pat)) if cdir.is_dir() else 0
        counts.append((pid, n))

    n_total = len(counts)
    n_zero  = sum(1 for _, n in counts if n == 0)
    n_lt    = sum(1 for _, n in counts if 0 < n < k_min)
    n_ge    = sum(1 for _, n in counts if n >= k_min)
    print(f"=== {args.resolution.upper()} resolution scan ({cdir_name}) ===")
    print(f"Total pairs scanned: {n_total}")
    print(f"  K = 0          : {n_zero}")
    print(f"  0 < K < {k_min}      : {n_lt}  (dropped)")
    print(f"  K >= {k_min}        : {n_ge}  (VIABLE)")
    print()

    ok = sorted([n for _, n in counts if n >= k_min])
    if ok:
        print("Cluster-count distribution (viable pairs):")
        print(f"  min  = {ok[0]}")
        print(f"  Q1   = {ok[len(ok) // 4]}")
        print(f"  med  = {ok[len(ok) // 2]}")
        print(f"  Q3   = {ok[3 * len(ok) // 4]}")
        print(f"  max  = {ok[-1]}")
        print(f"  mean = {sum(ok) / len(ok):.1f}")
        print()

    viable = [p for p, n in counts if n >= k_min]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(viable) + "\n")
    print(f"Wrote {len(viable)} viable pairs -> {out_path}")
    print()

    dropped = [(p, n) for p, n in counts if n < k_min]
    if dropped:
        print(f"--- Dropped pairs (K < {k_min}) ---")
        for pid, n in sorted(dropped, key=lambda x: x[1]):
            print(f"  K={n:<3d}  {pid}")


if __name__ == "__main__":
    main()
