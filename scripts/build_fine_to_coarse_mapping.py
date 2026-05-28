#!/usr/bin/env python3
"""Build the fine->coarse cluster mapping for one pair.

Two cluster sets (fine + coarse) live in parallel output directories
under a pair's Pipeline/FoldPairs/<pair>/ tree:

    output_msa_cluster/          ShallowMsa_NNN.a3m   (fine, ~10-50 seqs each)
    output_msa_cluster_coarse/   CoarseMsa_NNN.a3m    (coarse, ~100-500 seqs)

Both are tree-cuts on the SAME phylogenetic tree, just at different
depths. By construction every fine cluster is nested in exactly one
coarse cluster (modulo small-cluster merging at boundaries).

This script computes that nesting by sequence-membership majority:

    For each fine cluster F:
      find the coarse cluster C whose membership has the largest
      intersection with F. That is F's parent.

The output CSV columns:
    fine_id, coarse_id, n_fine, n_coarse, overlap, frac_overlap,
    medoid_fine_header

(frac_overlap = overlap / n_fine; should be ~1.0 when nesting is clean.)

Usage:
    python3 build_fine_to_coarse_mapping.py \
        --fine-dir output_msa_cluster \
        --coarse-dir output_msa_cluster_coarse \
        --output output_msa_cluster/fine_to_coarse.csv
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def read_a3m_headers(a3m_path: Path) -> set[str]:
    """Return the set of sequence ids (first token of each '>' header).
    Skips the first row if it looks like a query/seed marker."""
    ids: set[str] = set()
    with a3m_path.open() as f:
        for line in f:
            if line.startswith(">"):
                tok = line[1:].split()[0]
                if tok:
                    ids.add(tok)
    return ids


def first_header(a3m_path: Path) -> str:
    """Return the first record's id (used as 'medoid header' label)."""
    with a3m_path.open() as f:
        for line in f:
            if line.startswith(">"):
                return line[1:].split()[0]
    return ""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fine-dir", required=True,
                    help="Dir containing ShallowMsa_NNN.a3m")
    ap.add_argument("--coarse-dir", required=True,
                    help="Dir containing CoarseMsa_NNN.a3m")
    ap.add_argument("--output", required=True,
                    help="Output CSV path")
    ap.add_argument("--fine-pattern", default="ShallowMsa_*.a3m")
    ap.add_argument("--coarse-pattern", default="CoarseMsa_*.a3m")
    args = ap.parse_args()

    fine_dir = Path(args.fine_dir)
    coarse_dir = Path(args.coarse_dir)
    if not fine_dir.is_dir():
        sys.exit(f"[map] fine dir not found: {fine_dir}")
    if not coarse_dir.is_dir():
        sys.exit(f"[map] coarse dir not found: {coarse_dir}")

    fine_files = sorted(fine_dir.glob(args.fine_pattern))
    coarse_files = sorted(coarse_dir.glob(args.coarse_pattern))
    if not fine_files:
        sys.exit(f"[map] no fine files matching {args.fine_pattern} in {fine_dir}")
    if not coarse_files:
        sys.exit(f"[map] no coarse files matching {args.coarse_pattern} in {coarse_dir}")

    print(f"[map] loading {len(coarse_files)} coarse clusters...", flush=True)
    coarse_members: dict[str, set[str]] = {}
    for c in coarse_files:
        coarse_members[c.stem] = read_a3m_headers(c)

    print(f"[map] loading {len(fine_files)} fine clusters...", flush=True)
    rows = []
    for f in fine_files:
        fine_ids = read_a3m_headers(f)
        n_fine = len(fine_ids)
        if n_fine == 0:
            print(f"[map] WARN fine cluster {f.stem} is empty; skipping")
            continue
        best_coarse, best_overlap = None, 0
        for c_stem, c_ids in coarse_members.items():
            overlap = len(fine_ids & c_ids)
            if overlap > best_overlap:
                best_overlap = overlap
                best_coarse = c_stem
        if best_coarse is None:
            print(f"[map] WARN fine {f.stem}: no overlap with any coarse cluster")
            continue
        n_coarse = len(coarse_members[best_coarse])
        rows.append({
            "fine_id":          f.stem,
            "coarse_id":        best_coarse,
            "n_fine":           n_fine,
            "n_coarse":         n_coarse,
            "overlap":          best_overlap,
            "frac_overlap":     round(best_overlap / n_fine, 4),
            "medoid_fine_header": first_header(f),
        })

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"[map] wrote {len(rows)} rows -> {out_path}")

    # Quality summary
    if rows:
        clean = sum(1 for r in rows if r["frac_overlap"] >= 0.95)
        any_overlap = sum(1 for r in rows if r["frac_overlap"] >= 0.5)
        print(f"[map] cleanly nested (frac>=0.95): {clean}/{len(rows)} "
              f"({100*clean/len(rows):.0f}%)")
        print(f"[map] partial overlap (frac>=0.50): {any_overlap}/{len(rows)} "
              f"({100*any_overlap/len(rows):.0f}%)")
        bad = [r for r in rows if r["frac_overlap"] < 0.5]
        if bad:
            print(f"[map] WARN {len(bad)} fine clusters with frac_overlap<0.5 "
                  f"(non-nested due to small-cluster merging at boundaries):")
            for r in bad[:5]:
                print(f"    {r['fine_id']} -> {r['coarse_id']}  "
                      f"({r['overlap']}/{r['n_fine']} = {r['frac_overlap']})")


if __name__ == "__main__":
    main()
