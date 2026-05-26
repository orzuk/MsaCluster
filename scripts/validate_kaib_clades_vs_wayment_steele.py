#!/usr/bin/env python3
"""
Validate that our 5jytA_2qkeE (KaiB) F1/F2 clade split corresponds to
the ground-state / FS-state clade split that Wayment-Steele et al. 2024
mapped at per-variant resolution.

Approach: the Wayment-Steele paper names five experimentally validated
KaiB variants:
  ground state -- R. sphaeroides, T. elongatus, S. elongatus
  FS state     -- L. pneumonia, T. elongatus vestitus (KaiB-TV4)

We search the per-pair Deep MSA (jackhmmer/MMseqs2 .a3m output) for
sequences whose FASTA-style header contains the organism name string,
then look up which ShallowMsa cluster each match belongs to. If the
3 ground-state variants land predominantly in our "F1-leaning" clade
(clusters with red on DDG/Boltz2 in the heat-map: typically S04, S09,
S21, S12, S13, S01, S16, S05) and the 2 FS-state variants land in our
"F2-leaning" clade (S08, S10, S19, S20), the clade split is biologically
the same signal Wayment-Steele saw.

Usage:
    python3 scripts/validate_kaib_clades_vs_wayment_steele.py
        --pair 5jytA_2qkeE
        [--deep-msa Pipeline/FoldPairs/5jytA_2qkeE/output_get_msa/DeepMsa.a3m]
        [--clusters-dir Pipeline/FoldPairs/5jytA_2qkeE/output_msa_cluster]
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

# Strings to grep for in header text. Case-insensitive substring match.
# Map: human label -> (organism keywords list, expected fold state)
ORGANISM_KEYS = {
    "R. sphaeroides KaiB":   (["sphaeroides"], "ground"),
    "T. elongatus KaiB":     (["thermosynechococcus elongatus",
                              "T.elongatus", "T. elongatus"], "ground"),
    "S. elongatus KaiB":     (["synechococcus elongatus",
                              "S.elongatus", "S. elongatus"], "ground"),
    "L. pneumonia KaiB":     (["legionella pneumonia",
                              "L.pneumonia", "L. pneumonia",
                              "legionella"], "FS"),
    "KaiB-TV4 (T. elongatus vestitus)": (["vestitus",
                                         "T.elongatus vestitus",
                                         "T. elongatus vestitus"], "FS"),
}


def parse_a3m_headers(a3m_path: Path) -> list[tuple[int, str]]:
    """Return list of (line_index_of_header, header_line) for a multi-FASTA
    file. Order matters because jackhmmer-style MSAs typically order hits.
    """
    headers = []
    with a3m_path.open() as f:
        for i, line in enumerate(f):
            if line.startswith(">"):
                headers.append((i, line.strip()))
    return headers


def find_organism_matches(headers: list[tuple[int, str]]
                          ) -> dict[str, list[tuple[int, str, str]]]:
    """For each organism label, return list of (header_index, full_header,
    matched_keyword) for headers matching any of its keywords."""
    out: dict[str, list[tuple[int, str, str]]] = {}
    for label, (keys, _state) in ORGANISM_KEYS.items():
        hits = []
        for idx, h in headers:
            h_lower = h.lower()
            for k in keys:
                if k.lower() in h_lower:
                    hits.append((idx, h, k))
                    break
        out[label] = hits
    return out


def parse_cluster_assignments(clusters_dir: Path) -> dict[str, str]:
    """Read all ShallowMsa_NNN.a3m in clusters_dir; return mapping of
    sequence-ID (FASTA header up to first whitespace) -> cluster tag.
    Headers are normalised by stripping leading '>'.
    """
    out: dict[str, str] = {}
    for f in sorted(clusters_dir.glob("ShallowMsa_*.a3m")):
        tag_match = re.search(r"ShallowMsa_0*(\d+)\.a3m", f.name)
        if not tag_match:
            continue
        cluster_tag = f"S{int(tag_match.group(1)):02d}"
        with f.open() as fh:
            for line in fh:
                if line.startswith(">"):
                    seq_id = line.lstrip(">").strip().split()[0]
                    out[seq_id] = cluster_tag
    return out


def match_headers_to_clusters(
    header_hits: dict[str, list[tuple[int, str, str]]],
    cluster_map: dict[str, str],
) -> dict[str, list[tuple[str, str, str]]]:
    """For each organism label, map its matching headers to a cluster tag
    by looking up the first whitespace-token of the header in cluster_map.
    Returns: {label: [(short_header_id, cluster_tag, matched_keyword), ...]}
    """
    out: dict[str, list[tuple[str, str, str]]] = {}
    for label, hits in header_hits.items():
        rows = []
        for _, h, keyword in hits:
            seq_id = h.lstrip(">").split()[0]
            cluster = cluster_map.get(seq_id, "(not in any cluster file)")
            rows.append((seq_id, cluster, keyword))
        out[label] = rows
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pair", default="5jytA_2qkeE")
    ap.add_argument("--deep-msa", default=None,
                    help="Path to Deep MSA .a3m (default: "
                         "Pipeline/FoldPairs/<pair>/output_get_msa/DeepMsa.a3m)")
    ap.add_argument("--clusters-dir", default=None,
                    help="Path to ShallowMsa_*.a3m directory (default: "
                         "Pipeline/FoldPairs/<pair>/output_msa_cluster)")
    args = ap.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    pair_dir = repo_root / "Pipeline" / "FoldPairs" / args.pair
    deep_msa = Path(args.deep_msa) if args.deep_msa else \
        pair_dir / "output_get_msa" / "DeepMsa.a3m"
    clusters_dir = Path(args.clusters_dir) if args.clusters_dir else \
        pair_dir / "output_msa_cluster"

    if not deep_msa.is_file():
        print(f"ERROR: Deep MSA not found at {deep_msa}", file=sys.stderr)
        sys.exit(1)
    if not clusters_dir.is_dir():
        print(f"ERROR: Clusters dir not found at {clusters_dir}",
              file=sys.stderr)
        sys.exit(1)

    print(f"[validate_kaib] pair: {args.pair}")
    print(f"[validate_kaib] Deep MSA: {deep_msa}")
    print(f"[validate_kaib] Clusters dir: {clusters_dir}")
    print()

    headers = parse_a3m_headers(deep_msa)
    print(f"[validate_kaib] Deep MSA has {len(headers)} headers")
    header_hits = find_organism_matches(headers)
    n_total_hits = sum(len(v) for v in header_hits.values())
    print(f"[validate_kaib] Found {n_total_hits} organism-keyword matches "
          f"across {sum(1 for v in header_hits.values() if v)} of "
          f"{len(ORGANISM_KEYS)} expected organism labels")
    print()

    print(f"[validate_kaib] Loading cluster assignments...")
    cluster_map = parse_cluster_assignments(clusters_dir)
    print(f"[validate_kaib] {len(cluster_map)} sequence -> cluster mappings "
          f"across {len(set(cluster_map.values()))} clusters")
    print()

    matched = match_headers_to_clusters(header_hits, cluster_map)

    print("=" * 78)
    print("Wayment-Steele KaiB variants in our Deep MSA + ShallowMsa clusters")
    print("=" * 78)
    for label, rows in matched.items():
        expected_state = ORGANISM_KEYS[label][1]
        print(f"\n{label}   [expected: {expected_state} state]")
        if not rows:
            print("    (no matching headers in Deep MSA)")
            continue
        cluster_counts: dict[str, int] = defaultdict(int)
        for seq_id, cluster, keyword in rows[:20]:
            print(f"    cluster={cluster:<10s} keyword='{keyword:<20s}' "
                  f"seq_id={seq_id[:60]}")
            cluster_counts[cluster] += 1
        if len(rows) > 20:
            print(f"    ... and {len(rows) - 20} more (total {len(rows)} hits)")
            for seq_id, cluster, keyword in rows[20:]:
                cluster_counts[cluster] += 1
        print(f"    cluster distribution: "
              f"{dict(sorted(cluster_counts.items()))}")

    print()
    print("=" * 78)
    print("Validation summary")
    print("=" * 78)
    # In our 5jytA_2qkeE phytree, the F1-leaning clade contains:
    f1_clusters = {"S04", "S09", "S21", "S12", "S13", "S01", "S16", "S05"}
    f2_clusters = {"S08", "S10", "S19", "S20"}

    expected_F1 = {label for label, (_keys, state) in ORGANISM_KEYS.items()
                   if state == "ground"}
    expected_F2 = {label for label, (_keys, state) in ORGANISM_KEYS.items()
                   if state == "FS"}

    print(f"F1-leaning clade (ground state per WS):  {sorted(f1_clusters)}")
    print(f"F2-leaning clade (FS state per WS):      {sorted(f2_clusters)}")
    print()
    n_correct = 0
    n_evaluated = 0
    for label, rows in matched.items():
        if not rows:
            print(f"  {label}: NO MATCHES -- cannot validate")
            continue
        # Take majority cluster across the organism's matches
        cluster_counts: dict[str, int] = defaultdict(int)
        for _seq_id, cluster, _keyword in rows:
            cluster_counts[cluster] += 1
        if not cluster_counts:
            continue
        top_cluster = max(cluster_counts.items(), key=lambda kv: kv[1])[0]
        expected_state = ORGANISM_KEYS[label][1]
        if expected_state == "ground":
            ok = top_cluster in f1_clusters
            expected_clade = "F1 (ground)"
        else:
            ok = top_cluster in f2_clusters
            expected_clade = "F2 (FS)"
        marker = "OK " if ok else "MISS"
        print(f"  [{marker}] {label}: majority cluster {top_cluster} "
              f"-- expected to be in {expected_clade} clade")
        n_evaluated += 1
        if ok:
            n_correct += 1

    if n_evaluated > 0:
        print(f"\nOverall: {n_correct}/{n_evaluated} variants land in the "
              f"expected clade. Random baseline ~50% if cluster assignment "
              f"is uninformative.")


if __name__ == "__main__":
    main()
