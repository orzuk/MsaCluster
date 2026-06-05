#!/usr/bin/env python3
"""Fetch the molecule name + source organism for each fold-switch pair's two
chains from RCSB, and write a single mapping table the rest of the project can
join on:  docs/pair_protein_names.csv

Columns: pair_id, pdbid1, chain1, name1, organism1, pdbid2, chain2, name2,
         organism2, display_name

`display_name` is a short human label: the fold-1 molecule name, plus the fold-2
name in parentheses if it differs. Reuses the cached RCSB GraphQL fetcher in
enrich_pdb_mining_candidates.py (data/_rcsb_meta_cache/), so reruns are cheap.

Usage (from repo root):
    python3 scripts/fetch_pair_protein_names.py
    python3 scripts/fetch_pair_protein_names.py --pairs-csv docs/triggers_from_pdb.csv
"""
from __future__ import annotations
import argparse
import os
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts.enrich_pdb_mining_candidates import fetch_entry, chain_name


def display_name(n1, n2):
    if n1 and n2 and n1.strip().lower() != n2.strip().lower():
        return f"{n1} / {n2}"
    return n1 or n2 or ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs-csv", default="docs/triggers_from_pdb.csv",
                    help="CSV with pair_id, pdbid1, chain1_used, pdbid2, chain2_used")
    ap.add_argument("--out", default="docs/pair_protein_names.csv")
    ap.add_argument("--cache-dir", default="data/_rcsb_meta_cache")
    ap.add_argument("--sleep", type=float, default=0.1,
                    help="polite delay between uncached fetches (s)")
    args = ap.parse_args()

    cache = Path(args.cache_dir); cache.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.pairs_csv)

    entries = {}  # pdb_id -> entry payload (fetched once)

    def get_entry(pdb):
        pdb = str(pdb).lower()
        if pdb not in entries:
            cached = (cache / f"{pdb}.json").is_file()
            entries[pdb] = fetch_entry(pdb, cache)
            if not cached:
                time.sleep(args.sleep)
        return entries[pdb]

    rows = []
    for _, r in df.iterrows():
        p1, c1 = str(r["pdbid1"]), str(r["chain1_used"])
        p2, c2 = str(r["pdbid2"]), str(r["chain2_used"])
        n1, o1 = chain_name(get_entry(p1), c1)
        n2, o2 = chain_name(get_entry(p2), c2)
        rows.append(dict(pair_id=r["pair_id"], pdbid1=p1, chain1=c1, name1=n1,
                         organism1=o1, pdbid2=p2, chain2=c2, name2=n2,
                         organism2=o2, display_name=display_name(n1, n2)))
        print(f"  {r['pair_id']}: {display_name(n1, n2)}  [{o1 or '?'} / {o2 or '?'}]")

    out = pd.DataFrame(rows)
    Path(os.path.dirname(args.out)).mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    n_named = int(out["name1"].notna().sum())
    print(f"\n[ok] wrote {args.out}: {len(out)} pairs, {n_named} with a fold-1 name")


if __name__ == "__main__":
    main()
