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
import re
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts.enrich_pdb_mining_candidates import fetch_entry, chain_name

# Short common-name aliases for the well-known fold-switchers, so the canonical
# table carries the names scientists actually use (RCSB returns "Lymphotactin",
# we also want "XCL1"). Kept here so re-running the fetch preserves them.
KNOWN_ALIASES = {
    "2n54B_2hdmA": "XCL1",
    "1zk9A_3jv6A": "RelB",
    "2namA_1uxmK": "SOD1",
    "5jytA_2qkeE": "KaiB",
    "6c6sD_2ougC": "RfaH",
    "4qhhA_4qhfA": "Selecase",
}


def clean_name(s, maxlen=46):
    """Normalize an RCSB molecule description for display: title-case the
    legacy ALL-CAPS entries (keeping short acronyms/IDs like L1, GP2, RNA),
    and truncate very long chimera blurbs at a word boundary."""
    if s is None or (isinstance(s, float)):
        return None
    s = str(s).strip()
    if not s:
        return None
    letters = [c for c in s if c.isalpha()]
    if letters and all(c.isupper() for c in letters):          # legacy all-caps
        def tc(w):
            core = re.sub(r"[^A-Za-z0-9]", "", w)
            if core and len(core) <= 4 and not any(c.islower() for c in core):
                return w                                        # keep L1, GP2, RNA...
            return w.title()
        s = " ".join(tc(w) for w in s.split())
    if len(s) > maxlen:
        s = s[:maxlen].rsplit(" ", 1)[0].rstrip(" ,;:-") + "…"
    return s


def display_name(n1, n2, alias=None):
    if alias:
        return alias
    c1, c2 = clean_name(n1), clean_name(n2)
    if c1 and c2:
        # compare on alphanumerics only, so "Alpha1-Antitrypsin" == "Alpha-1-antitrypsin"
        a = re.sub(r"[^a-z0-9]", "", c1.lower().rstrip("…"))
        b = re.sub(r"[^a-z0-9]", "", c2.lower().rstrip("…"))
        if a == b or a in b or b in a:                         # same protein, two folds
            return c1 if len(c1) <= len(c2) else c2
        return f"{c1} / {c2}"
    return c1 or c2 or ""


def combine_organism(o1, o2):
    o1 = (o1 or "").strip() or None
    o2 = (o2 or "").strip() or None
    if o1 and o2 and o1 != o2:
        return f"{o1} / {o2}"
    return o1 or o2


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
        alias = KNOWN_ALIASES.get(str(r["pair_id"]))
        disp = display_name(n1, n2, alias)
        org = combine_organism(o1, o2)
        rows.append(dict(pair_id=r["pair_id"], display_name=disp, organism=org,
                         alias=alias, pdbid1=p1, chain1=c1, name1=n1, organism1=o1,
                         pdbid2=p2, chain2=c2, name2=n2, organism2=o2))
        print(f"  {r['pair_id']}: {disp}  [{org or '?'}]")

    out = pd.DataFrame(rows)
    Path(os.path.dirname(args.out)).mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    n_named = int(out["name1"].notna().sum())
    print(f"\n[ok] wrote {args.out}: {len(out)} pairs, {n_named} with a fold-1 name")


if __name__ == "__main__":
    main()
