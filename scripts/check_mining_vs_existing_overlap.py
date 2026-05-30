#!/usr/bin/env python3
"""check_mining_vs_existing_overlap.py -- diagnose overlap between
mining candidates and the existing 93-pair set.

The mining script's --existing-pairs filter is an EXACT canonical-pair-ID
match: {pdb1}{chain1}_{pdb2}{chain2} in either ordering. That's narrow.
If the per-cluster dedup picks a DIFFERENT chain from the same RCSB
cluster as one of our existing-93 pairs, the candidate gets tagged
is_existing=False even though it's biologically the same protein.

This script computes overlap at three resolutions:
  1. EXACT pair-ID match (same chains in either order)
  2. PDB-pair match (same two PDB IDs, any chain choice)
  3. Single-PDB overlap (mining pair shares at least one PDB with existing)

And reports which existing-93 PDB IDs the mining never even saw, so we
can diagnose why (resolution filter? length filter? NMR-only?).

Usage:
  python3 scripts/check_mining_vs_existing_overlap.py
  python3 scripts/check_mining_vs_existing_overlap.py \\
      --mining data/pdb_mining_candidates.csv \\
      --existing data/foldswitch_PDB_IDs_full.txt
"""

import argparse
from pathlib import Path
import pandas as pd


def load_existing_pairs(path):
    """Return (exact_pair_ids, pdb_pairs, all_pdbs)."""
    exact = set()
    pdb_pairs = set()
    all_pdbs = set()
    for ln in Path(path).read_text().splitlines():
        toks = ln.strip().split()
        if len(toks) < 2:
            continue
        a, b = toks[0], toks[1]
        if len(a) < 5 or len(b) < 5:
            continue
        pa = f"{a[:4].lower()}{a[4].upper()}"
        pb = f"{b[:4].lower()}{b[4].upper()}"
        exact.add(f"{pa}_{pb}")
        exact.add(f"{pb}_{pa}")
        pdb_a, pdb_b = a[:4].lower(), b[:4].lower()
        pdb_pairs.add(tuple(sorted([pdb_a, pdb_b])))
        all_pdbs.add(pdb_a)
        all_pdbs.add(pdb_b)
    return exact, pdb_pairs, all_pdbs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mining", default="data/pdb_mining_candidates.csv")
    ap.add_argument("--existing", default="data/foldswitch_PDB_IDs_full.txt")
    args = ap.parse_args()

    exact_pairs, pdb_pairs, ex_pdbs = load_existing_pairs(args.existing)
    n_existing = len(exact_pairs) // 2  # each pair stored twice (both orders)
    print(f"=== Existing reference set ({args.existing}) ===")
    print(f"  pairs:           {n_existing}")
    print(f"  PDB-pair tuples: {len(pdb_pairs)}")
    print(f"  unique PDB IDs:  {len(ex_pdbs)}")
    print()

    df = pd.read_csv(args.mining)
    mine_pdbs = set(df["pdb1"].str.lower()) | set(df["pdb2"].str.lower())
    print(f"=== Mining candidates ({args.mining}) ===")
    print(f"  pairs:          {len(df)}")
    print(f"  unique PDB IDs: {len(mine_pdbs)}")
    print()

    # 1. Exact pair-ID match
    exact_hits = df[df["pair_id"].isin(exact_pairs)]
    print(f"=== Overlap @ Level 1: EXACT pair-ID (same chains, either order) ===")
    print(f"  hits: {len(exact_hits)} / {len(df)} mining pairs")
    print(f"  (this should equal what is_existing=True gives)")
    is_existing_count = int(df["is_existing"].sum()) if "is_existing" in df.columns else "n/a"
    print(f"  is_existing column True count: {is_existing_count}")
    print()

    # 2. PDB-pair match (same two PDB IDs, any chain choice)
    def pdb_pair_key(r):
        return tuple(sorted([str(r["pdb1"]).lower(), str(r["pdb2"]).lower()]))
    df_pdb_pair = df.apply(pdb_pair_key, axis=1)
    pdb_pair_hits = df[df_pdb_pair.isin(pdb_pairs)]
    print(f"=== Overlap @ Level 2: PDB-pair (same two PDB IDs, ANY chain) ===")
    print(f"  hits: {len(pdb_pair_hits)} / {len(df)} mining pairs")
    print(f"  -> these are biologically the same proteins as existing-93 entries,")
    print(f"     but the mining picked different chain choices")
    if len(pdb_pair_hits):
        print()
        print(f"  Sample of PDB-pair-overlap hits (NOT in exact-match):")
        not_exact = pdb_pair_hits[~pdb_pair_hits["pair_id"].isin(exact_pairs)]
        cols = ["pair_id", "tm_max", "ss_switch_fraction", "cluster_size"]
        cols = [c for c in cols if c in not_exact.columns]
        print(not_exact[cols].head(15).to_string(index=False))
    print()

    # 3. Single-PDB overlap
    single_pdb_hits = df[
        df["pdb1"].str.lower().isin(ex_pdbs) | df["pdb2"].str.lower().isin(ex_pdbs)
    ]
    print(f"=== Overlap @ Level 3: SINGLE-PDB (mining hit shares >=1 PDB with existing) ===")
    print(f"  hits: {len(single_pdb_hits)} / {len(df)} mining pairs")
    print(f"  -> these mining pairs partially overlap the existing set on one side")
    print()

    # 4. Which existing PDB IDs are NOT in mining at all?
    missing = ex_pdbs - mine_pdbs
    print(f"=== Existing-93 PDB IDs INVISIBLE to mining: {len(missing)} / {len(ex_pdbs)} ===")
    if missing:
        print(f"  (likely filtered by --max-resolution / length / NMR-exclusion)")
        print(f"  First 20: {sorted(missing)[:20]}")
        print()
        # Identify which existing PAIRS contain only invisible PDB IDs
        invisible_pairs = []
        for pa, pb in pdb_pairs:
            if pa in missing or pb in missing:
                invisible_pairs.append((pa, pb))
        print(f"  Existing-93 PDB-pairs with >=1 invisible PDB: {len(invisible_pairs)}")
    print()

    # 5. Summary table
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  Existing pairs:                                {n_existing}")
    print(f"  Mining candidates:                             {len(df)}")
    print()
    print(f"  Exact pair-ID overlap (Level 1):               {len(exact_hits)}")
    print(f"  PDB-pair overlap (Level 2, the real one):      {len(pdb_pair_hits)}")
    print(f"  Single-PDB overlap (Level 3):                  {len(single_pdb_hits)}")
    print()
    print(f"  Existing-93 PDB IDs not seen by mining:        {len(missing)} / {len(ex_pdbs)}")
    print()
    print("Interpretation:")
    print(f"  - If Level 2 >> Level 1: mining is re-finding existing biology")
    print(f"    via different chain choices. The is_existing filter is too narrow.")
    print(f"  - If Level 2 ~ Level 1: mining is finding genuinely new proteins.")
    print(f"  - High 'invisible PDB IDs' count: mining missed entries due to")
    print(f"    metadata filters (resolution, length, NMR-exclusion).")


if __name__ == "__main__":
    main()
