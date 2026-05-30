#!/usr/bin/env python3
"""compute_ss_switch_for_pairs.py -- DSSP-based SS-switch fraction per pair.

For each pair in data/foldswitch_PDB_IDs_full.txt (or a custom pair
list), download both chains' PDB files and compute the Porter-style
secondary-structure switch fraction: fraction of aligned positions
where DSSP class transitions between H (helix) and E (strand).

The actual computation reuses compute_ss_switch_fraction from
scripts/pdb_mine_foldswitch_candidates.py.

Writes docs/ss_switch_per_pair.csv (columns: pair_id, ss_switch_fraction).
This is then auto-merged into protein_comparison_table.html by
TableResults/gen_html_table.py as the SS_SWITCH% column.

Usage:
  python3 scripts/compute_ss_switch_for_pairs.py
  python3 scripts/compute_ss_switch_for_pairs.py \\
      --pairs data/my_pairs.txt --output docs/my_ss_switch.csv
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

# Reuse the mining script's DSSP-based SS-switch implementation.
sys.path.insert(0, str(Path(__file__).parent))
from pdb_mine_foldswitch_candidates import (
    compute_ss_switch_fraction,
    download_pdb_file,
)


def parse_pair_file(path):
    pairs = []
    for ln in Path(path).read_text().splitlines():
        toks = ln.strip().split()
        if len(toks) >= 2 and len(toks[0]) >= 5 and len(toks[1]) >= 5:
            a, b = toks[0], toks[1]
            pairs.append((a[:4].lower(), a[4].upper(),
                          b[:4].lower(), b[4].upper()))
    return pairs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", default="data/foldswitch_PDB_IDs_full.txt")
    ap.add_argument("--output", default="docs/ss_switch_per_pair.csv")
    ap.add_argument("--cache-dir", default="data/pdb_cache")
    args = ap.parse_args()

    pdb_cache = Path(args.cache_dir)
    pdb_cache.mkdir(parents=True, exist_ok=True)

    pairs = parse_pair_file(args.pairs)
    print(f"Loaded {len(pairs)} pairs from {args.pairs}")

    rows = []
    n_ok = n_fail = 0
    for i, (p1, c1, p2, c2) in enumerate(pairs, 1):
        pair_id = f"{p1}{c1}_{p2}{c2}"
        pdb1 = download_pdb_file(p1, pdb_cache)
        pdb2 = download_pdb_file(p2, pdb_cache)
        if pdb1 is None or pdb2 is None:
            rows.append({"pair_id": pair_id, "ss_switch_fraction": None,
                         "error": "pdb_download_failed"})
            n_fail += 1
        else:
            try:
                v = compute_ss_switch_fraction(pdb1, c1, pdb2, c2)
            except Exception as e:
                rows.append({"pair_id": pair_id, "ss_switch_fraction": None,
                             "error": f"compute_failed: {e}"})
                n_fail += 1
            else:
                if v is None:
                    rows.append({"pair_id": pair_id, "ss_switch_fraction": None,
                                 "error": "dssp_returned_none"})
                    n_fail += 1
                else:
                    rows.append({"pair_id": pair_id, "ss_switch_fraction": v,
                                 "error": ""})
                    n_ok += 1
        if i % 10 == 0:
            print(f"  [{i}/{len(pairs)}]  ok={n_ok}  fail={n_fail}")

    df = pd.DataFrame(rows)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)
    print(f"\nWrote {args.output}  (ok={n_ok}, fail={n_fail})")

    v = df["ss_switch_fraction"].dropna()
    if len(v):
        print("\nDistribution:")
        print(f"  min={v.min():.3f}, 5%={v.quantile(0.05):.3f}, "
              f"25%={v.quantile(0.25):.3f}, median={v.median():.3f}, "
              f"75%={v.quantile(0.75):.3f}, 95%={v.quantile(0.95):.3f}, "
              f"max={v.max():.3f}")
        for thr in [0.05, 0.10, 0.15, 0.20, 0.30]:
            print(f"  ss_switch >= {thr}: {(v >= thr).sum()} / {len(v)}")


if __name__ == "__main__":
    main()
