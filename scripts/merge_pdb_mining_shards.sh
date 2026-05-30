#!/bin/bash
# ============================================================================
# merge_pdb_mining_shards.sh -- concatenate sharded mining outputs.
#
# Each shard wrote data/pdb_mining_candidates.shard_XofN.csv. This script:
#   1. Concatenates all non-empty shards into data/pdb_mining_candidates.csv
#   2. Deduplicates rows by pair_id (in case the same pair was found across
#      shards via different cluster paths -- shouldn't happen, but safety)
#   3. Re-applies dedup against existing 93-pair dataset
#   4. Writes data/pdb_mining_candidates.new_only.csv
#
# Usage:
#   bash scripts/merge_pdb_mining_shards.sh         # auto-detect N
#   bash scripts/merge_pdb_mining_shards.sh 8       # explicit N=8
# ============================================================================
set -eo pipefail

# Auto-detect N from shard filenames if not given
if [[ -n "$1" ]]; then
    N="$1"
else
    N=$(ls data/pdb_mining_candidates.shard_*of*.csv 2>/dev/null \
        | head -1 \
        | sed -E 's:.*shard_[0-9]+of([0-9]+)\.csv$:\1:' \
        || echo "")
fi

if [[ -z "$N" ]]; then
    echo "[merge] No shard CSVs found at data/pdb_mining_candidates.shard_*of*.csv"
    exit 1
fi

echo "[merge] Merging $N shards..."

python3 << EOF
import pandas as pd
from pathlib import Path
import sys

n = $N
shard_dfs = []
empty_shards = []
for x in range(1, n + 1):
    p = Path(f"data/pdb_mining_candidates.shard_{x}of{n}.csv")
    if not p.is_file() or p.stat().st_size == 0:
        empty_shards.append(x)
        continue
    try:
        df = pd.read_csv(p)
        if len(df) > 0:
            shard_dfs.append(df)
            print(f"[merge] shard {x}/{n}: {len(df)} candidates")
        else:
            empty_shards.append(x)
    except Exception as e:
        print(f"[merge] WARN: shard {x}/{n} failed to read: {e}")
        empty_shards.append(x)

if empty_shards:
    print(f"[merge] Empty shards: {empty_shards}")

if not shard_dfs:
    print("[merge] No candidates in any shard")
    sys.exit(0)

merged = pd.concat(shard_dfs, ignore_index=True)
# Dedup by pair_id (canonicalize by sorted PDB IDs so AB == BA)
def canon(r):
    a, b = f"{r.pdb1}{r.chain1}", f"{r.pdb2}{r.chain2}"
    return tuple(sorted([a, b]))
merged["_canon"] = merged.apply(canon, axis=1)
merged = merged.drop_duplicates(subset=["_canon"]).drop(columns=["_canon"])
print(f"[merge] {len(merged)} unique candidates after dedup")

merged.to_csv("data/pdb_mining_candidates.csv", index=False)
new_only = merged[~merged.get("is_existing", False).astype(bool)] \
    if "is_existing" in merged.columns else merged
new_only.to_csv("data/pdb_mining_candidates.new_only.csv", index=False)
print(f"[merge] wrote data/pdb_mining_candidates.csv ({len(merged)} rows)")
print(f"[merge] wrote data/pdb_mining_candidates.new_only.csv ({len(new_only)} rows)")
EOF

echo
echo "[merge] Done. View results:"
echo "  wc -l data/pdb_mining_candidates*.csv"
echo "  head -10 data/pdb_mining_candidates.new_only.csv"
