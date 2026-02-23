#!/usr/bin/env python3
"""Generate _seed_both.a3m for all pairs that don't have one yet."""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils.msa_utils import build_pair_seed_a3m_from_pair
from config import DATA_DIR

pairs = [d for d in sorted(os.listdir(DATA_DIR)) if re.match(r"\w{5}_\w{5}", d)]
ok = skip = fail = 0

for p in pairs:
    out = os.path.join(DATA_DIR, p, "_seed_both.a3m")
    if os.path.isfile(out):
        print(f"SKIP {p}")
        skip += 1
        continue
    try:
        build_pair_seed_a3m_from_pair(p)
        print(f"OK   {p}")
        ok += 1
    except Exception as e:
        print(f"FAIL {p}: {e}")
        fail += 1

print(f"\nDone: {ok} created, {skip} skipped, {fail} failed")
