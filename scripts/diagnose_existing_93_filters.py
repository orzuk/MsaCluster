#!/usr/bin/env python3
"""diagnose_existing_93_filters.py -- for each existing-93 pair, identify
which mining filter would have killed it.

For each pair (chain1, chain2), fetch RCSB metadata and check:
  - experimental method (X-ray / NMR / EM)
  - resolution
  - chain length
  - presence in the cached 95% cluster file (would the cluster pass size>=2?)

Reports a table: pair_id | method_1 | res_1 | len_1 | method_2 | res_2 |
len_2 | killed_by

Aggregates which filter(s) account for the most lost pairs, so we know
which defaults to relax for a re-mine that recovers Porter's set.

Usage:
  python3 scripts/diagnose_existing_93_filters.py
"""

import argparse
import json
import sys
import time
from collections import Counter
from pathlib import Path
import urllib.request
import urllib.error

import pandas as pd


GRAPHQL_URL = "https://data.rcsb.org/graphql"

QUERY = """
{
  entry(entry_id: "%s") {
    rcsb_entry_info { resolution_combined }
    exptl { method }
    polymer_entities {
      rcsb_polymer_entity_container_identifiers { auth_asym_ids }
      entity_poly { pdbx_seq_one_letter_code_can }
    }
  }
}
"""


def fetch(pdb_id, cache_dir, max_retries=3):
    cf = cache_dir / f"{pdb_id.lower()}.json"
    if cf.is_file():
        try:
            payload = json.loads(cf.read_text())
            if payload == "null":
                return None
            return payload
        except Exception:
            cf.unlink()
    q = QUERY % pdb_id.upper()
    data = json.dumps({"query": q}).encode()
    req = urllib.request.Request(GRAPHQL_URL, data=data,
        headers={"Content-Type": "application/json"})
    for attempt in range(max_retries):
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                result = json.loads(r.read().decode())
            entry = (result.get("data") or {}).get("entry")
            cf.write_text(json.dumps(entry) if entry else "null")
            return entry
        except urllib.error.HTTPError as e:
            if e.code in (429, 502, 503, 504):
                time.sleep(2 * (attempt + 1))
                continue
            return None
        except Exception:
            if attempt < max_retries - 1:
                time.sleep(2 * (attempt + 1))
                continue
            return None
    return None


def chain_info(entry, chain_id):
    if not entry:
        return None
    for pe in (entry.get("polymer_entities") or []):
        cid = pe.get("rcsb_polymer_entity_container_identifiers") or {}
        if chain_id in (cid.get("auth_asym_ids") or []):
            seq = (pe.get("entity_poly") or {}).get("pdbx_seq_one_letter_code_can")
            return {"seq_len": len(seq) if seq else None}
    return None


def entry_meta(entry):
    if not entry:
        return None, None
    info = entry.get("rcsb_entry_info") or {}
    res_list = info.get("resolution_combined") or []
    res = res_list[0] if res_list else None
    exptl = entry.get("exptl") or []
    method = exptl[0].get("method") if exptl else None
    return method, res


def diagnose_filter(method, res, length,
                    max_resolution=3.0, min_length=50, max_length=500,
                    allow_nmr=False):
    """Return list of filter names that would kill this chain."""
    killed = []
    if method is None:
        killed.append("metadata_missing")
        return killed
    m_upper = method.upper()
    is_nmr = "NMR" in m_upper or "NUCLEAR" in m_upper
    is_em = ("ELECTRON" in m_upper) or m_upper == "EM"

    if is_nmr and not allow_nmr:
        killed.append("nmr_excluded")
    # Resolution check (only meaningful for X-ray / EM)
    if not is_nmr and res is not None and res > max_resolution:
        killed.append(f"resolution>{max_resolution}")
    # EM with missing resolution would be killed by the same check if res is None
    if is_em and res is None:
        killed.append("em_no_resolution")
    if length is None:
        killed.append("length_unknown")
    elif length < min_length:
        killed.append("length<min")
    elif length > max_length:
        killed.append("length>max")
    return killed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--existing", default="data/foldswitch_PDB_IDs_full.txt")
    ap.add_argument("--cache-dir", default="data/_rcsb_meta_cache")
    ap.add_argument("--output", default="data/existing_93_filter_diagnosis.csv")
    args = ap.parse_args()

    cache_dir = Path(args.cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    pairs = []
    for ln in Path(args.existing).read_text().splitlines():
        toks = ln.strip().split()
        if len(toks) >= 2 and len(toks[0]) >= 5 and len(toks[1]) >= 5:
            a, b = toks[0], toks[1]
            pairs.append((a[:4].lower(), a[4].upper(),
                          b[:4].lower(), b[4].upper()))
    print(f"Loaded {len(pairs)} pairs")

    rows = []
    for i, (p1, c1, p2, c2) in enumerate(pairs, 1):
        e1 = fetch(p1, cache_dir)
        e2 = fetch(p2, cache_dir)
        m1, r1 = entry_meta(e1)
        m2, r2 = entry_meta(e2)
        ci1 = chain_info(e1, c1)
        ci2 = chain_info(e2, c2)
        l1 = ci1["seq_len"] if ci1 else None
        l2 = ci2["seq_len"] if ci2 else None
        kill1 = diagnose_filter(m1, r1, l1)
        kill2 = diagnose_filter(m2, r2, l2)

        pair_killed = []
        if kill1:
            pair_killed += [f"chain1:{k}" for k in kill1]
        if kill2:
            pair_killed += [f"chain2:{k}" for k in kill2]

        rows.append({
            "pair_id": f"{p1}{c1}_{p2}{c2}",
            "method_1": m1, "res_1": r1, "len_1": l1,
            "method_2": m2, "res_2": r2, "len_2": l2,
            "killed_by": "; ".join(pair_killed) if pair_killed else "passes_filters",
            "pair_passes": len(pair_killed) == 0,
        })
        if i % 20 == 0:
            print(f"  [{i}/{len(pairs)}]")

    df = pd.DataFrame(rows)
    df.to_csv(args.output, index=False)
    print(f"\nWrote {args.output}")

    print("\n=== Pair-level summary ===")
    n = len(df)
    passing = int(df["pair_passes"].sum())
    print(f"  Pairs passing all default filters: {passing} / {n}")
    print(f"  Pairs killed by at least one filter: {n - passing} / {n}")

    print("\n=== Method distribution (chain-level, n=186 chains) ===")
    methods = pd.concat([df["method_1"], df["method_2"]]).dropna()
    print(methods.value_counts().to_string())

    print("\n=== Why pairs were killed (multi-kill counts each cause) ===")
    counter = Counter()
    for row in df.itertuples():
        if not row.pair_passes and row.killed_by:
            for cause in row.killed_by.split("; "):
                # strip chain prefix
                bare = cause.split(":", 1)[1] if ":" in cause else cause
                counter[bare] += 1
    for cause, count in counter.most_common():
        print(f"  {cause}: {count} chain-level hits")

    # Targeted: would relaxing each filter rescue many pairs?
    print("\n=== Hypothetical re-mine: how many pairs survive if we relax... ===")
    def recheck(row, allow_nmr=False, max_res=3.0, max_len=500):
        chains = [
            (row.method_1, row.res_1, row.len_1),
            (row.method_2, row.res_2, row.len_2),
        ]
        for m, r, l in chains:
            killed = diagnose_filter(m, r, l,
                max_resolution=max_res, max_length=max_len, allow_nmr=allow_nmr)
            if killed:
                return False
        return True

    scenarios = [
        ("default (NMR excluded, res<=3.0)",       dict(allow_nmr=False, max_res=3.0, max_len=500)),
        ("+ allow NMR",                            dict(allow_nmr=True,  max_res=3.0, max_len=500)),
        ("+ allow NMR + res<=4.0",                 dict(allow_nmr=True,  max_res=4.0, max_len=500)),
        ("+ allow NMR + res<=4.0 + len<=800",      dict(allow_nmr=True,  max_res=4.0, max_len=800)),
        ("+ allow NMR + res<=10.0 + len<=1000",    dict(allow_nmr=True,  max_res=10.0, max_len=1000)),
    ]
    for name, params in scenarios:
        n_passing = sum(recheck(r, **params) for r in df.itertuples())
        print(f"  {name}: {n_passing} / {n}")

    print("\nTop 10 killed pairs (sample, to spot patterns):")
    print(df[~df["pair_passes"]][
        ["pair_id", "method_1", "res_1", "len_1",
         "method_2", "res_2", "len_2", "killed_by"]
    ].head(10).to_string(index=False))


if __name__ == "__main__":
    main()
