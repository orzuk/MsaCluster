#!/usr/bin/env python3
"""diagnose_existing_93_seq_id.py -- check whether existing-93 chain pairs
share enough sequence identity to be in the same 95% RCSB cluster.

Reads data/foldswitch_PDB_IDs_full.txt, fetches canonical sequences for
each chain pair via the RCSB GraphQL API (reusing the same cache as
enrich_pdb_mining_candidates.py), computes pairwise identity, and
reports the distribution.

If most pairs are <95% identical, the mining script's
--min-seq-identity 0.95 is the wrong setting -- it will systematically
miss most real fold-switchers. We'd need to rerun at 70% or 50%.

Usage:
  python3 scripts/diagnose_existing_93_seq_id.py
"""

import argparse
import json
import sys
import time
from pathlib import Path
import urllib.request
import urllib.error

import pandas as pd

try:
    from Bio.Align import PairwiseAligner
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False


GRAPHQL_URL = "https://data.rcsb.org/graphql"

QUERY_TMPL = """
{
  entry(entry_id: "%s") {
    polymer_entities {
      rcsb_polymer_entity_container_identifiers { auth_asym_ids }
      entity_poly { pdbx_seq_one_letter_code_can }
    }
  }
}
"""


def fetch_entry(pdb_id, cache_dir, max_retries=3):
    cache_file = cache_dir / f"{pdb_id.lower()}.json"
    if cache_file.is_file():
        try:
            return json.loads(cache_file.read_text())
        except json.JSONDecodeError:
            cache_file.unlink()
    query = QUERY_TMPL % pdb_id.upper()
    data = json.dumps({"query": query}).encode("utf-8")
    req = urllib.request.Request(
        GRAPHQL_URL, data=data,
        headers={"Content-Type": "application/json"},
    )
    for attempt in range(max_retries):
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                result = json.loads(r.read().decode("utf-8"))
            entry = (result.get("data") or {}).get("entry")
            cache_file.write_text(json.dumps(entry) if entry else "null")
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


def get_chain_seq(entry, chain_id):
    if not entry:
        return None
    for pe in (entry.get("polymer_entities") or []):
        cid = pe.get("rcsb_polymer_entity_container_identifiers") or {}
        if chain_id in (cid.get("auth_asym_ids") or []):
            return (pe.get("entity_poly") or {}).get("pdbx_seq_one_letter_code_can")
    return None


def compute_seq_id(s1, s2, aligner):
    if not s1 or not s2 or aligner is None:
        return None
    try:
        aln = aligner.align(s1, s2)[0]
        a1 = str(aln[0])
        a2 = str(aln[1])
        matches = sum(1 for c1, c2 in zip(a1, a2) if c1 == c2 and c1 != "-")
        return matches / max(len(s1), len(s2))
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--existing", default="data/foldswitch_PDB_IDs_full.txt")
    ap.add_argument("--cache-dir", default="data/_rcsb_meta_cache")
    ap.add_argument("--output", default="data/existing_93_seq_id.csv")
    args = ap.parse_args()

    if not HAS_BIOPYTHON:
        print("ERROR: biopython is required.", file=sys.stderr)
        sys.exit(1)

    cache_dir = Path(args.cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    pairs = []
    for ln in Path(args.existing).read_text().splitlines():
        toks = ln.strip().split()
        if len(toks) < 2:
            continue
        a, b = toks[0], toks[1]
        if len(a) >= 5 and len(b) >= 5:
            pairs.append((a[:4].lower(), a[4].upper(), b[:4].lower(), b[4].upper()))
    print(f"Loaded {len(pairs)} pairs from {args.existing}")

    rows = []
    for i, (pdb1, c1, pdb2, c2) in enumerate(pairs, 1):
        e1 = fetch_entry(pdb1, cache_dir)
        e2 = fetch_entry(pdb2, cache_dir)
        s1 = get_chain_seq(e1, c1) if e1 else None
        s2 = get_chain_seq(e2, c2) if e2 else None
        sid = compute_seq_id(s1, s2, aligner)
        rows.append({
            "pair_id": f"{pdb1}{c1}_{pdb2}{c2}",
            "pdb1": pdb1, "chain1": c1, "pdb2": pdb2, "chain2": c2,
            "seq_len_1": len(s1) if s1 else None,
            "seq_len_2": len(s2) if s2 else None,
            "seq_id": sid,
        })
        if i % 20 == 0:
            print(f"  [{i}/{len(pairs)}]")

    df = pd.DataFrame(rows)
    df.to_csv(args.output, index=False)
    print(f"\nWrote {args.output}")

    sid = df["seq_id"].dropna()
    print(f"\n=== Pairwise seq-id distribution (n={len(sid)}) ===")
    print(f"  min:    {sid.min():.3f}")
    print(f"  5%:     {sid.quantile(0.05):.3f}")
    print(f"  25%:    {sid.quantile(0.25):.3f}")
    print(f"  median: {sid.median():.3f}")
    print(f"  75%:    {sid.quantile(0.75):.3f}")
    print(f"  95%:    {sid.quantile(0.95):.3f}")
    print(f"  max:    {sid.max():.3f}")
    print()
    print(f"  pairs >= 95% identity: {(sid >= 0.95).sum()} / {len(sid)}")
    print(f"  pairs >= 90% identity: {(sid >= 0.90).sum()} / {len(sid)}")
    print(f"  pairs >= 70% identity: {(sid >= 0.70).sum()} / {len(sid)}")
    print(f"  pairs >= 50% identity: {(sid >= 0.50).sum()} / {len(sid)}")
    print(f"  pairs >= 30% identity: {(sid >= 0.30).sum()} / {len(sid)}")
    print(f"  pairs <  30% identity: {(sid < 0.30).sum()} / {len(sid)}")
    print()
    print("DIAGNOSIS:")
    pct95 = (sid >= 0.95).sum() / len(sid)
    if pct95 < 0.20:
        print(f"  *** {(1-pct95)*100:.0f}% of existing-93 pairs are NOT in the same 95% RCSB cluster.")
        print(f"  *** Mining at --min-seq-identity 0.95 systematically misses real fold-switchers.")
        # Suggest the right threshold
        for thr in (0.90, 0.70, 0.50, 0.30):
            covered = (sid >= thr).sum() / len(sid)
            if covered >= 0.90:
                print(f"  -> At {int(thr*100)}% identity threshold, {covered*100:.0f}% of existing pairs would be in the same cluster.")
                print(f"  -> RECOMMEND rerun mining with --min-seq-identity {thr}")
                break
    else:
        print(f"  {pct95*100:.0f}% of existing-93 are >=95% identical -- mining threshold OK.")
        print(f"  The miss must be due to other filters (NMR? resolution? length?).")

    print()
    print("Bottom 20 pairs by seq_id (most divergent in existing-93):")
    print(df.dropna(subset=["seq_id"]).nsmallest(20, "seq_id")[
        ["pair_id", "seq_len_1", "seq_len_2", "seq_id"]
    ].to_string(index=False))


if __name__ == "__main__":
    main()
