#!/usr/bin/env python3
"""enrich_pdb_mining_candidates.py -- add per-chain metadata + pairwise seq-id.

Reads data/pdb_mining_candidates.csv and writes <input>.enriched.csv with
these additional columns:
  seq_len_1, seq_len_2            : canonical-sequence length
  resolution_1, resolution_2      : X-ray/EM resolution in Angstrom
  method_1, method_2              : experimental method
  deposition_year_1, deposition_year_2
  organism_1, organism_2          : NCBI scientific name (first source organism)
  title_1, title_2                : one-line PDB entry title
  seq_id                          : pairwise identity between chain1 and chain2

All RCSB API calls are cached by PDB ID under data/_rcsb_meta_cache/, so
re-runs are essentially instant.

Usage:
  python3 scripts/enrich_pdb_mining_candidates.py
  python3 scripts/enrich_pdb_mining_candidates.py \
      --input data/pdb_mining_candidates.strong.csv
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
    rcsb_id
    rcsb_entry_info { resolution_combined }
    exptl { method }
    pdbx_database_status { recvd_initial_deposition_date }
    struct { title }
    polymer_entities {
      rcsb_polymer_entity_container_identifiers { auth_asym_ids }
      entity_poly { pdbx_seq_one_letter_code_can }
      rcsb_entity_source_organism { ncbi_scientific_name }
    }
  }
}
"""


def fetch_entry(pdb_id, cache_dir, max_retries=3):
    """Fetch RCSB entry metadata via GraphQL, with disk cache."""
    cache_file = cache_dir / f"{pdb_id.lower()}.json"
    if cache_file.is_file():
        try:
            payload = json.loads(cache_file.read_text())
            return payload
        except json.JSONDecodeError:
            cache_file.unlink()

    query = QUERY_TMPL % pdb_id.upper()
    data = json.dumps({"query": query}).encode("utf-8")
    req = urllib.request.Request(
        GRAPHQL_URL,
        data=data,
        headers={"Content-Type": "application/json"},
    )
    for attempt in range(max_retries):
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                result = json.loads(r.read().decode("utf-8"))
            entry = (result.get("data") or {}).get("entry")
            if entry is None:
                cache_file.write_text("null")
                return None
            cache_file.write_text(json.dumps(entry))
            return entry
        except urllib.error.HTTPError as e:
            if e.code in (429, 502, 503, 504):
                time.sleep(2 * (attempt + 1))
                continue
            print(f"  HTTP {e.code} for {pdb_id}: {e.reason}", file=sys.stderr)
            return None
        except Exception as e:
            if attempt < max_retries - 1:
                time.sleep(2 * (attempt + 1))
                continue
            print(f"  fetch failed for {pdb_id}: {e}", file=sys.stderr)
            return None
    return None


def chain_info(entry, chain_id):
    """Given an entry payload and an auth chain ID, return (seq, organism)."""
    if not entry:
        return None, None
    for pe in (entry.get("polymer_entities") or []):
        cid = pe.get("rcsb_polymer_entity_container_identifiers") or {}
        auth_ids = cid.get("auth_asym_ids") or []
        if chain_id in auth_ids:
            seq = (pe.get("entity_poly") or {}).get("pdbx_seq_one_letter_code_can")
            orgs = pe.get("rcsb_entity_source_organism") or []
            organism = orgs[0].get("ncbi_scientific_name") if orgs else None
            return seq, organism
    return None, None


def entry_basic(entry):
    """Extract resolution / method / deposition_year / title."""
    if not entry:
        return None, None, None, None
    info = entry.get("rcsb_entry_info") or {}
    res_list = info.get("resolution_combined") or []
    resolution = res_list[0] if res_list else None
    exptl = entry.get("exptl") or []
    method = exptl[0].get("method") if exptl else None
    dep_date = (entry.get("pdbx_database_status") or {}).get("recvd_initial_deposition_date")
    year = int(dep_date[:4]) if dep_date else None
    title = (entry.get("struct") or {}).get("title")
    return resolution, method, year, title


def compute_seq_id(s1, s2, aligner):
    """Global alignment -> identity = matches / max(len1, len2)."""
    if not s1 or not s2 or aligner is None:
        return None
    try:
        alignments = aligner.align(s1, s2)
        aln = alignments[0]
        a1 = str(aln[0])
        a2 = str(aln[1])
        matches = sum(1 for c1, c2 in zip(a1, a2) if c1 == c2 and c1 != "-")
        return matches / max(len(s1), len(s2))
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="data/pdb_mining_candidates.csv")
    ap.add_argument("--output", default=None,
                    help="Default: <input>.enriched.csv")
    ap.add_argument("--cache-dir", default="data/_rcsb_meta_cache")
    ap.add_argument("--sleep", type=float, default=0.05,
                    help="Pause between fresh API calls (cached calls skip this)")
    args = ap.parse_args()

    if args.output is None:
        if args.input.endswith(".csv"):
            args.output = args.input[:-4] + ".enriched.csv"
        else:
            args.output = args.input + ".enriched.csv"

    cache_dir = Path(args.cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    if not HAS_BIOPYTHON:
        print("WARN: biopython not available; seq_id column will be NaN")

    df = pd.read_csv(args.input)
    print(f"[enrich] Read {len(df)} pairs from {args.input}")

    unique_pdbs = sorted(
        set(df["pdb1"].str.lower().tolist()) | set(df["pdb2"].str.lower().tolist())
    )
    print(f"[enrich] {len(unique_pdbs)} unique PDB IDs to fetch")

    entries = {}
    n_cached = n_fetched = n_failed = 0
    for i, pdb_id in enumerate(unique_pdbs):
        was_cached = (cache_dir / f"{pdb_id}.json").is_file()
        e = fetch_entry(pdb_id, cache_dir)
        if e is None:
            n_failed += 1
        elif was_cached:
            n_cached += 1
        else:
            n_fetched += 1
            time.sleep(args.sleep)
        entries[pdb_id] = e
        if (i + 1) % 100 == 0:
            print(f"  [{i+1}/{len(unique_pdbs)}] "
                  f"cached={n_cached} fetched={n_fetched} failed={n_failed}")
    print(f"[enrich] Metadata: cached={n_cached} fetched={n_fetched} failed={n_failed}")

    aligner = None
    if HAS_BIOPYTHON:
        aligner = PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 1
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -0.5

    out_rows = []
    for _, row in df.iterrows():
        e1 = entries.get(row.pdb1.lower())
        e2 = entries.get(row.pdb2.lower())

        seq1, org1 = chain_info(e1, row.chain1)
        seq2, org2 = chain_info(e2, row.chain2)
        res1, meth1, year1, title1 = entry_basic(e1)
        res2, meth2, year2, title2 = entry_basic(e2)

        new = dict(row)
        new.update({
            "seq_len_1": len(seq1) if seq1 else None,
            "seq_len_2": len(seq2) if seq2 else None,
            "resolution_1": res1,
            "resolution_2": res2,
            "method_1": meth1,
            "method_2": meth2,
            "deposition_year_1": year1,
            "deposition_year_2": year2,
            "organism_1": org1,
            "organism_2": org2,
            "title_1": title1,
            "title_2": title2,
            "seq_id": compute_seq_id(seq1, seq2, aligner),
        })
        out_rows.append(new)

    out_df = pd.DataFrame(out_rows)
    out_df.to_csv(args.output, index=False)
    print(f"[enrich] Wrote {len(out_df)} enriched rows -> {args.output}")

    print("\n[enrich] Quick summary:")
    if "seq_id" in out_df.columns:
        s = out_df["seq_id"].dropna()
        if len(s) > 0:
            print(f"  seq_id (n={len(s)}): "
                  f"min={s.min():.3f}, median={s.median():.3f}, "
                  f"mean={s.mean():.3f}, max={s.max():.3f}")
            below_95 = (s < 0.95).sum()
            print(f"  seq_id < 0.95: {below_95} / {len(s)}  "
                  f"({'should be ~0; nonzero means RCSB cluster broke' if below_95>0 else 'as expected'})")
    rs = pd.concat([out_df["resolution_1"], out_df["resolution_2"]]).dropna()
    if len(rs) > 0:
        print(f"  resolution (n={len(rs)}): "
              f"min={rs.min():.2f}, median={rs.median():.2f}, "
              f"max={rs.max():.2f} Angstrom")
    yrs = pd.concat([out_df["deposition_year_1"], out_df["deposition_year_2"]]).dropna()
    if len(yrs) > 0:
        post2021 = (yrs >= 2021).sum()
        print(f"  deposition_year: median={int(yrs.median())}, "
              f"post-2021 entries={post2021}/{len(yrs)} "
              f"(post-AF2-training-cutoff)")
    if "method_1" in out_df.columns:
        methods = pd.concat([out_df["method_1"], out_df["method_2"]]).dropna()
        print(f"  method distribution:")
        print(methods.value_counts().to_string())


if __name__ == "__main__":
    main()
