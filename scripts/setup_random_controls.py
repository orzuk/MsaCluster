#!/usr/bin/env python3
"""
Select ~100 random single-chain proteins as broad negative controls.

Queries the RCSB PDB Search API for high-resolution, single-chain X-ray
structures (resolution < 2.5 A, 50-500 AA, one protein entity).  Randomly
samples N entries (default 100, seed=42) after excluding all PDB IDs present
in the fold-switching dataset.

Creates directories under Pipeline/Controls/RandomSingleFold/<pdb_chain>/
with FASTA files, and saves data/random_singlefold_controls.csv.

Usage:
    python scripts/setup_random_controls.py                # default 100
    python scripts/setup_random_controls.py --n 50         # custom count
    python scripts/setup_random_controls.py --dry-run      # no dirs created
    python scripts/setup_random_controls.py --seed 123     # custom seed
"""

import argparse
import json
import os
import random
import sys
import time
from pathlib import Path

import pandas as pd
import requests

# ---------------------------------------------------------------------------
# Project imports
# ---------------------------------------------------------------------------
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT)

from config import CONTROLS_BASE_DIR

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
PAIRS_FILE = os.path.join(ROOT, "data", "foldswitch_PDB_IDs_full.txt")
OUTPUT_CSV = os.path.join(ROOT, "data", "random_singlefold_controls.csv")
CONTROL_TYPE = "RandomSingleFold"
CONTROL_DIR = os.path.join(CONTROLS_BASE_DIR, CONTROL_TYPE)
CACHE_DIR = os.path.join(ROOT, "data", "cache")

RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_DATA_URL = "https://data.rcsb.org/rest/v1/core"

API_DELAY = 0.15  # seconds between RCSB data API calls

# Query thresholds
MAX_RESOLUTION = 2.5
MIN_SEQ_LEN = 50
MAX_SEQ_LEN = 500

DEFAULT_N = 100
DEFAULT_SEED = 42


# ---------------------------------------------------------------------------
# Caching
# ---------------------------------------------------------------------------

def _cache_path(category: str, key: str) -> str:
    d = os.path.join(CACHE_DIR, category)
    os.makedirs(d, exist_ok=True)
    safe_key = key.replace("/", "_").replace(":", "_")
    return os.path.join(d, f"{safe_key}.json")


def _cached_get(url: str, category: str, key: str) -> dict | None:
    """GET with JSON file cache."""
    cp = _cache_path(category, key)
    if os.path.isfile(cp):
        with open(cp) as f:
            return json.load(f)

    time.sleep(API_DELAY)
    try:
        resp = requests.get(url, timeout=30)
        if resp.status_code == 404:
            with open(cp, "w") as f:
                json.dump({"_error": "404"}, f)
            return None
        resp.raise_for_status()
        data = resp.json()
        with open(cp, "w") as f:
            json.dump(data, f)
        return data
    except Exception as e:
        print(f"  [API ERROR] {url}: {e}", flush=True)
        return None


# ---------------------------------------------------------------------------
# RCSB Search API
# ---------------------------------------------------------------------------

def query_rcsb_single_chain_proteins() -> list[str]:
    """Query RCSB for single-chain, high-resolution X-ray protein structures.

    Returns list of 'PDBID.chain' identifiers (e.g. '1UBQ.A').
    """
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "exptl.method",
                        "operator": "exact_match",
                        "value": "X-RAY DIFFRACTION",
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.resolution_combined",
                        "operator": "less",
                        "value": MAX_RESOLUTION,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.polymer_entity_count_protein",
                        "operator": "equals",
                        "value": 1,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "entity_poly.rcsb_sample_sequence_length",
                        "operator": "range",
                        "value": {"from": MIN_SEQ_LEN, "to": MAX_SEQ_LEN},
                    },
                },
            ],
        },
        "return_type": "polymer_instance",
        "request_options": {
            "paginate": {"start": 0, "rows": 10000},
            "sort": [{"sort_by": "score", "direction": "desc"}],
        },
    }

    # Check cache first (large query result)
    cp = _cache_path("rcsb_search", "single_chain_xray")
    if os.path.isfile(cp):
        with open(cp) as f:
            return json.load(f)

    print("Querying RCSB PDB Search API...", flush=True)
    all_ids = []
    start = 0
    page_size = 10000

    while True:
        query["request_options"]["paginate"] = {"start": start, "rows": page_size}
        resp = requests.post(RCSB_SEARCH_URL, json=query, timeout=60)
        resp.raise_for_status()
        data = resp.json()

        total = data.get("total_count", 0)
        results = data.get("result_set", [])
        for r in results:
            all_ids.append(r["identifier"])

        print(f"  Fetched {len(all_ids)}/{total} hits", flush=True)

        if len(all_ids) >= total or not results:
            break
        start += page_size
        time.sleep(0.5)

    # Cache the full list
    with open(cp, "w") as f:
        json.dump(all_ids, f)

    print(f"Total candidates: {len(all_ids)}", flush=True)
    return all_ids


# ---------------------------------------------------------------------------
# RCSB Data API
# ---------------------------------------------------------------------------

def fetch_protein_info(pdb_id: str) -> dict | None:
    """Fetch sequence and resolution for a single-entity PDB entry.

    For single-chain proteins, entity_id is always 1.
    """
    # Get sequence from polymer_entity
    entity_url = f"{RCSB_DATA_URL}/polymer_entity/{pdb_id}/1"
    entity_data = _cached_get(entity_url, "rcsb_entity", f"{pdb_id}_1")
    if not entity_data or "_error" in entity_data:
        return None

    ep = entity_data.get("entity_poly", {})
    sequence = ep.get("pdbx_seq_one_letter_code_can", "")
    if not sequence:
        return None

    # Get resolution from entry
    entry_url = f"{RCSB_DATA_URL}/entry/{pdb_id}"
    entry_data = _cached_get(entry_url, "rcsb_entry", pdb_id)
    resolution = None
    if entry_data and "_error" not in entry_data:
        res_list = entry_data.get("rcsb_entry_info", {}).get(
            "resolution_combined", [])
        if res_list:
            resolution = res_list[0]

    return {
        "sequence": sequence,
        "length": len(sequence),
        "resolution": resolution,
    }


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_exclude_pdbs() -> set[str]:
    """Load PDB IDs from fold-switching pairs to exclude."""
    exclude = set()
    if not os.path.isfile(PAIRS_FILE):
        return exclude
    with open(PAIRS_FILE) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                exclude.add(parts[0][:4].upper())
                exclude.add(parts[1][:4].upper())
    return exclude


def create_control_dir(pdb_chain: str, sequence: str,
                       resolution: float | None) -> str:
    """Create Pipeline/Controls/RandomSingleFold/{pdb_chain}/ with FASTA."""
    pdir = os.path.join(CONTROL_DIR, pdb_chain)
    os.makedirs(pdir, exist_ok=True)

    fasta_path = os.path.join(pdir, f"{pdb_chain}.fasta")
    pdb_id = pdb_chain[:4]
    chain = pdb_chain[4:]
    header = f">{pdb_id}:{chain}"
    if resolution is not None:
        header += f" resolution={resolution:.2f}A"
    with open(fasta_path, "w") as fh:
        fh.write(f"{header}\n{sequence}\n")

    return pdir


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Set up random single-fold controls from RCSB PDB")
    parser.add_argument("--n", type=int, default=DEFAULT_N,
                        help=f"Number of proteins to sample (default {DEFAULT_N})")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED,
                        help=f"Random seed (default {DEFAULT_SEED})")
    parser.add_argument("--dry-run", action="store_true",
                        help="Query API but don't create directories")
    args = parser.parse_args()

    # 1. Query RCSB for candidates
    all_hits = query_rcsb_single_chain_proteins()

    # 2. Exclude fold-switching PDB IDs
    exclude = load_exclude_pdbs()
    print(f"Excluding {len(exclude)} fold-switching PDB IDs", flush=True)

    filtered = []
    for hit in all_hits:
        # hit format: "1UBQ.A"
        pdb_id = hit.split(".")[0].upper()
        if pdb_id not in exclude:
            filtered.append(hit)

    print(f"After exclusion: {len(filtered)} candidates", flush=True)

    # 3. Random sample
    random.seed(args.seed)
    n = min(args.n, len(filtered))
    sampled = random.sample(filtered, n)
    print(f"Sampled {n} proteins (seed={args.seed})", flush=True)

    # 4. Fetch sequences and create directories
    results = []
    failed = 0

    for i, hit in enumerate(sampled):
        parts = hit.split(".")
        pdb_id = parts[0].upper()
        chain = parts[1] if len(parts) > 1 else "A"
        pdb_chain = f"{pdb_id}{chain}"  # e.g. "1UBQA"

        print(f"  [{i+1}/{n}] {pdb_chain}", end="", flush=True)

        info = fetch_protein_info(pdb_id)
        if not info:
            print(" - FAILED (no data)", flush=True)
            failed += 1
            continue

        if not args.dry_run:
            create_control_dir(pdb_chain, info["sequence"], info["resolution"])

        results.append({
            "pdb_chain": pdb_chain,
            "pdb_id": pdb_id,
            "chain": chain,
            "length": info["length"],
            "resolution": info["resolution"],
        })
        print(f" - OK (len={info['length']}, res={info['resolution']})", flush=True)

    # 5. Save CSV
    if results and not args.dry_run:
        os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
        df = pd.DataFrame(results)
        df.to_csv(OUTPUT_CSV, index=False)
        print(f"\nSaved {len(df)} controls to {OUTPUT_CSV}", flush=True)

    # Summary
    print(f"\n{'='*60}")
    print(f"Results: {len(results)} selected, {failed} failed")
    if results:
        lengths = [r["length"] for r in results]
        print(f"Length range: {min(lengths)}-{max(lengths)} AA")
        resolutions = [r["resolution"] for r in results if r["resolution"]]
        if resolutions:
            print(f"Resolution range: {min(resolutions):.2f}-{max(resolutions):.2f} A")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
