#!/usr/bin/env python3
"""
Find matched single-fold controls for each of the 93 fold-switching pairs.

For each pair, queries PDBe SIFTS for the Pfam family, then finds a reviewed
UniProt protein from the same family that:
  - Has similar length (within ±30%)
  - Has exactly 1 PDB structure (proxy for "not a fold-switcher")
  - Shares >30% sequence identity with the fold-switching protein

Creates directories under Pipeline/Controls/MatchedSingleFold/<pdb_chain>/
with FASTA files, and saves data/matched_controls.csv.

Usage:
    python scripts/setup_matched_controls.py                  # all 93 pairs
    python scripts/setup_matched_controls.py --pair-id 1dzlA_5keqF  # one pair
    python scripts/setup_matched_controls.py --dry-run        # no dirs created
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path

import requests

# ---------------------------------------------------------------------------
# Project imports
# ---------------------------------------------------------------------------
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT)

from config import CONTROLS_BASE_DIR, DATA_DIR

# Try BioPython for sequence identity; fall back to simple alignment
try:
    from Bio.Align import PairwiseAligner
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
PAIRS_FILE = os.path.join(ROOT, "data", "foldswitch_PDB_IDs_full.txt")
OUTPUT_CSV = os.path.join(ROOT, "data", "matched_controls.csv")
CONTROL_TYPE = "MatchedSingleFold"
CONTROL_DIR = os.path.join(CONTROLS_BASE_DIR, CONTROL_TYPE)
CACHE_DIR = os.path.join(ROOT, "data", "cache")

API_DELAY = 0.35  # seconds between API calls (EBI rate limit)

# Filtering thresholds
LENGTH_TOLERANCE = 0.30  # ±30% of query length
MIN_SEQ_IDENTITY = 0.30  # 30%
MAX_PDB_ENTRIES = 1      # exactly 1 PDB structure

# InterPro page size
INTERPRO_PAGE_SIZE = 200

# ---------------------------------------------------------------------------
# Caching helpers
# ---------------------------------------------------------------------------

def _cache_path(category: str, key: str) -> str:
    d = os.path.join(CACHE_DIR, category)
    os.makedirs(d, exist_ok=True)
    safe_key = key.replace("/", "_").replace(":", "_")
    return os.path.join(d, f"{safe_key}.json")


def _cached_get(url: str, category: str, key: str,
                headers: dict | None = None) -> dict | None:
    """GET with JSON file cache.  Returns parsed JSON or None on error."""
    cp = _cache_path(category, key)
    if os.path.isfile(cp):
        with open(cp) as f:
            return json.load(f)

    time.sleep(API_DELAY)
    try:
        resp = requests.get(url, headers=headers or {}, timeout=60)
        if resp.status_code == 404:
            # Cache the 404 so we don't re-query
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
# API wrappers
# ---------------------------------------------------------------------------

def get_pfam_for_pdb_chain(pdb_id: str, chain: str) -> list[str]:
    """Query PDBe SIFTS for Pfam families of a PDB chain.
    Returns list of Pfam accessions (e.g. ['PF00500']).
    """
    pdb_lower = pdb_id.lower()
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/pfam/{pdb_lower}"
    data = _cached_get(url, "sifts_pfam", pdb_lower)
    if not data or "_error" in data:
        return []

    pfam_ids = []
    pdb_data = data.get(pdb_lower, {})
    pfam_section = pdb_data.get("Pfam", {})
    for pfam_id, info in pfam_section.items():
        mappings = info.get("mappings", [])
        for m in mappings:
            if m.get("chain_id", "").upper() == chain.upper():
                if pfam_id not in pfam_ids:
                    pfam_ids.append(pfam_id)
    return pfam_ids


def get_pfam_members(pfam_id: str) -> list[dict]:
    """Query InterPro for reviewed (Swiss-Prot) proteins in a Pfam family.
    Returns list of {accession, name, length}.
    """
    members = []
    page = 1
    url = (f"https://www.ebi.ac.uk/interpro/api/protein/reviewed/"
           f"entry/pfam/{pfam_id}/?page_size={INTERPRO_PAGE_SIZE}&format=json")

    while url:
        cache_key = f"{pfam_id}_p{page}"
        data = _cached_get(url, "interpro_members", cache_key)
        if not data or "_error" in data:
            break

        for r in data.get("results", []):
            meta = r.get("metadata", r)  # InterPro nests under metadata
            members.append({
                "accession": meta.get("accession", ""),
                "name": meta.get("name", ""),
                "length": meta.get("length", 0),
            })

        url = data.get("next")
        page += 1

    return members


def get_uniprot_details(accession: str) -> dict | None:
    """Fetch UniProt entry: sequence + PDB cross-references.
    Returns {accession, length, sequence, pdb_ids: [str]}.
    """
    url = (f"https://rest.uniprot.org/uniprotkb/{accession}"
           f"?format=json&fields=accession,length,sequence,xref_pdb")
    headers = {"Accept": "application/json"}
    data = _cached_get(url, "uniprot", accession, headers=headers)
    if not data or "_error" in data:
        return None

    seq_obj = data.get("sequence", {})
    sequence = seq_obj.get("value", "")
    length = seq_obj.get("length", len(sequence))

    pdb_ids = []
    for xref in data.get("uniProtKBCrossReferences", []):
        if xref.get("database") == "PDB":
            pdb_ids.append(xref.get("id", ""))

    return {
        "accession": accession,
        "length": length,
        "sequence": sequence,
        "pdb_ids": pdb_ids,
    }


# ---------------------------------------------------------------------------
# Sequence identity
# ---------------------------------------------------------------------------

def compute_seq_identity(seq1: str, seq2: str) -> float:
    """Compute fractional sequence identity via global alignment."""
    if not seq1 or not seq2:
        return 0.0

    if HAS_BIOPYTHON:
        aligner = PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 1
        aligner.mismatch_score = 0
        aligner.open_gap_score = 0
        aligner.extend_gap_score = 0
        score = aligner.score(seq1, seq2)
        return float(score) / max(len(seq1), len(seq2))

    # Simple k-mer based estimate (fallback)
    k = 3
    set1 = set(seq1[i:i+k] for i in range(len(seq1) - k + 1))
    set2 = set(seq2[i:i+k] for i in range(len(seq2) - k + 1))
    if not set1 or not set2:
        return 0.0
    return len(set1 & set2) / len(set1 | set2)


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def _fetch_pdb_sequence(pdb_id: str, chain: str) -> str:
    """Fetch sequence from RCSB Data API for a PDB chain."""
    pdb_lower = pdb_id.lower()
    # First get entity_id for the chain
    inst_url = (f"https://data.rcsb.org/rest/v1/core/"
                f"polymer_entity_instance/{pdb_lower}/{chain}")
    inst = _cached_get(inst_url, "rcsb_chain", f"{pdb_lower}_{chain}")
    if not inst or "_error" in inst:
        return ""
    entity_id = (inst.get("rcsb_polymer_entity_instance_container_identifiers", {})
                 .get("entity_id", "1"))
    # Then get the sequence
    ent_url = (f"https://data.rcsb.org/rest/v1/core/"
               f"polymer_entity/{pdb_lower}/{entity_id}")
    ent = _cached_get(ent_url, "rcsb_entity", f"{pdb_lower}_{entity_id}")
    if not ent or "_error" in ent:
        return ""
    return ent.get("entity_poly", {}).get("pdbx_seq_one_letter_code_can", "")


def read_fold_sequence(fold_name: str, pair_id: str) -> str:
    """Read sequence from local FASTA or fetch from RCSB if missing."""
    fasta = os.path.join(DATA_DIR, pair_id, f"{fold_name}.fasta")
    if os.path.isfile(fasta):
        seq_lines = []
        with open(fasta) as f:
            for line in f:
                line = line.strip()
                if not line.startswith(">"):
                    seq_lines.append(line)
        seq = "".join(seq_lines)
        if seq:
            return seq
    # Fallback: fetch from RCSB (works on cluster without local data)
    pdb_id, chain = fold_name[:4], fold_name[4:]
    return _fetch_pdb_sequence(pdb_id, chain)


def find_control_for_pair(fold1: str, fold2: str, exclude_pdbs: set,
                          dry_run: bool = False) -> dict | None:
    """Find the best matched single-fold control for a pair.

    Tries Pfam families from fold1 first, then fold2.
    Returns dict with control info, or None if no match found.
    """
    pair_id = f"{fold1}_{fold2}"
    pdb1, chain1 = fold1[:4], fold1[4:]
    pdb2, chain2 = fold2[:4], fold2[4:]

    # Read query sequences
    seq1 = read_fold_sequence(fold1, pair_id)
    seq2 = read_fold_sequence(fold2, pair_id)
    query_len = max(len(seq1), len(seq2)) if (seq1 or seq2) else 0
    query_seq = seq1 if len(seq1) >= len(seq2) else seq2

    if not query_seq:
        print(f"  [WARN] No FASTA found for {pair_id}", flush=True)
        return None

    # Get Pfam families from both folds
    pfam_ids = []
    for pdb, chain in [(pdb1, chain1), (pdb2, chain2)]:
        for pf in get_pfam_for_pdb_chain(pdb, chain):
            if pf not in pfam_ids:
                pfam_ids.append(pf)

    if not pfam_ids:
        print(f"  [WARN] No Pfam found for {pair_id}", flush=True)
        return None

    print(f"  Pfam families: {pfam_ids}", flush=True)

    # Length bounds
    min_len = int(query_len * (1 - LENGTH_TOLERANCE))
    max_len = int(query_len * (1 + LENGTH_TOLERANCE))

    best = None

    for pfam_id in pfam_ids:
        members = get_pfam_members(pfam_id)
        print(f"  {pfam_id}: {len(members)} reviewed members", flush=True)

        # Pre-filter by length before making expensive UniProt API calls
        length_ok = [m for m in members if min_len <= m["length"] <= max_len]
        print(f"  Length filter ({min_len}-{max_len}): {len(length_ok)} candidates", flush=True)

        for m in length_ok:
            acc = m["accession"]

            # Fetch UniProt details (sequence + PDB xrefs)
            details = get_uniprot_details(acc)
            if not details:
                continue

            pdb_ids = details["pdb_ids"]

            # Must have exactly 1 PDB structure
            if len(pdb_ids) != MAX_PDB_ENTRIES:
                continue

            # Exclude known fold-switching proteins
            pdb_lower = pdb_ids[0].lower()
            if pdb_lower in exclude_pdbs:
                continue

            # Compute sequence identity
            identity = compute_seq_identity(query_seq, details["sequence"])
            if identity < MIN_SEQ_IDENTITY:
                continue

            # Track best candidate
            if best is None or identity > best["seq_identity"]:
                best = {
                    "pair_id": pair_id,
                    "fold1": fold1,
                    "fold2": fold2,
                    "pfam_id": pfam_id,
                    "control_uniprot": acc,
                    "control_pdb": pdb_ids[0],
                    "control_name": m.get("name", ""),
                    "seq_len": details["length"],
                    "seq_identity": round(identity, 4),
                    "sequence": details["sequence"],
                }

        # If we found a good match in this Pfam family, stop
        if best and best["seq_identity"] >= MIN_SEQ_IDENTITY:
            break

    return best


# ---------------------------------------------------------------------------
# Directory creation
# ---------------------------------------------------------------------------

def create_control_dir(ctrl: dict) -> str:
    """Create Pipeline/Controls/MatchedSingleFold/<pdb_id>/ with FASTA."""
    pdb_id = ctrl["control_pdb"]
    pdir = os.path.join(CONTROL_DIR, pdb_id)
    os.makedirs(pdir, exist_ok=True)

    fasta_path = os.path.join(pdir, f"{pdb_id}.fasta")
    with open(fasta_path, "w") as f:
        header = f">{pdb_id} UniProt:{ctrl['control_uniprot']}"
        if ctrl.get("control_name"):
            header += f" {ctrl['control_name']}"
        f.write(f"{header}\n{ctrl['sequence']}\n")

    return pdir


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def _save_csv(results: list[dict]) -> None:
    """Write results to CSV (merging with any existing file)."""
    import pandas as pd
    os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
    if os.path.isfile(OUTPUT_CSV):
        try:
            df_old = pd.read_csv(OUTPUT_CSV)
            df_new = pd.DataFrame(results)
            df_all = pd.concat([df_old, df_new], ignore_index=True)
            df_all = df_all.drop_duplicates(subset=["pair_id"], keep="last")
        except Exception:
            df_all = pd.DataFrame(results)
    else:
        df_all = pd.DataFrame(results)
    df_all.to_csv(OUTPUT_CSV, index=False)


def load_pairs(pair_id_filter: str | None = None) -> list[tuple[str, str]]:
    """Read fold-switching pairs from data/foldswitch_PDB_IDs_full.txt."""
    pairs = []
    with open(PAIRS_FILE) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                fold1, fold2 = parts
                pid = f"{fold1}_{fold2}"
                if pair_id_filter and pid != pair_id_filter:
                    continue
                pairs.append((fold1, fold2))
    return pairs


def build_exclude_set(pairs: list[tuple[str, str]]) -> set[str]:
    """Build set of lowercase PDB IDs that are known fold-switchers."""
    exclude = set()
    for f1, f2 in pairs:
        exclude.add(f1[:4].lower())
        exclude.add(f2[:4].lower())
    return exclude


def main():
    parser = argparse.ArgumentParser(description="Set up matched single-fold controls")
    parser.add_argument("--pair-id", default=None,
                        help="Process only this pair (e.g. 1dzlA_5keqF)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Query APIs but don't create directories")
    args = parser.parse_args()

    # Load all pairs (for exclude set) and filter if needed
    all_pairs = load_pairs()
    exclude_pdbs = build_exclude_set(all_pairs)

    pairs = load_pairs(args.pair_id)
    if not pairs:
        sys.exit(f"No pairs found (filter: {args.pair_id})")

    print(f"Processing {len(pairs)} pairs, excluding {len(exclude_pdbs)} PDB IDs",
          flush=True)

    # Load existing results to allow resuming
    import pandas as pd
    existing = set()
    if os.path.isfile(OUTPUT_CSV):
        try:
            df_existing = pd.read_csv(OUTPUT_CSV)
            existing = set(df_existing["pair_id"])
            print(f"Loaded {len(existing)} existing results from {OUTPUT_CSV}")
        except Exception:
            pass

    results = []
    skipped = 0
    failed = 0

    for i, (fold1, fold2) in enumerate(pairs):
        pair_id = f"{fold1}_{fold2}"
        print(f"\n[{i+1}/{len(pairs)}] {pair_id}", flush=True)

        if pair_id in existing and not args.pair_id:
            print(f"  [SKIP] already in CSV", flush=True)
            skipped += 1
            continue

        ctrl = find_control_for_pair(fold1, fold2, exclude_pdbs, args.dry_run)

        if ctrl:
            if not args.dry_run:
                pdir = create_control_dir(ctrl)
                ctrl["control_dir"] = pdir
            else:
                ctrl["control_dir"] = "(dry run)"

            # Remove sequence from CSV output (large)
            row = {k: v for k, v in ctrl.items() if k != "sequence"}
            results.append(row)
            existing.add(pair_id)
            print(f"  => MATCH: {ctrl['control_pdb']} "
                  f"(UniProt:{ctrl['control_uniprot']}, "
                  f"identity={ctrl['seq_identity']:.1%}, "
                  f"len={ctrl['seq_len']})", flush=True)

            # Save incrementally so progress survives crashes
            if not args.dry_run:
                _save_csv(results)
        else:
            failed += 1
            print(f"  => NO MATCH found", flush=True)

    # Summary
    print(f"\n{'='*60}")
    print(f"Results: {len(results)} matched, {failed} failed, {skipped} skipped")
    if results:
        ids = [r["seq_identity"] for r in results]
        print(f"Identity range: {min(ids):.1%} - {max(ids):.1%}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
