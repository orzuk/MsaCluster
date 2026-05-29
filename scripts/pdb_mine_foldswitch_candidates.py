#!/usr/bin/env python3
"""
PDB-mining for new fold-switching candidate pairs.

Follows the Porter & Looger 2018 methodology on the *current* PDB: identify
pairs of chains within the same sequence-identity cluster (default 95%) whose
TM-align TM-score is below a threshold (default 0.5). The 93-pair dataset
we currently use was mined from PDB ~2018; ~5 years of PDB growth should
yield 20-100 new candidate pairs once duplicates against the existing list
are removed.

Algorithm:
  1. Download the RCSB chain-level sequence-cluster file (bc-{N}.out) at the
     requested identity level (default 95%).
  2. Filter clusters to those with >= 2 chain members.
  3. For each cluster, fetch RCSB entry metadata (resolution, method, length)
     and drop members failing filters (default: X-ray + cryo-EM only,
     resolution <= 3.0 A, length in [50, 500]).
  4. For each surviving cluster, run pairwise TM-align across up to
     --max-chains-per-cluster representatives. Pairs with TM-score below
     --tm-threshold are recorded as candidates.
  5. Deduplicate the candidate list against the existing fold-switching
     dataset (data/foldswitch_PDB_IDs_full.txt).
  6. Write candidates CSV with pair_id, pdb1, chain1, pdb2, chain2,
     tm_by_1, tm_by_2, rmsd, aligned_length, cluster_idx.

CLI examples:
  # Full run (all clusters at 95% identity), parallel over 8 workers:
  python3 scripts/pdb_mine_foldswitch_candidates.py --workers 8

  # 100-cluster test run:
  python3 scripts/pdb_mine_foldswitch_candidates.py --max-clusters 100

  # Tight filter (99% identity, only most stringent fold-switch threshold):
  python3 scripts/pdb_mine_foldswitch_candidates.py \\
      --min-seq-identity 0.99 --tm-threshold 0.4

  # Resume an interrupted run (re-uses cached PDB files + cluster
  # CSV from previous attempt; --output remains the same):
  python3 scripts/pdb_mine_foldswitch_candidates.py --resume

Caches under data/cache/pdb_mining/ (created on first run). Set MSACLUSTER_CACHED_ONLY=1
in env to disable network I/O (errors out on cache misses).
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import os
import random
import shutil
import sys
import time
import urllib.request
import urllib.error
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import requests

# ---------------------------------------------------------------------------
# Project imports
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from config import TMALIGN_EXE
try:
    # tmalign_unified auto-falls back to tmtools (pure-Python) when the
    # TMalign C++ binary is missing, so the mining script works on any node.
    from utils.align_utils import tmalign_unified as run_tmalign_binary
except ImportError:
    run_tmalign_binary = None  # falls back to inline parsing below

# DSSP for secondary-structure-difference filter (Porter 2018 criterion)
try:
    import mdtraj as md  # mdtraj.compute_dssp is the cleanest API
    _HAS_MDTRAJ = True
except ImportError:
    _HAS_MDTRAJ = False

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
RCSB_CLUSTER_URL_TEMPLATES = [
    # New (post-2022, MMseqs2-based, entity-level) endpoint
    "https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-{identity}.txt",
    # Old (BLAST-based, chain-level) endpoint -- kept as fallback in case the
    # new endpoint is unavailable for the requested identity level.
    "https://cdn.rcsb.org/resources/sequence/clusters/bc-{identity}.out",
]
RCSB_DATA_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
RCSB_DATA_POLYMER_URL = "https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
RCSB_DOWNLOAD_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
RCSB_DOWNLOAD_CIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif.gz"

ALLOWED_METHODS = {"X-RAY DIFFRACTION", "ELECTRON MICROSCOPY"}

DEFAULT_CACHE_DIR = ROOT / "data" / "cache" / "pdb_mining"
DEFAULT_PDB_CACHE = DEFAULT_CACHE_DIR / "pdb_files"
DEFAULT_METADATA_CACHE = DEFAULT_CACHE_DIR / "metadata"
DEFAULT_EXISTING_PAIRS = ROOT / "data" / "foldswitch_PDB_IDs_full.txt"
DEFAULT_OUTPUT = ROOT / "data" / "pdb_mining_candidates.csv"


# ---------------------------------------------------------------------------
# Step 1: download + parse cluster file
# ---------------------------------------------------------------------------

def download_cluster_file(identity_pct: int, cache_dir: Path,
                          force: bool = False) -> Tuple[Path, str]:
    """Download the RCSB sequence-cluster file and cache.

    Tries the new MMseqs2-based endpoint
    (clusters-by-entity-{identity}.txt; entity-level) first, then falls
    back to the old BLAST-based endpoint (bc-{identity}.out; chain-level).
    Returns (local_path, format_tag) where format_tag is 'entity' (new) or
    'chain' (old).
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    # Check for either cached file first
    new_local = cache_dir / f"clusters-by-entity-{identity_pct}.txt"
    old_local = cache_dir / f"bc-{identity_pct}.out"
    if not force:
        if new_local.exists():
            return new_local, "entity"
        if old_local.exists():
            return old_local, "chain"
    # Try each URL in turn
    last_err = None
    for tmpl in RCSB_CLUSTER_URL_TEMPLATES:
        url = tmpl.format(identity=identity_pct)
        is_new_endpoint = "clusters-by-entity" in tmpl
        local = new_local if is_new_endpoint else old_local
        print(f"[mine] trying {url}")
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                data = r.read()
        except urllib.error.URLError as e:
            last_err = e
            print(f"[mine]   -> {e}")
            continue
        local.write_bytes(data)
        fmt = "entity" if is_new_endpoint else "chain"
        print(f"[mine] downloaded {local} ({len(data) / 1024:.1f} KB, "
              f"format={fmt})")
        return local, fmt
    raise RuntimeError(
        f"Failed to download RCSB cluster file at identity={identity_pct}; "
        f"last error: {last_err}. Supported identity levels are usually "
        f"{{30, 40, 50, 70, 90, 95, 100}} for the new MMseqs2-based "
        f"endpoint."
    )


def parse_clusters(cluster_file: Path, fmt: str) -> List[List[str]]:
    """Parse RCSB cluster file into list of clusters.

    fmt='entity' (new MMseqs2 format): members are 'PDBID_ENTITYID' (e.g.
    '4HHB_1'). These need a later entity->chain lookup pass before TM-align.
    fmt='chain' (old BLAST format): members are 'PDBID_CHAINID' (e.g.
    '4HHB_A'); usable directly for TM-align.
    Returns list of [member, ...] where member is the raw cluster token.
    """
    out: List[List[str]] = []
    with cluster_file.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            members = line.split()
            normalized = []
            for m in members:
                if "_" not in m:
                    continue
                p, suffix = m.split("_", 1)
                if fmt == "chain":
                    # chain ID -- upper-case is the convention
                    normalized.append(f"{p.lower()}_{suffix.upper()}")
                else:
                    # entity ID -- a numeric suffix, leave as-is
                    normalized.append(f"{p.lower()}_{suffix}")
            if normalized:
                out.append(normalized)
    return out


def fetch_entity_chain_id(pdb_id: str, entity_id: str,
                          cache_dir: Path) -> Optional[str]:
    """For an entity reference (PDBID_ENTITYID), return one representative
    author chain ID (the first one in the entity's auth_asym_ids list).
    Cached as {cache_dir}/{pdb_id}_{entity_id}.json.

    Required only for the new entity-level cluster format.
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / f"{pdb_id.lower()}_{entity_id}.json"
    if cache_path.exists():
        try:
            j = json.loads(cache_path.read_text())
            return j.get("chain_id")
        except Exception:
            cache_path.unlink(missing_ok=True)
    j = _http_get_json(RCSB_DATA_POLYMER_URL.format(
        pdb_id=pdb_id.lower(), entity_id=entity_id))
    if j is None:
        cache_path.write_text(json.dumps({"_unavailable": True}))
        return None
    container = j.get("rcsb_polymer_entity_container_identifiers", {}) or {}
    auth_asym_ids = container.get("auth_asym_ids", []) or []
    chain_id = auth_asym_ids[0] if auth_asym_ids else None
    cache_path.write_text(json.dumps({"chain_id": chain_id}))
    return chain_id


def entity_clusters_to_chain_clusters(
    clusters_entity: List[List[str]],
    cache_dir: Path,
    max_per_cluster: int = 20,
    n_threads: int = 16,
    seed: int = 42,
) -> List[List[str]]:
    """Convert entity-level cluster members (PDBID_ENTITYID) to chain-level
    (PDBID_CHAIN) by looking up the first author chain for each entity.

    Subsamples each cluster to at most max_per_cluster members BEFORE the
    lookup (which is the same cap we apply during TM-align, so doing the
    full lookup would just waste API calls on members we'll never use).
    HTTP lookups run in a ThreadPoolExecutor with n_threads workers --
    each call is short and network-bound, so threading scales well.
    Clusters that lose >=2 members after the lookup are dropped.
    """
    from concurrent.futures import ThreadPoolExecutor

    rng = random.Random(seed)
    # Subsample first
    sampled = []
    for cluster in clusters_entity:
        if len(cluster) > max_per_cluster:
            sampled.append(rng.sample(cluster, max_per_cluster))
        else:
            sampled.append(list(cluster))

    # Build the flat list of unique (pdb_id, entity_id) lookups to run
    unique_entities = sorted({
        e for cl in sampled for e in cl
    })
    n_total = len(unique_entities)
    print(f"[mine]   {n_total} unique entities to resolve "
          f"(across {len(sampled)} clusters, capped at "
          f"max_per_cluster={max_per_cluster})", flush=True)

    chain_lookup: Dict[str, Optional[str]] = {}

    def _lookup_one(entity_ref: str) -> Tuple[str, Optional[str]]:
        p, e = entity_ref.split("_", 1)
        return entity_ref, fetch_entity_chain_id(p, e, cache_dir)

    n_done = 0
    with ThreadPoolExecutor(max_workers=n_threads) as ex:
        for entity_ref, chain_id in ex.map(_lookup_one, unique_entities):
            chain_lookup[entity_ref] = chain_id
            n_done += 1
            if n_done % 200 == 0 or n_done == n_total:
                print(f"[mine]   entity->chain lookup {n_done}/{n_total}",
                      flush=True)

    # Build final chain-level clusters
    out: List[List[str]] = []
    for cluster in sampled:
        chain_members = []
        for entity_ref in cluster:
            chain_id = chain_lookup.get(entity_ref)
            if chain_id:
                p, _ = entity_ref.split("_", 1)
                chain_members.append(f"{p.lower()}_{chain_id.upper()}")
        if len(chain_members) >= 2:
            out.append(chain_members)
    return out


# ---------------------------------------------------------------------------
# Step 2: RCSB metadata
# ---------------------------------------------------------------------------

def _http_get_json(url: str, timeout: int = 30,
                   max_retries: int = 3) -> Optional[dict]:
    """GET a JSON URL with exponential backoff on transient errors."""
    for attempt in range(max_retries):
        try:
            r = requests.get(url, timeout=timeout)
        except requests.exceptions.RequestException as e:
            wait = 2 ** attempt
            print(f"[mine] HTTP error {e}; retry in {wait}s")
            time.sleep(wait)
            continue
        if r.status_code == 200:
            return r.json()
        if r.status_code == 404:
            return None
        # 429 / 5xx: backoff
        wait = 2 ** attempt
        print(f"[mine] HTTP {r.status_code} on {url}; retry in {wait}s")
        time.sleep(wait)
    return None


def fetch_entry_metadata(pdb_id: str, cache_dir: Path) -> Optional[Dict]:
    """Cache-aware fetch of RCSB entry metadata. Returns dict or None.

    Cached as {cache_dir}/{pdb_id_lower}.json.
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / f"{pdb_id.lower()}.json"
    if cache_path.exists():
        try:
            return json.loads(cache_path.read_text())
        except Exception:
            cache_path.unlink(missing_ok=True)
    j = _http_get_json(RCSB_DATA_ENTRY_URL.format(pdb_id=pdb_id.lower()))
    if j is None:
        # Cache a sentinel so we don't re-query missing entries forever
        cache_path.write_text(json.dumps({"_unavailable": True}))
        return None
    cache_path.write_text(json.dumps(j))
    return j


def extract_filterable_fields(entry_meta: Dict) -> Dict:
    """Pull (method, resolution, length) from RCSB entry JSON.

    Returns dict {method, resolution, length_max}. Missing fields -> None.
    Length is the maximum polymer_entity length (the longest chain
    in the entry).
    """
    if not entry_meta or "_unavailable" in entry_meta:
        return {"method": None, "resolution": None, "length_max": None}
    exptl = entry_meta.get("exptl", [{}])
    method = exptl[0].get("method") if exptl else None
    refine = entry_meta.get("refine", [{}])
    resolution = refine[0].get("ls_d_res_high") if refine else None
    if resolution is None:
        rsh = entry_meta.get("rcsb_entry_info", {}).get("resolution_combined", [])
        resolution = rsh[0] if rsh else None
    polymer_lengths = (entry_meta.get("rcsb_entry_info", {})
                       .get("polymer_molecular_weight_maximum") and
                       entry_meta.get("rcsb_entry_info", {}).get("polymer_monomer_count_maximum"))
    length_max = (entry_meta.get("rcsb_entry_info", {})
                  .get("polymer_monomer_count_maximum"))
    return {"method": method, "resolution": resolution, "length_max": length_max}


# ---------------------------------------------------------------------------
# Step 3: PDB file download
# ---------------------------------------------------------------------------

def download_pdb_file(pdb_id: str, cache_dir: Path) -> Optional[Path]:
    """Download a .pdb file (uncompressed) into cache. Falls back to .cif.gz
    if the .pdb endpoint 404s (very large structures).

    Returns local path or None.
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    pdb_id = pdb_id.lower()
    local = cache_dir / f"{pdb_id}.pdb"
    if local.exists() and local.stat().st_size > 0:
        return local
    # Try .pdb first
    url = RCSB_DOWNLOAD_URL.format(pdb_id=pdb_id)
    try:
        r = requests.get(url, timeout=60)
        if r.status_code == 200 and r.content:
            local.write_bytes(r.content)
            return local
    except requests.exceptions.RequestException:
        pass
    # Fall back to gzipped CIF
    url = RCSB_DOWNLOAD_CIF_URL.format(pdb_id=pdb_id)
    try:
        r = requests.get(url, timeout=120)
        if r.status_code == 200 and r.content:
            local_cif = cache_dir / f"{pdb_id}.cif"
            local_cif.write_bytes(gzip.decompress(r.content))
            return local_cif
    except (requests.exceptions.RequestException, OSError):
        pass
    return None


# ---------------------------------------------------------------------------
# Step 4: TM-align driver
# ---------------------------------------------------------------------------

def compute_ss_switch_fraction(pdb1_path: Path, chain1: str,
                                pdb2_path: Path, chain2: str
                                ) -> Optional[float]:
    """Fraction of aligned positions where DSSP secondary structure switches
    between H (helix) and E (strand) -- the Porter 2018 fold-switching
    criterion. Switches H<->C or E<->C are NOT counted (those are weak
    "ordering" changes, not topological refolds).

    Returns the fraction (in [0, 1]) or None if DSSP cannot be computed
    on either structure.

    Implementation. We DSSP both structures with mdtraj. For the aligned
    region we use the union of residue indices common to both chains
    (after stripping HETATM and water); this is a sequence-position-based
    alignment that misses register shifts but is robust to terminal
    truncations between two PDB constructs of the same protein. For the
    overlap, count positions where (ss1 in {H} and ss2 in {E}) or
    (ss1 in {E} and ss2 in {H}).
    """
    if not _HAS_MDTRAJ:
        return None
    try:
        t1 = md.load(str(pdb1_path))
        t2 = md.load(str(pdb2_path))
        # Subset to single chain by author chain ID
        sel1 = t1.topology.select(f"chainid {ord(chain1.upper()) - ord('A')} "
                                  f"and protein")
        sel2 = t2.topology.select(f"chainid {ord(chain2.upper()) - ord('A')} "
                                  f"and protein")
        if len(sel1) == 0 or len(sel2) == 0:
            # Fall back: use all protein atoms (single-chain structures)
            sel1 = t1.topology.select("protein")
            sel2 = t2.topology.select("protein")
        t1c = t1.atom_slice(sel1)
        t2c = t2.atom_slice(sel2)
        ss1 = md.compute_dssp(t1c, simplified=True)[0]
        ss2 = md.compute_dssp(t2c, simplified=True)[0]
    except Exception:
        return None
    # Match by residue position in each chain's sequence
    n = min(len(ss1), len(ss2))
    if n < 20:
        return None
    switches = sum(1 for i in range(n)
                   if (ss1[i] == "H" and ss2[i] == "E")
                   or (ss1[i] == "E" and ss2[i] == "H"))
    return switches / n


def run_pair_tmalign(member1: str, member2: str, pdb_cache: Path
                     ) -> Optional[Dict]:
    """Run TM-align on two 'PDBID_CHAIN' members. Returns dict with
    tm_by_1, tm_by_2, rmsd, aligned_length, or None on failure.
    """
    p1, c1 = member1.split("_", 1)
    p2, c2 = member2.split("_", 1)
    fp1 = download_pdb_file(p1, pdb_cache)
    fp2 = download_pdb_file(p2, pdb_cache)
    if fp1 is None or fp2 is None:
        return None
    if run_tmalign_binary is None:
        return None
    try:
        res = run_tmalign_binary(str(fp1), str(fp2), c1, c2)
    except Exception:
        return None
    return {
        "tm_by_1": res.get("tm_by_1"),
        "tm_by_2": res.get("tm_by_2"),
        "rmsd": res.get("rmsd"),
        "aligned_length": res.get("aligned_length"),
    }


def process_cluster(cluster_idx: int, members: List[str],
                    max_chains_per_cluster: int,
                    tm_threshold: float,
                    min_tm_threshold: float,
                    min_ss_switch_fraction: float,
                    pdb_cache: Path,
                    seed: int) -> List[Dict]:
    """Run pairwise TM-align within one cluster. Returns list of candidate
    pair dicts (min_tm_threshold <= TM <= tm_threshold AND
    SS-switch fraction >= threshold).

    Two-stage filter (Porter 2018 fold-switching criteria):
      1. TM-score in [min_tm_threshold, tm_threshold] -- the upper bound
         ensures topological rearrangement, the lower bound excludes
         unrelated proteins / degenerate alignments where the proteins
         are too divergent for the comparison to be meaningful.
      2. Fraction of residues where DSSP secondary structure switches
         between H (helix) and E (strand) is at least
         min_ss_switch_fraction. Set to 0.0 to disable the SS filter
         (keeping just the TM-score test).
    Caps compute by sampling up to max_chains_per_cluster representatives
    when the cluster is large.
    """
    if len(members) < 2:
        return []
    if len(members) > max_chains_per_cluster:
        rng = random.Random(seed + cluster_idx)
        members = rng.sample(members, max_chains_per_cluster)
    out: List[Dict] = []
    for i in range(len(members)):
        for j in range(i + 1, len(members)):
            r = run_pair_tmalign(members[i], members[j], pdb_cache)
            if r is None:
                continue
            tm_max = max(r["tm_by_1"] or 0.0, r["tm_by_2"] or 0.0)
            if tm_max > tm_threshold or tm_max < min_tm_threshold:
                continue
            # Compute the SS-switch fraction for surviving TM hits.
            p1, c1 = members[i].split("_", 1)
            p2, c2 = members[j].split("_", 1)
            fp1 = download_pdb_file(p1, pdb_cache)
            fp2 = download_pdb_file(p2, pdb_cache)
            ss_frac = None
            if fp1 and fp2 and min_ss_switch_fraction > 0.0:
                ss_frac = compute_ss_switch_fraction(fp1, c1, fp2, c2)
                if ss_frac is None or ss_frac < min_ss_switch_fraction:
                    continue
            pair_id = f"{p1.lower()}{c1.upper()}_{p2.lower()}{c2.upper()}"
            out.append({
                "pair_id": pair_id,
                "pdb1": p1.lower(), "chain1": c1.upper(),
                "pdb2": p2.lower(), "chain2": c2.upper(),
                "tm_by_1": r["tm_by_1"],
                "tm_by_2": r["tm_by_2"],
                "tm_max": tm_max,
                "rmsd": r["rmsd"],
                "aligned_length": r["aligned_length"],
                "ss_switch_fraction": ss_frac,
                "cluster_idx": cluster_idx,
                "cluster_size": len(members),
            })
    return out


# ---------------------------------------------------------------------------
# Step 5: dedup vs existing pairs
# ---------------------------------------------------------------------------

def load_existing_pairs(path: Path) -> set[str]:
    """Load the existing 93-pair list and return set of canonical pair IDs.

    A pair {AB, CD} is canonical-equivalent to {CD, AB}; both orderings are
    added to the set so a candidate hit on either ordering is filtered.
    """
    existing: set[str] = set()
    if not path.exists():
        return existing
    for line in path.read_text().splitlines():
        toks = line.strip().split()
        if len(toks) < 2:
            continue
        a, b = toks[0], toks[1]
        # Each token is like 1abcA (4-char PDB ID + 1-char chain)
        if len(a) >= 5 and len(b) >= 5:
            n_a = f"{a[:4].lower()}{a[4].upper()}"
            n_b = f"{b[:4].lower()}{b[4].upper()}"
            existing.add(f"{n_a}_{n_b}")
            existing.add(f"{n_b}_{n_a}")
    return existing


# ---------------------------------------------------------------------------
# Step 6: metadata filtering helpers
# ---------------------------------------------------------------------------

def filter_cluster_members_by_metadata(
    members: List[str],
    metadata_cache: Path,
    allowed_methods: set[str],
    max_resolution: float,
    min_length: int,
    max_length: int,
) -> List[str]:
    """Filter members of a single cluster by entry-level RCSB metadata.

    Drops any chain whose PDB entry (a) has missing metadata, (b) was solved
    by a method not in allowed_methods, (c) has resolution > max_resolution,
    or (d) has max polymer length outside [min_length, max_length].
    """
    kept = []
    seen_pdbs: dict[str, dict] = {}
    for m in members:
        pdb_id = m.split("_", 1)[0]
        if pdb_id in seen_pdbs:
            meta = seen_pdbs[pdb_id]
        else:
            meta = fetch_entry_metadata(pdb_id, metadata_cache)
            seen_pdbs[pdb_id] = meta or {}
        if not meta or "_unavailable" in meta:
            continue
        flds = extract_filterable_fields(meta)
        if flds["method"] not in allowed_methods:
            continue
        if flds["resolution"] is None or flds["resolution"] > max_resolution:
            continue
        L = flds["length_max"]
        if L is None or L < min_length or L > max_length:
            continue
        kept.append(m)
    return kept


# ---------------------------------------------------------------------------
# Worker entry point (top-level for picklability)
# ---------------------------------------------------------------------------

def _worker_process_cluster(args_tuple):
    """Wrapper for ProcessPoolExecutor — must be top-level for pickling."""
    (cluster_idx, members, max_per_cluster, tm_threshold,
     min_tm_threshold,
     min_ss_switch_fraction, pdb_cache_str, seed) = args_tuple
    return process_cluster(
        cluster_idx, members, max_per_cluster,
        tm_threshold, min_tm_threshold, min_ss_switch_fraction,
        Path(pdb_cache_str), seed,
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--min-seq-identity", type=float, default=0.95,
                    help="Sequence-identity threshold for RCSB clusters. "
                         "Must map to one of {0.30, 0.40, 0.50, 0.70, 0.90, "
                         "0.95, 1.00}. Default 0.95.")
    ap.add_argument("--max-clusters", type=int, default=None,
                    help="If set, process only the first N clusters with "
                         ">=2 members (after sorting by size descending). "
                         "Useful for test runs. Default: process all.")
    ap.add_argument("--tm-threshold", type=float, default=0.5,
                    help="TM-score upper bound: TM <= this -> pair is candidate "
                         "fold-switcher. TM uses max(tm_by_1, tm_by_2). Default 0.5.")
    ap.add_argument("--min-tm-threshold", type=float, default=0.3,
                    help="TM-score lower bound: TM < this -> pair is REJECTED "
                         "(unrelated proteins or degenerate alignment). "
                         "Porter 2018 standard is [0.3, 0.5]. Default 0.3.")
    ap.add_argument("--max-candidates-per-cluster", type=int, default=1,
                    help="Within a sequence-similarity cluster, keep at most "
                         "this many candidate pairs (the ones with highest "
                         "max(tm_by_1, tm_by_2), since those are the strongest "
                         "fold-switch candidates). Default 1, which avoids "
                         "redundant pairs from one fold-switcher appearing in "
                         "multiple PDB depositions. Set to 0 to disable.")
    ap.add_argument("--min-ss-switch-fraction", type=float, default=0.10,
                    help="Porter 2018 criterion: minimum fraction of residues "
                         "where DSSP secondary structure switches between H "
                         "(helix) and E (strand). Filters out domain swaps "
                         "and hinge motions whose SS is unchanged. Set to "
                         "0.0 to disable the SS filter (keep only the "
                         "TM-score test). Default 0.10.")
    ap.add_argument("--max-resolution", type=float, default=3.0,
                    help="Drop PDB entries with resolution worse than this "
                         "(Angstroms). Default 3.0.")
    ap.add_argument("--min-length", type=int, default=50,
                    help="Drop entries whose max polymer is shorter than "
                         "this many AA. Default 50.")
    ap.add_argument("--max-length", type=int, default=500,
                    help="Drop entries whose max polymer is longer than "
                         "this many AA. Default 500.")
    ap.add_argument("--max-chains-per-cluster", type=int, default=20,
                    help="Cap pairwise TM-aligns by sampling up to this many "
                         "chains per cluster (random with --seed). Default 20.")
    ap.add_argument("--workers", type=int, default=1,
                    help="Process-pool workers for parallel TM-align. "
                         "Default 1. Recommend 8-32 on cluster.")
    ap.add_argument("--cache-dir", default=str(DEFAULT_CACHE_DIR),
                    help=f"Root cache dir. Default {DEFAULT_CACHE_DIR}.")
    ap.add_argument("--output", default=str(DEFAULT_OUTPUT),
                    help=f"Output CSV path. Default {DEFAULT_OUTPUT}.")
    ap.add_argument("--existing-pairs", default=str(DEFAULT_EXISTING_PAIRS),
                    help="Path to existing fold-switching pair list (for "
                         "deduplication). Default data/foldswitch_PDB_IDs_full.txt.")
    ap.add_argument("--allow-nmr", action="store_true",
                    help="Include NMR structures (default: X-ray + cryo-EM only).")
    ap.add_argument("--seed", type=int, default=42,
                    help="Seed for cluster sub-sampling. Default 42.")
    ap.add_argument("--checkpoint-every", type=int, default=50,
                    help="Write partial CSV every N candidates. Default 50.")
    ap.add_argument("--resume", action="store_true",
                    help="If set, skip clusters whose chain pairs are already "
                         "in --output (resume an interrupted run).")
    ap.add_argument("--force-cluster-redownload", action="store_true",
                    help="Re-download the RCSB cluster file even if cached.")
    return ap.parse_args()


def main():
    args = parse_args()

    cache_dir = Path(args.cache_dir)
    pdb_cache = cache_dir / "pdb_files"
    metadata_cache = cache_dir / "metadata"
    cache_dir.mkdir(parents=True, exist_ok=True)
    pdb_cache.mkdir(parents=True, exist_ok=True)
    metadata_cache.mkdir(parents=True, exist_ok=True)

    # --- 1. Download + parse cluster file ---
    identity_pct = int(round(args.min_seq_identity * 100))
    if identity_pct not in (30, 40, 50, 70, 90, 95, 100):
        print(f"WARN: identity={identity_pct} is not a standard RCSB "
              f"clustering level; download may 404.")
    cluster_file, cluster_fmt = download_cluster_file(
        identity_pct, cache_dir,
        force=args.force_cluster_redownload,
    )
    clusters = parse_clusters(cluster_file, fmt=cluster_fmt)
    clusters_multi = [c for c in clusters if len(c) >= 2]
    # Sort by cluster size DESC so test runs hit the most informative clusters
    clusters_multi.sort(key=len, reverse=True)
    if args.max_clusters is not None:
        clusters_multi = clusters_multi[:args.max_clusters]

    print(f"[mine] {len(clusters)} total clusters at identity={identity_pct}%; "
          f"{len(clusters_multi)} have >=2 members; processing "
          f"{len(clusters_multi)} after --max-clusters")

    # If the new MMseqs2 endpoint gave us entity-level cluster members,
    # do the entity->chain lookup pass before metadata filtering. We
    # subsample to max_chains_per_cluster here (the same cap as TM-align)
    # so we don't waste API calls on members we'd never align.
    if cluster_fmt == "entity":
        entity_chain_cache = cache_dir / "entity_chains"
        print(f"[mine] resolving entity -> chain IDs (caching to "
              f"{entity_chain_cache})", flush=True)
        clusters_multi = entity_clusters_to_chain_clusters(
            clusters_multi, entity_chain_cache,
            max_per_cluster=args.max_chains_per_cluster,
            n_threads=max(args.workers, 8),
            seed=args.seed,
        )
        print(f"[mine] {len(clusters_multi)} clusters remain after "
              f"entity->chain", flush=True)

    # --- 2. Metadata filtering ---
    allowed = set(ALLOWED_METHODS)
    if args.allow_nmr:
        allowed.add("SOLUTION NMR")
        allowed.add("SOLID-STATE NMR")
    print(f"[mine] filtering by metadata: methods={sorted(allowed)}, "
          f"resolution<={args.max_resolution}, "
          f"length in [{args.min_length}, {args.max_length}]")
    # Prefetch metadata for all PDB IDs in parallel before the filter pass
    # (each entry call is cheap and cacheable; threading is the right
    # parallelism here since the calls are network-I/O-bound).
    from concurrent.futures import ThreadPoolExecutor
    unique_pdbs = sorted({m.split("_", 1)[0] for cl in clusters_multi
                          for m in cl})
    print(f"[mine] prefetching metadata for {len(unique_pdbs)} unique PDB "
          f"entries (threaded)", flush=True)
    n_done = 0
    with ThreadPoolExecutor(max_workers=max(args.workers, 8)) as ex:
        for _ in ex.map(lambda pid: fetch_entry_metadata(pid, metadata_cache),
                        unique_pdbs):
            n_done += 1
            if n_done % 500 == 0 or n_done == len(unique_pdbs):
                print(f"[mine]   metadata prefetch {n_done}/{len(unique_pdbs)}",
                      flush=True)

    filtered_clusters = []
    for i, cl in enumerate(clusters_multi):
        if i % 100 == 0 and i > 0:
            print(f"[mine]   metadata-filter cluster {i}/{len(clusters_multi)}",
                  flush=True)
        kept = filter_cluster_members_by_metadata(
            cl, metadata_cache, allowed,
            args.max_resolution, args.min_length, args.max_length,
        )
        if len(kept) >= 2:
            filtered_clusters.append(kept)
    print(f"[mine] {len(filtered_clusters)} clusters pass metadata filter",
          flush=True)

    # --- 3. Optional resume ---
    seen_pair_ids: set[str] = set()
    candidates: List[Dict] = []
    output_path = Path(args.output)
    if args.resume and output_path.exists():
        try:
            prev = pd.read_csv(output_path)
            seen_pair_ids = set(prev["pair_id"].astype(str).tolist())
            candidates.extend(prev.to_dict("records"))
            print(f"[mine] resume: skipping {len(seen_pair_ids)} pair_ids "
                  f"already in {output_path}")
        except Exception as e:
            print(f"[mine] WARN: --resume but failed to read {output_path}: {e}")

    # --- 4. Pairwise TM-align across clusters ---
    if run_tmalign_binary is None:
        print("ERROR: utils/align_utils.tmalign_unified not importable. "
              "TM-align driver is required.")
        sys.exit(1)
    if not os.path.isfile(TMALIGN_EXE) or not os.access(TMALIGN_EXE, os.X_OK):
        print(f"[mine] TM-align binary not found at {TMALIGN_EXE}; "
              "tmalign_unified will fall back to tmtools (Python).")

    if args.min_ss_switch_fraction > 0.0 and not _HAS_MDTRAJ:
        print("[mine] WARN: --min-ss-switch-fraction > 0 but mdtraj not "
              "installed. SS-switch filter will be skipped (all TM-passing "
              "pairs kept). Install mdtraj on the cluster to enable.")

    work = [
        (i, cl, args.max_chains_per_cluster, args.tm_threshold,
         args.min_tm_threshold,
         args.min_ss_switch_fraction, str(pdb_cache), args.seed)
        for i, cl in enumerate(filtered_clusters)
    ]

    n_clusters_done = 0
    if args.workers == 1:
        for w in work:
            results = _worker_process_cluster(w)
            for r in results:
                if r["pair_id"] not in seen_pair_ids:
                    seen_pair_ids.add(r["pair_id"])
                    candidates.append(r)
            n_clusters_done += 1
            if n_clusters_done % 10 == 0:
                print(f"[mine]   processed {n_clusters_done}/{len(work)} clusters; "
                      f"{len(candidates)} candidates so far", flush=True)
            if len(candidates) % args.checkpoint_every == 0 and candidates:
                pd.DataFrame(candidates).to_csv(output_path, index=False)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futures = {ex.submit(_worker_process_cluster, w): w[0]
                       for w in work}
            for fut in as_completed(futures):
                try:
                    results = fut.result()
                except Exception as e:
                    print(f"[mine] cluster {futures[fut]} failed: {e}")
                    results = []
                for r in results:
                    if r["pair_id"] not in seen_pair_ids:
                        seen_pair_ids.add(r["pair_id"])
                        candidates.append(r)
                n_clusters_done += 1
                if n_clusters_done % 10 == 0:
                    print(f"[mine]   processed {n_clusters_done}/{len(work)} "
                          f"clusters; {len(candidates)} candidates so far")
                if (len(candidates) and
                        len(candidates) % args.checkpoint_every == 0):
                    pd.DataFrame(candidates).to_csv(output_path, index=False)

    # --- 5. Dedup against existing 93 pairs ---
    existing = load_existing_pairs(Path(args.existing_pairs))
    print(f"[mine] loaded {len(existing) // 2} existing pairs for dedup")

    if not candidates:
        print("[mine] no candidates passed thresholds")
        pd.DataFrame(columns=["pair_id"]).to_csv(output_path, index=False)
        return

    df = pd.DataFrame(candidates)
    df["is_existing"] = df["pair_id"].isin(existing)
    # Also flag the reverse-ordering case
    df["is_existing"] |= df.apply(
        lambda r: (f"{r['pdb2']}{r['chain2']}_{r['pdb1']}{r['chain1']}"
                   in existing), axis=1,
    )
    # Per-cluster dedup: a fold-switcher with N PDB depositions of variants
    # generates ~N redundant pairs in one cluster. Keep at most
    # max_candidates_per_cluster per cluster, picking the strongest (highest tm_max,
    # i.e., closest to passing the fold-similarity floor). Set to 0 to keep all.
    n_before_dedup = len(df)
    if args.max_candidates_per_cluster > 0:
        df = (df.sort_values("tm_max", ascending=False)
                .groupby("cluster_idx", as_index=False)
                .head(args.max_candidates_per_cluster)
                .reset_index(drop=True))
        if len(df) < n_before_dedup:
            print(f"[mine] per-cluster dedup: {n_before_dedup} -> {len(df)} "
                  f"(kept top {args.max_candidates_per_cluster} per cluster)")
    new_pairs = df[~df["is_existing"]].sort_values("tm_max").reset_index(drop=True)
    print(f"[mine] {len(df)} total pairs with TM <= {args.tm_threshold} "
          f"(after per-cluster dedup), {len(new_pairs)} new (not in existing 93)")

    df.to_csv(output_path, index=False)
    new_pairs.to_csv(output_path.with_suffix(".new_only.csv"), index=False)
    print(f"[mine] wrote {output_path} and "
          f"{output_path.with_suffix('.new_only.csv')}")


if __name__ == "__main__":
    main()
