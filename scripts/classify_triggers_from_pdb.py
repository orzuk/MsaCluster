"""Classify fold-switching trigger types from PDB metadata.

For each pair (pdb1_chain1, pdb2_chain2) in data/foldswitch_PDB_IDs_full.txt
(or a user-supplied list), inspect the two PDB files and compare:

  - HETATM ligands (excluding crystallographic buffers/cofactors)
  - Chain count and chain-sequence multiplicity
  - SEQRES sequence (point mutations, length differences)
  - Heteromeric vs homomeric biological assemblies

and emit a trigger-type classification per pair with auditable evidence.

Trigger taxonomy (single-label, prioritised top-down):
  1. mutation              : pdb1.SEQRES != pdb2.SEQRES (non-trivial diff)
  2. protein_binding       : pair has a HETERO-COMPLEX partner in one entry,
                             monomer (or homo-oligomer of same chain) in the
                             other -- e.g. RfaH apo vs RfaH-RNAP
  3. ligand                : non-buffer HETATM present in one but not other
                             (or different non-buffer HETATM sets)
  4. oligomerization       : same sequence, same ligands, but chain count
                             differs (1 vs >=2 of the same chain)
  5. equilibrium_or_unknown: all else; needs paper inspection

Output: docs/triggers_from_pdb.csv with columns:
  pair_id, pdbid1, chain1, pdbid2, chain2, trigger_class,
  evidence, ligands_pdb1, ligands_pdb2, n_unique_chains_pdb1,
  n_unique_chains_pdb2, seq_identity, length_diff, notes

Usage:
    # From repo root, with Pipeline/FoldPairs/<pair_id>/chain_pdb_files/ populated:
    python3 scripts/classify_triggers_from_pdb.py
    python3 scripts/classify_triggers_from_pdb.py --pairs 2n54B_2hdmA,6c6sD_2ougC
    # Or download PDBs on the fly (slower; uses RCSB):
    python3 scripts/classify_triggers_from_pdb.py --fetch-rcsb
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import urllib.request
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_THIS_DIR))

from config import DATA_DIR, MAIN_DIR  # noqa: E402

# Buffers / common crystallization additives / cofactors that are not real
# biological "trigger ligands" -- exclude these from the ligand-trigger call.
# (Conservative: keep e.g. metal ions only if they're not the obvious crystal
# buffer ions. NAG/MAN sugars are kept as potential triggers; user can
# refine.)
CRYSTAL_BUFFERS: Set[str] = {
    # waters / cryoprotectants / detergents / buffer ions
    "HOH", "WAT", "DOD", "PEG", "PG4", "P6G", "PE4", "PEU",
    "GOL", "EDO", "MPD", "DMS", "BME", "DTT", "MES", "EPE",
    "TRS", "BCT", "FMT", "ACT", "ACY", "IPA", "IOH", "TFA",
    "SO4", "PO4", "CO3", "NO3", "CL", "BR", "IOD", "F",
    "NA", "K", "LI", "CS", "RB",
    "MG", "CA", "ZN", "MN", "FE", "NI", "CO", "CU",  # metals -- conservative
    # crystallization tags / linkers
    "GLY", "ALA",  # placeholder if used in linkers (rare)
}


def _is_real_ligand(resname: str) -> bool:
    """True if HETATM resname is plausibly a 'real' ligand, not a buffer."""
    rn = resname.strip().upper()
    if rn in CRYSTAL_BUFFERS:
        return False
    # Exclude single-atom ions whose codes start with the metal symbol
    # already handled by CRYSTAL_BUFFERS for the common ones.
    return True


def _pdb_path_for(pair_id: str, pdbid_full: str, fetch: bool) -> Optional[str]:
    """Find the PDB file for one entry of the pair.

    Looks under Pipeline/FoldPairs/<pair_id>/ for the canonical truth file
    (e.g. 4dxr.pdb if pdbid_full='4dxrA'). Falls back to RCSB download if
    --fetch-rcsb is set.
    """
    pdb4 = pdbid_full[:4]
    pair_dir = Path(DATA_DIR) / pair_id
    candidates = [
        pair_dir / f"{pdb4}.pdb",
        pair_dir / f"{pdb4}_cif.pdb",
        pair_dir / "chain_pdb_files" / f"{pdb4}.pdb",
    ]
    for c in candidates:
        if c.is_file():
            return str(c)
    if not fetch:
        return None
    # Download from RCSB
    url = f"https://files.rcsb.org/download/{pdb4.upper()}.pdb"
    dest = pair_dir / f"{pdb4}.pdb"
    pair_dir.mkdir(parents=True, exist_ok=True)
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            data = r.read()
        with open(dest, "wb") as f:
            f.write(data)
        return str(dest)
    except Exception as e:
        print(f"  [warn] could not fetch {pdb4} from RCSB: {e}", file=sys.stderr)
        return None


def _parse_pdb(path: str, target_chain: str) -> Dict:
    """Lightweight PDB parser (no biopython dependency).

    Reads HEADER/COMPND/SEQRES/ATOM/HETATM lines directly. Returns:
      {
        'chain_seqs': {chain_id: seqres_or_atom_derived_one_letter},
        'ligands_per_chain': {chain_id: set(of HETATM resnames seen near chain)},
        'all_ligands': set(of all non-buffer HETATM resnames in file),
        'compnd_chains': [list of chain ids declared in COMPND],
        'header': str,
        'title': str,
      }
    SEQRES is preferred; ATOM-derived sequence is used as fallback when
    SEQRES is missing for a chain (common in CIF-derived files).
    """
    aa3to1 = {
        'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
        'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
        'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
        'SEC':'U','PYL':'O','MSE':'M',
    }

    seqres_chain: Dict[str, List[str]] = {}
    atom_chain: Dict[str, List[Tuple[int, str]]] = {}  # (resnum, aa1) per chain
    ligands_per_chain: Dict[str, Set[str]] = {}
    all_ligands: Set[str] = set()
    header = ""
    title_parts: List[str] = []

    try:
        f = open(path, "r", errors="replace")
    except OSError:
        return {"chain_seqs": {}, "ligands_per_chain": {}, "all_ligands": set(),
                "compnd_chains": [], "header": "", "title": ""}

    with f:
        for line in f:
            rec = line[:6]
            if rec == "HEADER":
                header = line[10:50].strip()
            elif rec == "TITLE ":
                title_parts.append(line[10:80].rstrip())
            elif rec == "SEQRES":
                chain_id = line[11]
                resns = line[19:70].split()
                seq_chars = [aa3to1.get(r.upper(), "X") for r in resns]
                seqres_chain.setdefault(chain_id, []).extend(seq_chars)
            elif rec == "ATOM  " and line[12:16].strip() == "CA":
                # ATOM-record fallback for chains missing SEQRES
                resname = line[17:20].strip().upper()
                chain_id = line[21]
                try:
                    resnum = int(line[22:26])
                except ValueError:
                    continue
                aa1 = aa3to1.get(resname, "X")
                atom_chain.setdefault(chain_id, []).append((resnum, aa1))
            elif rec == "HETATM":
                resname = line[17:20].strip().upper()
                chain_id = line[21]
                if _is_real_ligand(resname):
                    all_ligands.add(resname)
                    ligands_per_chain.setdefault(chain_id, set()).add(resname)

    # Build final chain_seqs: prefer SEQRES, fall back to ATOM-derived
    chain_seqs: Dict[str, str] = {}
    for c, chars in seqres_chain.items():
        chain_seqs[c] = "".join(chars)
    for c, items in atom_chain.items():
        if c in chain_seqs and chain_seqs[c]:
            continue
        # Dedupe by residue number (alt-loc duplicates), keep monotone order
        seen = set()
        ordered = []
        for resnum, aa1 in sorted(items):
            if resnum in seen:
                continue
            seen.add(resnum)
            ordered.append(aa1)
        chain_seqs[c] = "".join(ordered)

    return {
        "chain_seqs": chain_seqs,
        "ligands_per_chain": ligands_per_chain,
        "all_ligands": all_ligands,
        "compnd_chains": sorted(chain_seqs.keys()),
        "header": header,
        "title": " ".join(title_parts),
    }


def _best_chain_match(meta1: Dict, meta2: Dict
                      ) -> Tuple[str, str, str, str, float]:
    """Across all chains in pdb1 and pdb2, return the (chain1, seq1, chain2,
    seq2, identity) tuple maximising sequence identity over the overlap.

    Robust to: chain-ID mismatches between deposited PDBs (the dataset's
    pair_id-encoded chain may not match the SEQRES chain ID); chains
    missing SEQRES (uses ATOM-record fallback); long terminal truncations
    in one entry but not the other (identity is computed over the overlap
    so a truncated terminus does not hurt the score).
    """
    best = ("", "", "", "", -1.0)
    for c1, s1 in meta1["chain_seqs"].items():
        if not s1 or len(s1) < 20:
            continue
        for c2, s2 in meta2["chain_seqs"].items():
            if not s2 or len(s2) < 20:
                continue
            iden = _seq_identity(s1, s2)
            if iden > best[4]:
                best = (c1, s1, c2, s2, iden)
    return best


def _seq_identity(a: str, b: str) -> float:
    """Fast ungapped identity over min-length window. Returns fraction."""
    if not a or not b:
        return 0.0
    L = min(len(a), len(b))
    same = sum(1 for x, y in zip(a[:L], b[:L]) if x == y and x != "X")
    return same / L if L else 0.0


def _classify_pair(meta1: Dict, meta2: Dict, chain1: str, chain2: str
                   ) -> Tuple[str, str]:
    """Return (trigger_class, evidence_string).

    v2: match chains by best-sequence-identity across all chains in both
    files (not by user-supplied chain IDs, which often don't match the
    SEQRES chain IDs in the deposited PDBs). Mutation is decided only on
    over-the-overlap identity; large length_diff alone (terminal
    truncations between the same protein) no longer triggers mutation.
    """
    lig1 = meta1["all_ligands"]
    lig2 = meta2["all_ligands"]

    # Best-matching chain pair (regardless of user-supplied IDs)
    bc1, s1, bc2, s2, id_frac = _best_chain_match(meta1, meta2)
    len_diff = abs(len(s1) - len(s2)) if s1 and s2 else 0

    # Count DISTINCT chain sequences (homomeric duplicates collapse)
    seqs1 = {s for s in meta1["chain_seqs"].values() if s and len(s) >= 20}
    seqs2 = {s for s in meta2["chain_seqs"].values() if s and len(s) >= 20}
    n_unique1 = len(seqs1)
    n_unique2 = len(seqs2)

    # Total chain instance count (multimers)
    n_chains1 = len([c for c, s in meta1["chain_seqs"].items() if s and len(s) >= 20])
    n_chains2 = len([c for c, s in meta2["chain_seqs"].items() if s and len(s) >= 20])

    # 1. MUTATION: best-matching chain pair shows < 95% identity over the
    #    overlap. This now correctly skips terminal truncations.
    if s1 and s2 and id_frac < 0.95:
        return ("mutation",
                f"best chain pair {bc1}/{bc2}: identity={id_frac:.3f} "
                f"over {min(len(s1), len(s2))} residues "
                f"(length_diff={len_diff})")

    # 2. PROTEIN_BINDING: heteromeric vs homomeric.
    #    One entry has >1 distinct chain sequence (a complex); the other
    #    is a homo-oligomer or monomer of a single chain.
    if (n_unique1 > 1) != (n_unique2 > 1):
        return ("protein_binding",
                f"pdb1 distinct_chains={n_unique1}; "
                f"pdb2 distinct_chains={n_unique2} "
                f"(best-match identity={id_frac:.3f})")

    # 3. LIGAND: same protein (identity >= 0.95), distinct non-buffer
    #    HETATM sets.
    if lig1 != lig2:
        diff = (lig1 ^ lig2)
        if diff:
            return ("ligand",
                    f"ligands_pdb1={sorted(lig1)}; ligands_pdb2={sorted(lig2)}; "
                    f"diff={sorted(diff)} (identity={id_frac:.3f})")

    # 4. OLIGOMERIZATION: same sequence, same ligands, different total
    #    chain INSTANCE count (monomer vs dimer/trimer/etc.).
    if abs(n_chains1 - n_chains2) > 0:
        return ("oligomerization",
                f"chain instances: pdb1={n_chains1}, pdb2={n_chains2} "
                f"(distinct seqs n1={n_unique1}, n2={n_unique2}; "
                f"identity={id_frac:.3f})")

    return ("equilibrium_or_unknown",
            f"identical sequence over {min(len(s1), len(s2))} residues "
            f"(identity={id_frac:.3f}), same ligands, same chain count "
            f"(n1={n_chains1}, n2={n_chains2}); paper inspection required")


def _load_pair_list(path: str) -> List[Tuple[str, str, str, str, str]]:
    """Parse data/foldswitch_PDB_IDs_full.txt -> [(pair_id, pdb1, ch1, pdb2, ch2), ...]"""
    out = []
    with open(path) as f:
        for line in f:
            toks = line.strip().split()
            if len(toks) < 2:
                continue
            e1, e2 = toks[0], toks[1]
            if len(e1) < 5 or len(e2) < 5:
                continue
            pdb1, ch1 = e1[:4].lower(), e1[4]
            pdb2, ch2 = e2[:4].lower(), e2[4]
            pair_id = f"{pdb1}{ch1}_{pdb2}{ch2}"
            out.append((pair_id, pdb1, ch1, pdb2, ch2))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", default=os.path.join(MAIN_DIR, "data", "foldswitch_PDB_IDs_full.txt"),
                    help="Path to fold-switching pair list (tab-separated).")
    ap.add_argument("--output", default=os.path.join(MAIN_DIR, "docs", "triggers_from_pdb.csv"),
                    help="Output CSV path.")
    ap.add_argument("--pairs", default="ALL",
                    help="Comma-separated pair_ids to classify; default ALL.")
    ap.add_argument("--fetch-rcsb", action="store_true",
                    help="Download missing PDB files from RCSB.")
    args = ap.parse_args()

    pair_list = _load_pair_list(args.input)
    if args.pairs.upper() != "ALL":
        wanted = {s.strip() for s in args.pairs.split(",") if s.strip()}
        pair_list = [p for p in pair_list if p[0] in wanted]

    print(f"Classifying {len(pair_list)} pairs")

    rows = []
    summary = Counter()
    for pair_id, pdb1, ch1, pdb2, ch2 in pair_list:
        p1 = _pdb_path_for(pair_id, f"{pdb1}{ch1}", fetch=args.fetch_rcsb)
        p2 = _pdb_path_for(pair_id, f"{pdb2}{ch2}", fetch=args.fetch_rcsb)
        if not p1 or not p2:
            rows.append({
                "pair_id": pair_id, "pdbid1": pdb1, "chain1": ch1,
                "pdbid2": pdb2, "chain2": ch2,
                "trigger_class": "unresolved_no_pdb",
                "evidence": f"missing PDB file(s): pdb1={p1}, pdb2={p2}",
                "ligands_pdb1": "", "ligands_pdb2": "",
                "n_unique_chains_pdb1": "", "n_unique_chains_pdb2": "",
                "seq_identity": "", "length_diff": "", "notes": "",
            })
            summary["unresolved_no_pdb"] += 1
            continue

        m1 = _parse_pdb(p1, ch1)
        m2 = _parse_pdb(p2, ch2)
        trig, evid = _classify_pair(m1, m2, ch1, ch2)
        bc1, s1, bc2, s2, id_frac = _best_chain_match(m1, m2)
        rows.append({
            "pair_id": pair_id, "pdbid1": pdb1, "chain1_requested": ch1,
            "pdbid2": pdb2, "chain2_requested": ch2,
            "chain1_best_match": bc1, "chain2_best_match": bc2,
            "trigger_class": trig,
            "evidence": evid,
            "ligands_pdb1": "|".join(sorted(m1["all_ligands"])),
            "ligands_pdb2": "|".join(sorted(m2["all_ligands"])),
            "n_unique_chains_pdb1": len({s for s in m1["chain_seqs"].values()
                                         if s and len(s) >= 20}),
            "n_unique_chains_pdb2": len({s for s in m2["chain_seqs"].values()
                                         if s and len(s) >= 20}),
            "seq_identity": f"{id_frac:.3f}" if id_frac >= 0 else "",
            "length_diff": abs(len(s1) - len(s2)) if s1 and s2 else "",
            "notes": f"{m1['header']} | {m2['header']}",
        })
        summary[trig] += 1
        print(f"  {pair_id:<20s}  {trig:<24s}  {evid[:80]}")

    # Write CSV
    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(args.output, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {args.output}: {len(rows)} rows")
    print("\nTrigger-class distribution:")
    for k, v in summary.most_common():
        print(f"  {k:<24s}  {v}")


if __name__ == "__main__":
    main()
