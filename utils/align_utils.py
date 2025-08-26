# Utilities for aligning two structures,contact maps and sequences
# utils/align_utils.py
from __future__ import annotations
from Bio import pairwise2
from Bio.Align import substitution_matrices

from tmtools import tm_align
from tmtools.io import get_structure, get_residue_data
from utils.msa_utils import *
from utils.protein_utils import *
import os, re, shutil, subprocess, tempfile, pathlib
from typing import Optional, Dict, Any


from config import * # uses TMALIGN_EXE and USE_TMALIGN_BINARY
import numpy as np
# ---------------- helpers (use YOUR loaders) ---------------- #

def _ensure_local_pdb(pdb_or_path: str) -> Tuple[str, bool]:
    """
    If given a 4-letter PDB ID, fetch to a temp .pdb via your util.
    Otherwise return the path as-is. Returns (path, is_temp).
    """
    if isinstance(pdb_or_path, str) and len(pdb_or_path) == 4 and pdb_or_path.isalnum():
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
        fetch_and_save_pdb_file(pdb_or_path, "pdb", tmp.name)
        return tmp.name, True
    return pdb_or_path, False

# Minimal CA+chain filter writer for the binary path (keeps TM-align’s CA convention)
from Bio.PDB import PDBIO, Select  # requires biopython

class _SelectCAChain(Select):
    def __init__(self, chain_id: Optional[str]): self.chain_id = chain_id
    def accept_chain(self, chain):  # keep only requested chain if provided
        return 1 if (self.chain_id is None or chain.id.strip() == str(self.chain_id)) else 0
    def accept_atom(self, atom):    # TM-align uses Cα
        return 1 if atom.get_id() == "CA" else 0

def _prep_pdb_for_binary(pdb_or_path: str, chain: Optional[str]) -> Tuple[str, bool]:
    """
    Return a PDB file path (maybe temp) suitable for the TM-align binary.
    If a chain is provided, write a CA-only, chain-filtered temp PDB using your loader.
    """
    local_pdb, is_temp = _ensure_local_pdb(pdb_or_path)
    if chain is None:
        return local_pdb, is_temp

    struct = load_pdb_structure(local_pdb)
    out = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb").name
    io = PDBIO()
    io.set_structure(struct)
    io.save(out, _SelectCAChain(chain))
    return out, True

# ---------------- engines ---------------- #

def _run_tmalign_binary(pdb1, pdb2, chain1, chain2) -> Dict[str, Any]:
    if not (os.path.isfile(TMALIGN_EXE) and os.access(TMALIGN_EXE, os.X_OK)):
        raise FileNotFoundError(f"TM-align binary not found or not executable: {TMALIGN_EXE}")

    p1, tmp1 = _prep_pdb_for_binary(pdb1, chain1)
    p2, tmp2 = _prep_pdb_for_binary(pdb2, chain2)

    try:
        res = subprocess.run([TMALIGN_EXE, p1, p2], check=True, capture_output=True, text=True)
        out = res.stdout
    finally:
        for path, is_tmp in [(p1, tmp1), (p2, tmp2)]:
            if is_tmp and os.path.exists(path):
                try: os.remove(path)
                except OSError: pass

    # Parse the two TM-scores (normalized by 1 and by 2), plus RMSD & aligned length
    scores = re.findall(r"TM-score\\s*=\\s*([0-9]*\\.?[0-9]+)", out)
    tm1 = float(scores[0]) if len(scores) >= 1 else None
    tm2 = float(scores[1]) if len(scores) >= 2 else None
    rmsd_m = re.search(r"RMSD\\s*=\\s*([0-9]*\\.?[0-9]+)", out)
    rmsd = float(rmsd_m.group(1)) if rmsd_m else None
    ali_m = re.search(r"Aligned length\\s*=\\s*(\\d+)", out)
    aligned_length = int(ali_m.group(1)) if ali_m else None

    if tm1 is None:
        raise RuntimeError("Could not parse TM-score from TMalign output:\\n" + out)

    print(f"[TM-align] Using binary: {TMALIGN_EXE}")
    return {
        "engine": "binary",
        "exe": TMALIGN_EXE,
        "tm_by_1": tm1,
        "tm_by_2": tm2,
        "rmsd": rmsd,
        "aligned_length": aligned_length,
        "raw_output": out,
    }


def _run_tmtools_python(pdb1, pdb2, chain1, chain2) -> dict:
    def _coords_seq_from_any(local_pdb: str, chain: str | None):
        """
        Get (CA_coords, seq) robustly:
        - If we get a true Biotite AtomArray -> use your read_seq_coord_contacts_from_pdb()
        - If we get a NumPy recarray -> derive (coords, seq) from its fields
        """
        atom_array = load_structure_to_atom_array(local_pdb)  # your loader

        # Case 1: true AtomArray (has array_length method) → use your function
        if hasattr(atom_array, "array_length"):
            _, _, seq, _, ca_coords = read_seq_coord_contacts_from_pdb(atom_array, chain=chain)
            return ca_coords, seq

        # Case 2: recarray (no array_length) → build (coords, seq) ourselves
        arr = atom_array  # expected fields: chain_id, res_name, atom_name, coord, res_id, ins_code
        mask = np.ones(len(arr), dtype=bool)
        if chain is not None:
            mask &= (arr.chain_id == chain)

        # keep amino acids we know how to map
        valid_res = np.isin(arr.res_name, np.array(list(aa_long_short.keys())))
        mask &= valid_res

        # CA-only (TM-align convention)
        mask &= (arr.atom_name == "CA")
        sub = arr[mask]

        if len(sub) == 0:
            raise ValueError("No CA atoms found after filtering (check chain id and input).")

        # ensure stable residue order
        order = np.lexsort((sub.ins_code, sub.res_id, sub.chain_id))
        sub = sub[order]

        coords = np.asarray(sub.coord, dtype=float)
        seq = "".join(aa_long_short.get(res, "X") for res in sub.res_name)
        return coords, seq

    def _coords_and_seq(pdb_or_path: str, chain: str | None):
        local_pdb, _ = _ensure_local_pdb(pdb_or_path)
        return _coords_seq_from_any(local_pdb, chain)

    coords1, seq1 = _coords_and_seq(pdb1, chain1)
    coords2, seq2 = _coords_and_seq(pdb2, chain2)

    res = tm_align(coords1, coords2, seq1, seq2)
    print("[TM-align] Using tmtools (Python bindings)")
    return {
        "engine": "tmtools",
        "tm_by_1": float(res.tm_norm_chain1),
        "tm_by_2": float(res.tm_norm_chain2),
        "rmsd": float(res.rmsd),
    }

# ---------------- public API ---------------- #

def tmalign_unified(
    pdb1: str, pdb2: str,
    chain1: Optional[str] = None,
    chain2: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Uses the binary if USE_TMALIGN_BINARY is True; otherwise falls back to tmtools.
    If the binary is selected but fails, falls back to tmtools automatically.
    """
    if USE_TMALIGN_BINARY:
        try:
            return _run_tmalign_binary(pdb1, pdb2, chain1, chain2)
        except Exception as e:
            print(f"[TM-align] Binary failed ({e}); falling back to tmtools...")
    return _run_tmtools_python(pdb1, pdb2, chain1, chain2)

def compute_tmscore_align(
    pdb1: str, pdb2: str,
    chain1: Optional[str] = None,
    chain2: Optional[str] = None,
    norm: str = "by_1"  # "by_1", "by_2", or "max"
) -> float:
    """
    Thin wrapper for pipeline use: returns a single float TM-score.
    """
    res = tmalign_unified(pdb1, pdb2, chain1=chain1, chain2=chain2)
    if norm == "by_1":
        return float(res["tm_by_1"])
    if norm == "by_2" and res.get("tm_by_2") is not None:
        return float(res["tm_by_2"])
    if norm == "max":
        vals = [v for k, v in res.items() if k in ("tm_by_1", "tm_by_2") and v is not None]
        return float(max(vals)) if vals else float(res["tm_by_1"])
    raise ValueError("norm must be one of: 'by_1', 'by_2', 'max'")



def align_cmaps_by_sequence(
    cmap1: np.ndarray, seq1: str,
    cmap2: np.ndarray, seq2: str,
    mode: str = "blosum",       # choose "blosum" by default; or "standard"
    gap_open: int = 11,
    gap_extend: int = 1) -> Tuple[np.ndarray, np.ndarray, Tuple[List[int], List[int]]]:
    """
    Build residue mapping from sequences, then slice both contact maps to matched positions.
    Returns (aligned_cmap1, aligned_cmap2, (idx1, idx2)).
    """
    idx1, idx2 = get_align_indexes(seq1, seq2, mode=mode, gap_open=gap_open, gap_extend=gap_extend)
    if not idx1:
        raise ValueError("No matched residues found from sequence alignment.")
    print("Finished pairwise-sequence alignment with " + mode + " mode; Now to cmaps")

    n1, n2 = cmap1.shape[0], cmap2.shape[0]
    paired = [(i, j) for i, j in zip(idx1, idx2) if i < n1 and j < n2]
    if not paired:
        raise ValueError("Matched indices exceed cmap sizes; ensure sequences match how maps were built (e.g., CA-only).")
    i_sel, j_sel = map(np.array, zip(*paired))

    cm1 = cmap1[np.ix_(i_sel, i_sel)]
    cm2 = cmap2[np.ix_(j_sel, j_sel)]
    return cm1, cm2, (i_sel.tolist(), j_sel.tolist())

