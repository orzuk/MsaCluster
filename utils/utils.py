from config import *
import argparse
import sys, warnings, re

if any(a in ("-h", "--help") for a in sys.argv) or ("--run_mode" in sys.argv and "help" in sys.argv):
    warnings.filterwarnings("ignore", category=DeprecationWarning, module=r"^Bio\.pairwise2$")

from Bio.PDB import PDBParser, PDBIO, Select
from Bio import SeqIO
import matplotlib.pyplot as plt
from pathlib import Path
from typing import List, Tuple, Union

import numpy as np
import pandas as pd
import time

def now():
    return time.strftime("%Y-%m-%d %H:%M:%S")


def str2bool(v):
    if isinstance(v, bool): return v
    s = str(v).strip().lower()
    if s in {"true", "t", "1", "yes", "y"}: return True
    if s in {"false", "f", "0", "no", "n"}: return False
    raise argparse.ArgumentTypeError("expected boolean (true/false)")


def load_csv_or_none(path, **kwargs):
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        print(f"[skip] missing CSV: {p}", flush=True)
        return None
    try:
        df = pd.read_csv(p, **kwargs)
        print(f"[ok] loaded: {p}  rows={len(df)}", flush=True)
        return df
    except Exception as e:
        print(f"[skip] failed to read {p}: {e}", flush=True)
        return None

def pdb_to_contact_map(
    pdb_file: str,
    chain: str | None = None,
    cutoff_A: float = 5.0,   # keep your current ~5 Å (0.5 nm) default
    include_diag: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """
    Return (binary_contact_map, distance_matrix) using Cα coordinates.

    - pdb_file : path to a PDB file
    - chain    : if given, use only this chain (e.g., 'A'); else use all chains in model 0
    - cutoff_A : threshold in Å for the binary contact map
    - include_diag : if True, distance[i,i]=0 and CM diagonal = 1

    Shapes: (N_res, N_res), dtype: CM=uint8, D=float64
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("x", pdb_file)

    # collect Cα coords (Å)
    coords = []
    model0 = next(iter(structure))  # first model
    for ch in model0:
        if chain is not None and ch.id != chain:
            continue
        for res in ch:
            if res.id[0] == ' ' and 'CA' in res:  # protein residue, has CA
                coords.append(res['CA'].coord)

    coords = np.asarray(coords, dtype=float)
    if coords.size == 0:
        raise ValueError(f"No Cα atoms found (file={pdb_file}, chain={chain!r}).")

    # pairwise distances (Å); vectorized
    diffs = coords[:, None, :] - coords[None, :, :]
    dist_mat = np.sqrt(np.einsum("ijk,ijk->ij", diffs, diffs, optimize=True))
    if include_diag:
        np.fill_diagonal(dist_mat, 0.0)
    cmap = (dist_mat <= float(cutoff_A)).astype(np.uint8)
    return cmap, dist_mat

_NUM = re.compile(r"^\s*([+-]?\d+(?:\.\d+)?)")

def save_org_cmaps(chain_pdb_file_path,fold):
    import mdtraj as md
    from contact_map import ContactFrequency
    traj = md.load(filename_or_filenames=f'{chain_pdb_file_path}/{fold}.pdb')
    frame_contacts = ContactFrequency(traj[0])
    np.save(f'./Pipeline/FoldPairs/org_cmaps/{fold}.npy',frame_contacts)


def find_max_keys(input_dict):
    """
    Given a dictionary of the format {key: {subkey1: value1, subkey2: value2, ...}},
    find the keys that maximize the values for each subkey.

    Parameters:
        input_dict (dict): Dictionary of dictionaries with numeric values.

    Returns:
        dict: A dictionary with subkeys as keys and tuples as values, where each tuple contains
              (key_with_max_value, max_value).
    """
    result = {}

    # Check the subkeys present in the first inner dictionary
    first_key = next(iter(input_dict))
    subkeys = input_dict[first_key].keys()

    for subkey in subkeys:
        max_key = None
        max_value = float('-inf')
        for key, subdict in input_dict.items():
            if subdict[subkey] > max_value:
                max_key = key
                max_value = subdict[subkey]
        result[subkey] = (max_key, max_value)

    return result


def pair_str_to_tuple(s: str) -> tuple[str, str]:
    """
    '1dzlA_5keqF' -> ('1dzlA','5keqF')
    """
    if isinstance(s, tuple) or isinstance(s, list): # ready correct format
        print("Return list/tuple!", s, type(s), len(s), flush=True)
        return s
    s = s.strip()
    if "_" not in s:
        raise ValueError(f"Pair must look like PDBidChain_PDBidChain, got: {s}")
    a, b = s.split("_", 1)
    return a.strip(), b.strip()


def list_protein_pairs(parsed: bool = True, sort_result: bool = True) -> (
        List)[Union[str, Tuple[str, str]]]:
    """
    Return all pair directories under `root`.
    A valid pair dir looks like: '4rmbB_4rmbA' (PDBid+chain)_(PDBid+chain).

    Parameters
    ----------
    parsed : bool
        If True -> return [('4rmbB','4rmbA'), ...]
        If False -> return ['4rmbB_4rmbA', ...]
    sort_result : bool
        If True, sort lexicographically by folder name.

    Returns
    -------
    list
        List of pairs (tuples or strings).
    """
    root = Path(DATA_DIR)
    pairs = []
    for child in root.iterdir():
        if not child.is_dir():
            continue
        m = PAIR_DIR_RE.match(child.name)
        if not m:
            continue
        if parsed:
            pairs.append((m.group(1), m.group(2)))
        else:
            pairs.append(child.name)

    if sort_result:
        pairs.sort(key=lambda x: "_".join(x) if isinstance(x, tuple) else x)
    return pairs


def list_control_proteins(control_type: str = "PorterIDP") -> List[str]:
    """
    Return sorted list of control protein directory names under
    Pipeline/Controls/<control_type>/.
    """
    root = Path(CONTROLS_BASE_DIR) / control_type
    if not root.is_dir():
        return []
    return sorted(d.name for d in root.iterdir() if d.is_dir())


def ensure_dir(p: str) -> None:
    Path(p).mkdir(parents=True, exist_ok=True)

def pick(d: pd.DataFrame, colnames):
    """Return the first present column as numeric series, else NaNs."""
    for c in colnames:
        if c in d.columns:
            return pd.to_numeric(d[c], errors="coerce")
    return pd.Series(index=d.index, dtype=float)
