from config import *
import argparse
from Bio.SeqUtils import seq1
from Bio.PDB import PDBParser, PDBIO, Select
from Bio import SeqIO
import mdtraj as md
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from contact_map import ContactFrequency, ContactDifference
from pathlib import Path
from typing import List, Tuple, Union
import sys

import numpy as np
import pandas as pd

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
    dist_mat = np.sqrt(np.einsum("ijk,ijk->ij", diffs, diffs, optimize=True), dtype=float)
    if include_diag:
        np.fill_diagonal(dist_mat, 0.0)
    cmap = (dist_mat <= float(cutoff_A)).astype(np.uint8)
    return cmap, dist_mat


class ChainSelect(Select):
    def __init__(self, chain_letter):
        self.chain_letter = chain_letter

    def accept_chain(self, chain):
        return chain.id == self.chain_letter


def create_chain_pdb_files(fold_1,fold_2, pdb_file_path, chain_pdb_file_path):
    chain_fold_1 = fold_1[-1]
    chain_fold_2 = fold_2[-1]
    # Load the original PDB file
    parser = PDBParser()
    structure_fold1 = parser.get_structure('PDB_structure' , f'{pdb_file_path}/{fold_1[:-1]}.pdb')
    structure_fold2 = parser.get_structure('PDB_structure' , f'{pdb_file_path}/{fold_2[:-1]}.pdb')
    io = PDBIO()
    # Set the structure for saving and use ChainSelect to filter the chain
    io.set_structure(structure_fold1)
    io.save(f'./{chain_pdb_file_path}/{fold_1}.pdb', ChainSelect(chain_fold_1))
    io.set_structure(structure_fold2)
    io.save(f'./{chain_pdb_file_path}/{fold_2}.pdb', ChainSelect(chain_fold_2))


# Extract the protein sequence from a pdb file
# Write the sequence into a fasta file
def get_fasta_chain_seq(pdb_file, fold_name, output_dir):
    # Create a PDB parser
    parser = PDBParser()

    print("Load parser for pdb file:", pdb_file)
    # Parse the structure
    structure = parser.get_structure('PDB_structure', pdb_file)
    # Extract the sequence
    for model in structure:
        for chain in model:
            sequence = ""
            for residue in chain:
                if residue.id[0] == ' ':
                    sequence += seq1(residue.resname)

    with open(f"./{output_dir}/fasta_chain_files/{fold_name}.fasta", "w") as output_handle:
        output_handle.write('>'+'\n'+sequence)



def save_org_cmaps(chain_pdb_file_path,fold):
    traj = md.load(filename_or_filenames=f'{chain_pdb_file_path}/{fold}.pdb')
    frame_contacts = ContactFrequency(traj[0])
    np.save(f'./Pipeline/org_cmaps/{fold}.npy',frame_contacts)


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


def ensure_dir(p: str) -> None:
    Path(p).mkdir(parents=True, exist_ok=True)


# For pipeline running:
def write_pair_pipeline_script(pair_id: str, args: argparse.Namespace) -> str:
    """
    Create Pipeline/<pair>/jobs/msaclust_pipeline_<pair>.sh that runs
    all steps sequentially in one Slurm job, logging stages.
    Uses absolute repo & python paths so Slurm CWD never matters.
    """
    # Absolute paths (no guessing inside the SBATCH job)
    repo_dir = os.path.abspath(os.path.dirname(__file__))     # repo root (dir of this file)
    py_exec  = sys.executable                                 # same python as caller

    pair_dir = f"Pipeline/{pair_id}"
    jobs_dir = os.path.join(pair_dir, "jobs")
    logs_dir = os.path.join(pair_dir, "logs")
    ensure_dir(jobs_dir); ensure_dir(logs_dir)

    script_path = os.path.join(jobs_dir, f"msaclust_pipeline_{pair_id}.sh")
    log_path    = os.path.join(logs_dir, f"msaclust_pipeline_{pair_id}.out")
    status_path = os.path.join(logs_dir, f"msaclust_pipeline_{pair_id}.status")

    cluster_sample_n = getattr(args, "cluster_sample_n", 10) or 10
    esm_device       = getattr(args, "esm_device", "auto") or "auto"

    sb = f"""#!/bin/bash
#SBATCH --job-name=pipe_{pair_id}
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --gres=gpu:a100:1
# #SBATCH --partition=a100
# #SBATCH --constraint=a100
#SBATCH -o {log_path}

set -euo pipefail

# --- absolute locations (filled in by Python) ---
REPO_DIR="{repo_dir}"
PY="{py_exec}"

cd "${{REPO_DIR}}"

PAIR_ID="{pair_id}"
PAIR_DIR="Pipeline/{pair_id}"
LOG="{log_path}"
STATUS="{status_path}"

mkdir -p "$(dirname "${{LOG}}")"
mkdir -p "${{PAIR_DIR}}/output_get_msa" "${{PAIR_DIR}}/output_msa_cluster" \\
         "${{PAIR_DIR}}/output_cmaps/msa_transformer" "${{PAIR_DIR}}/output_AF/AF2" \\
         "${{PAIR_DIR}}/output_esm_fold/esm2" "${{PAIR_DIR}}/output_esm_fold/esm3" \\
         "${{PAIR_DIR}}/fasta_chain_files"

# optional: activate env (export this before sbatch if you need it)
#   export COLABFOLD_ENV_ACTIVATE=/path/to/colabfold_env/bin/activate
if [[ -n "${{COLABFOLD_ENV_ACTIVATE:-}}" && -f "${{COLABFOLD_ENV_ACTIVATE}}" ]]; then
  # shellcheck disable=SC1090
  source "${{COLABFOLD_ENV_ACTIVATE}}"
fi

# capture all output to the log as well
exec > >(stdbuf -oL tee -a "${{LOG}}") 2>&1

# --- status helpers ---
STAGE=""
mark() {{  # $1=stage $2=STATUS (START|OK|FAIL)
  local now; now="$(date '+%Y-%m-%d %H:%M:%S')"
  echo "${{now}} [{pair_id}] $1 $2" | tee -a "${{STATUS}}"
  if [[ "$2" == "OK" ]]; then touch "${{PAIR_DIR}}/.stage_$1.ok"; fi
  if [[ "$2" == "FAIL" ]]; then touch "${{PAIR_DIR}}/.stage_$1.fail"; fi
}}
start() {{ STAGE="$1"; mark "$1" "START"; }}
ok()    {{ mark "$STAGE" "OK"; }}
fail()  {{ mark "$STAGE" "FAIL"; }}

trap 'fail; exit 1' ERR

CLUSTER_SAMPLE_N="{cluster_sample_n}"
ESM_DEVICE="{esm_device}"

# ----------------- STAGES -----------------

start load
"$PY" "{os.path.join(repo_dir, 'run_foldswitch_pipeline.py')}" --run_mode load --foldpair_ids "${{PAIR_ID}}" --run_job_mode inline
ok

start get_msa
"$PY" "{os.path.join(repo_dir, 'run_foldswitch_pipeline.py')}" --run_mode get_msa --foldpair_ids "${{PAIR_ID}}" --run_job_mode inline
ok

start cluster_msa
"$PY" "{os.path.join(repo_dir, 'run_foldswitch_pipeline.py')}" --run_mode cluster_msa --foldpair_ids "${{PAIR_ID}}" --run_job_mode inline
ok

start run_AF
"$PY" "{os.path.join(repo_dir, 'run_foldswitch_pipeline.py')}" --run_mode run_AF --foldpair_ids "${{PAIR_ID}}" --run_job_mode inline
ok

start run_cmap_esm
"$PY" "{os.path.join(repo_dir, 'run_foldswitch_pipeline.py')}" --run_mode run_cmap_esm --foldpair_ids "${{PAIR_ID}}" --run_job_mode inline
ok

start run_esmfold_esm2
"$PY" "{os.path.join(repo_dir, 'run_foldswitch_pipeline.py')}" --run_mode run_esmfold --foldpair_ids "${{PAIR_ID}}" --run_job_mode inline --esm_model esm2 --cluster_sample_n "${{CLUSTER_SAMPLE_N}}" --esm_device "${{ESM_DEVICE}}"
ok

start run_esmfold_esm3
"$PY" "{os.path.join(repo_dir, 'run_foldswitch_pipeline.py')}" --run_mode run_esmfold --foldpair_ids "${{PAIR_ID}}" --run_job_mode inline --esm_model esm3 --cluster_sample_n "${{CLUSTER_SAMPLE_N}}" --esm_device "${{ESM_DEVICE}}"
ok

start plot_pair
if "$PY" "{os.path.join(repo_dir, 'run_foldswitch_pipeline.py')}" --run_mode plot --foldpair_ids "${{PAIR_ID}}" --run_job_mode inline; then
  ok
else
  echo "[plot] skipped or failed (PyMOL not available?)"
  touch "${{PAIR_DIR}}/.stage_plot_pair.skipped"
fi

echo "ALL DONE for {pair_id}"
"""

    with open(script_path, "w") as f:
        f.write(sb)
    os.chmod(script_path, 0o755)
    return script_path

