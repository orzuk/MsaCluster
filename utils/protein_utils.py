# Some utilities for proteins and their mutations
# import esm
import string
import shutil

from config import *

# import pcmap
import torch
import torch.nn.functional as F  # for padding
from scipy.spatial.distance import squareform, pdist, cdist

from typing import List, Tuple, Optional, Dict, NamedTuple, Union, Callable
# import matplotlib.pyplot as plt
import Bio
import Bio.PDB
import Bio.SeqRecord
from Bio import SeqIO, PDB, AlignIO

import pickle
import os
import sys
# import urllib
import mdtraj as md
from Bio.PDB import *
# import biotite.structure as bs

from biotite.structure.io.pdb import PDBFile
from biotite.database import rcsb
# from biotite.structure.io.pdbx import get_structure
from biotite.structure import filter_amino_acids, distance, AtomArray
# from biotite.structure.residues import get_residues
from biotite.structure import get_residues

from Bio.PDB.MMCIFParser import MMCIFParser

# Biotite provides a mapping of residue codes to single-letter amino acid codes
# from biotite.sequence.residues import RESIDUE_CODES_3TO1
# from biotite.sequence import ProteinSequence
import numpy as np
import subprocess
import re

#from biotite.structure.io.pdbx import PDBxFile, get_structure


# import iminuit
# Helper function for loading

from Bio.PDB.MMCIFParser import MMCIFParser
from types import SimpleNamespace


# This is an efficient way to delete lowercase characters and insertion characters from a string
deletekeys = dict.fromkeys(string.ascii_lowercase)
deletekeys["."] = None
deletekeys["*"] = None
translation = str.maketrans(deletekeys)


def remove_insertions(sequence: str) -> str:
    """ Removes any insertions into the sequence. Needed to load aligned sequences in an MSA. """
    return sequence.translate(translation)


def read_msa(filename: str) -> List[Tuple[str, str]]:
    """ Reads the sequences from an MSA file, automatically removes insertions."""
    return [(record.description, remove_insertions(str(record.seq))) for record in SeqIO.parse(filename, "fasta")]


def genetic_code():
    """Return the standard genetic code as a dictionary."""
    code = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    return code


aa_long_short = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                 'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
                 'ASX': 'B', 'XLE': 'J', 'PYL': 'O', 'SEC': 'U', 'UNK': 'X', 'GLX': 'Z'}
aa_short_long = {y: x for x, y in aa_long_short.items()}

# Normalize common non-standard residue names to standard 3-letter codes
NONSTD_TO_STD = {"MSE":"MET","SEP":"SER","TPO":"THR","PTR":"TYR",
                 "HSD":"HIS","HSE":"HIS","HSP":"HIS",
                 "MLE":"LEU","MLY":"LYS","M3L":"LYS"}


def aa3_to_1(resname: str) -> str:
    """
    Convert a 3-letter residue name to 1-letter using our tables.
    Handles common non-standard names via NONSTD_TO_STD.
    Returns 'X' if unknown.
    """
    if not resname:
        return "X"
    r = resname.upper().strip()
    r = NONSTD_TO_STD.get(r, r)            # normalize non-standard (e.g., MSE->MET)
    return aa_long_short.get(r, "X")       # MET -> M, etc.

def seq3_to_1(resnames):
    """
    Convert an iterable of 3-letter residue names to a string of 1-letter codes.
    """
    return "".join(aa3_to_1(r) for r in resnames)




# Get all the possible amino acids that we get with a single point mutation
# for a specific codon
def point_mutation_to_aa(codon, genetic_code_dict):
    aa_mut_list = [genetic_code_dict.get(codon)]  # this is the current aa
    for pos in range(3):
        for s in ["A", "C", "G", "T"]:
            aa_mut_list.append(genetic_code_dict.get(codon[:pos] + s + codon[pos + 1:]))
    return list(set(aa_mut_list))


# Get all codons that code for a specific amino-acid
def aa_to_codon(aa, genetic_code_dict):
    return [x for x in genetic_code_dict if genetic_code_dict[x] == aa]


# Get all aa that can be obtained by a single point mutation from a given aa,
# where we don't know the codon
def aa_point_mutation_to_aa(aa, genetic_code_dict):
    return list(
        set([j for l in [point_mutation_to_aa(c, genetic_code_dict) for c in aa_to_codon(aa, genetic_code_dict)] for j
             in l]))



# Extract a sequence from a protein pdb file.
# Next do it by chain (?)

def process_sequence(seq):
    # Remove any whitespace and ensure uppercase
    seq = seq.strip().upper()
    # Replace gaps with the ESM mask token
    seq = seq.replace('-', '')  # ESMfold uses <mask> token for masked positions
    # Remove any other invalid characters
    seq = ''.join(c for c in seq if c in 'ACDEFGHIKLMNPQRSTVWY')
    return seq


def extract_protein_sequence(
    pdb_file: str,
    chain: Optional[str] = None,
    ca_only: bool = False) -> str:
    """
    Return a one-letter AA sequence from a PDB file.

    - If `chain` is given (e.g., 'A'), only that chain is used.
    - If `ca_only=True`, include only residues that have a Cα atom.
      (This makes the sequence length match a CA-based contact map built
       with the same chain and residue filtering.)

    Notes:
      * Non-standard AA names are normalized via NONSTD_TO_STD, then mapped
        with your `aa_long_short` dict; anything unknown becomes 'X'.
      * Only polymer residues are included (res.id[0] == ' ').
      * First model only (common for crystal/NMR single-model files).
    """
    structure = PDBParser(QUIET=True).get_structure("x", pdb_file)
    seq_chars = []
    try:
        model0 = next(iter(structure))
    except StopIteration:
        raise ValueError(f"No models found in PDB: {pdb_file}")

    for ch in model0:
        if chain is not None and ch.id != chain:
            continue
        for res in ch:
            # skip waters/hetero
            if res.id[0] != " ":
                continue
            # keep amino-acid residues (include nonstandard names too)
            if not is_aa(res, standard=False):
                continue
            # if CA-only, require a CA atom
            if ca_only and "CA" not in res:
                continue

            res3 = res.get_resname().upper()
            res3 = NONSTD_TO_STD.get(res3, res3)
            seq_chars.append(aa_long_short.get(res3, "X"))

        # if a specific chain was requested, stop after it
        if chain is not None:
            break

    seq = "".join(seq_chars)
    if not seq:
        raise ValueError(f"No amino-acid residues found for "
            f"{'chain '+chain if chain else 'protein'} in {pdb_file} "
            f"(ca_only={ca_only}).")

    # If you have a post-processor, use it; otherwise return as-is
    try:
        return process_sequence(seq)
    except Exception:
        return seq


def clean_sequence(residue_energies):
    """
    Clean residue energies to extract a valid amino acid sequence.

    Parameters:
    - residue_energies (list of dict): Residue data [{'residue_name': ..., 'residue_index': ..., 'energy': ...}, ...]

    Returns:
    - str: Cleaned amino acid sequence.
    """
    seq = ""
    for res in residue_energies:
        resname = res[0]  # res = ('ALA', 3, 2.5)
        if resname in aa_long_short:
            seq += aa_long_short[resname]
    return seq


def convert_atomarray_to_recarray(atom_array):
    """
    Convert a Biotite AtomArray (which may not be a recarray with a dtype)
    into a proper NumPy recarray with the following fields:
      - chain_id: Unicode string (length 1)
      - res_id: 32-bit integer
      - res_name: Unicode string (length 3)
      - atom_name: Unicode string (length 4)
      - ins_code: Unicode string (length 1)
      - coord: 3-element float32 array
    """
    n = len(atom_array)
    rec_dtype = np.dtype([
        ("chain_id", "U1"),
        ("res_id", "i4"),
        ("res_name", "U3"),
        ("atom_name", "U4"),
        ("ins_code", "U1"),
        ("coord", "f4", (3,))
    ])
    rec = np.empty(n, dtype=rec_dtype)
    rec["chain_id"] = atom_array.chain_id
    rec["res_id"] = atom_array.res_id
    rec["res_name"] = atom_array.res_name
    rec["atom_name"] = atom_array.atom_name
    rec["ins_code"] = atom_array.ins_code
    rec["coord"] = atom_array.coord
    # Convert to a recarray so fields can be accessed as attributes:
    return np.rec.array(rec)


# Taken from esm:
def extend(a, b, c, L, A, D):
    """
    input:  3 coords (a,b,c), (L)ength, (A)ngle, and (D)ihedral
    output: 4th coord
    """

    def normalize(x):
        return x / np.linalg.norm(x, ord=2, axis=-1, keepdims=True)

    bc = normalize(b - c)
    n = normalize(np.cross(b - a, bc))
    m = [bc, np.cross(n, bc), n]
    d = [L * np.cos(A), L * np.sin(A) * np.cos(D), -L * np.sin(A) * np.sin(D)]
    return c + sum([m * d for m, d in zip(m, d)])


def load_seq_and_struct(cur_family_dir, pdbids, pdbchains):
    """
    Load sequence and structure for given PDB IDs and chains using Biopython.
     Parameters:
    - cur_family_dir: directory
    - pdbids: PDB IDs
    - pdbchains: Chain IDs

    Returns:
    - Nothing, saves output to files
    """
    print("Loading seq and struct for " + pdbids[0] + " , " + pdbids[1])
    for fold in range(2):
        if not os.path.exists(cur_family_dir):
            print("Mkdir: " + cur_family_dir)
            os.mkdir(cur_family_dir)

        print(f"Get seq + struct for {pdbids[fold]}, out of {len(pdbids)}")
        cif_file_path = rcsb.fetch(pdbids[fold], "pdb")  # Fetch CIF file path
        print("Fetched pdbid=", pdbids[fold], " into file=", cif_file_path)

        # Fetch the CIF and PDB file and save it locally
        fetch_and_save_pdb_file(pdbids[fold], "pdb", cur_family_dir + "/" + pdbids[fold] + ".pdb")
        fetch_and_save_pdb_file(pdbids[fold], "cif", cur_family_dir + "/" + pdbids[fold] + "_cif.pdb")

        # Load the structure and convert it to a biotite-like AtomArray
        atom_array = load_structure_to_atom_array(cur_family_dir + "/" + pdbids[fold] + ".pdb")

#        struct = PDBFile.read(cur_family_dir + "/" + pdbids[fold] + ".pdb").get_structure(model=1)

        struct = PDBFile.read(cur_family_dir + "/" + pdbids[fold] + ".pdb").get_structure(model=1)
        print("Type of struct:", type(struct))
        try:
            print("struct.dtype:", struct.dtype)
            print("struct.dtype.names:", struct.dtype.names)
        except AttributeError as e:
            print("No dtype or dtype.names attribute:", e)
        print("Length/shape of struct:", len(struct))
        print("First record in struct:", struct[0])
        print("All keys in first record (if any):", getattr(struct[0], "dtype", "No dtype available"))

#        struct = convert_atomarray_to_recarray(struct)
#        print("Converted to recarray with fields:", struct.dtype.names)

        # Ensure that the structure is a recarray with named fields
        # If Biotite AtomArray, use directly
        # Check if Biotite AtomArray (modern version)
        # --- Determine if this is already an AtomArray ---
        print("==== DEBUG: Type of struct:", type(struct))

        if isinstance(struct, AtomArray):
            print("==== DEBUG: struct is a Biotite AtomArray")
            print("Length of AtomArray:", struct.array_length())
            print("First atom:", struct[0])
            print("Calling get_residues on AtomArray...")
            residue_starts, residues = get_residues(struct)
        else:
            print("==== DEBUG: struct is NOT an AtomArray")
            print("Falling back to legacy path")
            try:
                struct = np.array(struct).view(np.recarray)
                print("Successfully viewed as recarray")
                print("Struct dtype names:", struct.dtype.names)
                if struct.dtype.names is not None and all(
                        k in struct.dtype.names for k in ["chain_id", "res_id", "ins_code"]):
                    print("Sorting using lexsort")
                    sort_order = np.lexsort((struct.ins_code, struct.res_id, struct.chain_id))
                    struct = struct[sort_order]
                else:
                    print("Struct missing fields for sorting, skipping sort")
            except Exception as e:
                print("Failed to convert to recarray:", e)
                raise

            print("Calling get_residues on recarray")
            residue_starts, residues = get_residues(struct)
        print("Extracted residues outside loop:", residues, residue_starts)

        # If the extracted residues are numeric (i.e. not the expected three-letter codes),
        # then patch them using the res_name field.
        if residues.dtype.kind in "iuf":
            residue_starts = np.array(residue_starts, dtype=int)
            residues = struct.res_name[residue_starts]
            print("Patched residues (from res_name):", residues)

        # Now call your function that maps residues to one-letter codes
        pdb_seq = "".join(aa_long_short.get(res, "X") for res in residues)
        print("Extracted sequence:", pdb_seq, "len=", len(pdb_seq))

        # Process structure
        pdb_dists, pdb_contacts, pdb_seq, pdb_good_res_inds, cbeta_coord = \
            read_seq_coord_contacts_from_pdb(struct, chain=pdbchains[fold])  # structure
        print("pdb_seq=", pdb_seq, " len=", len(pdb_seq))

        # Save sequence to FASTA file
        fasta_file_name = os.path.join(
            cur_family_dir, f"{pdbids[fold]}{pdbchains[fold]}.fasta"
        )
        with open(fasta_file_name, "w") as text_file:
            text_file.writelines([
                f"> {pdbids[fold].upper()}:{pdbchains[fold].upper()}\n",
                pdb_seq
            ])

        # Save contacts to a binary file
        contact_file = os.path.join(
            cur_family_dir, f"{pdbids[fold]}{pdbchains[fold]}_pdb_contacts.npy"
        )
        np.save(contact_file, pdb_contacts)
        print(f"Saved contacts to: {contact_file}")


def read_seq_coord_contacts_from_pdb(
        structure: AtomArray,
        distance_threshold: float = 8.0,
        chain: Optional[str] = None
) -> Tuple[np.ndarray, np.ndarray, str, np.ndarray, np.ndarray]:
    """
    Extract distances, contacts, sequence, and coordinates from an AtomArray structure.

    Parameters:
    - structure (AtomArray): Biotite structure containing atomic data.
    - distance_threshold (float): Cutoff distance for contacts in Å.
    - chain (Optional[str]): Specific chain to process (if None, all chains are used).

    Returns:
    - dist (np.ndarray): Pairwise distance matrix of Cα atoms.
    - contacts (np.ndarray): Binary contact map (1 for contact, 0 otherwise).
    - pdb_seq (str): Protein sequence as a string.
    - good_res_ids (np.ndarray): Indices of valid residues.
    - CA_coords (np.ndarray): Coordinates of Cα atoms.
    """

    residue_starts, residues  = get_residues(structure)
    print("Extracted residues:", residues)

    # Filter by chain ID if specified
    print("BEFORE IF CHAIN IS", chain)
    if chain is not None:
        structure = structure[structure.chain_id == chain]
    print("Unique chain IDs before filtering:", np.unique(structure.chain_id))
    print("Number of atoms before filtering:", len(structure))

    # Filter amino acids only
    amino_acid_filter = filter_amino_acids(structure)
    structure = structure[amino_acid_filter]

    print("Unique chain IDs:", np.unique(structure.chain_id))
    print("Number of atoms after filtering:", len(structure))

    # Get residues and their starting indices
    residue_starts, residues  = get_residues(structure)

    # Map residues to single-letter codes (fallback to "X" for unknown residues)
    pdb_seq = "".join(aa_long_short.get(res, "X") for res in residues)

    # Extract Cα coordinates
    ca_mask = structure.atom_name == "CA"
    CA_coords = structure.coord[ca_mask]

    if len(CA_coords) == 0:
        raise ValueError("No Cα atoms found in structure.")

    # Extract valid residue indices
    good_res_ids = residue_starts

    # Calculate pairwise distances for Cα atoms
    dist = squareform(pdist(CA_coords))

    # Create binary contact map based on distance threshold
    contacts = (dist < distance_threshold).astype(int)

    return dist, contacts, pdb_seq, good_res_ids, CA_coords


class CustomAtomArray:
    def __init__(self, chain_id, atom_name, coord, res_id, res_name, ins_code):
        self.chain_id = chain_id  # numpy array of chain IDs (strings)
        self.atom_name = atom_name  # numpy array of atom names (strings)
        self.coord = coord  # numpy array of coordinates (n_atoms x 3)
        self.res_id = res_id  # numpy array of residue IDs
        self.res_name = res_name  # numpy array of residue names (strings)
        self.ins_code = ins_code  # numpy array of insertion codes (strings)

    def __getitem__(self, key):
        # Enable boolean indexing/filtering, e.g.:
        # filtered_array = atom_array[atom_array.chain_id == chain]
        return AtomArray(
            chain_id=self.chain_id[key],
            atom_name=self.atom_name[key],
            coord=self.coord[key],
            res_id=self.res_id[key],
            res_name=self.res_name[key],
            ins_code=self.ins_code[key]
        )

    def __len__(self):
        return len(self.coord)


###############################################################
# Function to convert a Biopython Structure to an AtomArray    #
###############################################################

# --- Function to fetch and save a file from RCSB ---
def fetch_and_save_pdb_file(pdb_id, file_format, local_filename):
    """
    Fetch the file for the given PDB ID and format ('pdb' or 'cif')
    using Biotite's fetch, and save it locally.

    This function handles both cases where the returned object is a
    StringIO (for cif) or a file path (for pdb).
    """
    file_obj = rcsb.fetch(pdb_id, file_format)
    try:
        # Try to get the content if file_obj is a StringIO-like object
        content = file_obj.getvalue()
    except AttributeError:
        # Otherwise assume it's a file path and read from it
        with open(file_obj, "r") as f:
            content = f.read()
    with open(local_filename, "w") as f:
        f.write(content)
    return local_filename


def load_structure_to_atom_array(local_pdb_file):
    # Use Biopython’s PDBParser (which usually retains the residue names correctly)
    parser = PDBParser(QUIET=True)
    bio_structure = parser.get_structure("structure", local_pdb_file)
    # Convert the Biopython structure into a Biotite AtomArray (structured NumPy array)
    atom_array = biopython_to_biotite_atom_array(bio_structure)
    return atom_array.view(np.recarray)


#def load_structure_to_atom_array(local_cif_file):
#    from Bio.PDB import MMCIFParser
#    parser = MMCIFParser(QUIET=True)
#    bio_structure = parser.get_structure("structure", local_cif_file)
#    atom_array = biopython_to_biotite_atom_array(bio_structure)
#    return atom_array

def biopython_to_biotite_atom_array(bio_structure):
    """
    Convert a Biopython Structure object to a proper Biotite AtomArray (a structured NumPy array)
    with the fields required by Biotite's functions, including get_residues.

    The returned AtomArray will have the following fields:
      - chain_id : a Unicode string of length 1 (e.g. "A")
      - res_id   : a 32-bit integer (residue sequence number)
      - res_name : a Unicode string of length 3 (three-letter residue code, e.g. "GLY")
      - atom_name: a Unicode string of length 4 (e.g. "CA")
      - ins_code : a Unicode string of length 1 (insertion code)
      - coord    : a 3-element float32 array (coordinates)
    """
    chain_ids = []
    atom_names = []
    coords = []
    res_ids = []
    res_names = []
    ins_codes = []

    for atom in bio_structure.get_atoms():
        # Get chain id: go from atom -> residue -> chain
        chain_ids.append(atom.get_parent().get_parent().get_id())
        # Get atom name (e.g. "CA", "CB", etc.)
        atom_names.append(atom.get_name())
        # Get coordinates (3-element vector)
        coords.append(atom.get_coord())
        # Get residue information
        res_id_tuple = atom.get_parent().get_id()  # typically (hetfield, resseq, icode)
        res_ids.append(res_id_tuple[1])
        # Normalize the residue name: strip whitespace and convert to uppercase
        res_names.append(atom.get_parent().get_resname().strip().upper())
        ins_codes.append(res_id_tuple[2])

    # Convert lists to numpy arrays with proper types
    chain_ids = np.array(chain_ids, dtype="U1")
    atom_names = np.array(atom_names, dtype="U4")
    coords = np.array(coords, dtype="f4")
    res_ids = np.array(res_ids, dtype="i4")
    res_names = np.array(res_names, dtype="U3")
    ins_codes = np.array(ins_codes, dtype="U1")

    # Define the structured dtype for a Biotite AtomArray
    atom_dtype = np.dtype([
        ("chain_id", "U1"),
        ("res_id", "i4"),
        ("res_name", "U3"),
        ("atom_name", "U4"),
        ("ins_code", "U1"),
        ("coord", "f4", (3,))
    ])

    n_atoms = len(chain_ids)
    atom_array = np.empty(n_atoms, dtype=atom_dtype)
    atom_array["chain_id"] = chain_ids
    atom_array["atom_name"] = atom_names
    atom_array["coord"] = coords
    atom_array["res_id"] = res_ids
    atom_array["res_name"] = res_names
    atom_array["ins_code"] = ins_codes

    return atom_array


def load_pdb_structure(pdb_file):
    """
    Load PDB file with specific handling for NMR structures
    """
    import warnings
    from Bio.PDB import PDBParser
    from Bio.PDB.PDBExceptions import PDBConstructionWarning

    # Suppress BiopythonWarning about PDB format
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', PDBConstructionWarning)

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('PDB', pdb_file)

        return structure


# Extract contact map from a pdb-file
# Also extract the distances themselves (more informative than the thresholded contacts)
# And the sequences
#
# Input:
# structure - atom array of biotite
# distance_threshold - threshold on distance for being considered as a contact
# chain - optional
#
# Output:
# dist - pairwise distance matrix between cBeta atoms
# contacts - pairwise binary contacts matrix for distances < threshold
# pdb_seq - sequence extracted from pdb file, after removing residues with missing atoms
# good_res_ids - indices of full good residues
# Cbeta - coordinates of Cbeta atoms (N*#chains numpy array)



# Evaluate precision of predicted contacts with respect to true contacts
# Input:
# predictions - matrix of predicted residue affinities
# targets - matrix of binary contacts
# Output :
# AUC - Area under the ROC curve for predicting contacts
# P@L - percent of top L contacts recovered among
def compute_precisions(
        predictions: torch.Tensor,
        targets: torch.Tensor,
        src_lengths: Optional[torch.Tensor] = None,
        minsep: int = 6,
        maxsep: Optional[int] = None,
        override_length: Optional[int] = None,  # for casp
):
    if isinstance(predictions, np.ndarray):
        predictions = torch.from_numpy(predictions)
    if isinstance(targets, np.ndarray):
        targets = torch.from_numpy(targets)
    if predictions.dim() == 2:
        predictions = predictions.unsqueeze(0)
    if targets.dim() == 2:
        targets = targets.unsqueeze(0)
    override_length = (targets[0, 0] >= 0).sum()

    # Check sizes
    if predictions.size() != targets.size():
        raise ValueError(
            f"Size mismatch. Received predictions of size {predictions.size()}, "
            f"targets of size {targets.size()}"
        )
    device = predictions.device

    batch_size, seqlen, _ = predictions.size()
    seqlen_range = torch.arange(seqlen, device=device)

    sep = seqlen_range.unsqueeze(0) - seqlen_range.unsqueeze(1)
    sep = sep.unsqueeze(0)
    valid_mask = sep >= minsep
    valid_mask = valid_mask & (targets >= 0)  # negative targets are invalid

    if maxsep is not None:
        valid_mask &= sep < maxsep

    if src_lengths is not None:
        valid = seqlen_range.unsqueeze(0) < src_lengths.unsqueeze(1)
        valid_mask &= valid.unsqueeze(1) & valid.unsqueeze(2)
    else:
        src_lengths = torch.full([batch_size], seqlen, device=device, dtype=torch.long)

    predictions = predictions.masked_fill(~valid_mask, float("-inf"))

    x_ind, y_ind = np.triu_indices(seqlen, minsep)
    predictions_upper = predictions[:, x_ind, y_ind]
    targets_upper = targets[:, x_ind, y_ind]

    topk = seqlen if override_length is None else max(seqlen, override_length)
    indices = predictions_upper.argsort(dim=-1, descending=True)[:, :topk]
    topk_targets = targets_upper[torch.arange(batch_size).unsqueeze(1), indices]

#    print("TOPK:")
#    print(topk_targets)
#    print(type(topk_targets))

#    print(topk_targets.size(1))
#    print(topk)
    if topk_targets.size(1) < topk:  # what is F???
        topk_targets = F.pad(topk_targets, [0, topk - topk_targets.size(1)])

    cumulative_dist = topk_targets.type_as(predictions).cumsum(-1)

    gather_lengths = src_lengths.unsqueeze(1)
    if override_length is not None:
        gather_lengths = override_length * torch.ones_like(
            gather_lengths, device=device
        )

    gather_indices = (
                             torch.arange(0.1, 1.1, 0.1, device=device).unsqueeze(0) * gather_lengths
                     ).type(torch.long) - 1

    binned_cumulative_dist = cumulative_dist.gather(1, gather_indices)
    binned_precisions = binned_cumulative_dist / (gather_indices + 1).type_as(
        binned_cumulative_dist)

    pl5 = binned_precisions[:, 1]
    pl2 = binned_precisions[:, 4]
    pl = binned_precisions[:, 9]
    auc = binned_precisions.mean(-1)

    return {"AUC": auc, "P@L": pl, "P@L2": pl2, "P@L5": pl5}


# General utility for dictionary of unique values
def unique_values_dict(original_dict):
    # Invert the dictionary. This will discard duplicate values.
    inverted_dict = {v: k for k, v in original_dict.items()}

    # Invert the dictionary again to get unique values.
    unique_dict = {v: k for k, v in inverted_dict.items()}

    return unique_dict


# Score the predictions
def evaluate_prediction(
        predictions: torch.Tensor,
        targets: torch.Tensor,
) -> Dict[str, float]:
    if isinstance(predictions, np.ndarray):
        predictions = torch.from_numpy(predictions)
    if isinstance(targets, np.ndarray):
        targets = torch.from_numpy(targets)
    contact_ranges = [
        ("local", 3, 6),
        ("short", 6, 12),
        ("medium", 12, 24),
        ("long", 24, None),
    ]
    metrics = {}
    targets = targets.to(predictions.device)
    for name, minsep, maxsep in contact_ranges:
        rangemetrics = compute_precisions(
            predictions,
            targets,
            minsep=minsep,
            maxsep=maxsep,
        )
        for key, val in rangemetrics.items():
            metrics[f"{name}_{key}"] = val.item()
    return metrics


# Select sequences from the MSA to maximize the hamming distance
# Alternatively, can use hhfilter
def greedy_select(msa: List[Tuple[str, str]], num_seqs: int, mode: str = "max") -> List[Tuple[str, str]]:
    assert mode in ("max", "min")
    if len(msa) <= num_seqs:
        return msa

    array = np.array([list(seq) for _, seq in msa], dtype=np.bytes_).view(np.uint8)

    optfunc = np.argmax if mode == "max" else np.argmin
    all_indices = np.arange(len(msa))
    indices = [0]
    pairwise_distances = np.zeros((0, len(msa)))
    for _ in range(num_seqs - 1):
        dist = cdist(array[indices[-1:]], array, "hamming")
        pairwise_distances = np.concatenate([pairwise_distances, dist])
        shifted_distance = np.delete(pairwise_distances, indices, axis=1).mean(0)
        shifted_index = optfunc(shifted_distance)
        index = np.delete(all_indices, indices)[shifted_index]
        indices.append(index)
    indices = sorted(indices)
    return [msa[idx] for idx in indices]



# Functions below are from András Aszódi:
# https://stackoverflow.com/questions/10324674/parsing-a-pdb-file-in-python
def read_pdb(pdbcode, pdbfilenm):
    """
    Read a PDB structure from a file.
    :param pdbcode: A PDB ID string
    :param pdbfilenm: The PDB file
    :return: a Bio.PDB.Structure object or None if something went wrong
    """
    try:
        pdbparser = Bio.PDB.PDBParser(QUIET=True)  # suppress PDBConstructionWarning
        struct = pdbparser.get_structure(pdbcode, pdbfilenm)
        return struct
    except Exception as err:
        print(str(err), file=sys.stderr)
        return None


# Match between cmaps, get only aligned indices
def get_matching_indices_two_cmaps(pairwise_alignment, true_cmap, pred_cmap):
    """
       Match between cmaps, get only aligned indices

       Parameters:
       - pairwise_alignment: Object with alignment of two sequences .
       - true_cmap : Matrix representing first contact map .
       - pred_cmap : Matrix representing second contact map .
       """
    #    n_true = len(true_cmap)  # always 2 !!
    #    n_pred = len(pred_cmap)  # variable number !!

    print("Pairwise alignment: ")
    print(pairwise_alignment)
#    print("Cmap sizes: ", true_cmap.shape, pred_cmap.shape)
    match_true_cmap = {}  # [None]*2
    match_pred_cmap = {}  # [None]*n_pred

    # good_inds = np.minimum(pairwise_alignment[0].indices[0], pairwise_alignment[0].indices[1])
    good_inds = np.where(np.minimum(pairwise_alignment[0].indices[0], pairwise_alignment[0].indices[1]) >= 0)[0]

    ctr = 0
    for fold in true_cmap.keys():  # get true (these are dictionaries !!)
        match_true_cmap[fold] = true_cmap[fold][np.ix_(pairwise_alignment[0].indices[ctr][good_inds],
                                                       pairwise_alignment[0].indices[ctr][good_inds])]
        ctr = ctr + 1
    #        cur_ind = pairwise_alignment[0].indices[i][pairwise_alignment[0].indices[i] >= 0]
    #        print(true_cmap[i])
    #        print(true_cmap[i].shape)
    #        print(true_cmap[i][cur_ind,cur_ind])

    ctr = 0
    for fold in pred_cmap.keys():  # range(n_pred):  # get predicted
        match_pred_cmap[fold] = pred_cmap[fold][np.ix_(pairwise_alignment[0].indices[ctr][good_inds],
                                                       pairwise_alignment[0].indices[ctr][good_inds])]
    return match_true_cmap, match_pred_cmap


def Calculate_RMSD(structure_1: str, structure_2: str, structure_1_index: List[int], structure_2_index: List[int]) -> int:
    """
    calculate the RMSD between two structures using MDtraj library
    this script will fail if mdtraj is not loaded in your python environment
    recommend python 3.10
    """

    #with warnings.catch_warnings(action="ignore"):
    #    turn_off_warnings()
    #load structure information in mdtraj
    pdb = md.load(structure_1)
    pdb_ca = pdb.atom_slice(structure_1_index) #select only CA atoms

    #load structure information in mdtraj
    reference = md.load(structure_2)
    reference_ca = reference.atom_slice(structure_2_index) #select only CA atoms

    # Calculate RMSD of CA atoms
    pdb_ca.superpose(reference_ca)
    return md.rmsd(pdb_ca, reference_ca, frame=0)



run_example = False
if run_example:
    # Example usage:
    print("Hello")
    genetic_code_dict = genetic_code()
    aa_list = list(set(genetic_code_dict.values()))  # all possible amino-acids
    aa_list.remove("*")
    codon_list = list(genetic_code_dict)
    codon = 'GGA'
    amino_acid = genetic_code_dict.get(codon, 'Unknown')
    print(f'The amino acid corresponding to {codon} is {amino_acid}')

    for c in codon_list:
        print(c, genetic_code_dict.get(c), point_mutation_to_aa(c, genetic_code_dict))

    for aa in aa_list:
        print(aa, aa_point_mutation_to_aa(aa, genetic_code_dict))

    # Methods for phylogentic reconstruction
    # constructor = DistanceTreeConstructor()
    # tree = constructor.nj(dm)

# ---------- helpers ----------
def normalize_pair(pair: Union[Tuple[str, str], str]) -> Tuple[str, str, str]:
    """
    Accept ('A','B') or 'A_B' and return (A, B, 'A_B').
    """
    if isinstance(pair, tuple):
        a, b = pair
    else:
        s = str(pair)
        if "_" not in s:
            raise ValueError(f"Pair string must look like 'A_B', got: {s}")
        a, b = s.split("_", 1)
    pair_a = str(a)
    pair_b = str(b)
    pair_id = f"{pair_a}_{pair_b}"
    return pair_a, pair_b, pair_id
