# Some utilities for proteins and their mutations
import pandas as pd
import esm
# import pcmap
import torch
from scipy.spatial.distance import squareform, pdist, cdist
import numpy as np
from typing import List, Tuple, Optional, Dict, NamedTuple, Union, Callable
import matplotlib as mpl
import Bio
import Bio.PDB
import Bio.SeqRecord
import os
import sys
import urllib

import biotite.structure as bs
from biotite.structure.io.pdbx import PDBxFile, get_structure
from biotite.database import rcsb

import iminuit
import tmscoring   # for comparing structures
from Bio import Phylo  # for phylogenetic trees
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
# from TreeConstruction import DistanceTreeConstructor


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


# Get all the possible amino acids that we get with a single point mutation
# for a specific codon
def point_mutation_to_aa(codon, genetic_code_dict):
    aa_mut_list = [genetic_code_dict.get(codon)] # this is the current aa
    for pos in range(3):
        for s in ["A", "C", "G", "T"]:
            aa_mut_list.append(genetic_code_dict.get(codon[:pos] + s + codon[pos+1:]))
    return list(set(aa_mut_list))


# Get all codons that code for a specific amino-acid
def aa_to_codon(aa, genetic_code_dict):
    return [x for x in genetic_code_dict if genetic_code_dict[x] == aa]


# Get all aa that can be obtained by a single point mutation from a given aa,
# where we don't know the codon
def aa_point_mutation_to_aa(aa, genetic_code_dict):
    return list(set([j for l in [point_mutation_to_aa(c, genetic_code_dict) for c in aa_to_codon(aa, genetic_code_dict)] for j in l]))


# Run design on every position, collect
def design_every_position(S, pdbID):
    n = len(S) # get number of amino-acids
    S_des = ['']*n # The output designed sequence
    for i in range(n):  # loop on positions
        cur_S = run_design(S, pdbID, i)  # run design function using ProteinMPNN, keeping all positions fixed except i
        S_des[i] = cur_S[i]  # change the i-th letter of the i-th protein
    S_des = "".join(S_des) # convert to string
    return S_des


# Run design to both A and B, then fold both of them with alpha-fold. See if we get something closer
def compare_designs(S, pdbID1, pdbID2):
    S1 = design_every_position(S, pdbID1)
    S2 = design_every_position(S, pdbID2)

    # Compare sequences:
    print("\nS: ", S, "\nS1: ", S1, "\nS2: ", S2)

    AF = Alphafold(S)  # Fold the natural sequence
    AF1 = Alphafold(S1)  # Fold the first one
    AF2 = Alphafold(S2)  # Fold the second one

    TM = ['TM1', 'TM2']
    SEQ = ['S', 'S1', 'S2']
    df_tm = pd.DataFrame(columns=TM, index=SEQ)
    S_true_list = [pdbID1, pdbID2]
    S_pred_list = [AF, AF1, AF2]
    for i_true in range(2):   # loop on the two true structures
        for j_pred in range(3):  # S_predicted in [AF, AF1, AF2]:  # loop on the three predicted structures
#            df_tm[TM[i_true]][SEQ[j_pred]] = TMScore(S_true_list[i_true], S_pred_list[j_pred])  # Compute TMScore similarity
            alignment = tmscoring.TMscoring(S_true_list[i_true], S_pred_list[j_pred]) #  'structure1.pdb', 'structure2.pdb')  # from installed tmscoring
            df_tm[TM[i_true]][SEQ[j_pred]] = alignment.tmscore(**alignment.get_current_values())

    print(df_tm)
    return df_tm, S, S1, S2, AF, AF1, AF2  # Return all sequences, structures and their similarity


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


def contacts_from_pdb(
        structure: bs.AtomArray,
        distance_threshold: float = 8.0,
        chain: Optional[str] = None,
) -> np.ndarray:
    mask = ~structure.hetero
    if chain is not None:
        mask &= structure.chain_id == chain

    N = structure.coord[mask & (structure.atom_name == "N")]
    CA = structure.coord[mask & (structure.atom_name == "CA")]
    C = structure.coord[mask & (structure.atom_name == "C")]

    Cbeta = extend(C, N, CA, 1.522, 1.927, -2.143)
    dist = squareform(pdist(Cbeta))

    contacts = dist < distance_threshold
    contacts = contacts.astype(np.int64)
    contacts[np.isnan(dist)] = -1
    return contacts


def evaluate_prediction(
    predictions: torch.Tensor,
    targets: torch.Tensor,
) -> Dict[str, float]:
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


"""Adapted from: https://github.com/rmrao/evo/blob/main/evo/visualize.py"""
def plot_contacts_and_predictions(
    predictions: Union[torch.Tensor, np.ndarray],
    contacts: Union[torch.Tensor, np.ndarray],
    ax: Optional[mpl.axes.Axes] = None,
    # artists: Optional[ContactAndPredictionArtists] = None,
    cmap: str = "Blues",
    ms: float = 1,
    title: Union[bool, str, Callable[[float], str]] = True,
    animated: bool = False,
) -> None:

    if isinstance(predictions, torch.Tensor):
        predictions = predictions.detach().cpu().numpy()
    if isinstance(contacts, torch.Tensor):
        contacts = contacts.detach().cpu().numpy()
    if ax is None:
        ax = plt.gca()

    seqlen = contacts.shape[0]
    relative_distance = np.add.outer(-np.arange(seqlen), np.arange(seqlen))
    bottom_mask = relative_distance < 0
    masked_image = np.ma.masked_where(bottom_mask, predictions)
    invalid_mask = np.abs(np.add.outer(np.arange(seqlen), -np.arange(seqlen))) < 6
    predictions = predictions.copy()
    predictions[invalid_mask] = float("-inf")

    topl_val = np.sort(predictions.reshape(-1))[-seqlen]
    pred_contacts = predictions >= topl_val
    true_positives = contacts & pred_contacts & ~bottom_mask
    false_positives = ~contacts & pred_contacts & ~bottom_mask
    other_contacts = contacts & ~pred_contacts & ~bottom_mask

    if isinstance(title, str):
        title_text: Optional[str] = title
    elif title:
        long_range_pl = compute_precisions(predictions, contacts, minsep=24)[
            "P@L"
        ].item()
        if callable(title):
            title_text = title(long_range_pl)
        else:
            title_text = f"Long Range P@L: {100 * long_range_pl:0.1f}"
    else:
        title_text = None

    img = ax.imshow(masked_image, cmap=cmap, animated=animated)
    oc = ax.plot(*np.where(other_contacts), "o", c="grey", ms=ms)[0]
    fn = ax.plot(*np.where(false_positives), "o", c="r", ms=ms)[0]
    tp = ax.plot(*np.where(true_positives), "o", c="b", ms=ms)[0]
    ti = ax.set_title(title_text) if title_text is not None else None
    # artists = ContactAndPredictionArtists(img, oc, fn, tp, ti)

    ax.axis("square")
    ax.set_xlim([0, seqlen])
    ax.set_ylim([0, seqlen])



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
def download_read_pdb(pdbcode, datadir, keepfile=True):
    """
    Downloads a PDB file from the Internet and saves it in a data directory.
    Then it reads and returns the structure inside.
    :param pdbcode: The standard PDB ID e.g. '3ICB'
    :param datadir: The directory where the downloaded file will be saved
    :param keepfile: if False, then the downloaded file will be deleted (default: keep the downloaded file)
    :return: a Bio.PDB Structure object or None if something went wrong
    """
    pdbfilenm = download_pdb(pdbcode, datadir)
    if pdbfilenm is None:
        return None
    struct = read_pdb(pdbcode, pdbfilenm)
    if not keepfile:
        os.remove(pdbfilenm)
    return struct


def download_pdb(pdbcode, datadir, downloadurl="http://files.rcsb.org/download/"):
    """
    Downloads a PDB file from the Internet and saves it in a data directory.
    :param pdbcode: The standard PDB ID e.g. '3ICB' or '3icb'
    :param datadir: The directory where the downloaded file will be saved
    :param downloadurl: The base PDB download URL, cf.
        `https://www.rcsb.org/pages/download/http#structures` for details
        Note that the unencrypted HTTP protocol is used by default
        to avoid spurious OpenSSL errors...
    :return: the full path to the downloaded PDB file or None if something went wrong
    """
    pdbfn = pdbcode + ".pdb"
    url = downloadurl + pdbfn
    outfnm = os.path.join(datadir, pdbfn)
    try:
        urllib.request.urlretrieve(url, outfnm)
        return outfnm
    except Exception as err:
        # all sorts of things could have gone wrong...
        print(str(err), file=sys.stderr)
        return None


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


def extract_seqrecords(pdbcode, struct):
    """
    Extracts the sequence records from a Bio.PDB structure.
    :param pdbcode: the PDB ID of the structure, needed to add a sequence ID to the result
    :param struct: a Bio.PDB.Structure object
    :return: a list of Bio.SeqRecord objects
    """
    ppb = Bio.PDB.PPBuilder()
    seqrecords = []
    for i, chain in enumerate(struct.get_chains()):
        print(i)
        print(chain)
        # extract and store sequences as list of SeqRecord objects
        pps = ppb.build_peptides(chain)  # polypeptides
        if len(pps) == 0:  # empty chain !! skip
            continue
        seq = pps[0].get_sequence()  # just take the first, hope there's no chain break
        seqid = pdbcode + chain.id
        seqrec = Bio.SeqRecord.SeqRecord(seq, id=seqid,
                                         description="Sequence #{}, {}".format(i + 1, seqid))
        seqrecords.append(seqrec)
    return seqrecords


def get_calphas(struct):
    """
    Extracts the C-alpha atoms from a PDB structure.
    :param struct: A Bio.PDB.Structure object.
    :return: A list of Bio.PDB.Atom objects representing the C-alpha atoms in `struct`.
    """
    calphas = [atom for atom in struct.get_atoms() if atom.get_fullname() == " CA "]
    return calphas


run_example = False
if run_example:
    # Example usage:
    print("Hello")
    genetic_code_dict = genetic_code()
    aa_list = list(set(genetic_code_dict.values())) # all possible amino-acids
    aa_list.remove("*")
    codon_list = list(genetic_code_dict)
    codon = 'GGA'
    amino_acid = genetic_code_dict.get(codon, 'Unknown')
    print(f'The amino acid corresponding to {codon} is {amino_acid}')

    for c in codon_list:
        print(c, genetic_code_dict.get(c), point_mutation_to_aa(c, genetic_code_dict))

    for aa in aa_list:
        print(aa,  aa_point_mutation_to_aa(aa, genetic_code_dict))


    # Methods for phylogentic reconstruction
    # constructor = DistanceTreeConstructor()
    # tree = constructor.nj(dm)
