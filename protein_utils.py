# Some utilities for proteins and their mutations
import copy

import pandas as pd
# import esm
import string

# import pcmap
import torch
import torch.nn.functional as F  # for padding

from scipy.spatial.distance import squareform, pdist, cdist
import numpy as np
from typing import List, Tuple, Optional, Dict, NamedTuple, Union, Callable
import matplotlib as mpl
import matplotlib.pyplot as plt
import Bio
import Bio.PDB
import Bio.SeqRecord
from Bio import SeqIO
from Bio.PDB import PDBParser


import pickle
import os
import sys
import urllib
# import math

import biotite.structure as bs
from biotite.structure.io.pdbx import PDBxFile, get_structure
from biotite.database import rcsb

from tmtools import tm_align
from tmtools.io import get_structure, get_residue_data

import iminuit
import tmscoring  # for comparing structures


from Bio import AlignIO

# from TreeConstruction import DistanceTreeConstructor


# This is an efficient way to delete lowercase characters and insertion characters from a string
deletekeys = dict.fromkeys(string.ascii_lowercase)
deletekeys["."] = None
deletekeys["*"] = None
translation = str.maketrans(deletekeys)


def read_sequence(filename: str) -> Tuple[str, str]:
    """ Reads the first (reference) sequences from a fasta or MSA file."""
    record = next(SeqIO.parse(filename, "fasta"))
    return record.description, str(record.seq)


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
                 'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
aa_short_long = {y: x for x, y in aa_long_short.items()}


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


# Run design on every position, collect
def design_every_position(S, pdbID):
    n = len(S)  # get number of amino-acids
    S_des = [''] * n  # The output designed sequence
    for i in range(n):  # loop on positions
        cur_S = run_design(S, pdbID, i)  # run design function using ProteinMPNN, keeping all positions fixed except i
        S_des[i] = cur_S[i]  # change the i-th letter of the i-th protein
    S_des = "".join(S_des)  # convert to string
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
    for i_true in range(2):  # loop on the two true structures
        for j_pred in range(3):  # S_predicted in [AF, AF1, AF2]:  # loop on the three predicted structures
            #            df_tm[TM[i_true]][SEQ[j_pred]] = TMScore(S_true_list[i_true], S_pred_list[j_pred])  # Compute TMScore similarity
            alignment = tmscoring.TMscoring(S_true_list[i_true], S_pred_list[j_pred])  # 'structure1.pdb', 'structure2.pdb')  # from installed tmscoring
            df_tm[TM[i_true]][SEQ[j_pred]] = alignment.tmscore(**alignment.get_current_values())

    print(df_tm)
    return df_tm, S, S1, S2, AF, AF1, AF2  # Return all sequences, structures and their similarity


def extract_protein_sequence(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure("protein", pdb_file)

    sequences = []
    for model in structure:
        for chain in model:
            seq = []
            for residue in chain:
                if residue.get_id()[0] == ' ':
                    seq.append(residue.get_resname())
            sequences.append(''.join(seq))

    return sequences


# Compute tmscores of two structures, interface to tmscore module
def compute_tmscore(pdb_file1, pdb_file2, chain1=[], chain2=[]):
    print("Start compute tmscore:")
    print(pdb_file1)
    print(pdb_file2)
    print(chain1)

    s1 = get_structure(pdb_file1)
    s2 = get_structure(pdb_file2)
    chain1 = next(s1.get_chains())
    coords1, seq1 = get_residue_data(chain1)
    chain2 = next(s2.get_chains())
    coords2, seq2 = get_residue_data(chain2)

    print("read sequences:")
    chain1 = []
    seq1 = ''
    for c in s1.get_chains():
        chain1.append(c)
        coords1, cur_seq1 = get_residue_data(c)
        seq1 += cur_seq1

    seq1_alt = extract_protein_sequence(pdb_file1)  # Alternative reading:

    print("Now align")
    print(coords1.shape)
    print(coords2.shape)
    print(seq1)
    print(seq1_alt)
    print(seq2)
    print(len(seq1))
    print(len(seq2))

    with open("temp_tm_align.pkl", "wb") as f:
        pickle.dump([coords1, coords2, seq1, seq2], f)

    res = tm_align(coords1, coords2, seq1, seq2)

    print("TM RES: ")
    print(res)
    print(res.tm_norm_chain1)
    return res.tm_norm_chain1  # return results (both scores, and matrices)


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



# Extract contact map from a pdb-file
# Also extract the distances themselves (more informative than the thresholded contacts)
# Output:
# dist - pairwise distance matrix between cBeta atoms
# contacts - pairwise binary contacts matrix for distances < threshold
# pdb_seq - sequence extracted from pdb file, after removing residues with missing atoms
# good_res_ids - indices of full good residues
def contacts_from_pdb(
        structure: bs.AtomArray,
        distance_threshold: float = 8.0,
        chain: Optional[str] = None,
) -> np.ndarray:
    mask = ~structure.hetero
    if chain is not None:
        mask &= structure.chain_id == chain
    # Problem: what if they're not compatible?
    N = structure.coord[mask & (structure.atom_name == "N")]
    CA = structure.coord[mask & (structure.atom_name == "CA")]
    C = structure.coord[mask & (structure.atom_name == "C")]
    pdb_seq = "".join([aa_long_short[a] for a in structure.res_name[np.where(mask & (structure.atom_name == "N"))[0]]])

    good_res_ids = structure.res_id[mask & (structure.atom_name == "N")]  # all residues indices
    if len(N) != len(CA) or len(N) != len(C) or len(C) != len(CA):  # missing atoms
        print("Missing atoms in PDB! remove residues!")
        print([len(N), len(CA), len(C)])
        good_res_ids = np.intersect1d(np.intersect1d(structure.res_id[mask & (structure.atom_name == "N")],
                                                     structure.res_id[mask & (structure.atom_name == "CA")]),
                                      structure.res_id[mask & (structure.atom_name == "C")])
        N = structure.coord[mask & (structure.atom_name == "N") & np.in1d(structure.res_id, good_res_ids)]
        CA = structure.coord[mask & (structure.atom_name == "CA") & np.in1d(structure.res_id, good_res_ids)]
        C = structure.coord[mask & (structure.atom_name == "C") & np.in1d(structure.res_id, good_res_ids)]
        pdb_seq = "".join([aa_long_short[a] for a in structure.res_name[
            np.where(mask & (structure.atom_name == "N") & np.in1d(structure.res_id, good_res_ids))[0]]])

        # return []
    Cbeta = extend(C, N, CA, 1.522, 1.927, -2.143)
    dist = squareform(pdist(Cbeta))

    contacts = dist < distance_threshold
    contacts = contacts.astype(np.int64)
    contacts[np.isnan(dist)] = -1
    return dist, contacts, pdb_seq, good_res_ids   # [aa_long_short[aa] for aa in structure.res_name[good_res_ids]]


# Evaluate precision of predicted contacts with respect to true contacts
# Input:
# predictions - matrix of predicted residue affinities
# targets - matrix of binary contacts
# Output :
# AUC - Area under the ROC curve for preditcing contacts
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
        binned_cumulative_dist
    )

    pl5 = binned_precisions[:, 1]
    pl2 = binned_precisions[:, 4]
    pl = binned_precisions[:, 9]
    auc = binned_precisions.mean(-1)

    return {"AUC": auc, "P@L": pl, "P@L2": pl2, "P@L5": pl5}


# General utilitiy for dictionary of unique values
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
        #        print(i)
        #        print(chain)
        # extract and store sequences as list of SeqRecord objects
        pps = ppb.build_peptides(chain)  # polypeptides
        if len(pps) == 0:  # empty chain !! skip
            continue
        seq = pps[0].get_sequence()  # just take the first, hope there's no chain break
        for i in range(1, len(pps)):  # New: add all parts !!!
            seq += pps[i].get_sequence()
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


# Match between cmaps, get only aligned indices
# Match between cmaps, get only aligned indices
def get_matching_indices_two_maps(pairwise_alignment, true_cmap, pred_cmap):
    #    n_true = len(true_cmap)  # always 2 !!
    #    n_pred = len(pred_cmap)  # variable number !!

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
