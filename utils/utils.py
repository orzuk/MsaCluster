from Bio.SeqUtils import seq1
from Bio.PDB import PDBParser, PDBIO, Select
from Bio import SeqIO
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from contact_map import ContactFrequency, ContactDifference

user = 'zuk_laptop'  # "steve_moriah", "steve_laptop", "zuk_moriah"

# Set here all directories:
if user == "steve_moriah":
    pdf_files_path = f'/Users/steveabecassis/Desktop/pdb_file'
if user == "zuk_laptop":
    pdf_files_path = f'/mnt/c/Data/pdb_file'

pdb_file = f'{pdf_files_path}/1jfk.pdb'


class ChainSelect(Select):
    def __init__(self, chain_letter):
        self.chain_letter = chain_letter

    def accept_chain(self, chain):
        return chain.id == self.chain_letter


def create_chain_pdb_files(fold_1,fold_2,pdb_file_path,chain_pdb_file_path):
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
def get_fasta_chain_seq(pdb_file, fold_name, output_dir):
    # Create a PDB parser
    parser = PDBParser()
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

# pdb_file = '/Users/steveabecassis/Desktop/1eboE.pdb'
# with open(f"/Users/steveabecassis/Desktop/test.fasta", "w") as output_handle:
#     output_handle.write('>'+'\n>'+sequence)


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