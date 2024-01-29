from typing import List
from Bio import PDB
import mdtraj as md

def Alpha_Carbon_Indices(pdb: str) -> List:
    """
    take in a pdb file and identify the index of every alpha carbon
    """
    structure = PDB.PDBParser(QUIET=True).get_structure('protein', pdb)

    alpha_carbons = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    resid = residue.resname
                    alpha_carbons.append([resid, residue['CA'].get_serial_number() - 1])
    return alpha_carbons


def Match_Alpha_Carbons(pdb_1: str, pdb_2: str) -> List[int]:
    """
    Take in two pdb structure files and search through them for matching alpha carbons
    This should identify positions correctly even if sequences are not identical
    """
    alpha_c_1 = Alpha_Carbon_Indices(pdb_1)
    alpha_c_2 = Alpha_Carbon_Indices(pdb_2)

    matching_alpha_carbons1 = []
    matching_alpha_carbons2 = []

    for i, (resname_1, ca_index1) in enumerate(alpha_c_1):
        for j, (resname_2, ca_index2) in enumerate(alpha_c_2):
            if resname_2 == resname_1 and ca_index1 not in [_[1] for _ in matching_alpha_carbons1] and ca_index2 not in [_[1] for _ in matching_alpha_carbons2]:
                #prevent erroneous match at NTD
                if i > 0 and j > 0:
                    if alpha_c_1[i-1][0] != alpha_c_2[j-1][0]: #check previous matches
                        continue
                # prevent erroneous backtracking
                if len(matching_alpha_carbons1) > 2 and len(matching_alpha_carbons2) > 2:
                    if ca_index2 < matching_alpha_carbons2[-1][-1]:
                        continue
                #prevent erroneous match at CTD
                if i < len(alpha_c_1) - 1 and j < len(alpha_c_2) - 1:
                    if alpha_c_1[i+1][0] != alpha_c_2[j+1][0]: #check next matches
                        continue

                matching_alpha_carbons1.append([resname_1, ca_index1])
                matching_alpha_carbons2.append([resname_2, ca_index2])
                break
    #skip first residue to avoid erroneous glycine match
    return matching_alpha_carbons1[1:], matching_alpha_carbons2[1:]



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


pdb1 = '/Users/steveabecassis/Desktop/1eboE_5fhcJ/chain_pdb_files/1eboE.pdb'
pdb2 = '/Users/steveabecassis/Desktop/1eboE_5fhcJ/chain_pdb_files/5fhcJ.pdb'

structure_index,reference_index = Match_Alpha_Carbons(pdb1, pdb2)
Calculate_RMSD(pdb1, pdb2, [_[1] for _ in structure_index],[_[1] for _ in reference_index])[0]



import subprocess

# Define the terminal command you want to run
command = "./TMalign ./chain_pdb_files/1eboE.pdb ./chain_pdb_files/5fhcJ.pdb"
command = "./TMalign"
command = "pwd"
command = "gfortran  -O3 -ffast-math -lm -o TMalign /Users/steveabecassis/Desktop/TMalign.f"
command = "./TMalign   ./chain_pdb_files/1eboE.pdb ./chain_pdb_files/5fhcJ.pdb"
command = "./TMalign   ./chain_pdb_files/1eboE.pdb ./chain_pdb_files/5fhcJ.pdb"

# Run the command and capture its output
output = subprocess.check_output(command, shell=True)
# Print the output
print(output.decode())
import re
match = re.search(r"TM-score=\s+(\d+\.\d+)", str(output))

if match:
    result = match.group(1)
    print(result)
float(result)


