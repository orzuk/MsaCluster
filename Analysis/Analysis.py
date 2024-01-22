from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from argparse import  ArgumentParser
from Bio.PDB import PDBParser, PDBIO, Select
from Bio import SeqIO
import mdtraj as md
from contact_map import ContactFrequency, ContactDifference
import numpy as np
import mdtraj as md

pdf_files_path = f'/Users/steveabecassis/Desktop/pdb_file'
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
    structure_fold1 = parser.get_structure('PDB_structure' , f'{pdb_file_path}/{fold_1}.pdb')
    structure_fold2 = parser.get_structure('PDB_structure' , f'{pdb_file_path}/{fold_2}.pdb')
    io = PDBIO()
    # Set the structure for saving and use ChainSelect to filter the chain
    io.set_structure(structure_fold1)
    io.save(f'./{chain_pdb_file_path}/{fold_1}.pdb', ChainSelect(chain_fold_1))
    io.set_structure(structure_fold2)
    io.save(f'./{chain_pdb_file_path}/{fold_2}.pdb', ChainSelect(chain_fold_2))


def chain_seq(pdb_file):
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
    return sequence


def save_org_cmaps(chain_pdb_file_path,fold):
    traj = md.load(filename_or_filenames=f'{chain_pdb_file_path}/{fold}.pdb')
    frame_contacts = ContactFrequency(traj[0])
    np.save(f'./Pipeline/org_cmaps/{fold}.npy',frame_contacts)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument("--fold_1", help="should be pdb codes _ and the chain. ex:1ebo_J")
    parser.add_argument("--fold_2", help="should be pdb codes _ and the chain. ex:1ebo_J")
    # parser.add_argument("--output", help="Directory to write the results to")
    # parser.add_argument("--pdb_file_path", default='./pdb_files', type=int,help='files with all the pdb')
    # parser.add_argument("--chain_pdb_file_path", default='./chain_pdb_file_path', type=int, help='files with all the pdb')
    args = parser.parse_args()

    fold_1              = args.fold_1
    fold_2              = args.fold_2
    pdb_file_path       = f'./Pipeline/{fold_1}_{fold_2}/pdb_file_path'
    chain_pdb_file_path = f'./Pipeline/{fold_1}_{fold_2}/chain_pdb_file_path'

    create_chain_pdb_files(fold_1,fold_2,pdb_file_path,chain_pdb_file_path)

    sequence = chain_seq(f'{chain_pdb_file_path}/{fold_1}.pdb')
    with open(f"./Pipeline/fasta_files/{fold_1}.fasta", "w") as output_handle:
        SeqIO.write(sequence, output_handle, "fasta")

    save_org_cmaps(chain_pdb_file_path, fold_1)
    save_org_cmaps(chain_pdb_file_path, fold_2)







    # sequence = chain_seq(pdb_file,'A')
    #
    # traj = md.load(filename_or_filenames='/Users/steveabecassis/Desktop/pdb_file/5fhc_A.pdb')
    # topology = traj.topology
    # frame_contacts = ContactFrequency(traj[0])
    # fig, ax = frame_contacts.residue_contacts.plot()
    # plt.xlabel("Residue")
    # _ = plt.ylabel("Residue")
    # a = frame_contacts.residue_contacts









# # Load the original PDB file
# parser = PDBParser()
# structure = parser.get_structure('PDB_structure', pdb_file)
#
# # Specify the chain you want to isolate
# chain_to_keep = 'J'  # Replace with your desired chain identifier
#
# # Create an object to save the PDB structure
# io = PDBIO()
#
# # Set the structure for saving and use ChainSelect to filter the chain
# io.set_structure(structure)
# io.save('/Users/steveabecassis/Desktop/pdb_file/5fhc_J.pdb', ChainSelect(chain_to_keep))



'''
1) Get the chain from the csv file              DONE
2) Create new pdb file of the csv file          DONE
3) Get sequence from the new pdb_file           DONE
4) Get contactmap from the new pdb file and save the cmaps in npy format DONE
4.a) Get stats about the cmap
5) Compare the concact map with the cmap from MSA transformers
6) Compare the pdb file with AF
7) Compare the pdb file with ESM
'''


