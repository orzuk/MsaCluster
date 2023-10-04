# Rub pipeline from a list of PDBs
from protein_utils import *
import subprocess


pdb_datadir = "Pipeline/pdb_files"  # where to store all PDB files
fasta_dir = "Pipeline"
pdbids_file = "foldswitch_ids.txt"   # file with all pdb ids

with open(pdbids_file) as file:  # read all pdb ids
    pdbids = [line.rstrip() for line in file]

n_fam = len(pdbids)  # number of families

for i in range(n_fam): # loop on families
    pdb_struct =  download_read_pdb(pdbids[i], pdb_datadir, keepfile=True)  # extract pdb file
    pdb_seq = extract_seqrecords(pdbids[i], pdb_struct)  # extract sequences
    fasta_file_name = fasta_dir + "/" + pdbids

    with open(fasta_file_name, "w") as text_file: # save to fasta file
        text_file.write(pdb_seq)

    pipeline_str = "./pipeline_get_params.sh " + fasta_file_name + " " + pdbids[i]   # ""./input/2qke.fasta 2qke
    subprocess.run(pipeline_str)     # run pipeline (should be a separate job1!!)

