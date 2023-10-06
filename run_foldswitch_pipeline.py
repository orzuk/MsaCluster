# Rub pipeline from a list of PDBs
# import copy
import os.path
# from os import *  # import os
import os
from protein_utils import *
import subprocess

run_pipeline = False  # run entire pipeline
load_seq_and_struct = True  # jsut get from pdb the sequence and 3D structure for each protein

# pdb_datadir = "Pipeline/pdb_files"  # where to store all PDB files
fasta_dir = "Pipeline"
pdbids_file = "data/foldswitch_PDB_IDs_full.txt"   # file with all pdb ids

with open(pdbids_file, "r") as file:  # read all pdb ids
    pdbids = [line.rstrip() for line in file]  # two per row
foldpair_ids = [s.replace("\t", "_") for s in pdbids]

pdbids = [s.split("\t") for s in pdbids]
pdbchains = [[s[0][-1], s[1][-1]] for s in pdbids]
pdbids = [[s[0][:-1], s[1][:-1]] for s in pdbids]

n_fam = len(pdbids)  # number of families

for i in range(0, n_fam):  # loop on families
    if load_seq_and_struct:
        for fold in range(2):
            cur_family_dir = fasta_dir + "/" + foldpair_ids[i]
            if not os.path.exists(cur_family_dir):
                print("Mkdir: " + cur_family_dir)
                os.mkdir(cur_family_dir)
            print("Get seq + struct for " + pdbids[i][fold] + ", " + str(i) + " out of " + str(n_fam))
            pdb_struct = download_read_pdb(pdbids[i][fold], cur_family_dir, keepfile=True)  # extract pdb file
            pdb_seq = extract_seqrecords(pdbids[i][fold], pdb_struct)  # extract sequences
            fasta_file_name = fasta_dir + "/" + foldpair_ids[i] + "/" + pdbids[i][fold] + '.fasta'
#            print( pdbids[i][fold] + ":" + pdbchains[i][fold])

#            print("Index: " + str([s.id[-1] for s in pdb_seq].index(pdbchains[i][fold])))

            with open(fasta_file_name, "w") as text_file:  # save to fasta file. Take the correct chain
                text_file.writelines([ "> " + pdbids[i][fold].upper() + ":" + pdbchains[i][fold].upper() + '\n', str(pdb_seq[[s.id[-1] for s in pdb_seq].index(pdbchains[i][fold])].seq) ])

            # Finally, make a contact map from each pdb file:
            # Read structure in slightly different format
#            print("Get struct for ontacts, Cbeta distance < 8Angstrom")
            pdbX_struct = get_structure(PDBxFile.read(rcsb.fetch(pdbids[i][fold], "cif")))[0]
#            print("Get contacts, Cbeta distance < 8Angstrom")
            pdb_contacts = contacts_from_pdb(pdbX_struct, chain=pdbchains[i][fold])  # =True)  # extract pdb file
            np.savetxt(cur_family_dir + "/" + pdbids[i][fold] + pdbchains[i][fold] + "_pdb_contacts.npy", pdb_contacts)  # save contacts

    if run_pipeline:
        fasta_file_name = fasta_dir + "/" + foldpair_ids[i] + "/" + pdbids[i][0] + '.fasta'
        pipeline_str = "./pipeline_get_params.sh " + fasta_file_name + " " + pdbids[i] [0]  # Take one of the two !!! # ""./input/2qke.fasta 2qke
        subprocess.run(pipeline_str)     # run pipeline (should be a separate job1!!)

