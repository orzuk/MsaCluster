# Rub pipeline from a list of PDBs
# import copy
import sys
import os.path
# from os import *  # import os
import os
import glob
from protein_utils import *
from MSA_Clust import *
import subprocess
import platform

if platform.system() == "Linux":
    print("Run on cluster command line")
    run_mode = sys.argv[1]
else:
    print("Run on windows")
    run_mode = "plot" # "run_esm"  # sys.argv[1]


run_pipeline = False  # run entire pipeline
load_seq_and_struct = False  # jsut get from pdb the sequence and 3D structure for each protein
plot_results = False
if run_mode == "load":
    load_seq_and_struct = True  # jsut get from pdb the sequence and 3D structure for each protein
    print("Load pdb files, contact maps and fasta sequences for all families")
if run_mode == "run_esm":
    run_pipeline = True  # run entire pipeline
    print("Run ESM transfoer for all families")
if run_mode == "plot": # here do analysis of the results
    plot_results = True

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

cmap_dists_vec = [None]*n_fam # Results arrays
seqs_dists_vec = [None]*n_fam

for i in range(3, n_fam):  # loop on families
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
        pipeline_str = "./pipeline_get_params.sh " + fasta_file_name + " " + pdbids[i][0]  # Take one of the two !!! # ""./input/2qke.fasta 2qke
        pipeline_str = "sbatch -o './Pipeline/" + foldpair_ids[i] + "/run_pipeline_for_" + foldpair_ids[i] + ".out' ./pipeline_get_params.sh " + fasta_file_name + " " + foldpair_ids[i]  # Take one of the two !!! # ""./input/2qke.fasta 2qke
        print(pipeline_str)
#        subprocess.run(pipeline_str)     # run pipeline (should be a separate job1!!)
        os.system(pipeline_str)  # run pipeline (should be a separate job!)

    if plot_results:
        # First load files
        fasta_file_name = fasta_dir + "/" + foldpair_ids[i] + "/" + pdbids[i][0] + '.fasta'
        msa_file = fasta_dir + "/" + foldpair_ids[i] + "/output_get_msa/DeepMsa.a3m"
        MSA = AlignIO.read(open(msa_file), "fasta")
        msa_pred_files = glob.glob(fasta_dir + "/" + foldpair_ids[i] + "/output_cmap_esm/*.npy")
        n_cmaps = len(msa_pred_files)
        n_cmaps = min(3, n_cmaps)  # temp for debug !!
#        msa_transformer_pred = [None]*n_cmaps
#        for file in msa_pred_files:
#            print(file)
#            msa_transformer_pred[i] = np.genfromtxt(file)
        msa_transformer_pred = { file.split("\\")[-1][14:-4] : np.genfromtxt(file) for file in msa_pred_files  }
        print("Predicted cmap sizes:")
        print([c.shape[0] for c in msa_transformer_pred.values()])

#        true_cmap = [None]*2
#        for fold in range(2):
#            true_cmap[fold] = np.genfromtxt(fasta_dir +
#                "/" + foldpair_ids[i] + "/" + pdbids[i][fold] + pdbchains[i][fold] + "_pdb_contacts.npy")  # load for pickle, genfromtxt for tab-delimited

        true_cmap = {pdbids[i][fold] :  np.genfromtxt(fasta_dir +  # problem with first !!
                    "/" + foldpair_ids[i] + "/" + pdbids[i][fold] + pdbchains[i][fold] + "_pdb_contacts.npy").astype(int) for fold in range(1)}
        print("True cmap sizes:")
        print([c.shape[0] for c in true_cmap.values()])
        print("plot")
        plot_array_contacts_and_predictions(msa_transformer_pred, true_cmap)
        cmap_dists_vec[i] = compute_cmap_distances(msa_transformer_pred)
        seqs_dists_vec[i] = np.mean(compute_seq_distances(MSA)) # can take the entire sequence !


        break
        # next plotP


#import numpy as np
#cmap_pred = np.genfromtxt("Pipeline/1dzlA_5keqF/output_cmap_esm/msa_t__ShallowMsa_000.npy")
#cmap_true = np.genfromtxt("Pipeline/1dzlA_5keqF/1dzlA_pdb_contacts.npy")
#plot_array_contacts_and_predictions(cmap_pred, cmap_true)
