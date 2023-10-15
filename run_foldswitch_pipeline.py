# Rub pipeline from a list of PDBs
# import copy
import sys
import os.path
import os
import glob
from protein_utils import *
from MSA_Clust import *
import subprocess
import platform

from Bio import Align

if platform.system() == "Linux":
    print("Run on cluster command line")
    run_mode = sys.argv[1]
else:
    print("Run on windows")
    run_mode = "run_pipeline"  # "load"  # "run_esm" # "plot" # "run_esm"  # sys.argv[1]

run_pipeline = False  # run entire pipeline
run_esm = False # run just esm contacts
load_seq_and_struct = False  # jsut get from pdb the sequence and 3D structure for each protein
plot_results = False

# can't use match (only from python 3.10)
# match run_mode:
if run_mode == "load":
    load_seq_and_struct = True  # just get from pdb the sequence and 3D structure for each protein
    print("Load pdb files, contact maps and fasta sequences for all families")
if run_mode == "run_esm":
    run_esm = True  # run entire pipeline
    print("Run ESM transformer for all families")
if run_mode == "run_pipeline":
    run_pipeline = True  # run entire pipeline
    print("Run Entire pipeline for all families")
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

for i in range(2, n_fam):  # loop on families
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

    if run_esm: # Just compute conctacts !!
        esm_str = "sbatch -o './Pipeline/" + foldpair_ids[i] + "/run_esm_for_" + foldpair_ids[i] + ".out' ./Pipeline/CmapESM_params.sh  " + foldpair_ids[i]  # Take one of the two !!! # ""./input/2qke.fasta 2qke
        print(esm_str)
#        os.system(esm_str)  # run pipeline (should be a separate job!)

    if run_pipeline:  # entire pipeline: get alignments, contacts, predicted structure
        fasta_file_name = fasta_dir + "/" + foldpair_ids[i] + "/" + pdbids[i][0] + '.fasta'
        pipeline_str = "sbatch -o './Pipeline/" + foldpair_ids[i] + "/run_pipeline_for_" + foldpair_ids[i] + ".out' ./pipeline_get_params.sh " + fasta_file_name + " " + foldpair_ids[i]  # Take one of the two !!! # ""./input/2qke.fasta 2qke
        print(pipeline_str)
        os.system(pipeline_str)  # run pipeline (should be a separate job!)

    if plot_results:
        fasta_file_names = {pdbids[i][fold] : fasta_dir + "/" + foldpair_ids[i] + "/" + pdbids[i][fold] + '.fasta' for fold in range(2)}
        # First load files
#        fasta_file_name = fasta_dir + "/" + foldpair_ids[i] + "/" + pdbids[i][0] + '.fasta'
#        fasta_file_name1 = fasta_dir + "/" + foldpair_ids[i] + "/" + pdbids[i][1] + '.fasta'  # Second file  !!!

#        print(fasta_file_name)
#        print(fasta_file_name1)
        seq = {}
        for fold in range(2):
            with open(fasta_file_names[pdbids[i][fold]], "r") as text_file:
                seq[pdbids[i][fold]]  = text_file.read().split("\n")[1]
#        with open(fasta_file_name1, "r") as text_file:
#            seq1 = text_file.read().split("\n")[1]

        print("Aligning:")
        print(seq[pdbids[i][0]])
        print(seq[pdbids[i][1]])
#        pairwise_alignment = pairwise2.align.globalxx(seq[pdbids[i][0]], seq[pdbids[i][1]])
#        print("Alignment:")
#        print(pairwise_alignment[0].seqA)
#        print(pairwise_alignment[0].seqB)
        # Alternative alg for alignment
        aligner = Align.PairwiseAligner()
        pairwise_alignment = aligner.align(seq[pdbids[i][0]], seq[pdbids[i][1]])

        msa_file = fasta_dir + "/" + foldpair_ids[i] + "/output_get_msa/DeepMsa.a3m"
        MSA = read_msa(msa_file)  # AlignIO.read(open(msa_file), "fasta")
        msa_pred_files = glob(fasta_dir + "/" + foldpair_ids[i] + "/output_cmap_esm/*.npy")
        n_cmaps = len(msa_pred_files)
        n_cmaps = min(3, n_cmaps)  # temp for debug !!
#        msa_transformer_pred = [None]*n_cmaps
#        for file in msa_pred_files:
#            print(file)
#            msa_transformer_pred[i] = np.genfromtxt(file)
        try:  # read in text format or python format
            msa_transformer_pred = {file.split("\\")[-1][14:-4] : np.genfromtxt(file) for file in msa_pred_files}
        except:
            msa_transformer_pred = {file.split("\\")[-1][14:-4]: np.load(file) for file in msa_pred_files}

        print("Predicted cmap sizes for:" + foldpair_ids[i])
        print([c.shape[0] for c in msa_transformer_pred.values()])

#        true_cmap = [None]*2
#        for fold in range(2):
#            true_cmap[fold] = np.genfromtxt(fasta_dir +
#                "/" + foldpair_ids[i] + "/" + pdbids[i][fold] + pdbchains[i][fold] + "_pdb_contacts.npy")  # load for pickle, genfromtxt for tab-delimited
        try:
            true_cmap = {pdbids[i][fold] : np.genfromtxt(fasta_dir +  # problem with first !!
                        "/" + foldpair_ids[i] + "/" + pdbids[i][fold] + pdbchains[i][fold] + "_pdb_contacts.npy").astype(int) for fold in range(1)}
            print("Loaded text true!")
        except:
            print("load binary true!")
            true_cmap = {pdbids[i][fold] : np.load(fasta_dir +  # problem with first !!
                        "/" + foldpair_ids[i] + "/" + pdbids[i][fold] + pdbchains[i][fold] + "_pdb_contacts.npy").astype(int) for fold in range(1)}

        print("True cmap sizes:")
        print([c.shape[0] for c in true_cmap.values()])
        print("plot")

        # Match indices of the two folds:
        match_true_cmap, match_predicted_cmaps = get_matching_indices_two_maps(pairwise_alignment, true_cmap, msa_transformer_pred)

#        plot_array_contacts_and_predictions(msa_transformer_pred, true_cmap)
        plot_array_contacts_and_predictions(match_predicted_cmaps, match_true_cmap)
        cmap_dists_vec[i] = compute_cmap_distances(msa_transformer_pred)
        seqs_dists_vec[i] = np.mean(compute_seq_distances(MSA)) # can take the entire sequence !

        print("Cmap dist: " + str(cmap_dists_vec[i]) + ", seq dist:" + str(seqs_dists_vec[i]))

        break
        # next plotP


#import numpy as np
#cmap_pred = np.genfromtxt("Pipeline/1dzlA_5keqF/output_cmap_esm/msa_t__ShallowMsa_000.npy")
#cmap_true = np.genfromtxt("Pipeline/1dzlA_5keqF/1dzlA_pdb_contacts.npy")
#plot_array_contacts_and_predictions(cmap_pred, cmap_true)
