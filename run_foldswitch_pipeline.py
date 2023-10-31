# Rub pipeline from a list of PDBs
# import copy
# import sys
# import os.path
# import os
# import glob
# from protein_utils import *
# import subprocess
import copy

import pandas as pd
from phytree_utils import *

from MSA_Clust import *
import platform

from Bio import Align

# problematic_families = ["1nqjB_1nqdA", "1qlnA_1h38D", "3l5nB_2a73B", "2bzyB_2lqwA", "4cmqB_4zt0C", "5tpnA_1g2cF"]  # second, third is large, no/one cmaps were generated !
# problematic_families = ["2bzyB_2lqwA", "4cmqB_4zt0C", "5tpnA_1g2cF"]  # 2bzyB_2lqwA bad (missing ) true cmap!!
foldpair_ids_to_run = 'ALL'  # '3j7vG_3j7wB' # '2vfxL_3gmhL' # '1xjuB_1xjtA'  # "ALL"
if platform.system() == "Linux":
    print("Run on cluster command line")
    run_mode = sys.argv[1]
    if len(sys.argv) > 2:
        foldpair_ids_to_run = sys.argv[2]  # enable running for a specific family (default is running on all of them)
else:
    print("Run on windows")
    run_mode = "load"   # "plot"  # "load"  # "run_esm" # "plot" # "run_esm"  # sys.argv[1]
    foldpair_ids_to_run = "1eboE_5fhcJ"  #  "4gqcB_4gqcC"  # problematic_families  # '1nqjB_1nqdA'  # Problem with pdb to contact  '2n54B_2hdmA'  #  '4yhdG_7ahlE' #  '5l35G_5l35D' # '1eboE_5fhcJ'

# print("Running on: " + foldpair_ids_to_run)

run_pipeline, get_msa, cluster_msa, tree_reconstruct, run_esm, load_seq_and_struct, plot_results = [False]*7  # run entire pipeline or parts of it/plot ...

# can't use match (only from python 3.10)
# match run_mode:
if run_mode == "load":
    load_seq_and_struct = True  # just get from pdb the sequence and 3D structure for each protein
    print("Load pdb files, contact maps and fasta sequences for all families")
if run_mode == "get_msa":  # here do analysis of the results
    get_msa = True
if run_mode == "cluster_msa":  # here do analysis of the results
    cluster_msa = True
if run_mode == "run_esm":
    run_esm = True  # run entire pipeline
    print("Run ESM transformer for all families")
if run_mode == "run_pipeline":
    run_pipeline = True  # run entire pipeline
    print("Run Entire pipeline for all families")
if run_mode == "plot":  # here do analysis of the results
    plot_results = True
if run_mode == "tree":  # here do analysis of the results
    tree_reconstruct = True

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
cmap_dists_vec = [None]*n_fam  # Results arrays
seqs_dists_vec = [None]*n_fam
num_seqs_msa_vec = [None]*n_fam
if foldpair_ids_to_run == "ALL":
    foldpair_ids_to_run = foldpair_ids
else:  # make a list
    if type(foldpair_ids_to_run) == str:
        foldpair_ids_to_run = [foldpair_ids_to_run]

for foldpair_id in foldpair_ids_to_run:   # for i in range(17, n_fam):  # loop on families
    i = foldpair_ids.index(foldpair_id)
    cur_family_dir = fasta_dir + "/" + foldpair_id
    print("Run: " + run_mode + " : " + str(i) + " : " + foldpair_id)

    fasta_file_name = fasta_dir + "/" + foldpair_id + "/" + pdbids[i][0] + pdbchains[i][0] + '.fasta'
    if load_seq_and_struct or run_pipeline:  # also for entire pipeline
        for fold in range(2):
            if not os.path.exists(cur_family_dir):
                print("Mkdir: " + cur_family_dir)
                os.mkdir(cur_family_dir)
            print("Get seq + struct for " + pdbids[i][fold] + ", " + str(i) + " out of " + str(n_fam-1) )
            # Old option for extracting sequence
#            pdb_struct = download_read_pdb(pdbids[i][fold], cur_family_dir, keepfile=True)  # extract pdb file
#            pdb_seq = extract_seqrecords(pdbids[i][fold], pdb_struct)  # extract sequences
#            with open(fasta_file_name, "w") as text_file:  # save to fasta file. Take the correct chain
#                text_file.writelines([ "> " + pdbids[i][fold].upper() + ":" + pdbchains[i][fold].upper() + '\n',
#                                       str(pdb_seq[[s.id[-1] for s in pdb_seq].index(pdbchains[i][fold])].seq) ])

            fasta_file_name = fasta_dir + "/" + foldpair_id + "/" + pdbids[i][fold] + pdbchains[i][fold] + '.fasta'  # added chain to file ID
#            print( pdbids[i][fold] + ":" + pdbchains[i][fold])

#            print("Index: " + str([s.id[-1] for s in pdb_seq].index(pdbchains[i][fold])))

            # Finally, make a contact map from each pdb file:
            # Read structure in slightly different format
#            print("Get struct for ontacts, Cbeta distance < 8Angstrom")
            pdbX_struct = get_structure(PDBxFile.read(rcsb.fetch(pdbids[i][fold], "cif")))[0]
#            print("Get contacts, Cbeta distance < 8Angstrom")
#            pdb_contacts = contacts_from_pdb(pdbX_struct, chain=pdbchains[i][fold])  # =True)  # extract pdb file
            # New option: extract sequence and structure togehter. Remove from sequence the residues without contacts
            pdb_dists, pdb_contacts, pdb_seq, pdb_good_res_inds = contacts_from_pdb(   # extract distances from pdb file
                get_structure(PDBxFile.read(rcsb.fetch(pdbids[i][fold], "cif")))[0], chain=pdbchains[i][fold])
            with open(fasta_file_name, "w") as text_file:  # save to fasta file. Take the correct chain
                text_file.writelines([ "> " + pdbids[i][fold].upper() + ":" + pdbchains[i][fold].upper() + '\n',
                                       pdb_seq ])
            print(cur_family_dir + "/" + pdbids[i][fold] + pdbchains[i][fold] + "_pdb_contacts.npy")
            np.save(cur_family_dir + "/" + pdbids[i][fold] + pdbchains[i][fold] + "_pdb_contacts.npy", pdb_contacts)  # save contacts (binary format)

    if get_msa:  # Cluster the Multiple Sequence Alignment for the family
        get_msa_str = "sbatch -o './Pipeline/" + foldpair_id + "/get_msa_for_" + foldpair_id + ".out' ./Pipeline/get_msa_params.sh " + \
                      fasta_file_name + " " + foldpair_id  # Take one of the two !!! # ""./input/2qke.fasta 2qke
        print(get_msa_str)
        os.system(get_msa_str)

    if cluster_msa:  # Cluster the Multiple Sequence Alignment for the family
        cluster_msa_str = "sbatch -o './Pipeline/" + foldpair_id + "/cluster_msa_for_" + foldpair_id + ".out' ./Pipeline/ClusterMSA_params.sh " + foldpair_id  # Take one of the two !!! # ""./input/2qke.fasta 2qke
        print(cluster_msa_str)
        os.system(cluster_msa_str)

    if tree_reconstruct:  # New: Add phylogenetic tree reconstruction of MSA sequences
        phytree_msa_str = "sbatch -o './Pipeline/" + foldpair_id + "/tree_reconstruct_for_" + foldpair_id + ".out' ./Pipeline/tree_reconstruct_params.sh " + foldpair_id  # Take one of the two !!! # ""./input/2qke.fasta 2qke
        print(phytree_msa_str)
        os.system(phytree_msa_str)

    if run_esm:  # Just compute contacts !!
        esm_str = "sbatch -o './Pipeline/" + foldpair_id + "/run_esm_for_" + foldpair_id + ".out' ./Pipeline/CmapESM_params.sh  " + foldpair_id  # Take one of the two !!! # ""./input/2qke.fasta 2qke
        print(esm_str)
        os.system(esm_str)  # run pipeline (should be a separate job!)

    if run_pipeline:  # entire pipeline: get alignments, contacts, predicted structure
        pipeline_str = "sbatch -o './Pipeline/" + foldpair_id + "/run_pipeline_for_" + foldpair_id + ".out' ./pipeline_get_params.sh " + \
                       fasta_file_name + " " + foldpair_id  # Take one of the two !!! # ""./input/2qke.fasta 2qke
        print(pipeline_str)
        os.system(pipeline_str)  # run pipeline (should be a separate job!)

    if plot_results:
        fasta_file_names = {pdbids[i][fold] + pdbchains[i][fold]: fasta_dir + "/" + foldpair_id + "/" + \
                                             pdbids[i][fold] + pdbchains[i][fold] + '.fasta' for fold in range(2)}  # Added chain to file ID
        msa_file = fasta_dir + "/" + foldpair_id + "/output_get_msa/DeepMsa.a3m"
        MSA = read_msa(msa_file)  # AlignIO.read(open(msa_file), "fasta")
        msa_pred_files = glob(fasta_dir + "/" + foldpair_id + "/output_cmap_esm/*.npy")
        n_cmaps = len(msa_pred_files)
        n_cmaps = min(3, n_cmaps)  # temp for debug !!
        msa_files = glob(fasta_dir + "/" + foldpair_id + "/output_msa_cluster/*.a3m")
        msa_clusters = {file.split("\\")[-1][:-4]: read_msa(file) for file in msa_files}

        # Filter 'bad' families: too shallow alignments (no clusters), same PDB ID, other reasons??
        if len(msa_files) == 0:
            print("Shallow alignment! No MSA Clusters! Skipping family")
            continue
        if pdbids[i][0] == pdbids[i][1]:
            print("Same PDB-ID for both folds! Might be buggy!" ) # Not supported yet!  Skipping family")
#            continue
#        if foldpair_id in ["1nqjB_1nqdA", "1qlnA_1h38D", "3l5nB_2a73B", "2bzyB_2lqwA", "4cmqB_4zt0C", "5tpnA_1g2cF"]:
#            print("Special problematic family! " + foldpair_id)
#            continue
#       if i < 15:
#            print("Ignore first!")
#            continue


        # First load files
        seqs = {}
        for fold in range(2):
            with open(fasta_file_names[pdbids[i][fold] + pdbchains[i][fold]], "r") as text_file:
                seqs[pdbids[i][fold] + pdbchains[i][fold]] = text_file.read().split("\n")[1]
        pairwise_alignment = Align.PairwiseAligner().align(seqs[pdbids[i][0] + pdbchains[i][0]],
                                                           seqs[pdbids[i][1] + pdbchains[i][1]])

        if max([len(seqs[fold]) for fold in seqs]) > 1024:
            print("Long sequence! Length = " + " > maximum supported length of 1024")
            continue

#        msa_transformer_pred = [None]*n_cmaps
#        for file in msa_pred_files:
#            print(file)
#            msa_transformer_pred[i] = np.genfromtxt(file)
        try:  # read in text format or python format
            msa_transformer_pred = {file.split("\\")[-1][14:-4] : np.genfromtxt(file) for file in msa_pred_files}
        except:
            msa_transformer_pred = {file.split("\\")[-1][14:-4]: np.load(file, allow_pickle=True) for file in msa_pred_files}

#        print("Predicted cmap sizes for:" + foldpair_id)
#        print([c.shape[0] for c in msa_transformer_pred.values()])

#        true_cmap = [None]*2
#        for fold in range(2):
#            true_cmap[fold] = np.genfromtxt(fasta_dir +
#                "/" + foldpair_id + "/" + pdbids[i][fold] + pdbchains[i][fold] + "_pdb_contacts.npy")  # load for pickle, genfromtxt for tab-delimited
        try:
            true_cmap = {pdbids[i][fold] + pdbchains[i][fold] : np.genfromtxt(fasta_dir +  # problem with first !!
                        "/" + foldpair_id + "/" + pdbids[i][fold] + pdbchains[i][fold] + "_pdb_contacts.npy").astype(int) for fold in range(2)}
#            print("Loaded text true!")
        except:
#            print("load binary true!")
            true_cmap = {pdbids[i][fold] + pdbchains[i][fold] : np.load(fasta_dir +  # problem with first !!
                        "/" + foldpair_id + "/" + pdbids[i][fold] + pdbchains[i][fold] + "_pdb_contacts.npy", allow_pickle=True).astype(int) for fold in range(2)}

#        print("seq lens:")
#        print( [len(seqs[fold]) for fold in seqs])
#        print("aligned lens:")
#        print(len((pairwise_alignment[0])))

#        print("True cmap sizes:")
#        print([c.shape[0] for c in true_cmap.values()])
#        print("plot")
#        print("Predicted cmap sizes:")
#        print([c.shape[0] for c in msa_transformer_pred.values()])


#        ttt = copy.deepcopy(true_cmap)
#        ttt['4rmbB'] = true_cmap['4rmbA']
#        ttt['4rmbA'] = true_cmap['4rmbB']

        # Match indices of the two folds:
#        print("Matching: ")
#        print(len(msa_transformer_pred))
#        print(true_cmap)
#        print(len(true_cmap))
#        match_true_cmap, match_predicted_cmaps = get_matching_indices_two_maps(pairwise_alignment, ttt, msa_transformer_pred)
        match_true_cmap, match_predicted_cmaps = get_matching_indices_two_maps(pairwise_alignment, true_cmap, msa_transformer_pred)

#        for fold in match_true_cmap.keys():
#            print(fold)
#            print(true_cmap[fold].shape)
#            print(match_true_cmap[fold].shape)
#        for fold in match_predicted_cmaps.keys():
#            print(fold)
#            print(match_predicted_cmaps[fold].shape)
#        plot_array_contacts_and_predictions(msa_transformer_pred, true_cmap)
        plot_array_contacts_and_predictions(match_predicted_cmaps, match_true_cmap, fasta_dir + "/Results/Figures/Cmap_MSA/" + foldpair_id)

        # load tree
        draw_tree_with_values(phytree, fasta_dir + "/Results/Figures/PhyTree/" + foldpair_id)

        # Collect :
        cmap_dists_vec[i] = compute_cmap_distances(match_predicted_cmaps)  # msa_transformer_pred)
        seqs_dists_vec[i] = np.mean(compute_seq_distances(msa_clusters))  # can take the entire sequence !
        num_seqs_msa_vec[i] = len(seqs)
#        print("Cmap dist: " + str(cmap_dists_vec[i]) + ", seq dist:" + str(seqs_dists_vec[i]))
#        break
        # next plotP

if plot_results:  # save results at the end:
    res_DF = pd.DataFrame(
        {'cmap_dists': cmap_dists_vec,
         'seq_dists': seqs_dists_vec,
         'n_msa': num_seqs_msa_vec
         })

    res_DF.to_csv(fasta_dir + "/Results/foldswitch_res.csv")

#import numpy as np
#cmap_pred = np.genfromtxt("Pipeline/1dzlA_5keqF/output_cmap_esm/msa_t__ShallowMsa_000.npy")
#cmap_true = np.genfromtxt("Pipeline/1dzlA_5keqF/1dzlA_pdb_contacts.npy")
#plot_array_contacts_and_predictions(cmap_pred, cmap_true)
