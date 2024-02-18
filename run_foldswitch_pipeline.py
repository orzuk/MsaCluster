# Rub pipeline from a list of PDBs
import copy

import pandas as pd
from phytree_utils import *
from protein_plot_utils import *
from utils import *
from MSA_Clust import *
import platform
import requests


def download_pdb(pdb_id):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)

    # Check if the request was successful (status code 200)
    if response.status_code == 200:
        filename = f"./pdb_files/{pdb_id}.pdb"
        with open(filename, 'w') as file:
            file.write(response.text)
        print(f"Downloaded {filename}")
    else:
        print(f"Failed to download PDB file for ID {pdb_id}. Status code: {response.status_code}")




# Run pipeline on a bunch of families (MSAs can be given, or read from file
# or generated on the fly)
# MSAs can be represented as a3m format







def run_fold_switch_pipeline(run_mode, foldpair_ids_to_run='ALL',output_dir ="Pipeline", pdbids_file="data/foldswitch_PDB_IDs_full.txt", run_job_mode="inline"):

    if type(pdbids_file) == str:                       # input as file
        with open(pdbids_file, "r") as file:           # read all pdb ids
            pdbids = [line.rstrip() for line in file]  # two per row
        foldpair_ids = [s.replace("\t", "_") for s in pdbids]

        pdbids    = [s.split("\t") for s in pdbids]
        pdbchains = [[s[0][-1], s[1][-1]] for s in pdbids]
        pdbids    = [[s[0][:-1], s[1][:-1]] for s in pdbids]




    n_fam = len(pdbids)  # number of families
    cmap_dists_vec, seqs_dists_vec, num_seqs_msa_vec = [None] * n_fam, [None] * n_fam, [None] * n_fam  # Results arrays
    if foldpair_ids_to_run == "ALL":
        foldpair_ids_to_run = foldpair_ids
    else:  # make a list
        if type(foldpair_ids_to_run) == str:
            foldpair_ids_to_run = [foldpair_ids_to_run]

    if run_mode == "plot":
        pymol.finish_launching(['pymol', '-c'])  # '-c' for command line (no GUI)


    # idx_test = 0
    # pred_vec = [0] * n_fam     # loop on MSAs
    for foldpair_id in foldpair_ids_to_run:
        try:
            # idx_test +=1
            # if idx_test>3:
            #     break
            output_pair_dir = f'{output_dir}/{foldpair_id}'
            if not os.path.exists(output_pair_dir):
                print("Mkdir: " + output_pair_dir)
                os.mkdir(output_pair_dir)
            if not os.path.exists(f'./{output_pair_dir}/chain_pdb_files'):
                print("Mkdir: " + f'./{output_pair_dir}/chain_pdb_files')
                os.mkdir(f'./{output_pair_dir}/chain_pdb_files')
            if not os.path.exists(f'./{output_pair_dir}/fasta_chain_files'):
                print("Mkdir: " + f'./{output_pair_dir}/fasta_chain_files')
                os.mkdir(f'./{output_pair_dir}/fasta_chain_files')


            i = foldpair_ids.index(foldpair_id)
            # download_pdb(pdbids[i][0])
            # download_pdb(pdbids[i][1])
            # create_chain_pdb_files(pdbids[i][0]+pdbchains[i][0], pdbids[i][1]+pdbchains[i][1], './pdb_files', f'./{output_pair_dir}/chain_pdb_files')

            get_fasta_chain_seq(f'./{output_pair_dir}/chain_pdb_files/{pdbids[i][0] + pdbchains[i][0]}.pdb', pdbids[i][0] + pdbchains[i][0],output_pair_dir)

    #        if i < 76:  # already done
    #            print("Already plotted")
    #            continue
            fasta_file_name = output_dir + "/" + foldpair_id + "/fasta_chain_files/" + pdbids[i][0]+pdbchains[i][0] + '.fasta'  # First file of two folds
    #        cur_family_dir = output_dir + "/" + foldpair_id
            print("Run: " + run_mode + " : " + foldpair_id + " : " + str(i) + " out of : " + str(len(foldpair_ids_to_run)))

            if run_job_mode == "inline":
                run_fold_switch_pipeline_one_family(run_mode, foldpair_id,pdbids[i],pdbchains[i],fasta_file_name)
            else:  # run job
                if run_mode == "get_msa":
                    run_str = "sbatch -o './Pipeline/" + foldpair_id + "/get_msa_for_" + foldpair_id + ".out' ./Pipeline/get_msa_params.sh " + fasta_file_name + " " + foldpair_id  # Take one of the two !!! # ""./input/2qke.fasta 2qke

                if run_mode == "cluster_msa":  # here do analysis of the results
                    run_str = "sbatch -o './Pipeline/" + foldpair_id + "/cluster_msa_for_" + foldpair_id + ".out' ./Pipeline/ClusterMSA_params.sh " + foldpair_id  # Take one of the two !!! # ""./input/2qke.fasta 2qke
                if run_mode == "run_cmap_esm":
                    run_str = "sbatch -o './Pipeline/" + foldpair_id + "/run_esm_for_" + foldpair_id + ".out' ./Pipeline/CmapESM_params.sh  " + output_dir+"/"+foldpair_id  # Take one of the two !!! # ""./input/2qke.fasta 2qke
                if run_mode == "run_esmfold":
                    run_str = f"sbatch /sci/labs/orzuk/steveabecassis/MsaCluster/Pipeline/RunEsmFold_params.sh  {foldpair_id}"
                if run_mode == "run_AF":  # run alpha-fold to predict structures
                    run_str = "sbatch -o './Pipeline/" + foldpair_id + "/run_AF_for_" + foldpair_id + ".out' ./Pipeline/RunAF_params.sh  " + foldpair_id  # Take one of the two !!! # ""./input/2qke.fasta 2qke
                if run_mode == "run_pipeline":
                    run_str = "sbatch   ./pipeline_get_params.sh " +  fasta_file_name + " " + output_dir+"/"+foldpair_id  # Take one of the two !!! # ""./input/2qke.fasta 2qke
                if run_mode == "plot":  # here do analysis of the results
                    cmap_dists_vec[i], seqs_dists_vec[i], num_seqs_msa_vec[i] = make_foldswitch_all_plots(pdbids[i], output_dir, foldpair_id, pdbchains[i], plot_tree_clusters)
                if run_mode == "tree":  # here do analysis of the results
                    run_str = "sbatch -o './Pipeline/" + foldpair_id + "/tree_reconstruct_for_" + foldpair_id + ".out' ./Pipeline/tree_reconstruct_params.sh " + foldpair_id  # Take one of the two !!! # ""./input/2qke.fasta 2qke
                if run_mode == 'Analysis':
                    run_str = f"sbatch  ./Analysis_params.sh  {foldpair_id}"

                print("Send job for " + run_mode + ":\n" + run_str)
                os.system(run_str)
        except:
            continue
    # for i in range(n_fam):
    #    cur_MSA = MSAs_dir[i] # This should be replaced by code generating/reading the MSA
    #    pred_vec[i] = predict_fold_switch_from_MSA_cluster(cur_MSA, clust_params)
    if run_mode == "plot":
        cmd.quit()

    res_DF = pd.DataFrame(
        {'cmap_dists': cmap_dists_vec,
         'seq_dists': seqs_dists_vec,
         'n_msa': num_seqs_msa_vec
         })
    return res_DF # return results


# Run inline one family. Shouldn't get command line arguments but use function input!!!
def run_fold_switch_pipeline_one_family(run_mode, foldpair_id, pdbids, pdbchains, fasta_file_name):
    cmap_dists_vec, seqs_dists_vec, num_seqs_msa_vec = [None]*3
    cur_family_dir = output_dir + "/" + foldpair_id
    if run_mode == "load_seq_and_struct":  #      if load_seq_and_struct or run_pipeline:  # also for entire pipeline
        run_str = ''  # no pipeline in this mode !!!!
        for fold in range(2):
            if not os.path.exists(cur_family_dir):
                print("Mkdir: " + cur_family_dir)
                os.mkdir(cur_family_dir)
            print("Get seq + struct for " + pdbids[fold] + ", " + " out of " + str(n_fam-1) )
            fasta_file_name = output_dir + "/" + foldpair_id + "/" + pdbids[fold] + pdbchains[fold] + '.fasta'  # added chain to file ID
    #             # Finally, make a contact map from each pdb file:
    #             # Read structure in slightly different format
                 # New option: extract sequence and structure togehter. Remove from sequence the residues without contacts

            print("STRUCTURE TYPE: ")
            print(type(PDBxFile.read(rcsb.fetch(pdbids[fold], "cif"))))
            print(type(get_structure(PDBxFile.read(rcsb.fetch(pdbids[fold], "cif")))[0]))

            with open("good_temp_tm_align_structs.pkl", "wb") as f:
                pickle.dump([get_structure(PDBxFile.read(rcsb.fetch(pdbids[fold], "cif")))[0], pdbchains[fold]], f)
            with open('good_temp_tm_align_structs.pkl', 'rb') as f:
                s_good, c_good = pickle.load(f)

            pdb_dists, pdb_contacts, pdb_seq, pdb_good_res_inds, cbeta_coord = read_seq_coord_contacts_from_pdb(   # extract distances from pdb file
                s_good, chain=c_good)

#            pdb_dists, pdb_contacts, pdb_seq, pdb_good_res_inds, cbeta_coord = read_seq_coord_contacts_from_pdb(   # extract distances from pdb file
#                get_structure(PDBxFile.read(rcsb.fetch(pdbids[fold], "cif")))[0], chain=pdbchains[fold])
            with open(fasta_file_name, "w") as text_file:  # save to fasta file. Take the correct chain
                text_file.writelines(["> " + pdbids[fold].upper() + ":" + pdbchains[fold].upper() + '\n',
                                            pdb_seq ])
            print(cur_family_dir + "/" + pdbids[fold] + pdbchains[fold] + "_pdb_contacts.npy")
            np.save(cur_family_dir + "/" + pdbids[fold] + pdbchains[fold] + "_pdb_contacts.npy", pdb_contacts)  # save true contacts (binary format)
            print("FINISHED read_seq_coord_contacts_from_pdb ALL IS GOOD!")
    if run_mode == "get_msa":
        run_str = "python3. / get_msa.py " + fasta_file_name + " ./Pipeline/" + foldpair_id + "/output_get_msa - name 'DeepMsa'"
    if run_mode == "cluster_msa":
        run_str = "python3  ./ClusterMSA_moriah.py --keyword ShallowMsa -i ./Pipeline/" + foldpair_id + \
                  "/output_get_msa/DeepMsa.a3m  -o ./Pipeline/" + foldpair_id + "/output_msa_cluster"
    if run_mode == "run_esm":
        run_str = 'python3  ./runESM.py  --input_msas ./Pipeline/' + foldpair_id + \
                  '/output_msa_cluster -o ./Pipeline/' + foldpair_id + '/output_cmap_esm'
        print(run_str)
        os.system(run_str) # Run the runESM function

    if run_mode == "run_AF":  # run alpha-fold to predict structures
        run_str = "python3 runAF.py -input ./Pipeline/" + foldpair_id + ' -o ./Pipeline/' + foldpair_id + '/output_AF'  # Take one of the two !!! # ""./input/2qke.fasta 2qke
        print(run_str)
        os.system(run_str) # Run the runESM function

    if run_mode == "run_esmfold":
        run_str = 'python3  ./ESMFoldHF.py  -input ./Pipeline/' + foldpair_id + \
                  '/output_msa_cluster -o ./Pipeline/' + foldpair_id + '/output_esm_fold'
        print(run_str)
        os.system(run_str) # Run the runESM function

    if run_mode == "tree":
        msa_file = 'Pipeline/' + foldpair_id + '/output_get_msa/DeepMsa.a3m'
        phytree_from_msa(msa_file, output_tree_file= 'Pipeline/' + foldpair_id + \
                    '/output_phytree/' + os.path.basename(msa_file).replace(".a3m", "_tree.nwk"))
        run_str = ''
#        run_str = 'python  ./run_tree_reconstruct.py  --method distance  --input_msas ./Pipeline/' + foldpair_id + \
#                  '/output_get_msa/DeepMsa.a3m -o ./Pipeline/' + foldpair_id + '/output_phytree'
    if run_mode == "ancestral": # perform ancestral reconstuction
        msa_file = 'Pipeline/' + foldpair_id + '/output_get_msa/DeepMsa.a3m'
        anc_output_file = 'Pipeline/' + foldpair_id + '/output_phytree/' + \
                          os.path.basename(msa_file).replace(".a3m", "_anc_seq.a3m")
        output_tree_file = 'Pipeline/' + foldpair_id + \
                    '/output_phytree/' + os.path.basename(msa_file).replace(".a3m", "_tree.nwk")
        reconstruct_ancestral_sequences(output_tree_file, msa_file, anc_output_file)
    if run_mode == "plot":
        cmap_dists_vec, seqs_dists_vec, num_seqs_msa_vec, concat_scores = \
            make_foldswitch_all_plots(pdbids, output_dir, foldpair_id, pdbchains, plot_tree_clusters)
        run_str = ''  # no plotting in this mode !!!!
    if run_mode == "run_pipeline":
        run_str = ''  # no pipeline in this mode !!!!

    print("Run command line for " + run_mode + ":\n" + run_str)
    return cmap_dists_vec, seqs_dists_vec, num_seqs_msa_vec
#    os.system(run_str)

from Bio import Align

# problematic_families = ["1nqjB_1nqdA", "1qlnA_1h38D", "3l5nB_2a73B", "2bzyB_2lqwA", "4cmqB_4zt0C", "5tpnA_1g2cF"]  # second, third is large, no/one cmaps were generated !
# problematic_families = ["2bzyB_2lqwA", "4cmqB_4zt0C", "5tpnA_1g2cF"]  # 2bzyB_2lqwA bad (missing ) true cmap!!
foldpair_ids_to_run = 'ALL'  # '3j7vG_3j7wB' # '2vfxL_3gmhL' # '1xjuB_1xjtA'  # "ALL"
if platform.system() == "Linux":
    print("Run on cluster command line")
    run_mode = sys.argv[1]
    if len(sys.argv) > 2:
        foldpair_ids_to_run = sys.argv[2]  # enable running for a specific family (default is running on all of them)
    run_job_mode = "job"
else:
    print("Run on windows")
    run_mode = "plot"  # "load_seq_and_struct" #    "plot" # "run_esmfold"   # "plot"  # "load"  # "run_esm" # "plot" # "run_esm"  # sys.argv[1]
    run_job_mode = "inline"
    foldpair_ids_to_run = "1dzlA_5keqF"  #   "4ydqB_4twaA"  #  "3t5oA_4a5wB" # "3meeA_4b3oB"  # "2kb8A_6vw2A"  #  problematic, needs padding !
    plot_tree_clusters = True

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
if run_mode == "run_esmfold":
    run_esmfold = True  # run entire pipeline
    print("Run ESM transformer Fold for all families")
if run_mode == "run_pipeline":
    run_pipeline = True  # run entire pipeline
    print("Run Entire pipeline for all families")
if run_mode == "plot":  # here do analysis of the results
    plot_results = True
if run_mode == "tree":  # here do analysis of the results
    tree_reconstruct = True

# pdb_datadir = "Pipeline/pdb_files"  # where to store all PDB files
output_dir = "Pipeline"
pdbids_file = "data/foldswitch_PDB_IDs_full.txt"   # file with all pdb ids

with open(pdbids_file, "r") as file:  # read all pdb ids
    pdbids = [line.rstrip() for line in file]  # two per row
foldpair_ids = [s.replace("\t", "_") for s in pdbids]

pdbids = [s.split("\t") for s in pdbids]
pdbchains = [[s[0][-1], s[1][-1]] for s in pdbids]
pdbids = [[s[0][:-1], s[1][:-1]] for s in pdbids]

n_fam = len(pdbids)  # number of families
cmap_dists_vec   = [None]*n_fam  # Results arrays
seqs_dists_vec   = [None]*n_fam
num_seqs_msa_vec = [None]*n_fam
if foldpair_ids_to_run == "ALL":
    foldpair_ids_to_run = foldpair_ids
else:  # make a list
    if type(foldpair_ids_to_run) == str:
        foldpair_ids_to_run = [foldpair_ids_to_run]

res_DF = run_fold_switch_pipeline(run_mode,foldpair_ids_to_run,output_dir="Pipeline",pdbids_file="data/foldswitch_PDB_IDs_full.txt",run_job_mode=run_job_mode)
# res_DF.to_csv(output_dir + "/Results/foldswitch_res.csv")


### TEMP CODE FOR TRYING STUFF
#with open('tree_draw.pkl', 'rb') as f:  # Python 3: open(..., 'rb')
#    phytree_file, tree_outfile, node_values = pickle.load(f)
#visualize_tree_with_heatmap(phytree_file, node_values, tree_outfile)