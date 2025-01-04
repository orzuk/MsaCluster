# Rub pipeline from a list of PDBs

from scripts.protein_plot_utils import *
from scripts.energy_utils import *
from utils import *
import platform
import requests
# from biotite.structure.io.pdbx import get_structure
from biotite.structure import filter_amino_acids
# from biotite.structure.io import PDBxFile
# import biotite.structure as bs
from biotite.structure.io.pdbx import get_structure

# from Bio.PDB import PDBxFile
from Bio.PDB.MMCIFParser import MMCIFParser


# Run pipeline on a bunch of families (MSAs can be given, or read from file
# or generated on the fly)
# MSAs can be represented as a3m format


def run_fold_switch_pipeline(run_mode, foldpair_ids_to_run='ALL',output_dir ="Pipeline", pdbids_file="data/foldswitch_PDB_IDs_full.txt", run_job_mode="inline"):
    if run_mode == "help":
        print("Available run modes for fold-switch pipeline analysis:")
        print("- load: Load PDB files, contact maps, and fasta sequences.")
        print("- get_msa: Generate MSAs for sequences.")
        print("- cluster_msa: Cluster MSAs.")
        print("- run_cmap_esm: Run ESM (MSA?-)transformer on sequences.")
        print("- run_esmfold: Run ESMFold for structure prediction.")
        print("- run_pipeline: Execute the entire pipeline.")
        print("- plot: Generate plots for fold-switch data.")
        print("- tree: Reconstruct phylogenetic trees.")
        print("- compute_deltaG: Compute free energy of the structure + sequence using Rosetta")
        print("- Analysis: Perform custom analyses.")
        return

    print("START RunFoldSwitch PIPELINE!!")
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
        print("Running plot!!")
        pymol.finish_launching(['pymol', '-cq'])  # '-c' for command line, q for quiet (no GUI)

        print("Running plot after pymol launch!!")
    # idx_test = 0
    # pred_vec = [0] * n_fam     # loop on MSAs
    for foldpair_id in foldpair_ids_to_run:
        print("Run on fold pairs ", foldpair_id)
        try: # Change to try catch !!!
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
            # create_chain_pdb_files(pdbids[i][0]+pdbchains[i][0], pdbids[i][1]+pdbchains[i][1], './pdb_files', f'./{output_pair_dir}/chain_pdb_files')
            print("Search chain: ", f'./{output_pair_dir}/{pdbids[i][0]}.pdb')
            get_fasta_chain_seq(f'./{output_pair_dir}/{pdbids[i][0]}.pdb', pdbids[i][0] + pdbchains[i][0], output_pair_dir)  # Why only first fold here?
            fasta_file_name = output_dir + "/" + foldpair_id + "/fasta_chain_files/" + pdbids[i][0]+pdbchains[i][0] + '.fasta'  # First file of two folds
            cur_family_dir = output_dir + "/" + foldpair_id
            print("Run: " + run_mode + " : " + foldpair_id + " : " + str(i) + " out of : " + str(len(foldpair_ids_to_run)))

            if run_job_mode == "inline":
                run_fold_switch_pipeline_one_family(run_mode, foldpair_id, pdbids[i], pdbchains[i],fasta_file_name)
            else:  # run job
                if run_mode == "load":
                    run_str = ''  # no pipeline in this mode !!!!
                    print("Enter load function:::")
                    print(i, cur_family_dir, foldpair_id, pdbids[i])
                    load_seq_and_struct(cur_family_dir, foldpair_id, pdbids[i], pdbchains[i])
                    print("Finished load function")
                if run_mode == "get_msa":
                    run_str = "sbatch -o './Pipeline/" + foldpair_id + "/get_msa_for_" + foldpair_id + ".out' ./Pipeline/get_msa_params.sh " + fasta_file_name + " " + foldpair_id  # Take one of the two !!! # ""./input/2qke.fasta 2qke

                if run_mode == "cluster_msa":  # here do analysis of the results
                    run_str = "sbatch -o './Pipeline/" + foldpair_id + "/cluster_msa_for_" + foldpair_id + ".out' ./Pipeline/ClusterMSA_params.sh " + foldpair_id  # Take one of the two !!! # ""./input/2qke.fasta 2qke
                if run_mode == "run_cmap_esm":  # MSA Transformer (??)
                    run_str = "sbatch -o './Pipeline/" + foldpair_id + "/run_esm_for_" + foldpair_id + ".out' ./Pipeline/CmapESM_params.sh  " + output_dir+"/"+foldpair_id  # Take one of the two !!! # ""./input/2qke.fasta 2qke
                if run_mode in ("run_cmap_esm", "run_esmfold"):
                    run_str = f"sbatch /sci/labs/orzuk/steveabecassis/MsaCluster/Pipeline/RunEsmFold_params.sh  {foldpair_id}"
                if run_mode == "run_AF":  # run alpha-fold to predict structures
                    run_str = "sbatch -o './Pipeline/" + foldpair_id + "/run_AF_for_" + foldpair_id + ".out' ./Pipeline/RunAF_params.sh  " + foldpair_id  # Take one of the two !!! # ""./input/2qke.fasta 2qke
                if run_mode == "compute_deltaG":  # Compute free energy !!! NEW !!
                    run_str = ''
                    pdb_pair_files = [("Pipeline/" + p + "/" + p[:4] + ".pdb" ,  "Pipeline/" + p + "/" + p[-5:-1] + ".pdb") for p in foldpair_ids_to_run]
                    deltaG_output_file = "Pipeline/output_deltaG/deltaG_results.txt"
#                    compute_global_and_residue_energies(pdb_pair_files, foldpair_ids_to_run, output_dir)
                    pdb_pair_files = [("Pipeline/" + foldpair_id + "/" + foldpair_id[:4] + ".pdb",
                                       "Pipeline/" + foldpair_id + "/" + foldpair_id[-5:-1] + ".pdb") ]
                    compute_global_and_residue_energies(pdb_pair_files, [foldpair_id], "Pipeline/output_deltaG")

#                    run_compute_deltaG_with_output(pdb_pair_files, foldpair_ids_to_run, deltaG_output_file)

                    '''
                    for foldpair_id in foldpair_ids_to_run:
                        cur_family_dir = os.path.join(output_dir, foldpair_id)
                        print("Family pair dir: ", cur_family_dir)
                        print("pdbids: ", pdbids, " pdbchains: ", pdbchains)
                        pdb_files = [
                            os.path.join(cur_family_dir, f"{pdbid}.pdb")
                            for pdbid, chain in zip(pdbids[i], pdbchains[i])
                        ]
                        print("pdbfiles: ", pdb_files)
                        for pdb_file in pdb_files:
                            print(f"Computing ΔG for {pdb_file} using PyRosetta...")
                            deltaG = compute_deltaG_with_pyrosetta(pdb_file)
                            print(f"Computed ΔG for {pdb_file}: {deltaG:.2f} kcal/mol")
                    '''
                if run_mode == "run_pipeline":
                    run_str = "sbatch   ./pipeline_get_params.sh " +  fasta_file_name + " " + output_dir+"/"+foldpair_id  # Take one of the two !!! # ""./input/2qke.fasta 2qke
                if run_mode == "plot":  # here do analysis of the results
                    run_str = ''
                    pymol.finish_launching(['pymol', '-cq'])
                    print("Running make_plot from outside33 again!")
                    print(pdbids[i])
                    print(output_dir)
                    print(foldpair_id)
                    print("pdbchjain", pdbchains[i])
                    print("plot tree clusters:", plot_tree_clusters)

                    cmap_dists_vec[i], seqs_dists_vec[i], num_seqs_msa_vec[i] = make_foldswitch_all_plots(pdbids[i], output_dir, foldpair_id, pdbchains[i], plot_tree_clusters)
                    print("Finished Running make_plot from outside333")
                if run_mode == "tree":  # here do analysis of the results
                    run_str = "sbatch -o './Pipeline/" + foldpair_id + "/tree_reconstruct_for_" + foldpair_id + ".out' ./Pipeline/tree_reconstruct_params.sh " + foldpair_id  # Take one of the two !!! # ""./input/2qke.fasta 2qke
                if run_mode == 'Analysis':
                    run_str = f"sbatch  /sci/labs/orzuk/steveabecassis/MsaCluster/Pipeline/Analysis.sh  {foldpair_id}"
                if run_str != '':
                    print("Send job for " + run_mode + ":\n" + run_str)
                    os.system(run_str)
        except:  # Uncomment !!
            print("Failed " + run_mode + " function")
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
    if run_mode == "load_seq_struct":  #      if load_seq_and_struct or run_pipeline:  # also for entire pipeline
        run_str = ''  # no pipeline in this mode !!!!
        load_seq_and_struct(cur_family_dir, foldpair_id)
    if run_mode == "get_msa":
        run_str = "python3. / get_msa.py " + fasta_file_name + " ./Pipeline/" + foldpair_id + "/output_get_msa - name 'DeepMsa'"
    if run_mode == "cluster_msa":
        run_str = "python3  ./ClusterMSA_moriah.py --keyword ShallowMsa -i ./Pipeline/" + foldpair_id + \
                  "/output_get_msa/DeepMsa.a3m  -o ./Pipeline/" + foldpair_id + "/output_msa_cluster"
    if run_mode == "run_cmap_esm":  # msa-transformer ??
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
    if run_mode == "ancestral": # perform ancestral reconstruction
        msa_file = 'Pipeline/' + foldpair_id + '/output_get_msa/DeepMsa.a3m'
        anc_output_file = 'Pipeline/' + foldpair_id + '/output_phytree/' + \
                          os.path.basename(msa_file).replace(".a3m", "_anc_seq.a3m")
        output_tree_file = 'Pipeline/' + foldpair_id + \
                    '/output_phytree/' + os.path.basename(msa_file).replace(".a3m", "_tree.nwk")
        reconstruct_ancestral_sequences(output_tree_file, msa_file, anc_output_file)
    if run_mode == "plot":  # Plotting, works only in windows?
        cmap_dists_vec, seqs_dists_vec, num_seqs_msa_vec, concat_scores = \
            make_foldswitch_all_plots(pdbids, output_dir, foldpair_id, pdbchains, plot_tree_clusters)
        run_str = ''  # no plotting in this mode !!!!
    if run_mode == "run_pipeline":
        run_str = ''  # no pipeline in this mode !!!!

    print("Run command line for " + run_mode + ":\n" + run_str)
    return cmap_dists_vec, seqs_dists_vec, num_seqs_msa_vec
#    os.system(run_str)

# problematic_families = ["1nqjB_1nqdA", "1qlnA_1h38D", "3l5nB_2a73B", "2bzyB_2lqwA", "4cmqB_4zt0C", "5tpnA_1g2cF"]  # second, third is large, no/one cmaps were generated !
# problematic_families = ["2bzyB_2lqwA", "4cmqB_4zt0C", "5tpnA_1g2cF"]  # 2bzyB_2lqwA bad (missing ) true cmap!!
foldpair_ids_to_run = 'ALL'  # '3j7vG_3j7wB' # '2vfxL_3gmhL' # '1xjuB_1xjtA'  # "ALL"
if platform.system() == "Linux":
    print("Run on cluster command line")
    run_mode = sys.argv[1]
    if len(sys.argv) > 2:
        foldpair_ids_to_run = sys.argv[2]  # enable running for a specific family (default is running on all of them)
    run_job_mode = "job"
    plot_tree_clusters = False
else:
    print("Run on windows")
    run_mode = "plot"  # "load_seq_and_struct" #    "plot" # "run_esmfold"   # "plot"  # "load"  # "run_esm" # "plot" # "run_esm"  # sys.argv[1]
    run_job_mode = "inline"
    foldpair_ids_to_run = "1dzlA_5keqF"  #   "4ydqB_4twaA"  #  "3t5oA_4a5wB" # "3meeA_4b3oB"  # "2kb8A_6vw2A"  #  problematic, needs padding !
    plot_tree_clusters = True

# print("Running on: " + foldpair_ids_to_run)
run_pipeline, get_msa, cluster_msa, tree_reconstruct, run_esm, load_seq_struct, plot_results = [False]*7  # run entire pipeline or parts of it/plot ...

# can't use match (only from python 3.10)
# match run_mode:
if run_mode == "load":
    load_seq_struct = True  # just get from pdb the sequence and 3D structure for each protein
    print("Load pdb files, contact maps and fasta sequences for all families")
if run_mode == "get_msa":  # here do analysis of the results
    get_msa = True
if run_mode == "cluster_msa":  # here do analysis of the results
    cluster_msa = True
if run_mode == "run_cmap_esm":
    run_cmap_esm = True  # run esm
    print("Run ESM transformer for all families")
if run_mode == "run_esmfold":
    run_esmfold = True  # run esm-fold
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

res_DF = run_fold_switch_pipeline(run_mode, foldpair_ids_to_run,
                                  output_dir="Pipeline", pdbids_file="data/foldswitch_PDB_IDs_full.txt", run_job_mode=run_job_mode)
# res_DF.to_csv(output_dir + "/Results/foldswitch_res.csv")


### TEMP CODE FOR TRYING STUFF
#with open('tree_draw.pkl', 'rb') as f:  # Python 3: open(..., 'rb')
#    phytree_file, tree_outfile, node_values = pickle.load(f)
#visualize_tree_with_heatmap(phytree_file, node_values, tree_outfile)