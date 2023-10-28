# Pipeline for clustering alignemnts, predicting contacts for each cluster and comparing them
# Goal is to identify evidence for difference in co-evolution in the different clusters,
# possible indicative of differences in the structure
from polyleven import levenshtein
import random

import Levenshtein as pylev

import subprocess

import seaborn as sns
from protein_utils import *

sys.path.append('alphafold')


# Run clustering of the Multiple Sequence alignemnts,
# Default is 2 clusters (?), but can be more
# USe the method from AF-clust
def cluster_MSAs(MSA, clust_params):
    # Run ClusterMSA script (??)

    # USe AF-cluster script
    if clust_params == True:  # Run clustering
        AF_cluster_str = 'python ../AF_cluster/scripts/ClusterMSA.py EX -i ' + \
        '../AF_cluster/data_sep2022/00_KaiB/2QKEE_colabfold.a3m -o subsampled_MSAs'  # change later to input MSA
        subprocess.call(AF_cluster_str)

    clusters_dir = "../AF_Cluster/subsampled_MSAs"
    clusters_file_names = glob(clusters_dir + "/EX_*a3m")
    n_clusters = len(clusters_file_names)
    print(n_clusters)
    MSA_clusters = [None]*n_clusters  # replace by reading from output file
    for i in range(n_clusters):
        MSA_clusters[i] = AlignIO.read(open(clusters_file_names[i]), "fasta")



#    subprocess.call('python hello.py') # small check
    return MSA_clusters, clusters_file_names  # an object containing the partition of the MSA into multiple clusters


# Compute pairwise distances between the different contact maps
def compute_cmap_distances(cmaps, cmap_main = []):
    D = 0
    n_maps = len(cmaps)  # number of clusters/contact maps

    keys_list = list(cmaps.keys())
    if len(cmap_main) == 0:  # no centroid, calculate it
        cmap_main = cmaps[keys_list[0]]
        for i in range(1, n_maps):
            cmap_main += cmaps[keys_list[i]]
        cmap_main = cmap_main / n_maps

    for cur_cmap in cmaps.keys():  # Compute distance to centroid
        D += sum(cmaps[cur_cmap] - cmap_main)**2
    return D/n_maps  # normalize


# Compute pairwise distances between the different clusters
def compute_seq_distances(MSA_clust):
    n_clusters = len(MSA_clust)
    D = np.zeros((n_clusters, n_clusters))

#    for i in range(n_clusters):
#        summary_align = AlignInfo.SummaryInfo(MSA_clust[i])
#        PSSM = summary_align.MSA_clust.pos_specific_score_matrix()  # Code here
#
#    avg_dist_to_query = np.mean([1-levenshtein(x, query_['sequence'].iloc[0])/L for x in df.loc[df.dbscan_label==-1]['sequence'].tolist()])
#    lprint('avg identity to query of unclustered: %.2f' % avg_dist_to_query,f)

    keys_list = list(MSA_clust.keys())
    max_seq_per_cluster = 10  # maximum number of sequences per cluster
    for i in range(n_clusters):  # loop on pairs of clusters
        for j in range(i, n_clusters):
            n_i = len(MSA_clust[keys_list[i]])
            n_j = len(MSA_clust[keys_list[j]])
            II = random.sample(range(n_i), min(n_i, max_seq_per_cluster))
            JJ = random.sample(range(n_j), min(n_j, max_seq_per_cluster))

            for seq_i in II:
                for seq_j in JJ:
                    D[i, j] += levenshtein(str(MSA_clust[keys_list[i]][seq_i][1]), str(MSA_clust[keys_list[j]][seq_j][1]))  # compare sequences
            D[i, j] = D[i, j] / (len(II)*len(JJ))  # normalize
            D[j, i] = D[i, j]  # make symmetric

    return D  # average sequence distance between


# Predict if a family of proteins is fold-switching or is having a single structure,
# from the co-evolutionary patterns of the MSA
def predict_fold_switch_from_MSA_cluster(MSA, clust_params):
    MSA_clust, MSA_names = cluster_MSAs(MSA, clust_params)  # first cluster the MSAs
    n_clust = len(MSA_clust)  # number of clusters (can be set as a parameter)
    print("Compute sequence similarity:")
    seq_dist = compute_seq_distances(MSA_clust)  # sequence similarity of clusters

    # Plot distance matrix
    df_seq_distances = pd.DataFrame(seq_dist).sort_index().sort_index(axis=1)
    sns.heatmap(df_seq_distances)
    plt.show()  # display the plot

    from collections import defaultdict
    print("Compute cmap and their similarities:")
    cmaps = MSA_transformer(MSA_clust, MSA_names)
#    cmap = [None]*n_clust
#    for i in range(n_clust):  # loop on clusters
#        cmap[i] = MSA_transformer(MSA_clust[i])  # compute pairwise attention map for cluster
    print("Compute cmap distances:")
    cmap_dist = compute_cmap_distances(cmap)  # similarity of cluster contact maps

    fold_switch_pred_y = (cmap_dist - clust_params@beta * seq_dist > 0)  # replace by a learned function

    return fold_switch_pred_y, cmap_dist, seq_dist


# Convert MSAs to string format, where each sequence is a tuple of two,
# to serve as input for MSA transformer
def MSA_to_str_format(MSAs, MSAs_names):
    num_msas = len(MSAs)
    MSAs_str = {}  # [None]*num_msas

    for i in range(num_msas):
        MSAs_str[MSAs_names[i]] = [None] * len(MSAs[i])
        for j in range(len(MSAs[i])):
            MSAs_str[MSAs_names[i]][j] = (MSAs[2][1].name, str(MSAs[i][j].seq))

    return MSAs_str


# PDB_IDS = ["1a3a", "5ahw", "1xcr"]
#
# structures = {
#    name.lower(): get_structure(PDBxFile.read(rcsb.fetch(name, "cif")))[0]
#    for name in PDB_IDS
# }

#contacts = {
#    name: contacts_from_pdb(structure, chain="A")
#    for name, structure in structures.items()
#}


# Compute attention map using MSA transformer
def MSA_transformer(MSAs, MSAs_names, true_contacts = {}):
    print("Load model:")
    #
    msa_transformer, msa_transformer_alphabet = esm.pretrained.esm_msa1b_t12_100M_UR50S()  # MSA-transformer model (1Billion parameter). Problem: Memory!!!
#    msa_transformer = msa_transformer.eval().cuda()  # if cude is available
#    msa_transformer, msa_transformer_alphabet = esm.pretrained.esm2_t33_650M_UR50D()  # smaller newer model?
    print("Model eval:")
    msa_transformer = msa_transformer.eval()
    print("Batch convert:")
    msa_transformer_batch_converter = msa_transformer_alphabet.get_batch_converter()

#    batch_converter = msa_transformer_alphabet.get_batch_converter()
#    msa_transformer = msa_transformer.eval().cuda()
#    model.eval()  # disables dropout for deterministic results

    print("Get tokens:")
#    msa_transformer_batch_converter = msa_transformer_alphabet.get_batch_converter()

#    batch_labels, batch_strs, batch_tokens = batch_converter(MSA)
#    print("Get tokens:")
#    batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

    print("Convert MSAs:")
    MSAs_str = MSA_to_str_format(MSAs, MSAs_names)

    print("Get results:")
    msa_transformer_predictions = {}
    msa_transformer_results = []
    for name, inputs in MSAs_str.items():  # Get list of MSAs
        print("Compute cmap for MSA:" + name + " out of : " + str(len(MSAs_names)))
        inputs = greedy_select(inputs, num_seqs=min(len(inputs), 128))  # can change this to pass more/fewer sequences
        msa_transformer_batch_labels, msa_transformer_batch_strs, msa_transformer_batch_tokens = msa_transformer_batch_converter(
            [inputs])  # should accept list of tuples ???
        msa_transformer_batch_tokens = msa_transformer_batch_tokens.to(next(msa_transformer.parameters()).device)
        msa_transformer_predictions[name] = msa_transformer.predict_contacts(msa_transformer_batch_tokens)[0].cpu()
        metrics = {"id": name, "model": "MSA Transformer (Unsupervised)"}
        if(len(true_contacts)>0): # if there are true contacts
            metrics.update(evaluate_prediction(msa_transformer_predictions[name], true_contacts[name]))
        msa_transformer_results.append(metrics)
    msa_transformer_results = pd.DataFrame(msa_transformer_results)
    print("Token representations:")
    print(msa_transformer_results.head())  # Instead of display


#    with torch.no_grad():
#        results = model(batch_tokens, repr_layers=[33], return_contacts=True)
#    token_representations = results["representations"][33]
#    return token_representations

    return msa_transformer_predictions, msa_transformer_results  # return compact and detailed results

# load numpy array:
# np.loadtxt('/Users/steveabecassis/Desktop/PipelineTest/output_pipeline_1jfk/esm_cmap_output/msa_t__cluster_000.npy')



# Plot the true vs. predicted contact map, predict for each contact if it is:
# 1. Present in both
# 2. Present in first
# 3. Present in second
# 4. Absent
def match_predicted_and_true_contact_maps(cmap_clusters, cmap_true):
    return []

# Run pipeline on a bunch of families (MSAs can be given, or read from file
# or generated on the fly)
# MSAs can be represented as a3m format
def predict_fold_switch_pipeline(MSAs_dir, clust_params):
    # loop on MSAs
    n_fam = len(MSA)

    pred_vec = [0] * n_fam
    for i in range(n_fam):
        cur_MSA = MSAs_dir[i] # This should be replaced by code generating/reading the MSA
        pred_vec[i] = predict_fold_switch_from_MSA_cluster(cur_MSA, clust_params)




# Taken from here:
# https://stackoverflow.com/questions/76682974/similarity-score-sequence-alignment
def lev_distance_matrix(seqs):
    """Calculate Levenshtein distance and ratio metrics
       on input pair of strings.
    """
    seqs = sorted(seqs)

    return {
        seqs[0]: {
            seqs[1]: {
                "distance": pylev.distance(*seqs),
                "ratio": pylev.ratio(*seqs),
            }
        }
    }


# Main: test
# t_init = time.time()

# MSA_clust, MSA_names = cluster_MSAs([], False)

# predict_fold_switch_from_MSA_cluster([], False)