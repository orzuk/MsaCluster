from utils.utils import *
from utils.protein_utils import *
from utils.phytree_utils import *
import re
import colabfold
print("colab dir: " , colabfold.__file__)


# compute_tmscore('Pipeline/1dzlA_5keqF/1dzl.pdb', 'Pipeline/1dzlA_5keqF/AF_preds/ShallowMsa_000_unrelaxed_rank_005_alphafold2_ptm_model_1_seed_000.pdb', 'A', 'A')


with (open('tree_clusters.pkl', 'rb') as f):  # Python 3: open(..., 'rb')
     clusters_subtree, concat_scores, out_tree_file, tmscores_df, phytree_file, representative_cluster_leaves, \
     ete_leaves_node_values, ete_leaves_cluster_ids = pickle.load(f)

msa_file = "Pipeline/1dzlA_5keqF/output_get_msa/DeepMsa.a3m"
seqs_IDs, seqs = load_fasta(msa_file)
seqs = [''.join([x for x in s if x.isupper() or x == '-']) for s in seqs]  # remove lowercase letters in alignment
print("MSA seqs_IDS=", seqs_IDs, "\n Len: ", len(seqs_IDs))

print("ete_leaves_node_values: ", ete_leaves_node_values, ete_leaves_node_values.shape)
print("ete_leaves_cluster_ids: ", ete_leaves_cluster_ids)
print("ete_leaves_cluster_ids_unique: ", set(ete_leaves_cluster_ids.values()))
print("len: ", len(ete_leaves_cluster_ids.values()), " len unique: ", len(set(ete_leaves_cluster_ids.values())))
print("concat_scores: ", concat_scores)
print("tmscores_df", tmscores_df)
print("phytree_file=", phytree_file)
bio_tree = Phylo.read(phytree_file, "newick")  # This is different from write_newick_with_quotes !!!!
print("Convert to ete3 tree:")
ete_tree = convert_biopython_to_ete3(bio_tree)
# print("ete_tree: ", ete_tree)
print("representative_cluster_leaves", representative_cluster_leaves)
ll = list(representative_cluster_leaves.values())
ll.sort()
print("representative_cluster_leaves values sorted: ", ll)


visualize_tree_with_heatmap(clusters_subtree, concat_scores, 'temp_local_tree.png')

# esmf_df = pd.read_csv(ESMF_MODEL_FILE)
# print("ESM DF: ", esmf_df)

# msa_trans_df = pd.read_csv(MSA_TRANS_MODEL_FILE)
# print("MSA-Trans DF: ", msa_trans_df)


# node_names = ['1dzlA_5keqF/output_cmap_esm/msa_t__MSA_deep',
#       '1dzlA_5keqF/output_cmap_esm/msa_t__ShallowMsa_000',
#       '1dzlA_5keqF/output_cmap_esm/msa_t__ShallowMsa_001',
#       '1dzlA_5keqF/output_cmap_esm/msa_t__ShallowMsa_002',
#       '1dzlA_5keqF/output_cmap_esm/msa_t__ShallowMsa_003',
#       '1dzlA_5keqF/output_cmap_esm/msa_t__ShallowMsa_004',
#       '1dzlA_5keqF/output_cmap_esm/msa_t__ShallowMsa_005',
#       '1dzlA_5keqF/output_cmap_esm/msa_t__ShallowMsa_006',
#       '1dzlA_5keqF/output_cmap_esm/msa_t__ShallowMsa_007',
#       '1dzlA_5keqF/output_cmap_esm/msa_t__ShallowMsa_008',
#       '1dzlA_5keqF/output_cmap_esm/msa_t__ShallowMsa_009',
#       '1dzlA_5keqF/output_cmap_esm/msa_t__ShallowMsa_010',
#       '1dzlA_5keqF/output_cmap_esm/msa_t__ShallowMsa_011',
#       '1dzlA_5keqF/output_cmap_esm/msa_t__ShallowMsa_012',
#       '1dzlA_5keqF/output_cmap_esm/msa_t__ShallowMsa_013',
#       '1dzlA_5keqF/output_cmap_esm/msa_t__ShallowMsa_014',
#       '1dzlA_5keqF/output_cmap_esm/msa_t__ShallowMsa_015',
#       '1dzlA_5keqF/output_cmap_esm/msa_t__ShallowMsa_016',
#       '1dzlA_5keqF/output_cmap_esm/msa_t__ShallowMsa_017',
#       '1dzlA_5keqF/output_cmap_esm/msa_t__ShallowMsa_018',
#       '1dzlA_5keqF/output_cmap_esm/msa_t__ShallowMsa_019',
#       '1dzlA_5keqF/output_cmap_esm/msa_t__ShallowMsa_020',
#       '1dzlA_5keqF/output_cmap_esm/msa_t__ShallowMsa_021',
#       '1dzlA_5keqF/output_cmap_esm/msa_t__ShallowMsa_022',
#       '1dzlA_5keqF/output_cmap_esm/msa_t__ShallowMsa_023']

# print("Parsed names: ", [match.group() for s in node_names for match in re.finditer(r'M[sS][aA][a-zA-Z0-9_].*', s)])

