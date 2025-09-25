# count_deep_vs_tree.py
from Bio import SeqIO
from ete3 import Tree
import os

pair = "2qqjA_4qdsA"
pair_dir = f"Pipeline/{pair}"

msa = os.path.join(pair_dir, "output_get_msa", "DeepMsa.a3m")
nwk = os.path.join(pair_dir, "output_phytree", "DeepMsa_tree.nwk")

n_msa = sum(1 for _ in SeqIO.parse(msa, "fasta"))
n_leaves = len(Tree(nwk, format=1).get_leaves())

print(f"DeepMsa.a3m sequences: {n_msa}")
print(f"DeepMsa_tree.nwk leaves: {n_leaves}")

from ete3 import Tree
import os
from utils.msa_utils import load_fasta

pair_id = "2qqjA_4qdsA"
pair_dir = f"Pipeline/{pair_id}"

# --- load Deep tree ---
tree_file = os.path.join(pair_dir, "output_phytree", "DeepMsa_tree.nwk")
ete_tree = Tree(tree_file, format=1)
leaf_names = [lf.name for lf in ete_tree.get_leaves()]
print(f"[INFO] Deep tree leaves: {len(leaf_names)}")

# --- load cluster membership from all .a3m files ---
msa_dir = os.path.join(pair_dir, "output_msa_cluster")
cluster_map = {}  # seq_id -> cluster_name
for fn in sorted(os.listdir(msa_dir)):
    if not fn.endswith(".a3m"):
        continue
    cluster = fn.replace(".a3m", "")
    ids, seqs = load_fasta(os.path.join(msa_dir, fn))
    for sid in ids:
        cluster_map[sid] = cluster
print(f"[INFO] Loaded {len(cluster_map)} sequence IDs from {len(os.listdir(msa_dir))} clusters")

# --- map tree leaves to clusters ---
leaf_to_cluster = {lf: cluster_map.get(lf, None) for lf in leaf_names}
clusters_in_tree = {c for c in leaf_to_cluster.values() if c is not None}
print(f"[RESULT] Clusters represented in tree: {len(clusters_in_tree)} / {len(set(cluster_map.values()))}")

# Debug print (first 10 leaves)
for lf, cl in list(leaf_to_cluster.items())[:10]:
    print(f"  {lf} -> {cl}")


exit(8888)


from utils.utils import *
from utils.protein_plot_utils import *
from utils.phytree_utils import *
import re
# import colabfold
# print("colab dir: " , colabfold.__file__)


# compute_tmscore_align('Pipeline/1dzlA_5keqF/1dzl.pdb', 'Pipeline/1dzlA_5keqF/AF_preds/ShallowMsa_000_unrelaxed_rank_005_alphafold2_ptm_model_1_seed_000.pdb', 'A', 'A')

# Load MSA:

# msa = download_and_parse_pfam_msa("PF00069", alignment_type="seed")
# print(msa)

# exit(99)

pair_id = "2qqjA_4qdsA" #  "1jfkA_2nxqB" ## "2qqjA_4qdsA" # "1jfkA_2nxqB"
pdbids = pair_str_to_tuple(pair_id)
pdbchains = [pdbids[0][-1], pdbids[1][-1]]
pdbids = [pdbids[0][:-1], pdbids[1][:-1]]

print("Tuple: ", pair_str_to_tuple(pair_id), " no chain : ", pdbids)
# Example for a single pair (tree only):
make_foldswitch_all_plots(
    pdbids=pdbids, #["1jfk","2nxq"],
    fasta_dir="Pipeline",
    foldpair_id=pair_id,
    pdbchains=pdbchains,
    plot_tree_clusters=True,   # <-- IMPORTANT (builds 7-col table)
    plot_contacts=False,       # focus on the tree for now
    global_plots=False)


exit(99)
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


col_groups = [
    ["TM-AF1", "TM-AF2"],
    ["TM-ESM1", "TM-ESM2"],
    ["RE-MSAT-COM", "RE-MSAT1", "RE-MSAT2"],
]
group_titles = ["AF", "ESMF", "MSAT"]


visualize_tree_with_heatmap(clusters_subtree, concat_scores, col_groups, group_titles, 'temp_local_tree.png')

# esmf_df = pd.read_csv(ESMF_MODEL_FILE)
# print("ESM DF: ", esmf_df)

# msa_trans_df = pd.read_csv(MSA_TRANS_MODEL_FILE)
# print("MSA-Trans DF: ", msa_trans_df)



