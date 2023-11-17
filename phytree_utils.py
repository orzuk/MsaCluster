from Bio import Phylo  # for phylogenetic trees
from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import ParsimonyScorer, NNITreeSearcher, ParsimonyTreeConstructor


import matplotlib.pyplot as plt
import pickle
from pylab import *

from matplotlib.colors import Normalize, to_hex

# Use ete3 package for visualization
from ete3 import *
from msa_utils import *
import random


def resolve_duplicated_ids(ids_list):
        """
        Takes a list of strings and adds a counter to duplicates to make each string unique.
        Unique strings remain unchanged.

        Args:
        strings (list): A list of strings, potentially with duplicates.

        Returns:
        list: A list of strings where duplicates are made unique with a counter.
        """
        count = {}
        result = []

        for s in ids_list:
            if s in count:
                count[s] += 1
                result.append(f"{s}_{count[s]}")
            else:
                count[s] = 0
                result.append(s)

        return result


# Reconstruct a phylogenetic tree
# Input:
# msa_file - file with MSA in a3m format
# output_tree_file - where to save
# max_seqs - sample sequences if too many
def phytree_from_msa(msa_file, output_tree_file=[], max_seqs = 100):
    # Load the multiple sequence alignment from a file

    seqs_IDs, seqs = load_fasta(msa_file)
    seqs = [''.join([x for x in s if x.isupper() or x == '-']) for s in seqs]  # remove lowercase letters in alignment

    num_seqs = len(seqs)
#    print(seqs)
#    print(num_seqs)
#    print(max_seqs)
    if num_seqs > max_seqs:  # too many sequences! sample!!!
        rand_inds = random.sample(range(num_seqs), max_seqs)
        seqs = [seqs[i] for i in rand_inds]  # random.sample(seqs, max_seqs)  # [1:max_seqs]  # sample randomly
        seqs_IDs = [seqs_IDs[i] for i in rand_inds]
        # Need to match also seq IDs!!!

    seqs_IDs = resolve_duplicated_ids(seqs_IDs)  # resolve duplicate ids:


#    alignment = AlignIO.read(msa_file, "fasta")
#    print(seqs)
#    print([len(s) for s in seqs])

#    seq_records = [SeqRecord(Seq(s), id=f"Sequence_{i}") for i, s in enumerate(seqs)]  # Here must give correct names to sequences!
    seq_records = [SeqRecord(Seq(seqs[i]), id=seqs_IDs[i]) for i in range(len(seqs))]  # Here must give correct names to sequences!

    # Create a MultipleSeqAlignment object from the SeqRecord objects
    alignment = MultipleSeqAlignment(seq_records)

    # Calculate a distance matrix from the alignment
    calculator = DistanceCalculator('identity')
#    print("Alignment: ")
#    print(alignment)
#    print([type(s) for s in alignment])
#    print(type(alignment))
#    with open('bad_msa.pkl', 'wb') as f:  # Python 3: open(..., 'rb')
#        pickle.dump([alignment, calculator], f)
#    with open('bad_msa.pkl', 'rb') as f:  # Python 3: open(..., 'rb')
#        alignment, calculator = pickle.load(f)
#    all_ids = [s.id for s in alignment]

    distance_matrix = calculator.get_distance(alignment)

    # Build a phylogenetic tree using the UPGMA (Unweighted Pair Group Method with Arithmetic Mean) method
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(distance_matrix) # , names = seqs_IDs)  # New: Add leave names !!!

    # Print or save the resulting tree
    if len(output_tree_file) == 0:
        output_tree_file = msa_file.replace(".a3m", "_tree.nwk")
    if not os.path.exists(os.path.dirname(output_tree_file)):
        os.makedirs(os.path.dirname(output_tree_file))
    #    tree_file = "phylogenetic_tree.nwk"  # Replace with your desired output file
#    print("Output tree file: " + output_tree_file)

    # Add quotes to node names (needed?):
##    for clade in tree.find_clades():
#        clade.name = "'" + clade.name + "'"
##        clade.name = "\"" + clade.name + "\""

#    ctr = 0
#    for node in tree.find_clades():
#        print(node.name)
#        node.name = "leaf" + str(ctr)
#        ctr += 1
    # write_newick_with_quotes(tree, output_tree_file)
    Phylo.write(tree, output_tree_file, "newick")  # This is different from write_newick_with_quotes !!!!

    return tree


def convert_biopython_to_ete3(biopy_tree, parent_ete_node=None):
    # TRY using files: need a dictionary for the leaf names
    node_dict = {}
    ctr = 0
    for node in biopy_tree.find_clades():
        #        clade.name = "'" + clade.name + "'"
        node_dict["node" + str(ctr)] = node.name
        node.name = "node" + str(ctr)
        ctr += 1

    Phylo.write(biopy_tree, "temp_tree.nwk", "newick")  # write to file

    ete_tree = Tree("temp_tree.nwk", 1) # read
    for node in ete_tree.traverse():  # modify name
        node.name = node_dict[node.name]

    return ete_tree



# Read tree and strip quotes from node names
def read_tree_ete(phytree_file):
    ete_tree = Tree(phytree_file, format=1)
    for node in ete_tree.traverse():
        node.name = node.name.strip("'")
        node.name = node.name.strip("\"")
    return ete_tree


# Define a dictionary to map numerical values to colors
value_to_color = {
    "Leaf_1": (0.2, 0.4, 0.6),  # Replace "Leaf_1" with the actual leaf name
    "Leaf_2": (0.8, 0.1, 0.3),  # Replace "Leaf_2" with the actual leaf name
    # Add more entries for other leaf nodes as needed
}


def write_newick_with_quotes(tree, file_path):
    def format_clade(clade):
        if clade.is_terminal():
            return f"'{clade.name}':{clade.branch_length}" if clade.name else ''
        else:
            subtrees = ','.join(format_clade(sub_clade) for sub_clade in clade)
            return f"({subtrees}){'' if clade.name is None else f'{clade.name}'}:{clade.branch_length}"

    newick_str = format_clade(tree.root) + ';'
    with open(file_path, 'w') as f:
        f.write(newick_str)


# Function to assign colors to tree nodes (should be depending on values!)
def set_node_color(node):
    if node.is_terminal():
        color = value_to_color.get(node.name, (0.5, 0.5, 0.5))  # Default to gray
    else:
        # You can set a color for internal nodes here if needed
        color = (0.5, 0.5, 0.5)  # Default to gray
    node.color = color


# Induced subtree for an ete3 tree
def extract_induced_subtree(tree, leaf_names):
    """
    Extracts the induced subtree that includes only the specified leaves and their ancestral nodes.

    Args:
    tree (ete3.Tree): The original tree from which to extract the subtree.
    leaf_names (list): A list of leaf names to include in the subtree.

    Returns:
    ete3.Tree: The induced subtree containing only the specified leaves and their ancestors.
    """
    # Copy the tree to avoid modifying the original
    subtree = tree.copy(method='deepcopy')

    # Prune the subtree to keep only the leaves specified
    leaves_to_keep = set(leaf_names)
    for leaf in subtree.iter_leaves():
        if leaf.name not in leaves_to_keep:
            leaf.delete()

    return subtree


# Draw phylogenetic tree, with values assigned to each leaf
# Input:
# tree - a phylogenetic tree object
# output_file - where to save image
# node_values - vector/matrix of values representing each node
def visualize_tree_with_heatmap(phylo_tree, node_values_matrix, output_file=None):
    # Ensure the data matrix is a numpy array
    node_names = node_values_matrix.index.tolist()
    node_values_matrix = np.array(node_values_matrix)

    # Load the phylogenetic tree
    bio_tree = Phylo.read(phylo_tree, "newick")  # This is different from write_newick_with_quotes !!!!
    print("Convert to ete3 tree:")
    tree = convert_biopython_to_ete3(bio_tree)

    tree = extract_induced_subtree(tree, node_names)

    # get only subtree based on node_values_matrix

#    tree = read_tree_ete(phylo_tree)

    # Create a Normalize object to scale values between 0 and 1
    flat_data = node_values_matrix.flatten()
    epsilon = 0.00000001
    norm = Normalize(vmin=flat_data.min()-epsilon, vmax=flat_data.max()+epsilon)

    # Create a colormap for continuous data
    cmap = plt.cm.viridis

    # Create a TreeStyle for the phylogenetic tree
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = True
    ts.show_scale = False

    def layout(node):
        if node.is_leaf():
            # Access the corresponding rows in the data_matrix using leaf names
            leaf_name = node.name
            index = tree.get_leaf_names().index(leaf_name)
            values = node_values_matrix[index]

            # Define NodeStyle for the node
            node_style = NodeStyle()
            node_style["size"] = 0  # Set the size to 0 to hide the default node circle
            node.set_style(node_style)

            # Create RectFace instances with dynamically assigned colors
            for i, value in enumerate(values):
                hex_color = to_hex(cmap(norm(value)))
                rect_face = RectFace(width=20, height=20, fgcolor='black', bgcolor=hex_color)
                faces.add_face_to_node(rect_face, node, column=i, position="aligned")

    # Render the tree
    print("Saving tree image to : " + output_file)
    if output_file:
        if '.' not in output_file:
            output_file = output_file + ".png"
        tree.render(output_file, w=800, units="px", tree_style=ts, layout=layout)
    else:
        tree.show(tree_style=ts, layout=layout)


# Draw phylogenetic tree, with values assigned to each leaf
# Input:
# tree - a phylogenetic tree object
# output_file - where to save image
# node_values - vector/matrix of values representing each node
def draw_tree_with_values(tree, output_file= '', node_values= []):

#    tree = Phylo.read("apaf.xml", "phyloxml")
    # Try the epe package :
    ete_tree = read_tree_ete(tree)

    # Basic tree style
    ts = TreeStyle()
    ts.show_leaf_name = True

    if len(node_values) == 0:
        node_values = {n.name: 1 for n in ete_tree}

#    cmap = cm.get_cmap('seismic', 5)    # PiYG
    cmap = plt.colormaps['seismic']

#    cmap = plt.get_cmap('viridis')
    epsilon = 0.00000001

#    normalize_column = lambda col: (col - col.min()) / (col.max() - col.min())
    norm_cmap = plt.Normalize(min(node_values.shared)-epsilon, max(node_values.shared)+epsilon)

#    color = cmap(norm(200.))

#    for i in range(cmap.N):
#        rgba = cmap(i)
#        print(matplotlib.colors.rgb2hex(rgba))
    # rgb2hex accepts rgb or rgba
#    print(matplotlib.colors.rgb2hex(rgba))

    ete_tree.add_face(CircleFace(), column=0, position = "branch-right")

    for n in ete_tree.traverse():
        if n.is_leaf():
            nstyle = NodeStyle()
            print(n.name + " " + str(node_values.at[n.name, 'shared']) + " " + str(norm_cmap(node_values.at[n.name, 'shared'])))
            print(cmap(norm_cmap(node_values.at[n.name, 'shared'])))
            nstyle["fgcolor"] = matplotlib.colors.rgb2hex(cmap(norm_cmap(node_values.at[n.name, 'shared'])))  # "red"  # color based on scale
            nstyle["size"] = 15
            n.set_style(nstyle)
            # nface = NodeFace()

    # Let's now modify the aspect of the root node
    ete_tree.img_style["size"] = 30
    ete_tree.img_style["fgcolor"] = "blue"

    # Draws nodes as small red spheres of diameter equal to 10 pixels

#    if type(tree) == str:
#        tree = Phylo.read(tree, "newick")
#    tree.ladderize()  # Flip branches so deeper clades are displayed at top

    # Plot the tree without colors
#    Phylo.draw(tree, branch_labels=lambda c: c.branch_length, do_show=False)

    if len(output_file) > 0:  # save and close plot (enable automatic saving of multiple plots)
        if '.' not in output_file:
            output_file = output_file + ".png"
        print("Save tree fig: " + output_file)
        ete_tree.render(output_file)
#        plt.savefig(output_file + '.png')
    else:
        ete_tree.show()

#    Phylo.draw(tree)
    return 0


# Perform ancestral reconstruction of sequences in the phylogenetic tree
# Input:
# tree_file - file with the phylogenetic tree
# alignment_file - file with the multiple sequence alignment
# output_file - where to save the ancestral sequences
def reconstruct_ancestral_sequences(tree_file, alignment_file, output_file):
    # Read the tree
    tree = Phylo.read(tree_file, 'newick')

    # Read the multiple sequence alignment
    alignment = AlignIO.read(alignment_file, 'fasta')

    # Perform the parsimony reconstruction
    scorer = ParsimonyScorer()
    searcher = NNITreeSearcher(scorer)
    constructor = ParsimonyTreeConstructor(searcher, alignment)
    parsimony_tree = constructor.build_tree(tree)

    # Print to file the reconstructed sequences for the internal nodes
    with open(output_file, 'w') as f:
        for clade in parsimony_tree.get_nonterminals():
            f.write(f"{clade.name}: {clade.confidences[0]}\n")

# Example usage
# reconstruct_ancestral_sequences('path/to/tree_file.newick', 'path/to/alignment_file.fasta', 'output.txt')

