from Bio import Phylo  # for phylogenetic trees
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import random
import matplotlib.pyplot as plt

from pylab import *

# Use ete3 package for visualization
from ete3 import *

from msa_utils import *


# Reconstruct a phylogenetic tree
def phytree_from_msa(msa_file, output_tree_file=[], max_seqs = 100):
    # Load the multiple sequence alignment from a file

    IDs, seqs = load_fasta(msa_file)
    seqs = [''.join([x for x in s if x.isupper() or x == '-']) for s in seqs]  # remove lowercase letters in alignment

    num_seqs = len(seqs)
    if num_seqs > max_seqs: # too many sequences! sample!!!
        seqs = random.sample(seqs, max_seqs) # [1:max_seqs]  # sample randomly

#    alignment = AlignIO.read(msa_file, "fasta")
#    print(seqs)
#    print([len(s) for s in seqs])

    seq_records = [SeqRecord(Seq(s), id=f"Sequence_{i}") for i, s in enumerate(seqs)]

    # Create a MultipleSeqAlignment object from the SeqRecord objects
    alignment = MultipleSeqAlignment(seq_records)


    # Calculate a distance matrix from the alignment
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)

    # Build a phylogenetic tree using the UPGMA (Unweighted Pair Group Method with Arithmetic Mean) method
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(distance_matrix)

    # Print or save the resulting tree
    if len(output_tree_file) == 0:
        output_tree_file = msa_file.replace(".a3m", "_tree.nwk")
#    tree_file = "phylogenetic_tree.nwk"  # Replace with your desired output file
    Phylo.write(tree, output_tree_file, "newick")

    return tree



# Define a dictionary to map numerical values to colors
value_to_color = {
    "Leaf_1": (0.2, 0.4, 0.6),  # Replace "Leaf_1" with the actual leaf name
    "Leaf_2": (0.8, 0.1, 0.3),  # Replace "Leaf_2" with the actual leaf name
    # Add more entries for other leaf nodes as needed
}


# Function to assign colors to tree nodes (should be depending on values!)
def set_node_color(node):
    if node.is_terminal():
        color = value_to_color.get(node.name, (0.5, 0.5, 0.5))  # Default to gray
    else:
        # You can set a color for internal nodes here if needed
        color = (0.5, 0.5, 0.5)  # Default to gray
    node.color = color


# Draw phylogenetic tree, with values assigned to each leaf
def draw_tree_with_values(tree, output_file = '', node_values = []):

#    tree = Phylo.read("apaf.xml", "phyloxml")
    # Try the epe package :
    ete_tree = Tree(tree, format=1)

    # Basic tree style
    ts = TreeStyle()
    ts.show_leaf_name = True

    if len(node_values) == 0:
        node_values = {n.name: 1 for n in ete_tree}

#    cmap = cm.get_cmap('seismic', 5)    # PiYG
    cmap = plt.colormaps['seismic']

#    cmap = plt.get_cmap('viridis')
    epsilon = 0.00000001
    norm_cmap = plt.Normalize(min(node_values.values())-epsilon, max(node_values.values())+epsilon)

#    color = cmap(norm(200.))

#    for i in range(cmap.N):
#        rgba = cmap(i)
#        print(matplotlib.colors.rgb2hex(rgba))
    # rgb2hex accepts rgb or rgba
#    print(matplotlib.colors.rgb2hex(rgba))

    for n in ete_tree.traverse():
        if n.is_leaf():
            nstyle = NodeStyle()
            print(n.name)
            print(node_values[n.name])
            print(norm_cmap(node_values[n.name]))
            print(cmap(norm_cmap(node_values[n.name])))
            nstyle["fgcolor"] = matplotlib.colors.rgb2hex(cmap(norm_cmap(node_values[n.name])))  # "red"  # color based on scale
            nstyle["size"] = 15
            n.set_style(nstyle)



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
