from Bio import Phylo  # for phylogenetic trees
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import random

# Use ete3 package for visualization
from ete3 import Tree

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
def draw_tree_with_values(tree, output_file = ''):
#    tree = Phylo.read("apaf.xml", "phyloxml")

    if type(tree) == str:
        tree = Phylo.read(tree, "newick")
    tree.ladderize()  # Flip branches so deeper clades are displayed at top

    # Apply color mapping to all nodes in the tree
#    tree.clade.traverse(set_node_color)

    # Plot the tree with colors
    Phylo.draw(tree, branch_labels=lambda c: c.branch_length)

    if len(output_file) > 0:  # save and close plot (enable automatic saving of multiple plots)
        print("Save tree fig: " + output_file + '.png')
        plt.savefig(output_file + '.png')

#    Phylo.draw(tree)
    return 0
