from Bio import Phylo  # for phylogenetic trees
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor


# Reconstruct a phylogenetic tree
def phytree_from_msa(msa_file, output_tree_file=[]):
    # Load the multiple sequence alignment from a file
    alignment = AlignIO.read(msa_file, "fasta")

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


# Draw phylogenetic tree, with values assigned to each leaf
def draw_tree_with_values():

    return 0