from Bio import Phylo  # for phylogenetic trees
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import random

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
    print(seqs)
    print([len(s) for s in seqs])

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


# Draw phylogenetic tree, with values assigned to each leaf
def draw_tree_with_values():

    return 0