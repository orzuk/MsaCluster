import copy


# for phylogenetic trees
from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor, ParsimonyScorer, NNITreeSearcher, ParsimonyTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

import pickle
from pylab import *

# Use ete3 package for visualization
from ete3 import *
from .msa_utils import *
import random
import os
os.environ["QT_QPA_PLATFORM"] = "offscreen"
os.environ["XDG_RUNTIME_DIR"] = "/tmp"

from ete3 import TreeStyle, TextFace, RectFace, NodeStyle
import numpy as np
from matplotlib.colors import Normalize, to_hex
from matplotlib import pyplot as plt
from matplotlib.colorbar import ColorbarBase
from matplotlib.transforms import Bbox



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

    seq_records = [SeqRecord(Seq(seqs[i]), id=seqs_IDs[i]) for i in range(len(seqs))]  # Here must give correct names to sequences!

    # Create a MultipleSeqAlignment object from the SeqRecord objects
    alignment = MultipleSeqAlignment(seq_records)

    # Calculate a distance matrix from the alignment
    calculator = DistanceCalculator('identity')

    distance_matrix = calculator.get_distance(alignment)

    # Build a phylogenetic tree using the UPGMA (Unweighted Pair Group Method with Arithmetic Mean) method
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(distance_matrix)  # , names = seqs_IDs)  # New: Add leave names !!!

    # Print or save the resulting tree
    if len(output_tree_file) == 0:
        output_tree_file = msa_file.replace(".a3m", "_tree.nwk")
    if not os.path.exists(os.path.dirname(output_tree_file)):
        os.makedirs(os.path.dirname(output_tree_file))
    #    tree_file = "phylogenetic_tree.nwk"  # Replace with your desired output file

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


# Induced subtree for an ete3 tree
# keep only the input leaves appearing in 'leaf_names'
def extract_induced_subtree(tree, leaf_names):
    """
    Extracts the induced subtree that includes only the specified leaves and their ancestral nodes.

    Args:
    tree (ete3.Tree): The original tree from which to extract the subtree.
    leaf_names (list): A list of leaf names to include in the subtree.

    Returns:
    ete3.Tree: The induced subtree containing only the specified leaves and their ancestors.
    """
    if type(tree) == str: # extract tree from file
        print("Read tree from file: " + tree)
        bio_tree = Phylo.read(tree, "newick")  # This is different from write_newick_with_quotes !!!!
        tree = convert_biopython_to_ete3(bio_tree)

    # Copy the tree to avoid modifying the original
    subtree = tree.copy(method='deepcopy')

    # Prune the subtree to keep only the leaves specified
    leaves_to_keep = set(leaf_names)
    for leaf in subtree.iter_leaves():
        if leaf.name not in leaves_to_keep:
            leaf.delete()

    return subtree


def visualize_tree_with_heatmap_new(phylo_tree, node_values_matrix, output_file=None):
    """
    Visualizes a phylogenetic tree with a heatmap next to each leaf and two vertical colorbars.

    Parameters:
    - phylo_tree: An ETE3 Tree object representing the full phylogenetic tree.
    - node_values_matrix: A pandas DataFrame indexed by leaf names, with 7 columns.
                         The first 4 columns are assumed to represent TM-Score values (Group 1),
                         and the last 3 columns represent Recall values (Group 2).
    - output_file: Optional string path to save the output image file (e.g., "tree_heatmap.png").
                   If not provided, the figure will not be saved to disk.
    """
    if output_file is not None:
        matplotlib.use('Agg')  # Use non-GUI backend when saving

    node_names = node_values_matrix.index.tolist()
    node_values_matrix = np.array(node_values_matrix)
    tree = extract_induced_subtree(phylo_tree, node_names)

    group1 = node_values_matrix[:, :4].flatten()
    group2 = node_values_matrix[:, 4:].flatten()
    norm1 = Normalize(vmin=group1.min(), vmax=group1.max())
    norm2 = Normalize(vmin=group2.min(), vmax=group2.max())
    cmap1 = plt.cm.viridis
    cmap2 = plt.cm.plasma

    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = True
    ts.show_scale = False

    def layout(node):
        if node.is_leaf():
            index = tree.get_leaf_names().index(node.name)
            node.name = str(index)
            values = node_values_matrix[index]

            node_style = NodeStyle()
            node_style["size"] = 0
            node.set_style(node_style)

            group_breaks = [2, 4]
            initial_spacer = RectFace(width=10, height=20, fgcolor="white", bgcolor="white")
            faces.add_face_to_node(initial_spacer, node, column=0, position="aligned")
            column_offset = 1

            for i, value in enumerate(values):
                hex_color = to_hex(cmap1(norm1(value))) if i < 4 else to_hex(cmap2(norm2(value)))
                rect_face = RectFace(width=20, height=20, fgcolor="black", bgcolor=hex_color)
                faces.add_face_to_node(rect_face, node, column=i + column_offset, position="aligned")

                if i + 1 in group_breaks:
                    spacer = RectFace(width=5, height=20, fgcolor="white", bgcolor="white")
                    column_offset += 1
                    faces.add_face_to_node(spacer, node, column=i + column_offset, position="aligned")

    ts.layout_fn = layout

    temp_tree_file = "temp_tree.png"
    tree.render(temp_tree_file, w=800, units="px", tree_style=ts, layout=layout)

    tree_img = plt.imread(temp_tree_file)
    fig, ax_tree = plt.subplots(figsize=(12, 8))
    fig.subplots_adjust(left=0.05, right=0.92, top=0.95, bottom=0.1)
    ax_tree.imshow(tree_img)
    ax_tree.axis("off")

    # --- Dynamic colorbar placement based on tree image dimensions ---
    img_height, img_width, _ = tree_img.shape
    bar_width = 0.015
    bar_height = 0.3
    bar_x = 0.95 - bar_width

    # Optionally use number of leaves to adjust vertical placement
    n_leaves = len(tree.get_leaf_names())
    top_y = 0.9
    bottom_y = 0.1
    available_height = top_y - bottom_y
    dynamic_bar_height = min(bar_height, available_height * 0.4)

    cbar_ax1 = fig.add_axes([bar_x, top_y - dynamic_bar_height, bar_width, dynamic_bar_height])
    cb1 = ColorbarBase(cbar_ax1, cmap=cmap1, norm=norm1, orientation="vertical")
    cbar_ax1.set_ylabel("TM-Score", fontsize=10, labelpad=10, rotation=90)

    cbar_ax2 = fig.add_axes([bar_x, bottom_y, bar_width, dynamic_bar_height])
    cb2 = ColorbarBase(cbar_ax2, cmap=cmap2, norm=norm2, orientation="vertical")
    cbar_ax2.set_ylabel("Recall", fontsize=10, labelpad=10, rotation=90)

    # Add label under heatmap
    fig.text(0.68, 0.08, "AF   ESMF  MSAT", ha='center', fontsize=10, fontweight='bold')

    if output_file is not None:
        if '.' not in output_file:
            output_file = output_file + ".png"
        fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)
    os.remove(temp_tree_file)

# Updated visualize_tree_with_heatmap function with vertical colorbars
def visualize_tree_with_heatmap(phylo_tree, node_values_matrix, output_file=None):
    if output_file is not None:
        matplotlib.use('Agg')  # Non-GUI backend (e.g., for saving figures)

    # Ensure the data matrix is a numpy array
    node_names = node_values_matrix.index.tolist()
    node_values_matrix = np.array(node_values_matrix)

    # Extract the induced subtree
    tree = extract_induced_subtree(phylo_tree, node_names)

    # Normalize the node values separately for the two groups
    group1 = node_values_matrix[:, :4].flatten()  # First 4 columns
    group2 = node_values_matrix[:, 4:].flatten()  # Last 3 columns

    norm1 = Normalize(vmin=group1.min(), vmax=group1.max())
    norm2 = Normalize(vmin=group2.min(), vmax=group2.max())

    cmap1 = plt.cm.viridis  # Colormap for first 4 columns
    cmap2 = plt.cm.plasma   # Colormap for last 3 columns

    # Create TreeStyle
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = True
    ts.show_scale = False

    # Define layout function for rendering heatmap
    def layout(node):
        if node.is_leaf():
            # Simplify the leaf names
            index = tree.get_leaf_names().index(node.name)
            node.name = str(index)  # Use numeric names for simplicity

            # Access the row in the data matrix
            values = node_values_matrix[index]

            # Set node style
            node_style = NodeStyle()
            node_style["size"] = 0  # Hide the default node circle
            node.set_style(node_style)

            # Add heatmap rectangles and spacers
            group_breaks = [2, 4]  # Define spacers after 2nd and 4th columns

            initial_spacer = RectFace(width=10, height=20, fgcolor="white", bgcolor="white")  # Increase width for more space
            faces.add_face_to_node(initial_spacer, node, column=0, position="aligned")
            column_offset = 1  # Shift all heatmap columns to the right


#            spacer = RectFace(width=50, height=20, fgcolor="white", bgcolor="white")  # Adjust width to control space


            for i, value in enumerate(values):
                # Determine the normalization and colormap
                if i < 4:  # First 4 columns
                    hex_color = to_hex(cmap1(norm1(value)))
                else:  # Last 3 columns
                    hex_color = to_hex(cmap2(norm2(value)))

                # Add the heatmap block
                rect_face = RectFace(width=20, height=20, fgcolor="black", bgcolor=hex_color)
                faces.add_face_to_node(rect_face, node, column=i + column_offset, position="aligned")

                # Add a spacer after the group ends
                if i + 1 in group_breaks:
                    spacer = RectFace(width=5, height=20, fgcolor="white", bgcolor="white")
                    column_offset += 1  # Increment column index for the spacer
                    faces.add_face_to_node(spacer, node, column=i + column_offset, position="aligned")

    # Assign the layout function to TreeStyle
    ts.layout_fn = layout

    # Render the tree temporarily to get its dimensions
    temp_tree_file = "temp_tree.png"
    tree.render(temp_tree_file, w=800, units="px", tree_style=ts, layout=layout)

    # Load the tree image to combine with colorbars
    tree_img = plt.imread(temp_tree_file)
    fig, ax_tree = plt.subplots(figsize=(12, 8))
    fig.subplots_adjust(left=0.05, right=0.92, top=0.95, bottom=0.1)  # Move right limit closer

    ax_tree.imshow(tree_img)
    ax_tree.axis("off")  # Hide axes

    # Colorbar for the first 4 columns (TM-Score)
#    cbar_ax1 = fig.add_axes([0.88, 0.55, 0.02, 0.35])  # Narrower and closer
    cbar_ax1 = fig.add_axes([0.75, 0.6, 0.015, 0.3])  # Move closer & narrower
    cb1 = ColorbarBase(cbar_ax1, cmap=cmap1, norm=norm1, orientation="vertical")
#    cbar_ax1.set_title("TM-Score", fontsize=10, pad=10, loc='left', rotation=90)
    cbar_ax1.set_ylabel("TM-Score", fontsize=10, labelpad=10, rotation=90)

    # Colorbar for the last 3 columns (Recall)
#    cbar_ax2 = fig.add_axes([0.88, 0.1, 0.02, 0.35])  # Narrower and closer
    cbar_ax2 = fig.add_axes([0.75, 0.2, 0.015, 0.3])  # Same for second colorbar
    cb2 = ColorbarBase(cbar_ax2, cmap=cmap2, norm=norm2, orientation="vertical")
#    cbar_ax2.set_title("Recall", fontsize=10, pad=10, loc='left', rotation=90)
    cbar_ax2.set_ylabel("Recall", fontsize=10, labelpad=10, rotation=90)

    # Add the major title dynamically (AF, ESMF, MSAT)
    tree_width = calculate_tree_width(temp_tree_file)  # Implement a helper to calculate width from the image
#    print("x offset=", int(tree_width * 0.5))
#    add_figure_title(temp_tree_file,
#        title="AF   ESMF   MSAT",  # The title for the heatmap
#        title_x_offset=int(tree_width * 0.5), # 0.775),  # Adjust this value to center the title
#        title_y_offset=-50,  # Negative to place below the heatmap
#        font_size=18,
#        extra_space_ratio=0.03)  # Adds space at the bottom for the title
    fig.text(0.68, 0.08, "AF   ESMF  MSAT", ha='center', fontsize=10, fontweight='bold')

    # Save the combined figure
    if output_file is not None:
        if '.' not in output_file:
            output_file = output_file + ".png"
        fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Clean up the temporary file
    os.remove(temp_tree_file)


def calculate_tree_width(rendered_tree_path):
    """
    Calculate the width of the rendered tree from the temporary file.
    """
    from PIL import Image
    with Image.open(rendered_tree_path) as img:
        return img.size[0]  # Returns the width of the image


def add_figure_title(image_path, title, title_x_offset=None, title_y_offset=20,
                     font_size=30, extra_space_ratio=0.05):
    from PIL import Image, ImageDraw, ImageFont

    # Open the image
    img = Image.open(image_path)
    img_width, img_height = img.size

    # Calculate extra space to add at the bottom
    extra_space = int(img_height * extra_space_ratio)

    # Create a new image with extra space
    new_img = Image.new("RGB", (img_width, img_height + extra_space), (255, 255, 255))
    new_img.paste(img, (0, 0))

    # Use a default font that supports resizing
    try:
        font = ImageFont.truetype("DejaVuSans-Bold.ttf", font_size)
    except IOError:
        print("DejaVu Sans font not found. Using PIL's default font (no size customization).")
        font = ImageFont.load_default()  # Fallback font

    # Calculate text dimensions using getbbox()
    draw = ImageDraw.Draw(new_img)
    text_bbox = draw.textbbox((0, 0), title, font=font)  # Returns (x0, y0, x1, y1)
    text_width = text_bbox[2] - text_bbox[0]
    text_height = text_bbox[3] - text_bbox[1]

    # Center the title horizontally if no offset provided
    if title_x_offset is None:
        title_x_offset = (img_width - text_width) // 2

    # Position the title in the added space
    text_x = title_x_offset
    text_y = img_height + (extra_space - text_height) // 2

    # Draw the title
    draw.text((text_x, text_y), title, fill="black", font=font)

    # Save the updated image
    new_img.save(image_path)


def add_figure_title_working(image_path, title,
                             title_x_offset=None, title_y_offset=20,
                             font_size=30):
    from PIL import Image, ImageDraw, ImageFont

    # Open the image
    img = Image.open(image_path)
    img_width, img_height = img.size

    # Use a default font that supports resizing
    try:
        font = ImageFont.truetype("DejaVuSans-Bold.ttf", font_size)
    except IOError:
        print("DejaVu Sans font not found. Using PIL's default font (no size customization).")
        font = ImageFont.load_default()  # Fallback font

    # Calculate text dimensions using getbbox()
    draw = ImageDraw.Draw(img)
    text_bbox = draw.textbbox((0, 0), title, font=font)  # Returns (x0, y0, x1, y1)
    text_width = text_bbox[2] - text_bbox[0]
    text_height = text_bbox[3] - text_bbox[1]

    # Center the title horizontally if no offset provided
    if title_x_offset is None:
        title_x_offset = (img_width - text_width) // 2

    text_x = title_x_offset
    text_y = img_height - title_y_offset - text_height  # Position near the bottom

    # Draw the title
    draw.text((text_x, text_y), title, fill="black", font=font)

    # Save the updated image
    img.save(image_path)


# Draw phylogenetic tree, with values assigned to each leaf
# Input:
# phylo_tree - a phylogenetic tree object
# node_values - vector/matrix of values representing each node
# output_file - where to save image

def visualize_tree_with_heatmap_old(phylo_tree, node_values_matrix, output_file=None):
    from copy import deepcopy

    # Ensure the data matrix is a numpy array
    node_names = node_values_matrix.index.tolist()
    node_values_matrix = np.array(node_values_matrix)

    if type(phylo_tree) == str:  # Load the phylogenetic tree
        bio_tree = Phylo.read(phylo_tree, "newick")  # This is different from write_newick_with_quotes !!!!
        print("Read from file and convert to ete3 tree:")
        tree = convert_biopython_to_ete3(bio_tree)
    else:
        print("input phylotree: ")
        print(phylo_tree)
        print("node names:", node_names)
        tree = deepcopy(phylo_tree)

    tree = extract_induced_subtree(tree, node_names)
    print("Induced sub-tree:")
    print(tree)

    # get only subtree based on node_values_matrix

    # Create a Normalized object to scale values between 0 and 1
    flat_data = node_values_matrix.flatten()
    epsilon = 0.00000001
#    print("flat_data=", flat_data)
#    print("min=", flat_data.min(), "max=", flat_data.max())
    norm = Normalize(vmin=flat_data.min()-epsilon, vmax=flat_data.max()+epsilon)

    # Create a colormap for continuous data
    cmap = plt.cm.viridis

    def layout(node):
        if node.is_leaf():
            index = tree.get_leaf_names().index(node.name)
            node.name = str(index)
            values = node_values_matrix[index]

            node_style = NodeStyle()
            node_style["size"] = 0
            node.set_style(node_style)

            group_breaks = [2, 4]
            column = 0
            title_columns = []  # Keep track of columns where titles should go

            for i, value in enumerate(values):
                hex_color = to_hex(cmap(norm(value)))
                rect_face = RectFace(width=20, height=20, fgcolor='black', bgcolor=hex_color)
                faces.add_face_to_node(rect_face, node, column=column, position="aligned")

                if i in [0, 2, 4]:  # First column of each group
                    title_columns.append(column)

                column += 1

                if i + 1 in group_breaks:
                    spacer = RectFace(width=5, height=20, fgcolor='white', bgcolor='white')
                    faces.add_face_to_node(spacer, node, column=column, position="aligned")
                    column += 1

            # Add titles to the last leaf node
            if node == tree.get_leaves()[-1]:
                titles = ["AF", "ESMF", "MSAT"]
                for col, title in zip(title_columns, titles):
                    title_face = TextFace(title, fsize=12, fgcolor="black", bold=True)
                    faces.add_face_to_node(title_face, node, column=col, position="aligned")


    # Create a TreeStyle for the phylogenetic tree
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = True
    ts.show_scale = False
    ts.branch_vertical_margin = 15  # Add some vertical margin for titles

    if output_file:
        print("Saving tree image to : " + output_file)
        if '.' not in output_file:
            output_file = output_file + ".png"
        tree.render(output_file, w=800, units="px", tree_style=ts, layout=layout)
        print("Did tree render!!! ")
    else:
        tree.show(tree_style=ts, layout=layout)


# Perform ancestral reconstruction of sequences in the phylogenetic tree
# Input:
# tree_file - file with the phylogenetic tree
# alignment_file - file with the multiple sequence alignment
# output_file - where to save the ancestral sequences
def reconstruct_ancestral_sequences(tree_file, alignment_file, output_file):
    # Read the tree
    tree = Phylo.read(tree_file, 'newick')

    # Read the multiple sequence alignment
    seqs_IDs, seqs = load_fasta(alignment_file)
    seqs = [''.join([x for x in s if x.isupper() or x == '-']) for s in seqs]  # remove lowercase letters in alignment

    seq_records = [SeqRecord(Seq(seqs[i]), id=seqs_IDs[i]) for i in range(len(seqs))]  # Here must give correct names to sequences!

    # Create a MultipleSeqAlignment object from the SeqRecord objects
    alignment = MultipleSeqAlignment(seq_records)
#    alignment = AlignIO.read(alignment_file, 'fasta') # a3m file not all are of the same length !!


    # Perform the parsimony reconstruction
    scorer = ParsimonyScorer()
    searcher = NNITreeSearcher(scorer)
    constructor = ParsimonyTreeConstructor(searcher, alignment)
    print("Alignment:")
    print(alignment)
    print("Tree:")
    print(tree)
    parsimony_tree = constructor.build_tree(tree)

    # Print to file the reconstructed sequences for the internal nodes
    with open(output_file, 'w') as f:
        for clade in parsimony_tree.get_nonterminals():
            f.write(f"{clade.name}: {clade.confidences[0]}\n")

# Example usage
# reconstruct_ancestral_sequences('path/to/tree_file.newick', 'path/to/alignment_file.fasta', 'output.txt')

