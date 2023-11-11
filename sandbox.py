from ete3 import Tree, TreeStyle, faces, AttrFace, NodeStyle, TextFace, RectFace

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, to_hex

# Create a simple example tree
newick_tree = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5)E;"
tree = Tree(newick_tree, format= 1)

# Create a matrix for the heatmap (replace this with your actual data)
data_matrix = np.array([
    [1.0, 2.0, 3.0, 4.0],
    [5.0, 6.0, 7.0, 8.0],
    [9.0, 10.0, 11.0, 12.0],
    [13.0, 14.0, 15.0, 16.0]
])

# Flatten the data matrix and get unique values
unique_values = np.unique(data_matrix.flatten())

# Create a colormap dynamically based on unique values
cmap = ListedColormap(plt.cm.tab10.colors[:len(unique_values)])


# cmap = ListedColormap(['blue', 'green', 'yellow', 'red'])

# Create a TreeStyle for the phylogenetic tree
ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = True
ts.show_scale = False



def layout(node):
    if node.is_leaf():
        # Access the corresponding rows in the data_matrix using leaf names
        leaf_name = node.get_leaf_names()[0]
        index = tree.get_leaf_names().index(leaf_name)
        values = data_matrix[index]

        # Define NodeStyle for the node
        node_style = NodeStyle()
        node_style["size"] = 0  # Set the size to 0 to hide the default node circle
        node.set_style(node_style)

        # Create RectFace instances with dynamically assigned colors
        print("values: ")
        print(values)
        for i, value in enumerate(values):
            print("value:")
            print(value)
            print("cmap-value:")
            print(cmap(value))
            hex_color = to_hex(cmap(value))
            print("hex_color:")
            print(hex_color)
            rect_face = RectFace(width=20, height=20, fgcolor='black', bgcolor=hex_color)
            faces.add_face_to_node(rect_face, node, column=i)


# Apply the layout function to the tree
tree.render("tree_with_heatmap.png", w=800, units="px", tree_style=ts , layout=layout)

# Display the figure
plt.show()

