from ete3 import Tree, TreeStyle, faces, AttrFace, NodeStyle, TextFace, RectFace

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
# from matplotlib.colors import ListedColormap, to_hex
from matplotlib.colors import ListedColormap, Normalize, to_hex

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


flat_data = data_matrix.flatten()

# Create a Normalize object to scale values between 0 and 1
norm = Normalize(vmin=flat_data.min(), vmax=flat_data.max())

# Create a colormap for continuous data
cmap = plt.cm.viridis

# Flatten the data matrix and get unique values
# unique_values = np.unique(data_matrix.flatten())

# Create a colormap dynamically based on unique values
# cmap = ListedColormap(plt.cm.tab10.colors[:len(unique_values)])


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

        # Calculate the starting x-position for the rectangles based on node distance
        x_start = node.dist
        print(x_start)

        # Create RectFace instances with dynamically assigned colors
        for i, value in enumerate(values):
            hex_color = to_hex(cmap(norm(value)))
            rect_face = RectFace(width=20, height=20, fgcolor='black', bgcolor=hex_color)
            # Adjust the position of the RectFace within the node
#            rect_face.x = 5*x_start + i * 20
#            if i == 0:
 #               save_y = rect_face.y  # Adjust the y-position for vertical alignment
 #           else:
#            rect_face.yzz = -10
            faces.add_face_to_node(rect_face, node, column=i, position="aligned")




# Apply the layout function to the tree
tree.render("tree_with_heatmap.png", w=800, units="px", tree_style=ts , layout=layout)

#print("many cmaps:")
#print([cmap(i) for i in range(17)])

#for i, value in enumerate([1,2,3,4]):
#    print("outside-value:")
#    print(value)
#    print("outside-cmap-value:")
#    print(cmap(value))


# Display the figure
plt.show()

