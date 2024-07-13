#
# '/Users/steveabecassis/Desktop/Pipeline_res/Pipeline/1eboE_5fhcJ/chain_pdb_files/1eboE.pdb')
# '/Users/steveabecassis/Desktop/Pipeline_res/Pipeline/1eboE_5fhcJ/chain_pdb_files/5fhcJ.pdb')


import pymol
from pymol import cmd

pymol.finish_launching(['pymol', '-c'])  # '-c' for command line (no GUI)

# Delete existing objects
cmd.delete('all')

pdb_file1 = '/Users/steveabecassis/Desktop/Pipeline_res/Pipeline/1eboE_5fhcJ/chain_pdb_files/1eboE.pdb'
pdb_file2 = '/Users/steveabecassis/Desktop/Pipeline_res/Pipeline/1eboE_5fhcJ/chain_pdb_files/5fhcJ.pdb'

# Load the PDB files
cmd.load(pdb_file1, 'protein1')
cmd.load(pdb_file2, 'protein2')

# Align the structures
cmd.align('protein2', 'protein1')

# Set different colors for visualization
cmd.color('red', 'protein1')
cmd.color('blue', 'protein2')

# Set the view (optional, you can customize this)
cmd.zoom('all', buffer=10)


cmd.show("cartoon", "protein1")
cmd.show("cartoon", "protein2")
cmd.hide("lines")
cmd.bg_color("white")
cmd.orient()


# Save the image
cmd.png('/Users/steveabecassis/Desktop/Pipeline_res/test.png', width=1920, height=1080, dpi=300, ray=1)

