import pymol
from pymol import cmd
import sys
from pathlib import Path

# --- repo imports (top-level only) ---
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from config import *            # expects MAIN_DIR, etc.

# pymol.cmd.set("quiet", 1)
pymol.finish_launching(['pymol', '-cq'])  # '-c' for command line (no GUI)

# Delete existing objects
cmd.delete('all')

pdb_file1 = TABLES_RES + '/Pipeline/1eboE_5fhcJ/chain_pdb_files/1eboE.pdb'
pdb_file2 = TABLES_RES + '/Pipeline/1eboE_5fhcJ/chain_pdb_files/5fhcJ.pdb'

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
cmd.png(TABLES_RES + '/test.png', width=1920, height=1080, dpi=300, ray=1)

