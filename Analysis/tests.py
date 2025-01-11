from typing import List
from utils.protein_utils import *

# Sandbox for trying and testing different things

pdb1 = '/Users/steveabecassis/Desktop/1eboE_5fhcJ/chain_pdb_files/1eboE.pdb'
pdb2 = '/Users/steveabecassis/Desktop/1eboE_5fhcJ/chain_pdb_files/5fhcJ.pdb'

structure_index,reference_index = Match_Alpha_Carbons(pdb1, pdb2)
Calculate_RMSD(pdb1, pdb2, [_[1] for _ in structure_index],[_[1] for _ in reference_index])[0]



import subprocess

# Define the terminal command you want to run
command = "./TMalign ./chain_pdb_files/1eboE.pdb ./chain_pdb_files/5fhcJ.pdb"
command = "./TMalign"
command = "pwd"
command = "gfortran  -O3 -ffast-math -lm -o TMalign /Users/steveabecassis/Desktop/TMalign.f"
command = "./TMalign   ./chain_pdb_files/1eboE.pdb ./chain_pdb_files/5fhcJ.pdb"
command = "./TMalign   ./chain_pdb_files/1eboE.pdb ./chain_pdb_files/5fhcJ.pdb"

# Run the command and capture its output
output = subprocess.check_output(command, shell=True)
# Print the output
print(output.decode())
import re
match = re.search(r"TM-score=\s+(\d+\.\d+)", str(output))

if match:
    result = match.group(1)
    print(result)
float(result)



# Debug cmap plot
with open('debug_match_index.pkl', 'rb') as f:  # Python 3: open(..., 'wb')
    cmap1, cmap2 = pickle.load(f)
