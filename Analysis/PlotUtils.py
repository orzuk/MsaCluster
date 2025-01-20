import prody as pr
import py3Dmol
import tempfile
from contact_map import OverrideTopologyContactDifference
import matplotlib.pyplot as plt
import mdtraj as md
from contact_map import ContactFrequency, ContactDifference
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from matplotlib.pyplot import *
import prody
import os
from scipy.ndimage import binary_dilation
import matplotlib.patches as mpatches
from IPython.display import display, HTML
from Bio import pairwise2, PDB, Align
from Bio.PDB import PDBParser, PDBIO, Superimposer
from Bio.PDB.Polypeptide import three_to_one
import ipywidgets as widgets
from io import StringIO
from utils.msa_utils import *


'''
In the PlotTool class define you own path
'''


def get_seq_from_structure(structure):
    return ''.join([PDB.Polypeptide.three_to_one(res.resname)
                    for res in structure.get_residues()
                    if PDB.Polypeptide.is_aa(res)])


def get_aligned_atoms(structure1, structure2):
    seq1 = get_seq_from_structure(structure1)
    seq2 = get_seq_from_structure(structure2)

    alignment = pairwise2.align.globalxx(seq1, seq2)[0]

    aligned_atoms1 = []
    aligned_atoms2 = []

    seq1_index = 0
    seq2_index = 0

    for a1, a2 in zip(alignment[0], alignment[1]):
        if a1 != '-' and a2 != '-':
            res1 = list(structure1.get_residues())[seq1_index]
            res2 = list(structure2.get_residues())[seq2_index]
            if 'CA' in res1 and 'CA' in res2:
                aligned_atoms1.append(res1['CA'])
                aligned_atoms2.append(res2['CA'])

        if a1 != '-':
            seq1_index += 1
        if a2 != '-':
            seq2_index += 1

    return aligned_atoms1, aligned_atoms2


def toggle_visibility(change):
    if change['name'] == 'value':
        if change['owner'] == toggle1:
            view.setStyle({'model': 0},
                          {'cartoon': {'color': 'red'}} if change['new'] else {'cartoon': {'hidden': True}})
        elif change['owner'] == toggle2:
            view.setStyle({'model': 1},
                          {'cartoon': {'color': 'blue'}} if change['new'] else {'cartoon': {'hidden': True}})
    view.render()


def visualize_alignement_structure(PATH_PDB1, PATH_PDB2):
    # Set up PDB parser
    parser = PDB.PDBParser()

    # Load structures
    structure1 = parser.get_structure("protein1", PATH_PDB1)
    structure2 = parser.get_structure("protein2", PATH_PDB2)

    # Get aligned atoms
    aligned_atoms1, aligned_atoms2 = get_aligned_atoms(structure1, structure2)

    # Align structures
    super_imposer = Superimposer()
    super_imposer.set_atoms(aligned_atoms1, aligned_atoms2)
    super_imposer.apply(structure2.get_atoms())

    # Save aligned structure
    io = PDB.PDBIO()
    io.set_structure(structure2)
    io.save("aligned_protein2.pdb")

    # Create py3Dmol view
    view = py3Dmol.view(width=800, height=600)

    # Add structures to the view
    with open(pdb_file2_, 'r') as f:
        view.addModel(f.read(), 'pdb')
    with open("aligned_protein2.pdb", 'r') as f:
        view.addModel(f.read(), 'pdb')

    # Set styles
    view.setStyle({'model': 0}, {'cartoon': {'color': 'red'}})
    view.setStyle({'model': 1}, {'cartoon': {'color': 'blue'}})

    # Zoom to fit
    view.zoomTo()

    # Create toggle buttons
    toggle1 = widgets.ToggleButton(description='Toggle Protein 1', value=True)
    toggle2 = widgets.ToggleButton(description='Toggle Protein 2', value=True)

    toggle1.observe(toggle_visibility)
    toggle2.observe(toggle_visibility)

    # Display the view and toggle buttons
    display(HTML(view._make_html()))
    display(widgets.HBox([toggle1, toggle2]))


def dilate_with_tolerance(array, tolerance):
    structure = np.ones((2 * tolerance + 1, 2 * tolerance + 1), dtype=int)
    return binary_dilation(array, structure=structure).astype(int)


def load_pred_cmap(fileName):
    cmap = np.load(f'{plot_tool.folder}/{plot_tool.fold_pair}/output_cmap_esm/{fileName}.npy')
    cmap[cmap > 0.4] = 1  # 0.3 or 0.4? should be input to function?
    cmap[cmap <= 0.4] = 0
    return cmap


def visualize_structure_alignment(pdb_file1, pdb_file2, chain1='A', chain2='A'):
    # Parse PDB files
    parser = PDBParser(QUIET=True)
    structure1 = parser.get_structure("protein1", pdb_file1)
    structure2 = parser.get_structure("protein2", pdb_file2)

    # Extract sequences
    seq1 = "".join([three_to_one(res.resname) for res in structure1[0][chain1] if res.id[0] == ' '])
    seq2 = "".join([three_to_one(res.resname) for res in structure2[0][chain2] if res.id[0] == ' '])

    # Perform sequence alignment
    aligner = Align.PairwiseAligner()
    alignments = aligner.align(seq1, seq2)
    best_alignment = next(alignments)

    # Get aligned positions
    aligned_positions = [(i, j) for i, j in zip(best_alignment.path[0], best_alignment.path[1])
                         if i != -1 and j != -1]

    # Extract aligned CA atoms
    ca_atoms1 = [atom for atom in structure1[0][chain1].get_atoms() if atom.name == 'CA']
    ca_atoms2 = [atom for atom in structure2[0][chain2].get_atoms() if atom.name == 'CA']

    aligned_atoms1 = [ca_atoms1[i] for i, _ in aligned_positions]
    aligned_atoms2 = [ca_atoms2[j] for _, j in aligned_positions]

    # Superimpose structures
    super_imposer = PDB.Superimposer()
    super_imposer.set_atoms(aligned_atoms1, aligned_atoms2)
    super_imposer.apply(structure2.get_atoms())

    # Function to get PDB string
    def get_pdb_string(structure):
        io = StringIO()
        pdb_io = PDBIO()
        pdb_io.set_structure(structure)
        pdb_io.save(io)
        return io.getvalue()

    # Create a py3Dmol view
    view = py3Dmol.view(width=800, height=600)

    # Add the first structure
    view.addModel(get_pdb_string(structure1), 'pdb')
    view.setStyle({'model': -1}, {'cartoon': {'color': 'blue'}})

    # Add the second (aligned) structure
    view.addModel(get_pdb_string(structure2), 'pdb')
    view.setStyle({'model': -1}, {'cartoon': {'color': 'red'}})

    # Zoom to fit the structures
    view.zoomTo()

    # Calculate RMSD
    rmsd = super_imposer.rms

    # Add RMSD and alignment info as labels
    view.addLabel(f"RMSD: {rmsd:.2f} Å",
                  {'position': {'x': -20, 'y': -20}, 'backgroundColor': 'white', 'fontColor': 'black'})
    view.addLabel(f"Aligned: {len(aligned_positions)} of {len(seq1)} (blue) and {len(seq2)} (red) residues",
                  {'position': {'x': -20, 'y': -40}, 'backgroundColor': 'white', 'fontColor': 'black'})

    return view


class PlotTool:
    def __init__(self, folder='/Users/steveabecassis/Desktop/Pipeline', fold_pair='3hdfA_3hdeA'):
        self.folder    = folder
        self.fold_pair = fold_pair
        self.fold1 = fold_pair.split('_')[0]
        self.fold2 = fold_pair.split('_')[1]
        self.structure1 = pr.parsePDB(f'{self.folder}/{self.fold_pair}/chain_pdb_files/{self.fold1}.pdb')
        self.structure2 = pr.parsePDB(f'{self.folder}/{self.fold_pair}/chain_pdb_files/{self.fold2}.pdb')

    def read_pdb_file(self,file_path):
        with open(file_path, 'r') as file:
            return file.read()

    def plot_fold_alignement(self):
        # Select CA atoms for alignment
        ca_atoms1 = self.structure1.select('name CA')
        ca_atoms2 = self.structure2.select('name CA')

        common_residues = set(ca_atoms1.getResnums()).intersection(ca_atoms2.getResnums())
        matched_ca1 = ca_atoms1.select('resnum ' + ' '.join(map(str, common_residues)))
        matched_ca2 = ca_atoms2.select('resnum ' + ' '.join(map(str, common_residues)))

        transformation = pr.calcTransformation(matched_ca2, matched_ca1)
        structure2_aligned = transformation.apply(self.structure2)

        # Write the aligned structure2 to a temporary file
        with tempfile.NamedTemporaryFile(delete=False, suffix='.pdb') as temp_file:
            pr.writePDB(temp_file.name, structure2_aligned)
            aligned_pdb2_path = temp_file.name
        # Read the aligned PDB files
        pdb1 = self.read_pdb_file(f'{self.folder}/{self.fold_pair}/chain_pdb_files/{self.fold1}.pdb')
        pdb2 = self.read_pdb_file(aligned_pdb2_path)
        # Create a py3Dmol viewer
        viewer = py3Dmol.view(width=800, height=600)
        # Add the first structure with cartoon style
        viewer.addModel(pdb1, 'pdb')
        viewer.setStyle({'model': 0}, {'cartoon': {'color': 'blue'}})
        # Add the aligned second structure with cartoon style
        viewer.addModel(pdb2, 'pdb')
        viewer.setStyle({'model': 1}, {'cartoon': {'color': 'red'}})

        # Align the second model to the first
        viewer.zoomTo()
        viewer.zoom(1.2)
        viewer.setBackgroundColor('white')
        viewer.show()
        os.remove(aligned_pdb2_path)


    def plot_fold_alignement_(self,path_pdb1,path_pdb2):
        # Select CA atoms for alignment
        structure1 = pr.parsePDB(path_pdb1)
        structure2 = pr.parsePDB(path_pdb2)
        ca_atoms1 = structure1.select('name CA')
        ca_atoms2 = structure2.select('name CA')

        common_residues = set(ca_atoms1.getResnums()).intersection(ca_atoms2.getResnums())
        matched_ca1 = ca_atoms1.select('resnum ' + ' '.join(map(str, common_residues)))
        matched_ca2 = ca_atoms2.select('resnum ' + ' '.join(map(str, common_residues)))

        transformation = pr.calcTransformation(matched_ca2, matched_ca1)
        structure2_aligned = transformation.apply(self.structure2)

        pref = path_pdb2.split('/')[-1][:-4] + '_TEMP'
        print(pref)
        # Write the aligned structure2 to a temporary file
        with tempfile.NamedTemporaryFile(delete=False, suffix='.pdb',prefix= pref) as temp_file:
            pr.writePDB(temp_file.name, structure2_aligned)
            aligned_pdb2_path = temp_file.name
        # Read the aligned PDB files
        pdb1 = self.read_pdb_file(path_pdb1)
        pdb2 = self.read_pdb_file(aligned_pdb2_path)
        # Create a py3Dmol viewer
        viewer = py3Dmol.view(width=800, height=600)
        # Add the first structure with cartoon style
        viewer.addModel(pdb1, 'pdb')
        # viewer.addLabel(f"{path_pdb1.split('/')[-1][:-4]}", {'fontOpacity': 0.3}, {'resi': 158})
        print(f"Blue:{path_pdb1.split('/')[-1][:-4]}")
        viewer.setStyle({'model': 0}, {'cartoon': {'color': 'blue'}})
        # Add the aligned second structure with cartoon style
        viewer.addModel(pdb2, 'pdb')
        print(f"Red:{path_pdb2.split('/')[-1][:-4]}")
        # viewer.addLabel(f"{path_pdb2.split('/')[-1][:-4]}", {'fontOpacity': 0.3}, {'resi': 158})
        viewer.setStyle({'model': 1}, {'cartoon': {'color': 'red'}})

        # Align the second model to the first
        viewer.zoomTo()
        viewer.zoom(1.2)
        viewer.setBackgroundColor('white')
        viewer.show()


    def plot_fold_alignment_1(self, path_pdb1, path_pdb2):
            try:
                # Select CA atoms for alignment
                structure1 = pr.parsePDB(path_pdb1)
                structure2 = pr.parsePDB(path_pdb2)
                ca_atoms1 = structure1.select('name CA')
                ca_atoms2 = structure2.select('name CA')

                # Ensure CA atoms are found
                if ca_atoms1 is None or ca_atoms2 is None:
                    raise ValueError("CA atoms not found in one or both structures")

                # Get common residues
                common_residues = set(ca_atoms1.getResnums()).intersection(ca_atoms2.getResnums())
                if not common_residues:
                    raise ValueError("No common residues found for alignment")

                matched_ca1 = ca_atoms1.select('resnum ' + ' '.join(map(str, common_residues)))
                matched_ca2 = ca_atoms2.select('resnum ' + ' '.join(map(str, common_residues)))

                # Ensure matched CA atoms are found
                if matched_ca1 is None or matched_ca2 is None:
                    raise ValueError("Matched CA atoms not found in one or both structures")

                # Calculate transformation and apply to structure2
                transformation = pr.calcTransformation(matched_ca2, matched_ca1)
                structure2_aligned = transformation.apply(structure2)

                pref = path_pdb2.split('/')[-1][:-4] + '_TEMP'
                print(pref)
                # Write the aligned structure2 to a temporary file
                with tempfile.NamedTemporaryFile(delete=False, suffix='.pdb', prefix=pref) as temp_file:
                    pr.writePDB(temp_file.name, structure2_aligned)
                    aligned_pdb2_path = temp_file.name

                # Read the aligned PDB files
                pdb1 = self.read_pdb_file(path_pdb1)
                pdb2 = self.read_pdb_file(aligned_pdb2_path)

                # Create a py3Dmol viewer
                viewer = py3Dmol.view(width=800, height=600)

                # Add the first structure with cartoon style
                viewer.addModel(pdb1, 'pdb')
                print(f"Blue: {path_pdb1.split('/')[-1][:-4]}")
                viewer.setStyle({'model': 0}, {'cartoon': {'color': 'blue'}})

                # Add the aligned second structure with cartoon style
                viewer.addModel(pdb2, 'pdb')
                print(f"Red: {path_pdb2.split('/')[-1][:-4]}")
                viewer.setStyle({'model': 1}, {'cartoon': {'color': 'red'}})

                # Align the second model to the first
                viewer.zoomTo()
                viewer.zoom(1.2)
                viewer.setBackgroundColor('white')
                viewer.show()

                # Delete the temporary file
                os.remove(aligned_pdb2_path)
            except Exception as e:
                print(f"An error occurred: {e}")



    def align_and_visualize_pdb(self,pdb_file1, pdb_file2,score=None):
        # Load the PDB files using ProDy
        structure1 = pr.parsePDB(pdb_file1)
        structure2 = pr.parsePDB(pdb_file2)

        # Select only the alpha carbons for alignment
        calphas1 = structure1.select('name CA')
        calphas2 = structure2.select('name CA')

        transformation = pr.calcTransformation(calphas2, calphas1)
        transformation.apply(structure2)

        # Save the aligned second structure to a new PDB file
        pr.writePDB('aligned2.pdb', structure2)

        # Read the PDB structures
        with open(pdb_file1, 'r') as file:
            pdb1 = file.read()
        with open('aligned2.pdb', 'r') as file:
            pdb2 = file.read()

        # Initialize the 3Dmol viewer
        viewer = py3Dmol.view(width=800, height=600)

        # Add the first PDB structure
        if score != None:
            viewer.addLabel(f"TM Score:{score}", {'fontOpacity': 0.3}, {'resi': 158})

        viewer.addModel(pdb1, 'pdb')
        print(f"Blue:{pdb_file1.split('/')[-1][:-4]}")
        viewer.setStyle({'model': 0}, {'cartoon': {'color': 'blue'}})

        # Add the aligned second PDB structure
        viewer.addModel(pdb2, 'pdb')
        print(f"Red:{pdb_file2.split('/')[-1][:-4]}")
        viewer.setStyle({'model': 1}, {'cartoon': {'color': 'red'}})

        # Set the background color and zoom
        viewer.setBackgroundColor('white')
        viewer.zoomTo()

        # Show the viewer
        viewer.show()



    def plot_single_fold(self,pdb,label=''):
        viewer = py3Dmol.view(width=400, height=300)
        viewer.addModel(pdb, 'pdb')
        viewer.setStyle({'model': 0}, {'stick': {}})
        viewer.setStyle({'model': 0}, {'cartoon': {'color': 'blue'}})
        viewer.addLabel(label, {'fontOpacity': 0.3}, {'resi': 158})
        viewer.zoomTo()
        viewer.zoom(1.2)
        viewer.setBackgroundColor('white')
        viewer.show()

    def get_cmap(self,pdp_file_path):
        traj = md.load(filename_or_filenames=pdp_file_path)
        # topology = traj.topology
        frame_contacts = ContactFrequency(traj[0])
        return frame_contacts
    
    def plot_cmap(self,frame_contacts):
        fig, ax = frame_contacts.residue_contacts.plot()
        plt.xlabel("Residue")
        xlim = frame_contacts.query_residue_range[-1]
        plt.xlim([0, xlim])
        plt.ylim([0, xlim])
        _ = plt.ylabel("Residue")
        
    def plot_diff_cmap(self,cmap1,cmap2):
        diff = OverrideTopologyContactDifference(cmap1,cmap2,topology=cmap1.topology)
        diff.residue_contacts.plot(figsize=(12, 8),cmap='seismic', vmin=-1, vmax=1)
        xlim = max(cmap1.query_residue_range[-1],cmap2.query_residue_range[-1])
        plt.xlim([0, xlim])
        plt.ylim([0, xlim])
        _ = plt.ylabel("Residue")


    def get_contact_map_from_pdb(self,pdb_file, size=None, start_pos=0, end_pos=-1):
        pdb_obj = md.load_pdb(pdb_file)
        distances, pairs = md.compute_contacts(pdb_obj)
        contacts = md.geometry.squareform(distances, pairs)[0]
        if size is None:
            size = contacts.shape[0]
        arr = np.zeros([size, size])
        arr[np.where(contacts[start_pos:end_pos, start_pos:end_pos] < 0.5)] = 1
        return arr

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch

def plot_viz_cmap(file,legend_plot):
    # Create a sample 2D numpy array (replace this with your actual data)
    visualization_map = np.load(file)
    size = 50

    # Create a custom colormap
    colors = ['grey','blue','lightblue','purple','blue','magenta']
    # values = [0, 1, 1.25, 1.5, 1.75]
    bounds = [0,0.49,0.99, 1.24, 1.49, 1.74,2]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(bounds, cmap.N)

    # Create the plot
    plt.figure(figsize=(10,8))
    im = plt.imshow(visualization_map, cmap=cmap, norm=norm, interpolation='nearest',origin='lower')

    # Create legend elements
    legend_elements = [
        Patch(facecolor='lightblue', edgecolor='black', label='Experiment contact'),
        Patch(facecolor='purple', edgecolor='black', label='Unique State Experiment contact'),
        Patch(facecolor='blue', edgecolor='black', label='Experiment contact predicted'),
        Patch(facecolor='magenta', edgecolor='black', label='Unique State Experiment contact predicted')
    ]

    # Add the legend
    plt.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5))

    # Set title and labels
    plt.title(legend_plot)
    plt.xlabel('Residue Index')
    plt.ylabel('Residue Index')

    # Adjust layout to make room for the legend
    plt.tight_layout()

    # Show the plot
    plt.show()

    
        
    

