import prody as pr
import py3Dmol
import tempfile

# Function to read PDB file
def read_pdb_file(file_path):
    with open(file_path, 'r') as file:
        return file.read()

class PlotTool():
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

    def plot_single_fold(self,fold,label=''):
        pdb = read_pdb_file(f'{self.folder}/{self.fold_pair}/chain_pdb_files/{fold}.pdb')
        viewer = py3Dmol.view(width=800, height=600)
        viewer.addModel(pdb, 'pdb')
        viewer.setStyle({'model': 0}, {'stick': {}})
        viewer.setStyle({'model': 0}, {'cartoon': {'color': 'blue'}})
        viewer.addLabel(label, {'fontOpacity': 0.3}, {'resi': 158})



