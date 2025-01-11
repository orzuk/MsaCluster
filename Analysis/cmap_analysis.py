import numpy as np
import pandas as pd
# import pcmap
import Bio.PDB

pdb_code = "2nxq"
pdb_filename = "/Users/steveabecassis/Desktop/PipelineTest/output_pipeline_1jfk/2nxq.pdb"


def calc_residue_dist(residue_one, residue_two) :
    """Returns the C-alpha distance between two residues"""
    diff_vector = residue_one["CA"].coord - residue_two["CA"].coord
    return np.sqrt(np.sum(diff_vector * diff_vector))


def calc_dist_matrix(chain_one, chain_two) :
    """Returns a matrix of C-alpha distances between two chains"""
    answer = np.zeros((len(chain_one), len(chain_two)), np.float)
    for row, residue_one in enumerate(chain_one) :
        for col, residue_two in enumerate(chain_two) :
            answer[row, col] = calc_residue_dist(residue_one, residue_two)
    return answer


a = np.loadtxt('/Users/steveabecassis/Desktop/PipelineTest/output_pipeline_1jfk/esm_cmap_output/msa_t__cluster_000.npy')

# c2 = pcmap.contactMap("/Users/steveabecassis/Desktop/PipelineTest/output_pipeline_1jfk/2nxq.pdb")  # can't install pcmap !!

structure = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_filename)
model = structure[0]

dist_matrix = calc_dist_matrix(model["A"], model["A"])
contact_map = dist_matrix < 12.0

import pylab
pylab.matshow(np.transpose(contact_map))
pylab.colorbar()
pylab.show()


for i,j in enumerate(model["A"]):
    print(j["CA"].coord)

["CA"].coord

aa= model["A"]
model["A"].keys()

import contact_map as cme

path = "/Users/steveabecassis/Desktop/PipelineTest/output_pipeline_1jfk/1jfk.pdb"
import matplotlib.pyplot as plt
import mdtraj as md
traj = md.load(filename_or_filenames=path)
topology = traj.topology
from contact_map import ContactFrequency, ContactDifference
frame_contacts = ContactFrequency(traj[0])
fig, ax = frame_contacts.residue_contacts.plot()
plt.xlabel("Residue")
_ = plt.ylabel("Residue")
a = frame_contacts.residue_contacts

z = pd.DataFrame(a.df.values)




