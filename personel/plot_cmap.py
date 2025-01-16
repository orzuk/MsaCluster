import numpy as np
import os
import mdtraj as md
import numpy as np
from matplotlib.pyplot import *
import matplotlib.pyplot as plt


if __name__ == '__main__':
    cmp_folder = '/Users/steveabecassis/Desktop/Kaib/esm/'
    cmp_files = [cmp for cmp in os.listdir(cmp_folder) if 'npy' in str(cmp)]
    for cmp in cmp_files:
        fig, ax = plt.subplots()
        a = np.loadtxt(cmp_folder+cmp)
        ax.imshow(np.clip(a,0,1),origin='lower', cmap='Blues')
        fig.savefig(cmp_folder+'img/'+ cmp +'.png')
        fig.close()
