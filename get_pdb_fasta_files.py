import pandas as pd
import os
from Bio.PDB import *
from Bio import SeqIO
from tqdm import tqdm
tqdm.pandas()

output_path = '/Users/steveabecassis/Desktop/Kaib/pdb_files'
fold_swiching_pair_df = pd.read_csv('/Users/steveabecassis/Desktop/Thesis/Papers Thesis/data/AlphaFold2 fails to predict protein fold switching/fold_swiching_list.csv')
fold_swiching_pair_df['Fold1_PdbId'] = fold_swiching_pair_df.Fold1.apply(lambda x: x[:-2].upper())
fold_swiching_pair_df['Fold1_chain'] = fold_swiching_pair_df.Fold1.apply(lambda x: x[-1:].upper())
fold_swiching_pair_df['Fold2_PdbId'] = fold_swiching_pair_df.Fold2.apply(lambda x: x[:-2].upper())
fold_swiching_pair_df['Fold2_chain'] = fold_swiching_pair_df.Fold2.apply(lambda x: x[-1:].upper())

fold_swiching_pair_df['Fold1_PdbId'].progress_apply(lambda x: pdb_list.retrieve_pdb_file(pdb_code=x,file_format='pdb',pdir=output_path))
fold_swiching_pair_df['Fold2_PdbId'].progress_apply(lambda x: pdb_list.retrieve_pdb_file(pdb_code=x,file_format='pdb',pdir=output_path))


pdb_id = fold_swiching_pair_df.Fold1.iloc[0][:-2].upper()
pdb_list = PDBList(server='ftp://ftp.wwpdb.org', pdb=None, obsolete_pdb=None, verbose=True)
pdb_file = pdb_list.retrieve_pdb_file(pdb_code=pdb_id,file_format='pdb',pdir=output_path)



pdb_files = os.listdir(output_path)
for pdb_f in pdb_files:
    SeqIO.convert(output_path + '/' + pdb_f, 'pdb-atom', '/Users/steveabecassis/Desktop/Kaib/fasta_files/'+ pdb_f[3:-4] + '.fasta', 'fasta')




import requests as r
from Bio import SeqIO
from io import StringIO
def get_seq_from_uniprotId(cID):
    # cID='A0A3N5IQ47'
    baseUrl="http://www.uniprot.org/uniprot/"
    currentUrl=baseUrl+cID+".fasta"
    response = r.post(currentUrl)
    cData=''.join(response.text)
    Seq=StringIO(cData)
    pSeq=list(SeqIO.parse(Seq,'fasta'))
    return str(pSeq[0].seq)



