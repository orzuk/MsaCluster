import os
import subprocess
from Bio import SeqIO

a3m_file = "/Users/steveabecassis/Desktop/FastTree/2QKEE_002.a3m"
a3m_file = '/Users/steveabecassis/Desktop/FastTree/2QKEE_000_FULL_1.a3m'# Replace with the path to your A3M alignment file
sequences = []
df = open(f'/Users/steveabecassis/Desktop/FastTree/msa1.a3m', 'w')
names = []
l = []
with open(a3m_file, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        # sequences.append(str(record.seq))
        if record.name not in names:
            df.write(f">{record.name}\n{record.seq}\n")
            names.append(record.name)
            l.append(len(record.seq[:388]))
        else:
            continue




from Bio import SeqIO
import numpy as np
import copy
from collections import Counter

msa_in_f =  "/Users/steveabecassis/Desktop/FastTree/2QKEE_000_FULL_1.a3m"

prefix_f, suffix_f = msa_in_f.rsplit('.')[0], msa_in_f.rsplit('.')[1]

msa_seqs = list(SeqIO.parse(msa_in_f, 'fasta'))

headers = [str(item.id) for item in msa_seqs]
seqs = [str(item.seq)[:388] for item in msa_seqs]
lenghts = [len(str(item.seq)) for item in msa_seqs]
seqs = [str(item.seq)[:min(lenghts)] for item in msa_seqs]

# th = np.quantile(lenghts,0.05)
# idxs = [index for index, num in enumerate(lenghts) if num < th]


# seqs = seqs[1:]
# headers = headers[1:]

np_headers = np.array(headers)

concat_seqs = ''.join(seqs)

matrix = np.array(list(concat_seqs), dtype=str).reshape(len(headers), len(seqs[0]))

clean_mat = copy.deepcopy(matrix)
n_rows, n_cols = matrix.shape

#remove columns having more than 50% gaps
remove_cols = list()
for i in range(n_cols):
    col = str(''.join(matrix[:, i]))
    if col.count('-') > round(n_rows/2):
        remove_cols.append(i)

clean_mat = np.delete(matrix, remove_cols, 1)

#remove columns where the consensus amino acid is less than 25%
remove_cols = list()
for i in range(clean_mat.shape[1]):
    col = str(''.join(clean_mat[:, i]))
    char_counter = Counter(col)
    del char_counter['-']
    if max(char_counter.values())/n_rows < 0.25:
        remove_cols.append(i)

clean_mat = np.delete(clean_mat, remove_cols, 1)

headers = [i.replace('.','').replace(':','') for i in   headers]
with open(prefix_f + '_clean_columns.'+ suffix_f, 'w') as out_f:
    for i in range(clean_mat.shape[0]):
        row = str(''.join(clean_mat[i, :]))
        if headers.count(headers[i])<2:
            out_f.write('>'+str(headers[i])+'\n')
            out_f.write(row+'\n')
        else:
            print(headers[i])
            continue

#remove rows where the the alignment is more than 50% gaps
new_n_rows, new_n_cols = clean_mat.shape
remove_rows = list()
for i in range(clean_mat.shape[0]):
    row = str(''.join(clean_mat[i, :]))
    if row.count('-') > new_n_cols/2:
        remove_rows.append(i)
clean_mat_reduced_rows = np.delete(clean_mat, remove_rows, 0)
np_headers_reduced_rows = np.delete(np_headers, remove_rows)
np_headers_reduced_rows = [i.replace('.','').replace(':','') for i in   np_headers_reduced_rows]

with open(prefix_f + '_clean_columns_clean_rows.'+ suffix_f, 'w') as out_f:
    for i in range(clean_mat_reduced_rows.shape[0]):
        row = str(''.join(clean_mat_reduced_rows[i, :]))
        out_f.write('>'+str(np_headers_reduced_rows[i])+'\n')
        out_f.write(row+'\n')
        if list(np_headers_reduced_rows).count(np_headers_reduced_rows[i])<2:
            row = str(''.join(clean_mat_reduced_rows[i, :]))
            out_f.write('>'+str(np_headers_reduced_rows[i])+'\n')
            out_f.write(row+'\n')
        else:
            print(headers[i])
            continue

print('original input matrix shape,', matrix.shape)
print('clean columns matrix shape', clean_mat.shape)
print('clean columns, clean rows matrix shape', clean_mat_reduced_rows.shape)

with open(prefix_f + 'reduction_stats.txt', 'w') as out_f:
    out_f.write('original input matrix shape,' + str( matrix.shape)+'\n')
    out_f.write('clean columns matrix shape' + str(clean_mat.shape)+'\n')
    out_f.write('clean columns, clean rows matrix shape' + str( clean_mat_reduced_rows.shape)+'\n')