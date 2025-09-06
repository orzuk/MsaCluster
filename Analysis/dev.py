import pandas as pd
import os
from Bio import SeqIO

df_esmfold_analysis = pd.read_csv('/Users/steveabecassis/Desktop/Pipeline_res/df_esmfold_analysis.csv')
df_esmfold_analysis.fold_pair.nunique()

summary_table = pd.read_parquet('/Users/steveabecassis/Desktop/Pipeline_res/final_res_df_2510.parq')
summary_table.fold_pair.nunique()

esm_err_fold_pair = [i for i in summary_table.fold_pair.unique() if i not in df_esmfold_analysis.fold_pair.unique()]

list_len = []
fold_pairs = os.listdir('/Users/steveabecassis/Desktop/Pipeline')
for fold_pair in fold_pairs:
    try:
        if '_' not in fold_pair:
            continue
        sequence_lengths = {}
        fold1 = fold_pair.split('_')[0]
        fasta_file = f'/Users/steveabecassis/Desktop/Pipeline/{fold_pair}/fasta_chain_files/{fold1}.fasta'
        with open(fasta_file, "r") as file:
            for record in SeqIO.parse(file, "fasta"):
                sequence_lengths[record.id] = len(record.seq)
        length = sequence_lengths.get('')
        list_len.append({'fold_pair':fold_pair,'L':length})
    except:
        print(fold_pair)
        continue


df_len = pd.DataFrame(list_len)
df_len.sort_values(by='L',inplace=True,ascending=False)
fold_pair = '1jfkA_2nxqB'


r = []
esm_err_fold_pair = []
for fold_pair in fold_pairs:
    try:
        L=len(os.listdir(f'/Users/steveabecassis/Desktop/Pipeline/{fold_pair}/output_esm_fold'))
        if L==0:
            r.append({'fold_pair':fold_pair,'esm_fold_len_pred':L})
    except:
        continue

esm_err_fold_pair_df = pd.DataFrame(r)

'2jmrA_4j3oF'

r = []
cmap_err_fold_pair = []
for fold_pair in fold_pairs:
    try:
        L=len(os.listdir(f'/Users/steveabecassis/Desktop/Pipeline/{fold_pair}/output_cmaps/msa_transformer'))
        if L<3:
            cmap_err_fold_pair.append({'fold_pair':fold_pair,'esm_fold_len_pred':L})
    except:
        continue

cmap_err_fold_pair_df = pd.DataFrame(cmap_err_fold_pair)
cmap_err_fold_pair_df.sort_values(by='esm_fold_len_pred',ascending=False,inplace=True)


'''
CONTACT MAP
4cmqB_4zt0C Cmap is not working

ESM FOLD
EsmFold not worked:
1. Too many lot of clusters
2. Too big input sequence (4cmqB_4zt0C)

    fold_pair    esm_fold_len_pred
   1xjuB_1xjtA          0
   2bzyB_2lqwA          0
   3v0tA_3uyiA          0
   4q79F_4uv2D          0
   2grmB_2axzA          0
   4cmqB_4zt0C          0
   2kkwA_2n0aD          0
   2lepA_4hddA          0
   2jmrA_4j3oF          0
   2pbkB_3njqA          0
   4xwsD_4y0mJ          0
'''
