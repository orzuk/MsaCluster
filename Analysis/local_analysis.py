import os
import pandas
import pandas as pd

pipeline_path = '/Users/steveabecassis/Desktop/Pipeline_res/Pipeline'
pipeline_files = os.listdir(pipeline_path)
protein_pairs_folders = [i for i in pipeline_files if (len(i)==11 and '.sh' not in i)]

# c=-1
# for i in protein_pairs_folders[79]:
#     c+=1
#     if i=='1jfkA_2nxqB':
#         print(c)
#         break
#

final_res = []
for i in range(len(protein_pairs_folders)):


    # i = 79

    len_msa_cluster = len(os.listdir(f'{pipeline_path}/{protein_pairs_folders[i]}/output_msa_cluster'))
    len_cmap_esm = len(os.listdir(f'{pipeline_path}/{protein_pairs_folders[i]}/output_cmap_esm'))

    df_AF_preds = pd.DataFrame(os.listdir(f'{pipeline_path}/{protein_pairs_folders[i]}/AF_preds'))
    if len(df_AF_preds) > 0:
        df_AF_preds = pd.DataFrame(os.listdir(f'{pipeline_path}/{protein_pairs_folders[i]}/AF_preds'))
        df_AF_preds.columns = ['name']
        df_AF_preds = df_AF_preds[df_AF_preds.name.str.contains('pdb')][df_AF_preds.name.str.contains('Sha')]
        df_AF_preds['cluster'] = df_AF_preds.name.apply(lambda x: int(x[11:14]))
        len_AF_preds = len(df_AF_preds['cluster'].drop_duplicates())
    else:
        len_AF_preds = 0

    df_esm_fold = pd.DataFrame(os.listdir(f'{pipeline_path}/{protein_pairs_folders[i]}/output_esm_fold'))
    if len(df_esm_fold) > 0:
        df_esm_fold.columns = ['name']
        df_esm_fold = df_esm_fold[df_esm_fold.name.str.contains('pdb')][df_esm_fold.name.str.contains('Sha')]
        df_esm_fold['cluster'] = df_esm_fold.name.apply(lambda x: int(x[11:14]))
        len_EsmFold_preds = len(df_esm_fold['cluster'].drop_duplicates())
    else:
        len_EsmFold_preds = 0
    final_res.append({'protein_pair':protein_pairs_folders[i],'len_msa_cluster':len_msa_cluster,'len_cmap_esm':len_cmap_esm,
                       'len_AF_preds':len_AF_preds,'len_EsmFold_preds':len_EsmFold_preds})


summary_df = pd.DataFrame(final_res)
summary_df['MsaCluster_SUCCESS'] = (summary_df['len_msa_cluster']>0).astype(int)
summary_df['CmapESM_SUCCESS']    = (summary_df['len_msa_cluster'] == (summary_df['len_cmap_esm']-1)).astype(int)
summary_df['AF_SUCCESS']         = (summary_df['len_msa_cluster'] == (summary_df['len_AF_preds'])).astype(int)
summary_df['EsmFold_SUCCESS']         = (summary_df['len_msa_cluster'] == (summary_df['len_EsmFold_preds'])).astype(int)

summary_df[summary_df['len_msa_cluster']>0]['MsaCluster_SUCCESS'].mean()
summary_df[summary_df['len_msa_cluster']>0]['CmapESM_SUCCESS'].mean()
summary_df[summary_df['len_msa_cluster']>0]['AF_SUCCESS'].mean()
summary_df[summary_df['len_msa_cluster']>0]['EsmFold_SUCCESS'].mean()

summary_df.sort_values(by='protein_pair',inplace=True)
# summary_df[summary_df['len_msa_cluster']>0][summary_df['len_msa_cluster'] != (summary_df['len_cmap_esm']-1)]
# summary_df[summary_df['len_msa_cluster']>0][summary_df['len_msa_cluster'] != (summary_df['len_AF_preds'])]

success_data = summary_df[summary_df['len_msa_cluster']>0][summary_df['len_msa_cluster'] == (summary_df['len_AF_preds'])]



#summary_df['EsmFold_SUCCESS'].sum()





