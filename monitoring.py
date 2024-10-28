import os
import pandas as pd

main_folder = '/Users/steveabecassis/Desktop/Pipeline'
#output_msa_cluster
fold_pairs = [i for i in os.listdir(main_folder) if ('_' in i) and ('.' not in i)]
fold_pair_nb = [i[:-5] for i in os.listdir('/Users/steveabecassis/Desktop/Pipeline_res/HTMLs_new_2110')]
res = []
for fold_pair in fold_pairs:
    cluster_count = len(os.listdir(f'{main_folder}/{fold_pair}/output_msa_cluster'))
    af_pred_txt = os.listdir(f'{main_folder}/{fold_pair}/AF_preds')
    AF_pred_count = max(len( [i for i in af_pred_txt if 'txt' in i]) - 2 , 0)
    cmap_count = len(os.listdir(f'{main_folder}/{fold_pair}/output_cmap_esm')) - 2
    cmap_viz_count = len([i for i in os.listdir(f'{main_folder}/{fold_pair}/output_cmap_esm/VizCmaps') if 'tol_0' in i])/2
    have_notebook = fold_pair in fold_pair_nb
    res.append({'fold_pair':fold_pair,'cluster_count':cluster_count,'AF_pred_count':AF_pred_count,'cmap_count':cmap_count,'cmap_viz_count':cmap_viz_count,'have_notebook':have_notebook})

df  = pd.DataFrame(res)
df_ = df[df.AF_pred_count != 0]
df_nb_check = df_[~df['have_notebook']]
len(df_nb_check)


'''
1) Add contact map in monitoring DONE
2) Add esm fold monitoring       DONE
3) Fix contact map               DONE
4) Fix notebooks                 DONE
5) Full runs get results         DONE
6) Get full satistics table      DONE
7) Push new notebook and new tables to github
8) Fix thesis text with new cmaps etc...
'''
