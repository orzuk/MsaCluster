import tmscoring
import os
import pandas as pd



if __name__ == '__main__':
    p_5jyt ='/Users/steveabecassis/Desktop/Kaib/PDB/5jyt.pdb'
    p_2qke ='/Users/steveabecassis/Desktop/Kaib/PDB/2qke.pdb'

    folders = '/Users/steveabecassis/Desktop/Kaib/AF_cluster/HDBSCAN15'
    af_outputs_folder = os.listdir(folders)

    res = []
    for folder in af_outputs_folder:
        if '73815' in str(folder):
            pdb_files = os.listdir(folders +'/'+ folder)
            for pdb_file in pdb_files:
                if 'pdb' in str(pdb_file):
                    pdb_file_path = folders +'/'+ folder+'/'+pdb_file
                    p_2qke_temp = tmscoring.get_tm(pdb_file_path, p_2qke)
                    p_5jyt_temp = tmscoring.get_tm(pdb_file_path, p_5jyt)
                    temp = {'folder':folder,'pdb_file':pdb_file,'score_5jyt':p_5jyt_temp,'score_2qke':p_2qke_temp}
                    res.append(temp)
    print('Finish !')


    df = pd.DataFrame(res)
    # df.to_parquet('/Users/steveabecassis/Desktop/Kaib/AF_cluster/HDBSCAN15/tm_scores.parq')
    df = pd.read_parquet('/Users/steveabecassis/Desktop/Kaib/AF_cluster/HDBSCAN15/tm_scores.parq')





