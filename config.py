import os
import platform

def linux_to_windows_path(linux_path: str) -> str:
    """
    Convert a Linux-style path to Windows-style path.
    Example: '/mnt/c/Code/Github/MsaCluster' -> 'C:/Code/Github/MsaCluster'
    """
    if linux_path.startswith('/mnt/'):
        drive_letter = linux_path[5].upper() + ':'
        return drive_letter + linux_path[6:].replace('/', '\\')
    return linux_path.replace('/', '\\')

def windows_to_linux_path(windows_path: str) -> str:
    """
    Convert a Windows-style path to Linux-style path.
    Example: 'C:/Code/Github/MsaCluster' -> '/mnt/c/Code/Github/MsaCluster'
    """
    if ':' in windows_path:
        drive_letter = windows_path[0].lower()
        return f"/mnt/{drive_letter}" + windows_path[2:].replace('\\', '/')
    return windows_path.replace('\\', '/')


"""
Configuration file for the project.

This file contains all paths and constants required for the project,
specific to different users ('zuk_laptop' and 'steve_moriah').

To add a new user, extend the logic below with the appropriate paths.
"""

# Define the user: change this to match the current user's environment
user =  os.getenv('USER') #  "zuk_laptop"  # Options: "steve_moriah", "zuk_laptop"


if user is None: # for local run
    print("Current directory: ", os.getcwd())
    if 'C:\\Code\\Github\\MsaCluster' in os.getcwd():
        user = "zuk_laptop"
    else:
        user = 'steve_laptop'


print("RUNNING FOLD-SWITCH PIEPLINE WITH USER: " +  user + "  ENVIRONMENT: " +  platform.system())

if user == 'steveabecassis':
    MAIN_DIR = '/Users/steveabecassis/Desktop'
    DATA_DIR = MAIN_DIR + '/Pipeline'
    TABLES_RES = MAIN_DIR + '/Pipeline_res'
    OUTPUT_PATH_NOTEBOOKS = TABLES_RES + 'HTMLs_new_3001'
    CMAP_RES_PATH = '/Users/steveabecassis/PycharmProjects/MsaCluster/data/df_cmap_all.csv'
    SUMMARY_RESULTS_TABLE = TABLES_RES + '/final_res_df_2510.parq'
    SIMILARITY_RESULTS_TABLE = TABLES_RES + '/fold1_fold2_sim.parq'
    LOCAL_RUN = True

    # File create by the script cmap_analysis.py (in the analysis folder)
    CMAP_ANALYSIS_FILE = DATA_DIR + '/cmap_exact_analysis_tol0_2510.parq'
    # File create by the script esmfold_analysis.py (in the analysis folder)
    ESMFOLD_ANALYSIS_FILE = TABLES_RES + '/df_esmfold_analysis.csv'

    # Set directories for each part of the pipeline
    PAIR_DIRS = {"AF":"AF_Preds", "ESMFold":"esm_fold_output", "Cmap":"Cmap_Preds", "MSA":"MSA_Preds", "Cluster":"Cluster" }


if user in ['zuk_laptop', 'orzuk']:
    if user == 'zuk_laptop':  # WINDOWS PC
        MAIN_DIR = 'C://Code//Github//MsaCluster'
    else:
        MAIN_DIR = '/mnt/c/Code/Github/MsaCluster'
    DATA_DIR = MAIN_DIR + '/Pipeline'
    TABLES_DIR = MAIN_DIR + '/data'
    OUTPUT_PATH_NOTEBOOKS = DATA_DIR + '/HTMLs_new_3001'
    FIGURE_RES_DIR = DATA_DIR + '/Results/Figures'
    CMAP_RES_PATH = DATA_DIR + '/data/df_cmap_all.csv'
    AF_MODEL_FILE = TABLES_DIR + '/df_af_all.csv'
    ESMF_MODEL_FILE = TABLES_DIR + '/df_esmfold_all.csv'
    MSA_TRANS_MODEL_FILE = TABLES_DIR + '/df_cmap_all.csv'

    TABLES_RES = MAIN_DIR + '/Pipeline_res'
    SUMMARY_RESULTS_TABLE = TABLES_RES + '/final_res_df_2510.parq'
    SIMILARITY_RESULTS_TABLE = TABLES_RES + '/fold1_fold2_sim.parq'
    ENERGY_DIR = DATA_DIR + '/output_deltaG'
    ENERGY_FILE = ENERGY_DIR + '/deltaG_results.txt'


    LOCAL_RUN = False



# Which of the five AF models to use for prediction
AF2_MODEL = 2
