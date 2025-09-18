import os, re
import platform

"""
Configuration file for the project.

This file contains all paths and constants required for the project,
specific to different users ('zuk_laptop' and 'steve_moriah').

To add a new user, extend the logic below with the appropriate paths.
"""


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


def is_pycharm() -> bool:
    # PyCharm typically sets one of these
    return os.environ.get("PYCHARM_HOSTED") == "1" or "PYCHARM_DISPLAY_PORT" in os.environ


# Define the user: change this to match the current user's environment
user =  os.getenv('USER') #  "zuk_laptop"  # Options: "steve_moriah", "zuk_laptop"
moriah = os.path.isdir('/sci/labs/orzuk')  # Check if running locally or in moriah cluster

if user is None: # for local run
    print("Current directory: ", os.getcwd())
    if 'C:\\Code\\Github\\MsaCluster' in os.getcwd():
        user = "zuk_laptop"
    else:
        user = 'steve_laptop'


print("CONFIG FOLD-SWITCH PIEPLINE WITH USER: " +  user + "  ENVIRONMENT: " +  platform.system())

if user == 'steveabecassis':
    if moriah:
        MAIN_DIR = '/sci/labs/orzuk/steveabecassis/MsaCluster'
    else:
        MAIN_DIR = '/Users/MICHAEL_JORDAN/Desktop' # '/Users/steveabecassis/Desktop'
    DATA_DIR = MAIN_DIR + '/Pipeline'
    CMAP_RES_PATH = '/Users/steveabecassis/PycharmProjects/MsaCluster/data/df_cmap_all.csv'
    LOCAL_RUN = True

    TMALIGN_EXE = MAIN_DIR + '/TMalign'

    # File create by the script cmap_analysis.py (in the analysis folder)
    CMAP_ANALYSIS_FILE = DATA_DIR + '/cmap_exact_analysis_tol0_2510.parq'
    # File create by the script esmfold_analysis.py (in the analysis folder)

    # Set subdirectories for each part of the pipeline
    PAIR_DIRS = {"AF":"output_AF", "ESMFold":"output_esm_fold", "Cmap":"output_cmaps", "MSA":"output_get_msa", "Cluster":"output_msa_cluster" }


if user in ['zuk_laptop', 'orzuk']:
    if user == 'zuk_laptop':  # WINDOWS PC
        MAIN_DIR = 'C://Code//Github//MsaCluster'
        ESM_PATH = 'C://Code//Github//esm'
    else:
        if moriah:
            MAIN_DIR = '/sci/labs/orzuk/orzuk/github/MsaCluster'
            ESM_PATH = '/sci/labs/orzuk/orzuk/github/esm'
        else:
            MAIN_DIR = '/mnt/c/Code/Github/MsaCluster'
            # Scripts for running different tools
            ESM_PATH = "/mnt/c/Code/Github/esm"

    DATA_DIR = MAIN_DIR + '/Pipeline'
    TABLES_DIR = MAIN_DIR + '/data'
    FIGURE_RES_DIR = DATA_DIR + '/Results/Figures'
    CMAP_RES_PATH = DATA_DIR + '/data/df_cmap_all.csv'
    AF_MODEL_FILE = TABLES_DIR + '/df_af_all.csv'
    ESMF_MODEL_FILE = TABLES_DIR + '/df_esmfold_all.csv'
    MSA_TRANS_MODEL_FILE = TABLES_DIR + '/df_cmap_all.csv'
    TMALIGN_EXE = MAIN_DIR + '/TMalign'

USE_TMALIGN_BINARY = (os.path.isfile(TMALIGN_EXE) and os.access(TMALIGN_EXE, os.X_OK))


# Relative paths common to all users
TABLES_RES = MAIN_DIR + '/docs'
GITHUB_URL_HTML = "https://steveabecassis.github.io/MsaCluster/docs/HTML/"

OUTPUT_PATH_NOTEBOOKS = TABLES_RES + 'HTML'
SUMMARY_RESULTS_TABLE = TABLES_RES + '/summary_final_res_all_pairs_df.csv' #  '/final_res_df_2510.parq'
DETAILED_RESULTS_TABLE = TABLES_RES + '/detailed_final_res_all_pairs_df.csv' #  '/final_res_df_2510.parq'
SIMILARITY_RESULTS_TABLE = TABLES_RES + '/fold1_fold2_sim.parq'
ESMFOLD_ANALYSIS_FILE = TABLES_RES + '/df_esmfold_analysis.csv'
ENERGY_DIR = DATA_DIR + '/output_deltaG'
ENERGY_FILE = ENERGY_DIR + '/deltaG_results.txt'


LOCAL_RUN = False

PAIR_DIR_RE = re.compile(r'^([0-9A-Za-z]{4}[A-Za-z0-9])_([0-9A-Za-z]{4}[A-Za-z0-9])$')

# Which of the five AF models to use for prediction
AF2_MODEL = 2
