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
        MAIN_DIR = '/Users/steveabecassis/Desktop'
    CMAP_RES_PATH = '/Users/steveabecassis/PycharmProjects/MsaCluster/data/df_cmap_all.csv'
    LOCAL_RUN = True

    # File create by the script cmap_analysis.py (in the analysis folder)
    CMAP_ANALYSIS_FILE = DATA_DIR + '/cmap_exact_analysis_tol0_2510.parq'
    # File create by the script esmfold_analysis.py (in the analysis folder)

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

# Set subdirectories for each part of the pipeline
PAIR_DIRS = {"AF":"output_AF", "ESMFold":"output_esm_fold", "Cmap":"output_cmaps", "MSA":"output_get_msa", "Cluster":"output_msa_cluster" }


DATA_DIR = MAIN_DIR + '/Pipeline'
TABLES_DIR = MAIN_DIR + '/data'
CMAP_RES_PATH = DATA_DIR + '/data/df_cmap_all.csv'
AF_MODEL_FILE = TABLES_DIR + '/df_af_all.csv'
ESMF_MODEL_FILE = TABLES_DIR + '/df_esmfold_all.csv'
MSA_TRANS_MODEL_FILE = TABLES_DIR + '/df_cmap_all.csv'


TMALIGN_EXE = MAIN_DIR + '/TMalign'
USE_TMALIGN_BINARY = (os.path.isfile(TMALIGN_EXE) and os.access(TMALIGN_EXE, os.X_OK))

CCMPRED_EXE = MAIN_DIR.replace("MsaCluster", "CCMpred/bin/ccmpred")

PYMOL_EXE = MAIN_DIR.replace("github/MsaCluster", "/sci/labs/orzuk/orzuk/pymol-venv/bin/pymol")

# Relative paths common to all users
TABLES_RES = MAIN_DIR + '/docs'
REPO_NAME = os.path.basename(MAIN_DIR.rstrip('/\\'))
SITE_BASE  = f"/{REPO_NAME}"
GITHUB_URL_HTML = f"{SITE_BASE}/HTML"


OUTPUT_PATH_NOTEBOOKS = TABLES_RES + '/HTML'
SUMMARY_RESULTS_TABLE = TABLES_RES + '/summary_final_res_all_pairs_df.csv' #  '/final_res_df_2510.parq'
DETAILED_RESULTS_TABLE = TABLES_RES + '/detailed_final_res_all_pairs_df.csv' #  '/final_res_df_2510.parq'
SIMILARITY_RESULTS_TABLE = TABLES_RES + '/fold1_fold2_sim.parq'
ESMFOLD_ANALYSIS_FILE = TABLES_RES + '/df_esmfold_analysis.csv'
ENERGY_DIR = DATA_DIR + '/output_deltaG'
ENERGY_FILE = ENERGY_DIR + '/deltaG_results.txt'

# Local figures dir
FIGURE_RES_DIR = DATA_DIR + '/Results/Figures'

# Where to publish images for GitHub Pages (under docs/...)
FIGURE_PUBLISH_DIR = os.path.join(TABLES_RES, 'figures', 'global')   # docs/figures/global
PER_PAIR_PUBLISH_SUBDIR = os.path.join('HTML', 'figs')               # docs/HTML/figs/<pair_id> (per-pair)


# --- add near the bottom of config.py (after existing constants) ---

# Inclusive alias map for per-pair outputs (used by targeted clean)
# Keys are common names/commands; values are RELATIVE dirs inside Pipeline/<pair_id>
PAIR_TARGET_DIRS = {
    # MSA acquisition
    "get_msa":            ["output_get_msa", "tmp_msa_files"],
    "msa":                ["output_get_msa", "tmp_msa_files"],
    "MSA":                ["output_get_msa", "tmp_msa_files"],
    "colabfold":          ["output_get_msa", "tmp_msa_files"],

    # MSA clustering
    "msa_clusters":       ["output_msa_cluster", "output_msa_cluster_pairtrim"],
    "cluster":            ["output_msa_cluster"],
    "Cluster":            ["output_msa_cluster"],
    "clusters":           ["output_msa_cluster"],
    "shallow_msa":        ["output_msa_cluster"],

    # Contact maps
    "cmap":               ["output_cmaps/msa_transformer", "output_cmaps/ccmpred"],
    "Cmap":               ["output_cmaps/msa_transformer", "output_cmaps/ccmpred"],
    "cmap_msa_transformer": ["output_cmaps/msa_transformer"],
    "msa_transformer":    ["output_cmaps/msa_transformer"],
    "cmap_ccmpred":       ["output_cmaps/ccmpred"],
    "ccmpred":            ["output_cmaps/ccmpred"],

    # Structure prediction
    "esmfold":            ["output_esm_fold", "tmp_esmfold"],
    "esm":                ["output_esm_fold", "tmp_esmfold"],
    "esm_fold":           ["output_esm_fold", "tmp_esmfold"],
    "ESMFold":            ["output_esm_fold", "tmp_esmfold"],

    "alphafold":          ["output_AF", "tmp_AF"],
    "alphafold2":         ["output_AF", "tmp_AF"],
    "alphafold3":         ["output_AF", "tmp_AF"],
    "af":                 ["output_AF", "tmp_AF"],
    "AF":                 ["output_AF", "tmp_AF"],

    # Energy / ΔG (if you produce per-pair dirs; keep if relevant)
    "deltag":             ["output_deltaG"],
    "deltaG":             ["output_deltaG"],
    "energy":             ["output_deltaG"],

    # Analysis / plots / trees
    "plots":              ["Analysis"],
    "analysis":           ["Analysis"],
    "phytree":            ["output_phytree"],
    "tree":               ["output_phytree"],

    # Logs / scratch
    "logs":               ["logs"],        # plus root *.out/*.log handled in cleaner
    "tmp":                ["jobs", "tmp_msa_files", "tmp_esmfold"],
}



# Whether to auto-copy figures from compute dirs to publish dirs when building HTML
PUBLISH_FIGURES = True
LOCAL_RUN = False

PAIR_DIR_RE = re.compile(r'^([0-9A-Za-z]{4}[A-Za-z0-9])_([0-9A-Za-z]{4}[A-Za-z0-9])$')

# Which of the five AF2 models to use for prediction
AF2_MODEL = 2

CONTACT_CUTOFF = 8.0  # Å; adjust if you use a different truth definition


# ------------ RENAME: USED IN SOME CODE PLACES ! -------------------- #
PIPELINE_DIR = DATA_DIR
DOCS_DIR = TABLES_RES
# ------------ END RENAME: USED IN SOME CODE PLACES ! -------------------- #
