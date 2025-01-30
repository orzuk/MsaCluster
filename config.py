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
user = "zuk_laptop"  # Options: "steve_moriah", "zuk_laptop"

# User-specific configurations
if user == "steve_moriah":
    PROJECT_DIR = "??" # Change

elif user == "zuk_laptop":
    PROJECT_DIR = "/mnt/c/Code/Github/MsaCluster"
else:
    raise ValueError(f"Unknown user: {user}")

# Adjust paths for running in pc
if not platform.system() == "Linux":
    PROJECT_DIR = linux_to_windows_path(PROJECT_DIR)

# Add any additional global constants here
# File paths for AlphaFold
DATA_DIR = PROJECT_DIR + "/Data"
AF_MODEL_FILE = DATA_DIR + "/df_af_all.csv"
MSA_TRANS_MODEL_FILE = DATA_DIR +  "/df_cmap_all.csv"
ESMF_MODEL_FILE = DATA_DIR + "/df_esmfold_all.csv"

# Directory paths
PDF_FILES_PATH = DATA_DIR + "/pdb_file"
