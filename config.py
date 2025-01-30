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
user = "zuk_laptop"  # Options: "steve_moriah", "zuk_laptop"


if os.getenv('USER') == 'steveabecassis':
    DATA_DIR = '/Users/steveabecassis/Desktop/Pipeline'
    OUTPUT_PATH_NOTEBOOKS = '/Users/steveabecassis/Desktop/Pipeline_res/HTMLs_new_3001'
    CMAP_RES_PATH = '/Users/steveabecassis/PycharmProjects/MsaCluster/data/df_cmap_all.csv'
    LOCAL_RUN = True



