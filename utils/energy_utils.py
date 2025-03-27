try:
    import pyrosetta
    PYROSETTA_AVAILABLE = True
except ImportError:
    print("Can't run pyrosetta in windows. Re-run in linux")
    PYROSETTA_AVAILABLE = False

import csv
import matplotlib.pyplot as plt
import os
import sys
import numpy as np
# sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
#current_dir = os.path.dirname(os.path.abspath(__file__))
#parent_dir = os.path.dirname(current_dir)
#if parent_dir not in sys.path:
#    sys.path.append(parent_dir)

from utils.protein_utils import *
from Bio.pairwise2 import format_alignment


# Direct the output of rosetta to log file
class RedirectStdStreams:
    """
    Context manager to redirect stdout and stderr.
    """
    def __init__(self, stdout=None, stderr=None):
        self._stdout = stdout or sys.stdout
        self._stderr = stderr or sys.stderr

    def __enter__(self):
        self.old_stdout, self.old_stderr = sys.stdout, sys.stderr
        sys.stdout, sys.stderr = self._stdout, self._stderr

    def __exit__(self, exc_type, exc_value, traceback):
        sys.stdout, sys.stderr = self.old_stdout, self.old_stderr


def compute_global_and_residue_energies(pdb_pairs, foldpair_ids, output_dir):
    """
    Pipeline to compute global and residue-specific energies for PDB pairs.

    Parameters:
    - pdb_pairs (list of tuples): List of PDB file pairs [(pdb1, chain1), (pdb2, chain2)].
    - output_dir (str): Directory to save results.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # File to save global ΔG results
    global_results_file = os.path.join(output_dir, "deltaG_results.txt")

    with open(global_results_file, "w", newline="") as global_file:
        writer = csv.writer(global_file, delimiter="\t")
        writer.writerow(["Pair", "PDB_ID", "Delta_G", "Residue_Count", "Chain_Count"])

        print("All pdb-pairs: ", pdb_pairs)
#        for pdb_pair_file, foldpair_id in zip(pdb_pair_files, foldpair_ids):
        for pair in pdb_pairs:
            pdb1_path, pdb2_path = pair # pdb1, pdb2
            pdb1 = os.path.basename(pdb1_path)[:-4]
            pdb2 = os.path.basename(pdb2_path)[:-4]

            print(f"Processing pair: {pdb1_path} and {pdb2_path}")

            # Compute global ΔG and residue-specific energies for both PDBs
            residue_energies = [[], []]
            ctr = 0
            for pdb in pair:
                pdb_file = os.path.basename(pdb) # .replace(".pdb", f"{chain}.pdb")
                log_file_path = "Pipeline/output_deltaG/pyrosetta_" + pdb_file[:-4] + ".log"
#                print("pdb_file: ", pdb_file, " log_file: ", log_file_path)
                deltaG, residue_energies[ctr] = compute_deltaG_and_residue_energies(pdb, log_file_path=log_file_path, silent=True)

                # Save global ΔG
                writer.writerow([pair[0].split("/")[-2], pdb_file, deltaG,
                                 len(residue_energies[ctr]), len(set(r[1] for r in residue_energies[ctr]))])

                # Save residue-specific energies
                residue_file = os.path.join(output_dir, f"deltaG_{pdb_file[:-4]}.txt")
                print("Saving to residue_file: ", residue_file)
                with open(residue_file, "w") as residue_out:
                    residue_out.write("Residue_Name\tIndex\tEnergy\n")
                    for res_name, res_index, res_energy in residue_energies[ctr]:
                        residue_out.write(f"{res_name}\t{res_index}\t{res_energy:.3f}\n")

                # Compare residue energies and visualize differences
#                residue_energies[ctr] = compute_deltaG_and_residue_energies(pdb1_path, log_file_path=log_file_path, silent=True)[1]
                ctr += 1

            # Save visualization
            image_file = os.path.join(output_dir, f"deltaG_diff_{pdb1}_{pdb2}.jpg")
#            print("residue_energies_1:", residue_energies[0])
#            print("residue_energies_2:", residue_energies[1])

            align_and_compare_residues(residue_energies[0], residue_energies[1], os.path.basename(pdb1_path), os.path.basename(pdb2_path), image_file)
#            plot_residue_energy_differences(comparison, os.path.basename(pdb1), os.path.basename(pdb2), top_n=5, save_path=image_file)

    print(f"Global results saved to: {global_results_file}")
    print(f"Results for each PDB pair saved to {output_dir}")



def read_energy_tuples(file_path):
    """
    Reads a file and extracts energy values as a list of tuples (residue, index, energy).

    Args:
        file_path (str): Path to the file containing energies.

    Returns:
        list: A list of tuples, each containing (residue_name, residue_index, energy_value).
    """
    energy_tuples = []
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) >= 3:
                try:
                    residue = parts[0]
                    index = int(parts[1])  # Convert the residue index to integer
                    energy = float(parts[2])  # Convert the energy value to float
                    energy_tuples.append((residue, index, energy))
                except ValueError:
                    pass  # Ignore lines that don't have valid data
    return energy_tuples

# Supporting Functions
def compute_deltaG_and_residue_energies(
    pdb_path: str, log_file_path: str = None, silent: bool = True
):
    """
    Compute global ΔG and per-residue energies using PyRosetta, with options to suppress or redirect logs.

    Parameters:
    - pdb_path (str): Path to the PDB file.
    - log_file_path (str, optional): Path to save PyRosetta and custom log outputs.
    - silent (bool): If True, suppress all PyRosetta output completely.

    Returns:
    - tuple: (Global ΔG, Per-residue energies [(residue_name, residue_index, residue_energy)]).
    """
    # Redirect stdout and stderr if silent mode is enabled
    if silent:
        devnull = open(os.devnull, 'w')
        sys.stdout = devnull
        sys.stderr = devnull

    pyrosetta.init("-mute all")

    # Restore stdout and stderr after initialization
    if silent:
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__

    # Redirect custom logs to a file if log_file_path is provided
    log_file = open(log_file_path, "w") if log_file_path else None

    try:
        # Load PDB and compute energies
        if log_file:
            log_file.write(f"Processing file: {pdb_path}\n")
        else:
            print(f"Processing file: {pdb_path}")

        print(f"Processing file for ΔG pyrosseta: {pdb_path}")
        pose = pyrosetta.pose_from_file(pdb_path)
        scorefxn = pyrosetta.get_fa_scorefxn()

        deltaG = scorefxn(pose)

        # Get per-residue energies
        energies = pose.energies()
        residue_energies = [
            (pose.residue(i).name(), i, energies.residue_total_energy(i))
            for i in range(1, pose.total_residue() + 1)
        ]

        # Write results to log file if provided
        if log_file:
            log_file.write(f"Global ΔG: {deltaG}\n")
            log_file.write("Per-residue energies:\n")
            for residue_name, residue_index, residue_energy in residue_energies:
                log_file.write(f"{residue_name} {residue_index}: {residue_energy}\n")
        else:
            print(f"Global ΔG: {deltaG}")
            for residue_name, residue_index, residue_energy in residue_energies:
                print(f"{residue_name} {residue_index}: {residue_energy}")

        return deltaG, residue_energies

    finally:
        # Close log file if opened
        if log_file:
            log_file.close()


def plot_residue_energy_differences(residue_comparison, pdb_id_1, pdb_id_2, top_n=5, save_path=None):
    """
    Plot and optionally save per-residue energy differences between two structures.

    Parameters:
    - residue_comparison (list of tuples): [(residue_name, residue_index, energy1, energy2, difference)].
    - pdb_id_1: PDB ID of the first structure.
    - pdb_id_2: PDB ID of the second structure.
    - top_n (int): Number of largest differences to highlight.
    - save_path (str): File path to save the plot. If None, the plot is only displayed.
    """
    indices = [x[1] for x in residue_comparison]
    diffs = [x[4] for x in residue_comparison]
    labels = [x[0] for x in residue_comparison]

    sorted_by_diff = sorted(residue_comparison, key=lambda x: abs(x[4]), reverse=True)
    top_residues = sorted_by_diff[:top_n]

    plt.figure(figsize=(12, 6))
    bars = plt.bar(indices, diffs, color="blue", edgecolor="black", alpha=0.7)

    for residue in top_residues:
        index = residue[1]
        diff = residue[4]
        label = f"{residue[0]}{index}"

        bars[indices.index(index)].set_color("red")
        bars[indices.index(index)].set_alpha(0.9)
        plt.text(index, diff, label, fontsize=9, ha="center", va="bottom" if diff > 0 else "top")

    plt.axhline(0, color="black", linewidth=0.8, linestyle="--")
    plt.xlabel("Residue Index", fontsize=12)
#    plt.ylabel("Energy Difference (kcal/mol)", fontsize=12)
#    plt.title("Per-Residue Energy Differences", fontsize=14)
    plt.ylabel("ΔΔG (kcal/mol)")  # Use Greek delta for the axis label
    plt.title(f"Per-Residue ΔG({pdb_id_1})- ΔG({pdb_id_2})")
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path)
        print(f"Plot saved to: {save_path}")
    else:
        plt.show()


def align_and_compare_residues(residue_energies_1, residue_energies_2, pdb_id_1, pdb_id_2, output_file=None, top_n=5):
    """
    Align sequences from two PDB files and compare residue energies.

    Parameters:
    - residue_energies_1: List of tuples [(residue_name, residue_index, residue_energy)] for structure 1.
    - residue_energies_2: List of tuples [(residue_name, residue_index, residue_energy)] for structure 2.
    - pdb_id_1: PDB ID of the first structure.
    - pdb_id_2: PDB ID of the second structure.
    - output_file: File path to save the plot.
    - top_n: Number of residues with the largest differences to highlight.
    """
    # Clean the sequences
    seq1 = clean_sequence(residue_energies_1)
    seq2 = clean_sequence(residue_energies_2)

    # Perform sequence alignment
    from Bio import pairwise2
    alignments = pairwise2.align.globalxx(seq1, seq2)
    alignment = alignments[0]
    aligned_seq1, aligned_seq2 = alignment.seqA, alignment.seqB

    # Map energies to aligned positions
    aligned_energies_1, aligned_energies_2 = [], []
    index_1 = index_2 = 0
    for res1, res2 in zip(aligned_seq1, aligned_seq2):
        if res1 == "-" and res2 != "-":
            aligned_energies_1.append(None)
            aligned_energies_2.append(residue_energies_2[index_2][2])
            index_2 += 1
        elif res1 != "-" and res2 == "-":
            aligned_energies_1.append(residue_energies_1[index_1][2])
            aligned_energies_2.append(None)
            index_1 += 1
        elif res1 != "-" and res2 != "-":
            aligned_energies_1.append(residue_energies_1[index_1][2])
            aligned_energies_2.append(residue_energies_2[index_2][2])
            index_1 += 1
            index_2 += 1

    # Compute energy differences
    delta_energies = [
        (e1 - e2 if e1 is not None and e2 is not None else None)
        for e1, e2 in zip(aligned_energies_1, aligned_energies_2)
    ]

    # Find top N differences
    indexed_delta_energies = [(i, d) for i, d in enumerate(delta_energies) if d is not None]
    top_indices = sorted(indexed_delta_energies, key=lambda x: abs(x[1]), reverse=True)[:top_n]

    # Prepare for plotting
    indices = np.arange(len(delta_energies))
    delta_energies_filtered = [d if d is not None else 0 for d in delta_energies]
    if output_file is None:
        return delta_energies, delta_energies_filtered

    bar_colors = ["blue"] * len(delta_energies)
    for idx, _ in top_indices:
        bar_colors[idx] = "red"

    # Plot energy differences
    plt.figure(figsize=(12, 6))
    plt.bar(indices, delta_energies_filtered, color=bar_colors)

    # Annotate top residues
    for idx, delta in top_indices:
        residue_name = aa_short_long[aligned_seq1[idx]] if aligned_seq1[idx] != "-" else aligned_seq2[idx]
        plt.text(idx, delta + (1 if delta > 0 else -3), f"{residue_name}\n{idx + 1}",
                 color="red", ha="center", fontsize=8)


    plt.title(f"ΔG({pdb_id_1[:-4]}) - ΔG({pdb_id_2[:-4]})")
    plt.xlabel("Aligned Residue Index")
    plt.ylabel("ΔΔG (kcal/mol)")
    plt.axhline(0, color="black", linewidth=0.8)
    plt.ylim(min(delta_energies_filtered) - 5, max(delta_energies_filtered) + 5)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

    print(f"Plot saved to {output_file}")

    return delta_energies, delta_energies_filtered


def compute_deltaG_with_pyrosetta(pdb_path: str, log_file_path: str = None):
    """
    Compute the free energy (ΔG) of a protein structure using PyRosetta.

    Parameters:
    - pdb_path (str): Path to the PDB file.
    - log_file_path (str, optional): Path to a log file for PyRosetta output.
                                     If None, PyRosetta output is not redirected.

    Returns:
    - tuple: (PDB file path, ΔG, additional_info_dict)
    """
    if log_file_path:
        with open(log_file_path, "w") as log_file:
            with RedirectStdStreams(stdout=log_file, stderr=log_file):
                pyrosetta.init()
                return _compute_deltaG_internal(pdb_path)
    else:
        pyrosetta.init()
        return _compute_deltaG_internal(pdb_path)


def _compute_deltaG_internal(pdb_path: str):
    """
    Internal helper to compute ΔG without managing output redirection.

    Parameters:
    - pdb_path (str): Path to the PDB file.

    Returns:
    - tuple: (PDB file path, ΔG, additional_info_dict)
    """
    try:
        pose = pyrosetta.pose_from_file(pdb_path)
        scorefxn = pyrosetta.get_fa_scorefxn()
        deltaG = scorefxn(pose)

        # Additional information
        additional_info = {
            "residue_count": pose.total_residue(),
            "chain_count": pose.num_chains(),
        }

        return pdb_path, deltaG, additional_info
    except Exception as e:
        print(f"Error processing {pdb_path}: {e}")
        return pdb_path, None, None


def run_compute_deltaG_with_output(pdb_pair_files, foldpair_ids, output_file="/deltaG_pairs.csv"):
    """
    Compute ΔG for a list of PDB files and save results to a file.

    Parameters:
    - pdb_files (list): List of PDB file paths.
    - output_file (str): Path to save the results.
    """
    results = []

    for pdb_pair_file, foldpair_id in zip(pdb_pair_files, foldpair_ids):
        print(f"Computing ΔG for {pdb_pair_file}  pair using PyRosetta...")
        for fold in range(2):
            log_file_path = "Pipeline/output_deltaG/pyrosetta_" + foldpair_ids[fold] + ".log"
            print("Rosetta Log file: ", log_file_path)
            pdb_path, deltaG, additional_info = compute_deltaG_with_pyrosetta(pdb_pair_file[fold], log_file_path)

            if deltaG is not None:
                print(f"Computed ΔG for {pdb_pair_file[fold]}: {deltaG:.2f} kcal/mol")
                results.append({
                    "Pair": foldpair_id,
                    "PDB_ID": pdb_pair_file[fold][-8:-4],
                    "Delta_G": deltaG,
                    **additional_info  # Add additional outputs as columns
                })
            else:
                print(f"Failed to compute ΔG for {pdb_pair_file[fold]}.")
                results.append({
                    "Pair": foldpair_id,
                    "PDB_ID": pdb_pair_file[fold][-8:-4],
                    "Delta_G": "ERROR",
                    "residue_count": "N/A",
                    "chain_count": "N/A"
                })

    # Write results to a file
    with open(output_file, "w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=results[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(results)

    print(f"Results saved to {output_file}")
