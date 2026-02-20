try:
    import pyrosetta
    PYROSETTA_AVAILABLE = True
except ImportError:
    print("PyRosetta not installed; skipping energy/ΔG panels")
    PYROSETTA_AVAILABLE = False

import csv
import contextlib
import matplotlib.pyplot as plt

from config import *
import os
import sys
import numpy as np
import pandas as pd

from utils.protein_utils import *
from Bio.pairwise2 import format_alignment
from Bio import pairwise2


def clean_sequence(residue_energies):
    """
    Clean residue energies to extract a valid amino acid sequence.

    Accepts either:
      - list of dicts: [{'residue_name': 'ALA', 'residue_index': 3, 'energy': 2.5}, ...]
      - list of tuples: ('ALA', 3, 2.5)  # legacy

    Returns:
      str: one-letter AA sequence
    """
    seq_chars = []
    for res in residue_energies:
        if isinstance(res, dict):
            res3 = res.get("residue_name")
        else:  # tuple/list fallback
            res3 = res[0] if len(res) > 0 else None
        if not res3:
            continue
        # normalize nonstandard, then map to 1-letter (unknown → 'X')
        res3 = NONSTD_TO_STD.get(res3, str(res3).upper())
        seq_chars.append(aa_long_short.get(res3, "X"))
    return "".join(seq_chars)


def compute_global_and_residue_energies(pdb_pairs, foldpair_ids, output_dir, plot_dir: str | None = None):

    """
    Pipeline to compute global and residue-specific energies for PDB pairs.

    Parameters:
    - pdb_pairs (list of tuples): List of PDB file pairs [(pdb1, chain1), (pdb2, chain2)].
    - output_dir (str): Directory to save results.
    """
    if not PYROSETTA_AVAILABLE:
        raise RuntimeError("PyRosetta is not installed in this environment. "
                           "Install it or run this step on the Linux cluster.")
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # File to save global ΔG results

    pyrosetta.init() # ("-ignore_unrecognized_res")  # Init pyrosetta

    # At the beginning of the function, add these lines to load existing data:
    global_data = {}
    if os.path.exists(ENERGY_FILE):
        with open(ENERGY_FILE, "r", newline="") as global_file:
            reader = csv.DictReader(global_file, delimiter="\t")
            for row in reader:
                key = (row["Pair"], row["PDB_ID"])
                global_data[key] = row

#    with open(ENERGY_FILE, "w", newline="") as global_file:
#        writer = csv.writer(global_file, delimiter="\t")
#        writer.writerow(["Pair", "PDB_ID", "Delta_G", "Residue_Count", "Chain_Count"])

    print("All pdb-pairs: ", pdb_pairs)
#        for pdb_pair_file, foldpair_id in zip(pdb_pair_files, foldpair_ids):
    for pair in pdb_pairs:
        pdb1_path, pdb2_path = pair # pdb1, pdb2
        foldpair_id = pdb1_path[9:20] # remove pipeline, extract pair id pair[0] + "_" + pair[1]
        pair_anal_dir = os.path.join(DATA_DIR, foldpair_id, f"Analysis")
        # Load PDBs
        print(f"Processing pair: {pdb1_path} and {pdb2_path} ; foldpair_id: {foldpair_id}")

        # Compute global ΔG and residue-specific energies for both PDBs
        residue_energies = [[], []]
        ctr = 0
        for pdb in pair:
            pdb_file = os.path.basename(pdb) # .replace(".pdb", f"{chain}.pdb")
            deltaG, residue_energies[ctr], chain_count = compute_deltaG_and_residue_energies(pdb) # , log_file_path=log_file_path, silent=True)

            # Instead of writing directly, update the global_data dict:
            key = (pair[0].split("/")[-2], pdb_file)
            global_data[key] = {
                "Pair": pair[0].split("/")[-2],
                "PDB_ID": pdb_file,
                "Delta_G": deltaG,
                "Residue_Count": len(residue_energies[ctr]),
                "Chain_Count": chain_count,
            }

            # Save residue-specific energies
            residue_file = os.path.join(output_dir, f"deltaG_{pdb_file[:-4]}.txt")
            print("Saving to residue_file: ", residue_file)
            with open(residue_file, "w") as residue_out:
                residue_out.write("Residue_Name\tIndex\tEnergy\n")
#                    print("All residue energies for ctr: ", residue_energies[ctr])
                for res in residue_energies[ctr]:
#                        print("res: ", res)
                    try:
                        energy_val = float(res["energy"])
                        residue_out.write(f"{res['residue_name']}\t{res['residue_index']}\t{energy_val:.3f}\n")
                    except ValueError:
                        print(
                            f"Warning: Skipping residue {res['residue_index']} with non-numeric energy: {res['energy']}")

            # NEW: write per-structure per-residue CSV (pair-local)
            perres_csv = os.path.join(pair_anal_dir, f"df_energy_{pdb_file[:-4]}.csv")
            print("perres_csv: ", perres_csv)
            pd.DataFrame(residue_energies[ctr]).to_csv(perres_csv, index=False)

            # Compare residue energies and visualize differences
#                residue_energies[ctr] = compute_deltaG_and_residue_energies(pdb1_path, log_file_path=log_file_path, silent=True)[1]
            ctr += 1

        # NEW: pair-local global ΔG table (2 rows)
        pair_global_csv = os.path.join(pair_anal_dir, "df_energy_global.csv")
        # Find the two rows we just put into global_data (by keys)
        rows_local = []
        for pb in (pdb1_path, pdb2_path):
            k = (pair[0].split("/")[-2], os.path.basename(pb))
            rows_local.append(global_data[k])
        pd.DataFrame(rows_local).to_csv(pair_global_csv, index=False)


        # (optional) make the ΔΔG image only if plot_dir was requested
        if plot_dir:
            os.makedirs(plot_dir, exist_ok=True)
            image_file = os.path.join(plot_dir,
                                      f"{foldpair_ids[0] if len(foldpair_ids) == 1 else 'pair'}_ddg_aligned.png")
            print("plot_dir, image file: ", plot_dir, image_file)
        else:
            image_file = None

        deltas_raw, deltas_filled, aligned_seq1, aligned_seq2, aligned_e1, aligned_e2 = align_and_compare_residues(
            residue_energies[0], residue_energies[1],
            os.path.basename(pdb1_path), os.path.basename(pdb2_path), image_file)

        # Build aligned residue table (pair-local CSV in Analysis/)
        aligned_csv = os.path.join(pair_anal_dir, "df_ddg_aligned.csv")

        # Reconstruct aligned AA symbols for convenience (same logic as in align function)
        seq1 = clean_sequence(residue_energies[0])
        seq2 = clean_sequence(residue_energies[1])
        _aln = pairwise2.align.globalxx(seq1, seq2)[0]
        a1, a2 = _aln.seqA, _aln.seqB

        # Map original per-res energies to aligned positions (mirror logic)
        aligned_e1, aligned_e2 = [], []
        i1 = i2 = 0
        for r1, r2 in zip(a1, a2):
            if r1 == "-" and r2 != "-":
                aligned_e1.append(None); aligned_e2.append(residue_energies[1][i2]["energy"]); i2 += 1
            elif r1 != "-" and r2 == "-":
                aligned_e1.append(residue_energies[0][i1]["energy"]); aligned_e2.append(None); i1 += 1
            else:  # both not gap
                aligned_e1.append(residue_energies[0][i1]["energy"])
                aligned_e2.append(residue_energies[1][i2]["energy"])
                i1 += 1; i2 += 1

        pd.DataFrame({
            "aln_pos": list(range(1, len(deltas_filled)+1)),
            "aa1": list(a1),
            "aa2": list(a2),
            "E1": aligned_e1,
            "E2": aligned_e2,
            "dEdiff": deltas_raw,       # None where gaps; keep for strictness
            "dEdiff_filled": deltas_filled  # zeros at gaps (useful for plotting)
        }).to_csv(aligned_csv, index=False)


    # --- After processing all pairs, add these lines to write the updated data ---
    with open(ENERGY_FILE, "w", newline="") as global_file:
        fieldnames = ["Pair", "PDB_ID", "Delta_G", "Residue_Count", "Chain_Count"]
        writer = csv.DictWriter(global_file, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in global_data.values():
            writer.writerow(row)

    print(f"Global results saved to: {ENERGY_FILE}")
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


def remove_ssbond_records(pdb_in, pdb_out):
    """
    Reads a PDB file and writes a new version without SSBOND records.
    """
    with open(pdb_in, 'r') as fin:
        lines = fin.readlines()
    # Filter out SSBOND records
    new_lines = [line for line in lines if not line.startswith("SSBOND")]
    with open(pdb_out, 'w') as fout:
        fout.writelines(new_lines)
    return pdb_out


def compute_deltaG_and_residue_energies(pdb_path: str):
    # Use a full-atom scorefunction (Ref2015 weights by default in recent PyRosetta)
    scorefxn = pyrosetta.get_fa_scorefxn()
    try:
        # Attempt to load the PDB and score it normally
        pose = pyrosetta.pose_from_file(pdb_path)
        total_energy = scorefxn(pose)  # compute total energy for the pose
    except Exception as e:
        error_message = str(e)
        if "disulfide" in error_message:
            # **Disable disulfide detection** for the retry
            try:
                pyrosetta.init(options="-detect_disulf false")  # turn off automatic disulfide formation
            except RuntimeError:
                # PyRosetta may already be initialized; if so, we proceed without re-init
                pass
            # Remove GSH (or other problematic thiol ligand) from the PDB content
            with open(pdb_path, 'r') as f:
                pdb_lines = f.readlines()
            filtered_lines = []
            for line in pdb_lines:
                # Skip any HETATM lines for GSH residue
                if line.startswith("HETATM") and line[17:20].strip() == "GSH":
                    continue
                # Skip LINK/CONECT records involving GSH
                if (line.startswith("LINK") or line.startswith("CONECT")) and "GSH" in line:
                    continue
                filtered_lines.append(line)
            # Recreate the pose from the filtered PDB string (GSH removed)
            pdb_str = "".join(filtered_lines)
            pose = pyrosetta.Pose()
            pyrosetta.rosetta.core.import_pose.pose_from_pdbstring(pose, pdb_str)
            total_energy = scorefxn(pose)
        elif "Unrecognized residue: GTP" in error_message:
            print("GTP residue found. Removing GTP entries from the PDB and retrying...")
            with open(pdb_path, 'r') as f:
                pdb_lines = f.readlines()
            filtered_lines = []
            for line in pdb_lines:
                # Skip any HETATM lines for GTP residue
                if line.startswith("HETATM") and line[17:20].strip() == "GTP":
                    continue
                # Skip LINK/CONECT records involving GTP
                if (line.startswith("LINK") or line.startswith("CONECT")) and "GTP" in line:
                    continue
                filtered_lines.append(line)
            pdb_str = "".join(filtered_lines)
            pose = pyrosetta.Pose()
            pyrosetta.rosetta.core.import_pose.pose_from_pdbstring(pose, pdb_str)
            total_energy = scorefxn(pose)
        else:
            # If it's a different error, re-raise to avoid masking unexpected issues
            raise

    # Extract per-residue energies now that scoring is done
    energies = pose.energies()  # Energies object populated by scorefxn
    residue_energies = []
    chain_count = pose.num_chains()
    for i in range(1, pose.total_residue() + 1):
        res = pose.residue(i)
        res_energy = energies.residue_total_energy(i)
        residue_energies.append({
            "residue_name": res.name3(),      # 3-letter residue code (e.g., "ALA", "CYS")
            "residue_index": i,               # index in the pose
            "energy": res_energy
        })
    return total_energy, residue_energies, chain_count


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


def align_and_compare_residues(
    residue_energies_1, residue_energies_2,
    pdb_id_1, pdb_id_2,
    output_file=None, top_n=5
):
    """
    Align sequences from two PDB files and compare residue energies.

    Parameters:
      - residue_energies_1: list of dicts or tuples for structure 1
      - residue_energies_2: list of dicts or tuples for structure 2
      - pdb_id_1, pdb_id_2: labels for plotting
      - output_file: optional path to save plot
      - top_n: annotate N largest |ΔΔG|

    Returns:
      (delta_energies, delta_energies_filtered,
       aligned_seq1, aligned_seq2,
       aligned_energies_1, aligned_energies_2)
    """
    # 1) sequences
    seq1 = clean_sequence(residue_energies_1)
    seq2 = clean_sequence(residue_energies_2)

    # 2) global alignment
    alignments = pairwise2.align.globalxx(seq1, seq2)
    alignment = alignments[0]
    aligned_seq1, aligned_seq2 = alignment.seqA, alignment.seqB

    # 3) map energies along the alignment
    def _energy(lst, idx):
        x = lst[idx]
        return x["energy"] if isinstance(x, dict) else x[2]

    aligned_energies_1, aligned_energies_2 = [], []
    i1 = i2 = 0
    for r1, r2 in zip(aligned_seq1, aligned_seq2):
        if r1 == "-" and r2 != "-":
            aligned_energies_1.append(None)
            aligned_energies_2.append(_energy(residue_energies_2, i2))
            i2 += 1
        elif r1 != "-" and r2 == "-":
            aligned_energies_1.append(_energy(residue_energies_1, i1))
            aligned_energies_2.append(None)
            i1 += 1
        else:  # both residues present
            aligned_energies_1.append(_energy(residue_energies_1, i1))
            aligned_energies_2.append(_energy(residue_energies_2, i2))
            i1 += 1
            i2 += 1

    # 4) ΔΔG per aligned position
    delta_energies = [
        (e1 - e2) if (e1 is not None and e2 is not None) else None
        for e1, e2 in zip(aligned_energies_1, aligned_energies_2)
    ]
    indexed = [(i, d) for i, d in enumerate(delta_energies) if d is not None]
    top_indices = sorted(indexed, key=lambda x: abs(x[1]), reverse=True)[:top_n]
    delta_energies_filtered = [d if d is not None else 0 for d in delta_energies]

    # 5) optional plot
    if output_file is not None:
        plt.figure(figsize=(12, 6))
        colors = ["blue"] * len(delta_energies_filtered)
        for idx, _ in top_indices:
            colors[idx] = "red"
        plt.bar(range(len(delta_energies_filtered)), delta_energies_filtered, color=colors)
        for idx, delta in top_indices:
            aa = aligned_seq1[idx] if aligned_seq1[idx] != "-" else aligned_seq2[idx]
            plt.text(idx, delta + (1 if delta > 0 else -3), f"{aa}\n{idx+1}", ha="center", fontsize=8)
        plt.title(f"ΔG({pdb_id_1[:-4]}) - ΔG({pdb_id_2[:-4]})")
        plt.xlabel("Aligned Residue Index")
        plt.ylabel("ΔΔG (kcal/mol)")
        plt.axhline(0, color="black", linewidth=0.8)
        plt.tight_layout()
        plt.savefig(output_file)
        plt.close()
        print(f"Plot saved to {output_file}")

    return (delta_energies, delta_energies_filtered,
            aligned_seq1, aligned_seq2,
            aligned_energies_1, aligned_energies_2)



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
    if not PYROSETTA_AVAILABLE:
        raise RuntimeError("PyRosetta is not installed in this environment. "
                           "Install it or run this step on the Linux cluster.")
    if log_file_path:
        with open(log_file_path, "w") as log_file:
            with contextlib.redirect_stdout(log_file), contextlib.redirect_stderr(log_file):
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
            log_file_path = os.path.join(ENERGY_DIR, "pyrosetta_" + foldpair_ids[fold] + ".log")
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
