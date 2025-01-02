import pyrosetta
import csv

def compute_deltaG_with_pyrosetta(pdb_path: str):
    """
    Compute the free energy (ΔG) of a protein structure using PyRosetta.

    Parameters:
    - pdb_path (str): Path to the PDB file.

    Returns:
    - tuple: (PDB file path, ΔG, additional_info_dict)
    """
    pyrosetta.init()

    try:
        pose = pyrosetta.pose_from_file(pdb_path)
        scorefxn = pyrosetta.get_fa_scorefxn()
        deltaG = scorefxn(pose)

        # Add any additional outputs if needed
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
            pdb_path, deltaG, additional_info = compute_deltaG_with_pyrosetta(pdb_pair_file[fold])

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
