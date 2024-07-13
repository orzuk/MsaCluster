import requests


def download_pdb(pdb_id):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)

    # Check if the request was successful (status code 200)
    if response.status_code == 200:
        filename = f"./pdb_files/{pdb_id}.pdb"
        with open(filename, 'w') as file:
            file.write(response.text)
        print(f"Downloaded {filename}")
    else:
        print(f"Failed to download PDB file for ID {pdb_id}. Status code: {response.status_code}")


# Example usage
pdb_id = "1A2B"  # Replace '1A2B' with your PDB ID
download_pdb(pdb_id)
