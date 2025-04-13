import requests

def download_pdb(pdb_id):
    """Download PDB and FASTA files from RCSB Protein Data Bank"""
    pdb_id = pdb_id.upper()
    print(f"Downloading {pdb_id} from RCSB PDB...")

    try:
        pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        pdb_response = requests.get(pdb_url)
        pdb_response.raise_for_status()  # Raises HTTPError for bad responses
        with open(f"{pdb_id}.pdb", "w") as pdb_file:
            pdb_file.write(pdb_response.text)
    except requests.RequestException as e:
        print(f"Error: Failed to download PDB file for {pdb_id}. Reason: {e}")
        return None

    try:
        fasta_url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
        fasta_response = requests.get(fasta_url)
        fasta_response.raise_for_status()
        with open(f"{pdb_id}.fasta", "w") as fasta_file:
            fasta_file.write(fasta_response.text)
    except requests.RequestException as e:
        print(f"Error: Failed to download FASTA file for {pdb_id}. Reason: {e}")
        return None

    print(f"Downloaded {pdb_id}.pdb and {pdb_id}.fasta successfully.")
    return pdb_id

def main():
    """Main function to download PDB and FASTA files for a list of PDB IDs"""
    pdb_ids = ["1V11", "4HHB", "6LU7"]  # Add more PDB IDs here
    for pdb_id in pdb_ids:
        download_pdb(pdb_id)

if __name__ == "__main__":
    main()