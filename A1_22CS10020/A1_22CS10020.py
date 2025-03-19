

import os
import sys
import requests
# source myenv/bin/activate
from Bio import SeqIO
from Bio.PDB import PDBParser, PPBuilder
import subprocess
import re
import argparse

# Amino acid molecular weights from Thermo Fischer Scientific
AA_WEIGHTS = {
    'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.15,
    'E': 147.13, 'Q': 146.15, 'G': 75.07, 'H': 155.16, 'I': 131.17,
    'L': 131.17, 'K': 146.19, 'M': 149.21, 'F': 165.19, 'P': 115.13,
    'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15
}

def download_pdb(pdb_id):
    """Download PDB and FASTA files from RCSB Protein Data Bank"""
    print(f"Downloading {pdb_id} from RCSB PDB...")
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    pdb_response = requests.get(pdb_url)
    if pdb_response.status_code != 200:
        print(f"Error: Failed to download PDB file for {pdb_id}. Status code: {pdb_response.status_code}")
        return None
    with open(f"{pdb_id}.pdb", "w") as pdb_file:
        pdb_file.write(pdb_response.text)
    fasta_url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
    fasta_response = requests.get(fasta_url)
    if fasta_response.status_code != 200:
        print(f"Error: Failed to download FASTA file for {pdb_id}. Status code: {fasta_response.status_code}")
        return None
    with open(f"{pdb_id}.fasta", "w") as fasta_file:
        fasta_file.write(fasta_response.text)
    print(f"Downloaded {pdb_id}.pdb and {pdb_id}.fasta successfully")
    return pdb_id

def extract_sequences(pdb_id):
    """Extract sequences from PDB and FASTA files"""
    print(f"Extracting sequences from {pdb_id}.pdb and {pdb_id}.fasta...")
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure(pdb_id, f"{pdb_id}.pdb")
    except Exception as e:
        print(f"Error parsing PDB file: {e}")
        return {}, {}
    ppb = PPBuilder()
    pdb_sequences = {}
    for model in structure:
        for chain in model:
            chain_id = chain.id
            if chain_id == " ":
                continue
            seq = ""
            for pp in ppb.build_peptides(chain):
                seq += str(pp.get_sequence())
            if seq and chain_id not in pdb_sequences:
                pdb_sequences[chain_id] = seq
    fasta_sequences = {}
    try:
        for record in SeqIO.parse(f"{pdb_id}.fasta", "fasta"):
            header = record.description
            chain_id = None
            if "|" in header:
                parts = header.split("|")
                if len(parts) >= 3:
                    chain_id = parts[2][0]
            elif "_" in header:
                parts = header.split("_")
                if len(parts) >= 2:
                    chain_id = parts[1][0]
            if not chain_id and "Chain " in header:
                match = re.search(r"Chain ([A-Za-z0-9])", header)
                if match:
                    chain_id = match.group(1)
            if not chain_id:
                match = re.search(rf"{pdb_id}[_|\s]([A-Za-z0-9])", header, re.IGNORECASE)
                if match:
                    chain_id = match.group(1)
            if chain_id and chain_id not in fasta_sequences:
                fasta_sequences[chain_id] = str(record.seq)
                print(f"Found FASTA sequence for chain {chain_id}: {len(str(record.seq))} residues")
    except Exception as e:
        print(f"Error parsing FASTA file: {e}")
    print(f"PDB chains found: {', '.join(pdb_sequences.keys())}")
    print(f"FASTA chains found: {', '.join(fasta_sequences.keys())}")
    return pdb_sequences, fasta_sequences

def identify_chain_breaks(pdb_sequences, fasta_sequences):
    """Identify breaks by comparing PDB and FASTA sequences"""
    print("Checking for chain breaks...")
    chain_breaks = {}
    for chain_id, pdb_seq in pdb_sequences.items():
        if chain_id not in fasta_sequences:
            print(f"Warning: Chain {chain_id} found in PDB but not in FASTA")
            continue
        fasta_seq = fasta_sequences[chain_id]
        breaks = []
        pdb_seq_str = str(pdb_seq)
        fasta_seq_str = str(fasta_seq)
        min_length = min(len(pdb_seq_str), len(fasta_seq_str))
        for i in range(min_length):
            if pdb_seq_str[i] != fasta_seq_str[i]:
                breaks.append((i+1, f"{fasta_seq_str[i]}→{pdb_seq_str[i]}"))
        if len(pdb_seq_str) != len(fasta_seq_str):
            if len(pdb_seq_str) < len(fasta_seq_str):
                breaks.append((len(pdb_seq_str)+1, "Missing residues in PDB"))
            else:
                breaks.append((len(fasta_seq_str)+1, "Extra residues in PDB"))
        if breaks:
            chain_breaks[chain_id] = breaks
    return chain_breaks

import os
import subprocess

def run_naccess(pdb_id):
    print("Running NACCESS for surface area calculations...")

    #directory for naccess and required files
    naccess_dir = os.path.expanduser(r"/mnt/c/ComputationalBiophysics/for_student")
    naccess_path = os.path.join(naccess_dir, "naccess")
    vdw_path = os.path.join(naccess_dir, "vdw.radii")
    std_path = os.path.join(naccess_dir, "standard.data")

    if not os.path.exists(naccess_path) or not os.access(naccess_path, os.X_OK):
        print("Error: NACCESS executable not found or not executable.")
        return False

    missing_files = []
    if not os.path.exists(vdw_path):
        missing_files.append("vdw.radii")
    if not os.path.exists(std_path):
        missing_files.append("standard.data")

    if missing_files:
        print(f"Error: Missing required NACCESS files: {', '.join(missing_files)}")
        return False

    print(f"Found NACCESS at: {naccess_path}")
    print(f"Found vdw.radii at: {vdw_path}")
    print(f"Found standard.data at: {std_path}")

    cmd = f'"{naccess_path}" "{pdb_id}.pdb" -r "{vdw_path}" -s "{std_path}"'

    try:
        print(f" Running NACCESS command: {cmd}")
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if result.returncode != 0:
            print(f" Warning: NACCESS failed with exit status {result.returncode}")
            print(f"Error output:\n{result.stderr}")
            return False

        print(" NACCESS completed successfully")
        return True

    except Exception as e:
        print(f" Error: Failed to execute NACCESS: {e}")
        return False


def compute_asa(pdb_id):
    """Extract ASA data from NACCESS output file (.rsa)"""
    asa_data = {}
    rsa_file_path = f"{pdb_id}.rsa"
    
    if os.path.exists(rsa_file_path):
        with open(rsa_file_path) as rsa_file:
            for line in rsa_file:
                if line.startswith("CHAIN"):
                    parts = line.split()
                    chain_id = parts[2]  
                    asa = float(parts[3]) 
                    asa_data[chain_id] = asa  
    
    return asa_data

def compute_molecular_weight(sequence):
    weight = 0.0
    for aa in sequence:
        if aa in AA_WEIGHTS:
            weight += AA_WEIGHTS[aa]
        else:
            print(f"Warning: Non-standard amino acid '{aa}' encountered. Using average weight.")
            weight += sum(AA_WEIGHTS.values()) / len(AA_WEIGHTS)
    weight += 18.02 
    return weight

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Analyze protein structures from PDB.')
    parser.add_argument('pdb_id', nargs='?', type=str, help='PDB ID (e.g., 1ABC)')
    parser.add_argument('--roll', type=str, help='Your roll number')
    return parser.parse_args()

def main():
    print("Protein Structure Analysis Program")
    print("Enter the PDB ID and roll number")
    print("Example: python3 A1_22CS10020.py 2UV0 --roll 22CS10020")
    default_pdb_id = input("Enter PDB ID (default: 2UV0): ") or "2UV0"
    default_roll_no = input("Enter Roll Number (default: 22CS10020): ") or "22CS10020"
    args = parse_arguments()
    pdb_id = args.pdb_id if args.pdb_id else default_pdb_id
    roll_no = args.roll if args.roll else default_roll_no
    print(f"Using PDB ID: {pdb_id}, Roll No: {roll_no}")
    downloaded_pdb_id = download_pdb(pdb_id)
    if not downloaded_pdb_id:
        print("Error: Could not download required files. Exiting.")
        return
    pdb_sequences, fasta_sequences = extract_sequences(pdb_id)
    chain_breaks = identify_chain_breaks(pdb_sequences, fasta_sequences)
    naccess_success = run_naccess(pdb_id)
    asa_data = {}
    if naccess_success:
        asa_data = compute_asa(pdb_id)
    output_filename = f"A1_{pdb_id}.txt"
    with open(output_filename, "w") as output_file:
        output_file.write(f"Analysis of PDB ID: {pdb_id}\n")
        output_file.write("=" * 50 + "\n\n")
        output_file.write(f"Number of chains identified: {len(pdb_sequences)}\n\n")
        output_file.write("Chain Information:\n")
        output_file.write("-" * 50 + "\n")
        for chain_id, sequence in pdb_sequences.items():
            num_aa = len(sequence)
            mol_weight = compute_molecular_weight(sequence)
            asa = asa_data.get(chain_id, 0.0)
            output_file.write(f"Chain {chain_id}: {num_aa} amino acids, Molecular Weight: {mol_weight:.2f} Da, ASA: {asa:.2f} Å²\n")
        output_file.write("\n")
        if chain_breaks:
            output_file.write("Chain breaks detected:\n")
            output_file.write("-" * 50 + "\n")
            for chain_id, breaks in chain_breaks.items():
                output_file.write(f"Chain {chain_id}: {len(breaks)} break(s)\n")
                for pos, desc in breaks:
                    output_file.write(f"  Position {pos}: {desc}\n")
        else:
            output_file.write("No chain breaks detected.\n")
    print(f"Analysis completed. Results written to {output_filename}")
    dirname = f"A1_{roll_no}"
    os.makedirs(dirname, exist_ok=True)
    for filename in [output_filename, f"{pdb_id}.pdb", f"{pdb_id}.fasta"]:
        if os.path.exists(filename):
            try:
                subprocess.run(["cp", filename, f"{dirname}/"], check=True)
            except Exception as e:
                print(f"Error copying {filename}: {e}")
                os.system(f"cp {filename} {dirname}/")
    readme_content = (
        "PDB Structure Analysis Program\n"
        "==============================\n\n"
        "This program analyzes protein structures downloaded from the RCSB Protein Data Bank.\n\n"
        "Usage:\n"
        "    python3 A1_22CS10020.py <pdb_id> --roll <roll_number>\n\n"
        "Example:\n"
        "    python3 A1_22CS10020.py 2UV0 --roll 22CS10020\n\n"
        "Requirements:\n"
        "    - Python 3.6+\n"
        "    - Biopython (pip install biopython)\n"
        "    - Requests (pip install requests)\n"
        "    - NACCESS (optional, for ASA calculations)\n\n"
        "Note:\n"
        "    Ensure that the NACCESS executable is located in one of the following paths:\n"
        "    - naccess\n"
        "    - /mnt/c/ComputationalBiophysics/for_student\n"
        "    - or any other directory specified in the candidates list in the script.\n"
    )
    with open(f"{dirname}/README.txt", "w") as readme_file:
        readme_file.write(readme_content)
    try:
        script_path = os.path.abspath(__file__)
        if os.path.exists(script_path):
            subprocess.run(["cp", script_path, f"{dirname}/A1_{roll_no}.py"], check=True)
        else:
            subprocess.run(["cp", "A1_22CS10020.py", f"{dirname}/A1_{roll_no}.py"], check=True)
    except Exception as e:
        print(f"Error copying script: {e}")
        os.system(f"cp A1_22CS10020.py {dirname}/A1_{roll_no}.py")
    try:
        subprocess.run(["zip", "-r", f"A1_{roll_no}.zip", dirname], check=True)
    except Exception as e:
        print(f"Error creating zip file: {e}")
        os.system(f"zip -r A1_{roll_no}.zip {dirname}")
    print(f"All files organized and zipped as A1_{roll_no}.zip")

if __name__ == "__main__":
    main()