# # # import os
# # # import sys
# # # import requests
# # # from Bio import SeqIO
# # # from Bio.PDB import PDBParser, PPBuilder
# # # import subprocess
# # # import re
# # # import argparse

# # # # Amino acid molecular weights from Thermo Fischer Scientific
# # # AA_WEIGHTS = {
# # #     'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.15,
# # #     'E': 147.13, 'Q': 146.15, 'G': 75.07, 'H': 155.16, 'I': 131.17,
# # #     'L': 131.17, 'K': 146.19, 'M': 149.21, 'F': 165.19, 'P': 115.13,
# # #     'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15
# # # }

# # # def download_pdb(pdb_id):
# # #     """Download PDB and FASTA files from RCSB Protein Data Bank"""
# # #     print(f"Downloading {pdb_id} from RCSB PDB...")
    
# # #     # Download PDB file
# # #     pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
# # #     pdb_response = requests.get(pdb_url)
    
# # #     if pdb_response.status_code != 200:
# # #         print(f"Error: Failed to download PDB file for {pdb_id}. Status code: {pdb_response.status_code}")
# # #         return None
    
# # #     with open(f"{pdb_id}.pdb", "w") as pdb_file:
# # #         pdb_file.write(pdb_response.text)
    
# # #     # Download FASTA file
# # #     fasta_url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
# # #     fasta_response = requests.get(fasta_url)
    
# # #     if fasta_response.status_code != 200:
# # #         print(f"Error: Failed to download FASTA file for {pdb_id}. Status code: {fasta_response.status_code}")
# # #         return None
    
# # #     with open(f"{pdb_id}.fasta", "w") as fasta_file:
# # #         fasta_file.write(fasta_response.text)
    
# # #     print(f"Downloaded {pdb_id}.pdb and {pdb_id}.fasta successfully")
# # #     return pdb_id

# # # def extract_sequences(pdb_id):
# # #     """Extract sequences from PDB and FASTA files"""
# # #     print(f"Extracting sequences from {pdb_id}.pdb and {pdb_id}.fasta...")
    
# # #     # Parse PDB file
# # #     parser = PDBParser(QUIET=True)
# # #     try:
# # #         structure = parser.get_structure(pdb_id, f"{pdb_id}.pdb")
# # #     except Exception as e:
# # #         print(f"Error parsing PDB file: {e}")
# # #         return {}, {}
    
# # #     # Extract sequences from PDB file
# # #     ppb = PPBuilder()
# # #     pdb_sequences = {}
# # #     for model in structure:
# # #         for chain in model:
# # #             chain_id = chain.id
# # #             if chain_id == " ":  # Skip chains with no ID
# # #                 continue
            
# # #             seq = ""
# # #             for pp in ppb.build_peptides(chain):
# # #                 seq += str(pp.get_sequence())
            
# # #             if seq and chain_id not in pdb_sequences:
# # #                 pdb_sequences[chain_id] = seq
    
# # #     # Extract sequences from FASTA file
# # #     fasta_sequences = {}
# # #     try:
# # #         for record in SeqIO.parse(f"{pdb_id}.fasta", "fasta"):
# # #             # Extract chain ID from the FASTA header
# # #             header = record.description
# # #             chain_id = None
            
# # #             # Try different parsing patterns
# # #             # Common format: >pdb|1ABC|A Chain A, Some protein
# # #             # or >1ABC_A mol:protein length:153 Some protein Chain A
# # #             if "|" in header:
# # #                 parts = header.split("|")
# # #                 if len(parts) >= 3:
# # #                     chain_id = parts[2][0]  # Take first character after second |
# # #             elif "_" in header:
# # #                 parts = header.split("_")
# # #                 if len(parts) >= 2:
# # #                     chain_id = parts[1][0]  # Take first character after first _
            
# # #             # Try to extract from "Chain X" in description
# # #             if not chain_id and "Chain " in header:
# # #                 match = re.search(r"Chain ([A-Za-z0-9])", header)
# # #                 if match:
# # #                     chain_id = match.group(1)
            
# # #             # Last resort: try first character after PDB ID
# # #             if not chain_id:
# # #                 match = re.search(rf"{pdb_id}[_|\s]([A-Za-z0-9])", header, re.IGNORECASE)
# # #                 if match:
# # #                     chain_id = match.group(1)
            
# # #             # If we found a chain ID and it's not already in our dict
# # #             if chain_id and chain_id not in fasta_sequences:
# # #                 fasta_sequences[chain_id] = str(record.seq)
# # #                 print(f"Found FASTA sequence for chain {chain_id}: {len(str(record.seq))} residues")
# # #     except Exception as e:
# # #         print(f"Error parsing FASTA file: {e}")
    
# # #     # Debug output
# # #     print(f"PDB chains found: {', '.join(pdb_sequences.keys())}")
# # #     print(f"FASTA chains found: {', '.join(fasta_sequences.keys())}")
    
# # #     return pdb_sequences, fasta_sequences

# # # def identify_chain_breaks(pdb_sequences, fasta_sequences):
# # #     """Identify breaks by comparing PDB and FASTA sequences"""
# # #     print("Checking for chain breaks...")
    
# # #     chain_breaks = {}
    
# # #     for chain_id, pdb_seq in pdb_sequences.items():
# # #         if chain_id not in fasta_sequences:
# # #             print(f"Warning: Chain {chain_id} found in PDB but not in FASTA")
# # #             continue
        
# # #         fasta_seq = fasta_sequences[chain_id]
        
# # #         # Find the best alignment between sequences
# # #         # This is a simple implementation; for complex cases, use alignment algorithms
# # #         breaks = []
        
# # #         # Convert sequences to strings to ensure proper comparison
# # #         pdb_seq_str = str(pdb_seq)
# # #         fasta_seq_str = str(fasta_seq)
        
# # #         # Compare sequences
# # #         min_length = min(len(pdb_seq_str), len(fasta_seq_str))
# # #         for i in range(min_length):
# # #             if pdb_seq_str[i] != fasta_seq_str[i]:
# # #                 breaks.append((i+1, f"{fasta_seq_str[i]}→{pdb_seq_str[i]}"))
        
# # #         # Check for length differences
# # #         if len(pdb_seq_str) != len(fasta_seq_str):
# # #             if len(pdb_seq_str) < len(fasta_seq_str):
# # #                 breaks.append((len(pdb_seq_str)+1, "Missing residues in PDB"))
# # #             else:
# # #                 breaks.append((len(fasta_seq_str)+1, "Extra residues in PDB"))
        
# # #         if breaks:
# # #             chain_breaks[chain_id] = breaks
    
# # #     return chain_breaks

# # # # def run_naccess(pdb_id):
# # # #     """Run NACCESS to calculate solvent accessible surface area"""
# # # #     print("Running NACCESS for surface area calculations...")
    
# # # #     # Try to locate NACCESS executable
# # # #     naccess_locations = [
# # # #         "naccess",  # If in PATH
# # # #         "/usr/local/bin/naccess",
# # # #         "/usr/bin/naccess",
# # # #         "./naccess",
# # # #         "../naccess",
# # # #         os.path.expanduser("~/Downloads/naccess/naccess"),
# # # #         "/home/dadi/Downloads/naccess/naccess/naccess",  # From the provided code
# # # #     ]
    
# # # #     naccess_path = None
# # # #     for loc in naccess_locations:
# # # #         if os.path.exists(loc):
# # # #             naccess_path = loc
# # # #             break
# # # #         elif os.system(f"which {loc} > /dev/null 2>&1") == 0:
# # # #             # Get the path from which command
# # # #             try:
# # # #                 naccess_path = subprocess.check_output(["which", loc.split("/")[-1]]).decode().strip()
# # # #                 break
# # # #             except:
# # # #                 pass
    
# # # #     if not naccess_path:
# # # #         print("Warning: NACCESS executable not found. ASA calculations will be skipped.")
# # # #         return False
    
# # # #     print(f"Found NACCESS at: {naccess_path}")
    
# # # #     # Try to find standard data files
# # # #     vdw_locations = [
# # # #         "vdw.radii",
# # # #         "/usr/local/share/naccess/vdw.radii",
# # # #         "/usr/share/naccess/vdw.radii",
# # # #         "./vdw.radii",
# # # #         "../vdw.radii",
# # # #         os.path.expanduser("~/Downloads/naccess/vdw.radii"),
# # # #         "/home/dadi/Downloads/naccess/naccess/vdw.radii",  # From the provided code
# # # #     ]
    
# # # #     std_locations = [
# # # #         "standard.data",
# # # #         "/usr/local/share/naccess/standard.data",
# # # #         "/usr/share/naccess/standard.data",
# # # #         "./standard.data",
# # # #         "../standard.data",
# # # #         os.path.expanduser("~/Downloads/naccess/standard.data"),
# # # #         "/home/dadi/Downloads/naccess/naccess/standard.data",  # From the provided code
# # # #     ]
    
# # # #     # Look in same directory as naccess executable
# # # #     naccess_dir = os.path.dirname(naccess_path)
# # # #     vdw_locations.append(os.path.join(naccess_dir, "vdw.radii"))
# # # #     std_locations.append(os.path.join(naccess_dir, "standard.data"))
    
# # # #     vdw_path = None
# # # #     for loc in vdw_locations:
# # # #         if os.path.exists(loc):
# # # #             vdw_path = loc
# # # #             break
    
# # # #     std_path = None
# # # #     for loc in std_locations:
# # # #         if os.path.exists(loc):
# # # #             std_path = loc
# # # #             break
    
# # # #     if not vdw_path:
# # # #         print("Warning: vdw.radii file not found. NACCESS may not work correctly.")
# # # #     else:
# # # #         print(f"Found vdw.radii at: {vdw_path}")
        
# # # #     if not std_path:
# # # #         print("Warning: standard.data file not found. NACCESS may not work correctly.")
# # # #     else:
# # # #         print(f"Found standard.data at: {std_path}")
    
# # # #     # Check if accall is in the same directory as naccess
# # # #     accall_path = os.path.join(naccess_dir, "accall")
# # # #     if not os.path.exists(accall_path):
# # # #         print(f"Warning: 'accall' not found at {accall_path}. NACCESS may not work.")
    
# # # #     # Create the command
# # # #     cmd = [naccess_path, f"{pdb_id}.pdb"]
# # # #     if vdw_path:
# # # #         cmd.extend(["-r", vdw_path])
# # # #     if std_path:
# # # #         cmd.extend(["-s", std_path])
    
# # # #     # Add path to environment
# # # #     env = os.environ.copy()
# # # #     env["PATH"] = f"{naccess_dir}:{env.get('PATH', '')}"
    
# # # #     try:
# # # #         print(f"Running NACCESS command: {' '.join(cmd)}")
# # # #         process = subprocess.Popen(cmd, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
# # # #         stdout, stderr = process.communicate()
# # # #         if process.returncode != 0:
# # # #             print(f"Warning: NACCESS failed with error: {stderr}")
# # # #             return False
# # # #         print("NACCESS completed successfully")
# # # #         return True
# # # #     except Exception as e:
# # # #         print(f"Warning: NACCESS failed with error: {e}")
# # # #         return False


# # # # def run_naccess(pdb_id):
# # # #     """Run NACCESS to calculate solvent accessible surface area"""
# # # #     print("Running NACCESS for surface area calculations...")
    
# # # #     # Try to locate NACCESS executable
# # # #     naccess_locations = [
# # # #         "naccess",  # If in PATH
# # # #         "/usr/local/bin/naccess",
# # # #         "/usr/bin/naccess",
# # # #         "./naccess",
# # # #         "../naccess",
# # # #         os.path.expanduser("~/Downloads/naccess/naccess"),
# # # #         "/home/dadi/Downloads/naccess/naccess/naccess",  # From the provided code
# # # #     ]
    
# # # #     naccess_path = None
# # # #     for loc in naccess_locations:
# # # #         if os.path.exists(loc):
# # # #             naccess_path = loc
# # # #             break
# # # #         elif os.system(f"which {loc} > /dev/null 2>&1") == 0:
# # # #             # Get the path from which command
# # # #             try:
# # # #                 naccess_path = subprocess.check_output(["which", loc.split("/")[-1]]).decode().strip()
# # # #                 break
# # # #             except:
# # # #                 pass
    
# # # #     if not naccess_path:
# # # #         print("Warning: NACCESS executable not found. ASA calculations will be skipped.")
# # # #         return False
    
# # # #     print(f"Found NACCESS at: {naccess_path}")
    
# # # #     # Try to find standard data files
# # # #     vdw_locations = [
# # # #         "vdw.radii",
# # # #         "/usr/local/share/naccess/vdw.radii",
# # # #         "/usr/share/naccess/vdw.radii",
# # # #         "./vdw.radii",
# # # #         "../vdw.radii",
# # # #         os.path.expanduser("~/Downloads/naccess/vdw.radii"),
# # # #         "/home/dadi/Downloads/naccess/naccess/vdw.radii",  # From the provided code
# # # #     ]
    
# # # #     std_locations = [
# # # #         "standard.data",
# # # #         "/usr/local/share/naccess/standard.data",
# # # #         "/usr/share/naccess/standard.data",
# # # #         "./standard.data",
# # # #         "../standard.data",
# # # #         os.path.expanduser("~/Downloads/naccess/standard.data"),
# # # #         "/home/dadi/Downloads/naccess/naccess/standard.data",  # From the provided code
# # # #     ]
    
# # # #     # Look in same directory as naccess executable
# # # #     naccess_dir = os.path.dirname(naccess_path)
# # # #     vdw_locations.append(os.path.join(naccess_dir, "vdw.radii"))
# # # #     std_locations.append(os.path.join(naccess_dir, "standard.data"))
    
# # # #     vdw_path = None
# # # #     for loc in vdw_locations:
# # # #         if os.path.exists(loc):
# # # #             vdw_path = loc
# # # #             break
    
# # # #     std_path = None
# # # #     for loc in std_locations:
# # # #         if os.path.exists(loc):
# # # #             std_path = loc
# # # #             break
    
# # # #     if not vdw_path:
# # # #         print("Warning: vdw.radii file not found. NACCESS may not work correctly.")
# # # #     else:
# # # #         print(f"Found vdw.radii at: {vdw_path}")
        
# # # #     if not std_path:
# # # #         print("Warning: standard.data file not found. NACCESS may not work correctly.")
# # # #     else:
# # # #         print(f"Found standard.data at: {std_path}")
    
# # # #     # Check if accall is in the same directory as naccess
# # # #     accall_path = os.path.join(naccess_dir, "accall")
# # # #     if not os.path.exists(accall_path):
# # # #         print(f"Warning: 'accall' not found at {accall_path}. NACCESS may not work.")
# # # #     else:
# # # #         print(f"Found accall at: {accall_path}")
    
# # # #     # Create the command
# # # #     cmd = [naccess_path, f"{pdb_id}.pdb"]
# # # #     if vdw_path:
# # # #         cmd.extend(["-r", vdw_path])
# # # #     if std_path:
# # # #         cmd.extend(["-s", std_path])
    
# # # #     # Add path to environment
# # # #     env = os.environ.copy()
# # # #     env["PATH"] = f"{naccess_dir}:{env.get('PATH', '')}"
    
# # # #     try:
# # # #         print(f"Running NACCESS command: {' '.join(cmd)}")
# # # #         process = subprocess.Popen(cmd, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
# # # #         stdout, stderr = process.communicate()
# # # #         if process.returncode != 0:
# # # #             print(f"Warning: NACCESS failed with error: {stderr}")
# # # #             return False
# # # #         print("NACCESS completed successfully")
# # # #         return True
# # # #     except Exception as e:
# # # #         print(f"Warning: NACCESS failed with error: {e}")
# # # #         return False


# # # def run_naccess(pdb_id):
# # #     """Run NACCESS to calculate solvent accessible surface area"""
# # #     print("Running NACCESS for surface area calculations...")
    
# # #     # Try to locate NACCESS executable
# # #     naccess_locations = [
# # #         "naccess",  # If in PATH
# # #         "/usr/local/bin/naccess",
# # #         "/usr/bin/naccess",
# # #         "./naccess",
# # #         "../naccess",
# # #         os.path.expanduser("~/Downloads/naccess/naccess"),
# # #         "/home/dadi/Downloads/naccess/naccess/naccess",  # From the provided code
# # #     ]
    
# # #     naccess_path = None
# # #     for loc in naccess_locations:
# # #         if os.path.exists(loc):
# # #             naccess_path = loc
# # #             break
# # #         elif os.system(f"which {loc} > /dev/null 2>&1") == 0:
# # #             # Get the path from which command
# # #             try:
# # #                 naccess_path = subprocess.check_output(["which", loc.split("/")[-1]]).decode().strip()
# # #                 break
# # #             except:
# # #                 pass
    
# # #     if not naccess_path:
# # #         print("Warning: NACCESS executable not found. ASA calculations will be skipped.")
# # #         return False
    
# # #     print(f"Found NACCESS at: {naccess_path}")
    
# # #     # Try to find standard data files
# # #     vdw_locations = [
# # #         "vdw.radii",
# # #         "/usr/local/share/naccess/vdw.radii",
# # #         "/usr/share/naccess/vdw.radii",
# # #         "./vdw.radii",
# # #         "../vdw.radii",
# # #         os.path.expanduser("~/Downloads/naccess/vdw.radii"),
# # #         "/home/dadi/Downloads/naccess/naccess/vdw.radii",  # From the provided code
# # #     ]
    
# # #     std_locations = [
# # #         "standard.data",
# # #         "/usr/local/share/naccess/standard.data",
# # #         "/usr/share/naccess/standard.data",
# # #         "./standard.data",
# # #         "../standard.data",
# # #         os.path.expanduser("~/Downloads/naccess/standard.data"),
# # #         "/home/dadi/Downloads/naccess/naccess/standard.data",  # From the provided code
# # #     ]
    
# # #     # Look in same directory as naccess executable
# # #     naccess_dir = os.path.dirname(naccess_path)
# # #     vdw_locations.append(os.path.join(naccess_dir, "vdw.radii"))
# # #     std_locations.append(os.path.join(naccess_dir, "standard.data"))
    
# # #     vdw_path = None
# # #     for loc in vdw_locations:
# # #         if os.path.exists(loc):
# # #             vdw_path = loc
# # #             break
    
# # #     std_path = None
# # #     for loc in std_locations:
# # #         if os.path.exists(loc):
# # #             std_path = loc
# # #             break
    
# # #     if not vdw_path:
# # #         print("Warning: vdw.radii file not found. NACCESS may not work correctly.")
# # #     else:
# # #         print(f"Found vdw.radii at: {vdw_path}")
        
# # #     if not std_path:
# # #         print("Warning: standard.data file not found. NACCESS may not work correctly.")
# # #     else:
# # #         print(f"Found standard.data at: {std_path}")
    
# # #     # Check if accall is in the same directory as naccess
# # #     accall_path = os.path.join(naccess_dir, "accall")
# # #     if not os.path.exists(accall_path):
# # #         print(f"Warning: 'accall' not found at {accall_path}. NACCESS may not work.")
# # #     else:
# # #         print(f"Found accall at: {accall_path}")
    
# # #     # Create the command
# # #     cmd = f"{naccess_path} {pdb_id}.pdb"
# # #     if vdw_path:
# # #         cmd += f" -r {vdw_path}"
# # #     if std_path:
# # #         cmd += f" -s {std_path}"
    
# # #     # Add path to environment
# # #     env = os.environ.copy()
# # #     env["PATH"] = f"{naccess_dir}:{env.get('PATH', '')}"
    
# # #     try:
# # #         print(f"Running NACCESS command: {cmd}")
# # #         result = subprocess.run(cmd, env=env, shell=True, check=True, capture_output=True, text=True)
# # #         if result.returncode != 0:
# # #             print(f"Warning: NACCESS failed with error: {result.stderr}")
# # #             return False
# # #         print("NACCESS completed successfully")
# # #         return True
# # #     except Exception as e:
# # #         print(f"Warning: NACCESS failed with error: {e}")
# # #         return False

# # # def compute_asa(pdb_id):
# # #     """Extract ASA data from NACCESS output"""
# # #     asa_data = {}
    
# # #     # Try different NACCESS output file extensions
# # #     for ext in ['.asa', '.rsa']:
# # #         asa_file_path = f"{pdb_id}{ext}"
# # #         if os.path.exists(asa_file_path):
# # #             try:
# # #                 with open(asa_file_path) as asa_file:
# # #                     for line in asa_file:
# # #                         if line.startswith("RES") or line.startswith("CHAIN"):
# # #                             parts = line.split()
# # #                             if len(parts) >= 5:
# # #                                 chain_id = parts[2]
# # #                                 # Handle both total and per-residue ASA
# # #                                 if parts[0] == "CHAIN":
# # #                                     asa = float(parts[4])
# # #                                     asa_data[chain_id] = asa
# # #                                 elif chain_id not in asa_data:
# # #                                     asa = float(parts[4])
# # #                                     asa_data[chain_id] = asa
# # #                                 else:
# # #                                     asa_data[chain_id] += float(parts[4])
# # #             except Exception as e:
# # #                 print(f"Warning: Error reading ASA file: {e}")
            
# # #             break  # Stop if we found and processed a file
    
# # #     return asa_data

# # # def compute_molecular_weight(sequence):
# # #     """Calculate molecular weight of a protein sequence"""
# # #     weight = 0.0
# # #     for aa in sequence:
# # #         if aa in AA_WEIGHTS:
# # #             weight += AA_WEIGHTS[aa]
# # #         else:
# # #             # Handle non-standard amino acids
# # #             print(f"Warning: Non-standard amino acid '{aa}' encountered. Using average weight.")
# # #             weight += sum(AA_WEIGHTS.values()) / len(AA_WEIGHTS)
    
# # #     # Add weight of water (18.02 Da) to account for the terminal -OH and H-
# # #     weight += 18.02
    
# # #     return weight

# # # def parse_arguments():
# # #     """Parse command line arguments"""
# # #     parser = argparse.ArgumentParser(description='Analyze protein structures from PDB.')
# # #     parser.add_argument('pdb_id', nargs='?', type=str, help='PDB ID (e.g., 2UV0)')
# # #     parser.add_argument('--roll', type=str, help='Your roll number')
# # #     return parser.parse_args()

# # # def main():
# # #     # Hardcoded values - use these if no command line arguments are provided
# # #     default_pdb_id = "2UV0"
# # #     default_roll_no = "22CS10020"  # Replace with your actual roll number
    
# # #     # Parse command line arguments (optional)
# # #     args = parse_arguments()
    
# # #     # Use command line arguments if provided, otherwise use defaults
# # #     pdb_id = args.pdb_id if args.pdb_id else default_pdb_id
# # #     roll_no = args.roll if args.roll else default_roll_no
    
# # #     print(f"Using PDB ID: {pdb_id}, Roll No: {roll_no}")
    
# # #     # Download files
# # #     downloaded_pdb_id = download_pdb(pdb_id)
# # #     if not downloaded_pdb_id:
# # #         print("Error: Could not download required files. Exiting.")
# # #         return
    
# # #     # Extract sequences
# # #     pdb_sequences, fasta_sequences = extract_sequences(pdb_id)
    
# # #     # Identify chain breaks
# # #     chain_breaks = identify_chain_breaks(pdb_sequences, fasta_sequences)
    
# # #     # Run NACCESS
# # #     naccess_success = run_naccess(pdb_id)
    
# # #     # Compute ASA if NACCESS succeeded
# # #     asa_data = {}
# # #     if naccess_success:
# # #         asa_data = compute_asa(pdb_id)
    
# # #     # Write output file
# # #     output_filename = f"A1_{pdb_id}.txt"
# # #     with open(output_filename, "w") as output_file:
# # #         # Write header
# # #         output_file.write(f"Analysis of PDB ID: {pdb_id}\n")
# # #         output_file.write("=" * 50 + "\n\n")
        
# # #         # Number of chains
# # #         output_file.write(f"Number of chains identified: {len(pdb_sequences)}\n\n")
        
# # #         # Chain information
# # #         output_file.write("Chain Information:\n")
# # #         output_file.write("-" * 50 + "\n")
# # #         for chain_id, sequence in pdb_sequences.items():
# # #             num_aa = len(sequence)
# # #             mol_weight = compute_molecular_weight(sequence)
# # #             asa = asa_data.get(chain_id, 0.0)
            
# # #             output_file.write(f"Chain {chain_id}: {num_aa} amino acids, Molecular Weight: {mol_weight:.2f} Da, ASA: {asa:.2f} Å²\n")
        
# # #         output_file.write("\n")
        
# # #         # Chain breaks
# # #         if chain_breaks:
# # #             output_file.write("Chain breaks detected:\n")
# # #             output_file.write("-" * 50 + "\n")
# # #             for chain_id, breaks in chain_breaks.items():
# # #                 output_file.write(f"Chain {chain_id}: {len(breaks)} break(s)\n")
# # #                 for pos, desc in breaks:
# # #                     output_file.write(f"  Position {pos}: {desc}\n")
# # #         else:
# # #             output_file.write("No chain breaks detected.\n")
    
# # #     print(f"Analysis completed. Results written to {output_filename}")
    
# # #     # Create directory and organize files
# # #     dirname = f"A1{roll_no}"
# # #     os.makedirs(dirname, exist_ok=True)
    
# # #     # Copy files to directory
# # #     for filename in [output_filename, f"{pdb_id}.pdb", f"{pdb_id}.fasta"]:
# # #         if os.path.exists(filename):
# # #             try:
# # #                 subprocess.run(["cp", filename, f"{dirname}/"], check=True)
# # #             except Exception as e:
# # #                 print(f"Error copying {filename}: {e}")
# # #                 # Fallback to os.system
# # #                 os.system(f"cp {filename} {dirname}/")
    
# # #     # Create README
# # #     with open(f"{dirname}/README.txt", "w") as readme_file:
# # #         readme_file.write("PDB Structure Analysis Program\n")
# # #         readme_file.write("=" * 30 + "\n\n")
# # #         readme_file.write("This program analyzes protein structures from the RCSB Protein Data Bank.\n\n")
# # #         readme_file.write("To run the program with default values (PDB ID: 2UV0):\n")
# # #         readme_file.write("python3 main.py\n\n")
# # #         readme_file.write("To run with custom PDB ID:\n")
# # #         readme_file.write("python3 main.py <pdb_id> --roll <your_roll_number>\n\n")
# # #         readme_file.write("Example:\n")
# # #         readme_file.write(f"python3 main.py 2UV0 --roll {roll_no}\n\n")
# # #         readme_file.write("Requirements:\n")
# # #         readme_file.write("- Python 3.6+\n")
# # #         readme_file.write("- Biopython\n")
# # #         readme_file.write("- Requests\n")
# # #         readme_file.write("- NACCESS (optional, for ASA calculations)\n\n")
# # #         readme_file.write("Install dependencies with:\n")
# # #         readme_file.write("pip install biopython requests\n")
    
# # #     # Copy the python script to the directory with the correct name
# # #     try:
# # #         script_path = os.path.abspath(__file__)
# # #         if os.path.exists(script_path):
# # #             # Copy with roll number in filename
# # #             subprocess.run(["cp", script_path, f"{dirname}/A1_{roll_no}.py"], check=True)
# # #             # Also copy as main.py for convenience
# # #             subprocess.run(["cp", script_path, f"{dirname}/main.py"], check=True)
# # #         else:
# # #             # If __file__ doesn't work, try to copy the current script
# # #             subprocess.run(["cp", "main.py", f"{dirname}/A1_{roll_no}.py"], check=True)
# # #             subprocess.run(["cp", "main.py", f"{dirname}/main.py"], check=True)
# # #     except Exception as e:
# # #         print(f"Error copying script: {e}")
# # #         # Fallback with os.system
# # #         os.system(f"cp main.py {dirname}/A1_{roll_no}.py")
# # #         os.system(f"cp main.py {dirname}/main.py")
    
# # #     # Create zip file
# # #     try:
# # #         subprocess.run(["zip", "-r", f"A1{roll_no}.zip", dirname], check=True)
# # #     except Exception as e:
# # #         print(f"Error creating zip file: {e}")
# # #         # Fallback to os.system
# # #         os.system(f"zip -r A1{roll_no}.zip {dirname}")
    
# # #     print(f"All files organized and zipped as A1{roll_no}.zip")

# # # if __name__ == "__main__":
# # #     main()



# # import os
# # import sys
# # import requests
# # from Bio import SeqIO
# # from Bio.PDB import PDBParser, PPBuilder
# # import subprocess
# # import re
# # import argparse

# # # Amino acid molecular weights from Thermo Fischer Scientific
# # AA_WEIGHTS = {
# #     'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.15,
# #     'E': 147.13, 'Q': 146.15, 'G': 75.07, 'H': 155.16, 'I': 131.17,
# #     'L': 131.17, 'K': 146.19, 'M': 149.21, 'F': 165.19, 'P': 115.13,
# #     'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15
# # }

# # def download_pdb(pdb_id):
# #     """Download PDB and FASTA files from RCSB Protein Data Bank"""
# #     print(f"Downloading {pdb_id} from RCSB PDB...")
    
# #     # Download PDB file
# #     pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
# #     pdb_response = requests.get(pdb_url)
    
# #     if pdb_response.status_code != 200:
# #         print(f"Error: Failed to download PDB file for {pdb_id}. Status code: {pdb_response.status_code}")
# #         return None
    
# #     with open(f"{pdb_id}.pdb", "w") as pdb_file:
# #         pdb_file.write(pdb_response.text)
    
# #     # Download FASTA file
# #     fasta_url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
# #     fasta_response = requests.get(fasta_url)
    
# #     if fasta_response.status_code != 200:
# #         print(f"Error: Failed to download FASTA file for {pdb_id}. Status code: {fasta_response.status_code}")
# #         return None
    
# #     with open(f"{pdb_id}.fasta", "w") as fasta_file:
# #         fasta_file.write(fasta_response.text)
    
# #     print(f"Downloaded {pdb_id}.pdb and {pdb_id}.fasta successfully")
# #     return pdb_id

# # def extract_sequences(pdb_id):
# #     """Extract sequences from PDB and FASTA files"""
# #     print(f"Extracting sequences from {pdb_id}.pdb and {pdb_id}.fasta...")
    
# #     # Parse PDB file
# #     parser = PDBParser(QUIET=True)
# #     try:
# #         structure = parser.get_structure(pdb_id, f"{pdb_id}.pdb")
# #     except Exception as e:
# #         print(f"Error parsing PDB file: {e}")
# #         return {}, {}
    
# #     # Extract sequences from PDB file
# #     ppb = PPBuilder()
# #     pdb_sequences = {}
# #     for model in structure:
# #         for chain in model:
# #             chain_id = chain.id
# #             if chain_id == " ":  # Skip chains with no ID
# #                 continue
            
# #             seq = ""
# #             for pp in ppb.build_peptides(chain):
# #                 seq += str(pp.get_sequence())
            
# #             if seq and chain_id not in pdb_sequences:
# #                 pdb_sequences[chain_id] = seq
    
# #     # Extract sequences from FASTA file
# #     fasta_sequences = {}
# #     try:
# #         for record in SeqIO.parse(f"{pdb_id}.fasta", "fasta"):
# #             header = record.description
# #             chain_id = None
            
# #             # Try different parsing patterns
# #             if "|" in header:
# #                 parts = header.split("|")
# #                 if len(parts) >= 3:
# #                     chain_id = parts[2][0]  # First character after second |
# #             elif "_" in header:
# #                 parts = header.split("_")
# #                 if len(parts) >= 2:
# #                     chain_id = parts[1][0]  # First character after first _
            
# #             # Try to extract from "Chain X" in description
# #             if not chain_id and "Chain " in header:
# #                 match = re.search(r"Chain ([A-Za-z0-9])", header)
# #                 if match:
# #                     chain_id = match.group(1)
            
# #             # Last resort: try first character after PDB ID in header
# #             if not chain_id:
# #                 match = re.search(rf"{pdb_id}[_|\s]([A-Za-z0-9])", header, re.IGNORECASE)
# #                 if match:
# #                     chain_id = match.group(1)
            
# #             # If we found a chain ID that isn't in our dict, add it.
# #             if chain_id and chain_id not in fasta_sequences:
# #                 fasta_sequences[chain_id] = str(record.seq)
# #                 print(f"Found FASTA sequence for chain {chain_id}: {len(str(record.seq))} residues")
# #     except Exception as e:
# #         print(f"Error parsing FASTA file: {e}")
    
# #     # Debug output for chain IDs
# #     print(f"PDB chains found: {', '.join(pdb_sequences.keys())}")
# #     print(f"FASTA chains found: {', '.join(fasta_sequences.keys())}")
    
# #     return pdb_sequences, fasta_sequences

# # def identify_chain_breaks(pdb_sequences, fasta_sequences):
# #     """Identify breaks by comparing PDB and FASTA sequences"""
# #     print("Checking for chain breaks...")
    
# #     chain_breaks = {}
    
# #     for chain_id, pdb_seq in pdb_sequences.items():
# #         if chain_id not in fasta_sequences:
# #             print(f"Warning: Chain {chain_id} found in PDB but not in FASTA")
# #             continue
        
# #         fasta_seq = fasta_sequences[chain_id]
# #         breaks = []
        
# #         # Convert sequences to strings
# #         pdb_seq_str = str(pdb_seq)
# #         fasta_seq_str = str(fasta_seq)
        
# #         # Compare each residue
# #         min_length = min(len(pdb_seq_str), len(fasta_seq_str))
# #         for i in range(min_length):
# #             if pdb_seq_str[i] != fasta_seq_str[i]:
# #                 breaks.append((i+1, f"{fasta_seq_str[i]}→{pdb_seq_str[i]}"))
        
# #         # Check for length differences
# #         if len(pdb_seq_str) != len(fasta_seq_str):
# #             if len(pdb_seq_str) < len(fasta_seq_str):
# #                 breaks.append((len(pdb_seq_str)+1, "Missing residues in PDB"))
# #             else:
# #                 breaks.append((len(fasta_seq_str)+1, "Extra residues in PDB"))
        
# #         if breaks:
# #             chain_breaks[chain_id] = breaks
    
# #     return chain_breaks

# # def run_naccess(pdb_id):
# #     """Run NACCESS to calculate solvent accessible surface area"""
# #     print("Running NACCESS for surface area calculations...")
    
# #     # Try to locate NACCESS executable
# #     naccess_locations = [
# #         "naccess",  # If in PATH
# #         "/usr/local/bin/naccess",
# #         "/usr/bin/naccess",
# #         "./naccess",
# #         "../naccess",
# #         os.path.expanduser("~/Downloads/naccess/naccess"),
# #         "/home/dadi/Downloads/naccess/naccess/naccess",  # Example location
# #     ]
    
# #     naccess_path = None
# #     for loc in naccess_locations:
# #         if os.path.exists(loc):
# #             naccess_path = loc
# #             break
# #         elif os.system(f"which {loc} > /dev/null 2>&1") == 0:
# #             try:
# #                 naccess_path = subprocess.check_output(["which", loc.split("/")[-1]]).decode().strip()
# #                 break
# #             except Exception:
# #                 pass
    
# #     if not naccess_path:
# #         print("Warning: NACCESS executable not found. ASA calculations will be skipped.")
# #         return False
    
# #     print(f"Found NACCESS at: {naccess_path}")
    
# #     # Try to find standard data files for NACCESS
# #     vdw_locations = [
# #         "vdw.radii",
# #         "/usr/local/share/naccess/vdw.radii",
# #         "/usr/share/naccess/vdw.radii",
# #         "./vdw.radii",
# #         "../vdw.radii",
# #         os.path.expanduser("~/Downloads/naccess/vdw.radii"),
# #         "/home/dadi/Downloads/naccess/naccess/vdw.radii",
# #     ]
    
# #     std_locations = [
# #         "standard.data",
# #         "/usr/local/share/naccess/standard.data",
# #         "/usr/share/naccess/standard.data",
# #         "./standard.data",
# #         "../standard.data",
# #         os.path.expanduser("~/Downloads/naccess/standard.data"),
# #         "/home/dadi/Downloads/naccess/naccess/standard.data",
# #     ]
    
# #     # Look in same directory as NACCESS executable
# #     naccess_dir = os.path.dirname(naccess_path)
# #     vdw_locations.append(os.path.join(naccess_dir, "vdw.radii"))
# #     std_locations.append(os.path.join(naccess_dir, "standard.data"))
    
# #     vdw_path = None
# #     for loc in vdw_locations:
# #         if os.path.exists(loc):
# #             vdw_path = loc
# #             break
    
# #     std_path = None
# #     for loc in std_locations:
# #         if os.path.exists(loc):
# #             std_path = loc
# #             break
    
# #     if not vdw_path:
# #         print("Warning: vdw.radii file not found. NACCESS may not work correctly.")
# #     else:
# #         print(f"Found vdw.radii at: {vdw_path}")
        
# #     if not std_path:
# #         print("Warning: standard.data file not found. NACCESS may not work correctly.")
# #     else:
# #         print(f"Found standard.data at: {std_path}")
    
# #     # Check if accall exists in same directory
# #     accall_path = os.path.join(naccess_dir, "accall")
# #     if not os.path.exists(accall_path):
# #         print(f"Warning: 'accall' not found at {accall_path}. NACCESS may not work.")
# #     else:
# #         print(f"Found accall at: {accall_path}")
    
# #     # Build the command for NACCESS
# #     cmd = f"{naccess_path} {pdb_id}.pdb"
# #     if vdw_path:
# #         cmd += f" -r {vdw_path}"
# #     if std_path:
# #         cmd += f" -s {std_path}"
    
# #     # Add NACCESS directory to PATH in environment
# #     env = os.environ.copy()
# #     env["PATH"] = f"{naccess_dir}:{env.get('PATH', '')}"
    
# #     try:
# #         print(f"Running NACCESS command: {cmd}")
# #         result = subprocess.run(cmd, env=env, shell=True, check=True, capture_output=True, text=True)
# #         if result.returncode != 0:
# #             print(f"Warning: NACCESS failed with error: {result.stderr}")
# #             return False
# #         print("NACCESS completed successfully")
# #         return True
# #     except Exception as e:
# #         print(f"Warning: NACCESS failed with error: {e}")
# #         return False

# # def compute_asa(pdb_id):
# #     """Extract ASA data from NACCESS output files (.asa or .rsa)"""
# #     asa_data = {}
    
# #     for ext in ['.asa', '.rsa']:
# #         asa_file_path = f"{pdb_id}{ext}"
# #         if os.path.exists(asa_file_path):
# #             try:
# #                 with open(asa_file_path) as asa_file:
# #                     for line in asa_file:
# #                         if line.startswith("RES") or line.startswith("CHAIN"):
# #                             parts = line.split()
# #                             if len(parts) >= 5:
# #                                 chain_id = parts[2]
# #                                 # Depending on file type, sum the ASA data per chain
# #                                 if parts[0] == "CHAIN":
# #                                     asa = float(parts[4])
# #                                     asa_data[chain_id] = asa
# #                                 elif chain_id not in asa_data:
# #                                     asa = float(parts[4])
# #                                     asa_data[chain_id] = asa
# #                                 else:
# #                                     asa_data[chain_id] += float(parts[4])
# #             except Exception as e:
# #                 print(f"Warning: Error reading ASA file: {e}")
# #             break  
# #     return asa_data

# # def compute_molecular_weight(sequence):
# #     """Calculate molecular weight of a protein sequence"""
# #     weight = 0.0
# #     for aa in sequence:
# #         if aa in AA_WEIGHTS:
# #             weight += AA_WEIGHTS[aa]
# #         else:
# #             print(f"Warning: Non-standard amino acid '{aa}' encountered. Using average weight.")
# #             weight += sum(AA_WEIGHTS.values()) / len(AA_WEIGHTS)
    
# #     # Add weight of water (18.02 Da) for the terminal groups
# #     weight += 18.02
# #     return weight

# # def parse_arguments():
# #     """Parse command line arguments"""
# #     parser = argparse.ArgumentParser(description='Analyze protein structures from PDB.')
# #     parser.add_argument('pdb_id', nargs='?', type=str, help='PDB ID (e.g., 1ABC)')
# #     parser.add_argument('--roll', type=str, help='Your roll number')
# #     return parser.parse_args()

# # def main():
# #     # Default values in case no command line arguments are provided
# #     default_pdb_id = "2UV0"
# #     default_roll_no = "22CS10020"  # Replace with your actual roll number
    
# #     args = parse_arguments()
# #     pdb_id = args.pdb_id if args.pdb_id else default_pdb_id
# #     roll_no = args.roll if args.roll else default_roll_no
    
# #     print(f"Using PDB ID: {pdb_id}, Roll No: {roll_no}")
    
# #     # Download PDB and FASTA files
# #     downloaded_pdb_id = download_pdb(pdb_id)
# #     if not downloaded_pdb_id:
# #         print("Error: Could not download required files. Exiting.")
# #         return
    
# #     # Extract sequences from downloaded files
# #     pdb_sequences, fasta_sequences = extract_sequences(pdb_id)
    
# #     # Identify chain breaks by comparing sequences
# #     chain_breaks = identify_chain_breaks(pdb_sequences, fasta_sequences)
    
# #     # Run NACCESS to calculate ASA; only works if NACCESS is installed properly
# #     naccess_success = run_naccess(pdb_id)
    
# #     # Compute ASA from NACCESS output if it succeeded
# #     asa_data = {}
# #     if naccess_success:
# #         asa_data = compute_asa(pdb_id)
    
# #     # Write analysis results to output file
# #     output_filename = f"A1_{pdb_id}.txt"
# #     with open(output_filename, "w") as output_file:
# #         output_file.write(f"Analysis of PDB ID: {pdb_id}\n")
# #         output_file.write("=" * 50 + "\n\n")
# #         output_file.write(f"Number of chains identified: {len(pdb_sequences)}\n\n")
# #         output_file.write("Chain Information:\n")
# #         output_file.write("-" * 50 + "\n")
# #         for chain_id, sequence in pdb_sequences.items():
# #             num_aa = len(sequence)
# #             mol_weight = compute_molecular_weight(sequence)
# #             asa = asa_data.get(chain_id, 0.0)
# #             output_file.write(f"Chain {chain_id}: {num_aa} amino acids, Molecular Weight: {mol_weight:.2f} Da, ASA: {asa:.2f} Å²\n")
        
# #         output_file.write("\n")
# #         if chain_breaks:
# #             output_file.write("Chain breaks detected:\n")
# #             output_file.write("-" * 50 + "\n")
# #             for chain_id, breaks in chain_breaks.items():
# #                 output_file.write(f"Chain {chain_id}: {len(breaks)} break(s)\n")
# #                 for pos, desc in breaks:
# #                     output_file.write(f"  Position {pos}: {desc}\n")
# #         else:
# #             output_file.write("No chain breaks detected.\n")
    
# #     print(f"Analysis completed. Results written to {output_filename}")
    
# #     # Create directory to organize outputs and code
# #     dirname = f"A1{roll_no}"
# #     os.makedirs(dirname, exist_ok=True)
    
# #     # Copy required files to the directory
# #     for filename in [output_filename, f"{pdb_id}.pdb", f"{pdb_id}.fasta"]:
# #         if os.path.exists(filename):
# #             try:
# #                 subprocess.run(["cp", filename, f"{dirname}/"], check=True)
# #             except Exception as e:
# #                 print(f"Error copying {filename}: {e}")
# #                 os.system(f"cp {filename} {dirname}/")
    
# #     # Create README file with instructions
# #     readme_content = (
# #         "PDB Structure Analysis Program\n"
# #         + "="*30 + "\n\n"
# #         "This program analyzes protein structures downloaded from the RCSB Protein Data Bank.\n\n"
# #         "Usage:\n"
# #         "    python3 A1_22CS10020.py <pdb_id> --roll <your_roll_number>\n\n"
# #         "Example:\n"
# #         "    python3 A1_22CS10020.py 2UV0 --roll 22CS10020\n\n"
# #         "Requirements:\n"
# #         "    - Python 3.6+\n"
# #         "    - Biopython (pip install biopython)\n"
# #         "    - Requests (pip install requests)\n"
# #         "    - NACCESS (optional, for ASA calculations)\n"
# #     )
# #     with open(f"{dirname}/README.txt", "w") as readme_file:
# #         readme_file.write(readme_content)
    
# #     # Copy the python script into the folder
# #     try:
# #         script_path = os.path.abspath(__file__)
# #         if os.path.exists(script_path):
# #             subprocess.run(["cp", script_path, f"{dirname}/A1_{roll_no}.py"], check=True)
# #             subprocess.run(["cp", script_path, f"{dirname}/main.py"], check=True)
# #         else:
# #             subprocess.run(["cp", "A1_22CS10020.py", f"{dirname}/A1_{roll_no}.py"], check=True)
# #             subprocess.run(["cp", "A1_22CS10020.py", f"{dirname}/main.py"], check=True)
# #     except Exception as e:
# #         print(f"Error copying script: {e}")
# #         os.system(f"cp A1_22CS10020.py {dirname}/A1_{roll_no}.py")
# #         os.system(f"cp A1_22CS10020.py {dirname}/main.py")
    
# #     # Create a ZIP archive of the directory
# #     try:
# #         subprocess.run(["zip", "-r", f"A1{roll_no}.zip", dirname], check=True)
# #     except Exception as e:
# #         print(f"Error creating zip file: {e}")
# #         os.system(f"zip -r A1{roll_no}.zip {dirname}")
    
# #     print(f"All files organized and zipped as A1{roll_no}.zip")

# # if __name__ == "__main__":
# #     main()


# import os
# import sys
# import requests
# from Bio import SeqIO
# from Bio.PDB import PDBParser, PPBuilder
# import subprocess
# import re
# import argparse

# # Amino acid molecular weights from Thermo Fischer Scientific
# AA_WEIGHTS = {
#     'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.15,
#     'E': 147.13, 'Q': 146.15, 'G': 75.07, 'H': 155.16, 'I': 131.17,
#     'L': 131.17, 'K': 146.19, 'M': 149.21, 'F': 165.19, 'P': 115.13,
#     'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15
# }

# def download_pdb(pdb_id):
#     """Download PDB and FASTA files from RCSB Protein Data Bank"""
#     print(f"Downloading {pdb_id} from RCSB PDB...")
#     # Download PDB file
#     pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
#     pdb_response = requests.get(pdb_url)
#     if pdb_response.status_code != 200:
#         print(f"Error: Failed to download PDB file for {pdb_id}. Status code: {pdb_response.status_code}")
#         return None
#     with open(f"{pdb_id}.pdb", "w") as pdb_file:
#         pdb_file.write(pdb_response.text)
#     # Download FASTA file
#     fasta_url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
#     fasta_response = requests.get(fasta_url)
#     if fasta_response.status_code != 200:
#         print(f"Error: Failed to download FASTA file for {pdb_id}. Status code: {fasta_response.status_code}")
#         return None
#     with open(f"{pdb_id}.fasta", "w") as fasta_file:
#         fasta_file.write(fasta_response.text)
#     print(f"Downloaded {pdb_id}.pdb and {pdb_id}.fasta successfully")
#     return pdb_id

# def extract_sequences(pdb_id):
#     """Extract sequences from PDB and FASTA files"""
#     print(f"Extracting sequences from {pdb_id}.pdb and {pdb_id}.fasta...")
#     # Parse PDB file
#     parser = PDBParser(QUIET=True)
#     try:
#         structure = parser.get_structure(pdb_id, f"{pdb_id}.pdb")
#     except Exception as e:
#         print(f"Error parsing PDB file: {e}")
#         return {}, {}
#     # Extract sequences from PDB file
#     ppb = PPBuilder()
#     pdb_sequences = {}
#     for model in structure:
#         for chain in model:
#             chain_id = chain.id
#             if chain_id == " ":  # Skip chains with no ID
#                 continue
#             seq = ""
#             for pp in ppb.build_peptides(chain):
#                 seq += str(pp.get_sequence())
#             if seq and chain_id not in pdb_sequences:
#                 pdb_sequences[chain_id] = seq
#     # Extract sequences from FASTA file
#     fasta_sequences = {}
#     try:
#         for record in SeqIO.parse(f"{pdb_id}.fasta", "fasta"):
#             header = record.description
#             chain_id = None
#             # Try different parsing patterns
#             if "|" in header:
#                 parts = header.split("|")
#                 if len(parts) >= 3:
#                     chain_id = parts[2][0]  # First character after second |
#             elif "_" in header:
#                 parts = header.split("_")
#                 if len(parts) >= 2:
#                     chain_id = parts[1][0]  # First character after first _
#             # Try to extract from "Chain X" in description
#             if not chain_id and "Chain " in header:
#                 match = re.search(r"Chain ([A-Za-z0-9])", header)
#                 if match:
#                     chain_id = match.group(1)
#             # Last resort: try first character after PDB ID in header
#             if not chain_id:
#                 match = re.search(rf"{pdb_id}[_|\s]([A-Za-z0-9])", header, re.IGNORECASE)
#                 if match:
#                     chain_id = match.group(1)
#             if chain_id and chain_id not in fasta_sequences:
#                 fasta_sequences[chain_id] = str(record.seq)
#                 print(f"Found FASTA sequence for chain {chain_id}: {len(str(record.seq))} residues")
#     except Exception as e:
#         print(f"Error parsing FASTA file: {e}")
#     print(f"PDB chains found: {', '.join(pdb_sequences.keys())}")
#     print(f"FASTA chains found: {', '.join(fasta_sequences.keys())}")
#     return pdb_sequences, fasta_sequences

# def identify_chain_breaks(pdb_sequences, fasta_sequences):
#     """Identify breaks by comparing PDB and FASTA sequences"""
#     print("Checking for chain breaks...")
#     chain_breaks = {}
#     for chain_id, pdb_seq in pdb_sequences.items():
#         if chain_id not in fasta_sequences:
#             print(f"Warning: Chain {chain_id} found in PDB but not in FASTA")
#             continue
#         fasta_seq = fasta_sequences[chain_id]
#         breaks = []
#         pdb_seq_str = str(pdb_seq)
#         fasta_seq_str = str(fasta_seq)
#         min_length = min(len(pdb_seq_str), len(fasta_seq_str))
#         for i in range(min_length):
#             if pdb_seq_str[i] != fasta_seq_str[i]:
#                 breaks.append((i+1, f"{fasta_seq_str[i]}→{pdb_seq_str[i]}"))
#         if len(pdb_seq_str) != len(fasta_seq_str):
#             if len(pdb_seq_str) < len(fasta_seq_str):
#                 breaks.append((len(pdb_seq_str)+1, "Missing residues in PDB"))
#             else:
#                 breaks.append((len(fasta_seq_str)+1, "Extra residues in PDB"))
#         if breaks:
#             chain_breaks[chain_id] = breaks
#     return chain_breaks

# def run_naccess(pdb_id):
#     """Run NACCESS to calculate solvent accessible surface area"""
#     print("Running NACCESS for surface area calculations...")
#     # Candidate locations for the NACCESS executable.
#     candidates = [
#         "naccess",  # if present in PATH
#         "/usr/local/bin/naccess",
#         "/usr/bin/naccess",
#         "./naccess",
#         "../naccess",
#         os.path.expanduser("~/Downloads/naccess/naccess/naccess"),
#         os.path.expanduser("~/Downloads/naccess/naccess")
#     ]
#     naccess_path = None
#     for loc in candidates:
#         candidate = loc
#         # If the candidate is a directory, look for 'naccess' inside it.
#         if os.path.isdir(candidate):
#             candidate_exe = os.path.join(candidate, "naccess")
#             if os.path.exists(candidate_exe) and os.access(candidate_exe, os.X_OK):
#                 candidate = candidate_exe
#         if os.path.exists(candidate) and os.access(candidate, os.X_OK):
#             naccess_path = candidate
#             break
#         else:
#             # Try using "which" command if not an absolute path
#             try:
#                 which_output = subprocess.check_output(["which", loc], universal_newlines=True).strip()
#                 if which_output:
#                     naccess_path = which_output
#                     break
#             except Exception:
#                 continue
#     if not naccess_path:
#         print("Warning: NACCESS executable not found. ASA calculations will be skipped.")
#         return False
#     print(f"Found NACCESS at: {naccess_path}")
#     # Set naccess_dir as the directory containing the executable.
#     naccess_dir = os.path.dirname(naccess_path)
#     # Locate standard data files.
#     vdw_locations = [
#         "vdw.radii",
#         "/usr/local/share/naccess/vdw.radii",
#         "/usr/share/naccess/vdw.radii",
#         "./vdw.radii",
#         "../vdw.radii",
#         os.path.expanduser("~/Downloads/naccess/naccess/vdw.radii")
#     ]
#     std_locations = [
#         "standard.data",
#         "/usr/local/share/naccess/standard.data",
#         "/usr/share/naccess/standard.data",
#         "./standard.data",
#         "../standard.data",
#         os.path.expanduser("~/Downloads/naccess/naccess/standard.data")
#     ]
#     # Add files from the same directory as the executable.
#     vdw_locations.append(os.path.join(naccess_dir, "vdw.radii"))
#     std_locations.append(os.path.join(naccess_dir, "standard.data"))
#     vdw_path = None
#     for loc in vdw_locations:
#         if os.path.exists(loc):
#             vdw_path = loc
#             break
#     std_path = None
#     for loc in std_locations:
#         if os.path.exists(loc):
#             std_path = loc
#             break
#     if not vdw_path:
#         print("Warning: vdw.radii file not found. NACCESS may not work correctly.")
#     else:
#         print(f"Found vdw.radii at: {vdw_path}")
#     if not std_path:
#         print("Warning: standard.data file not found. NACCESS may not work correctly.")
#     else:
#         print(f"Found standard.data at: {std_path}")
#     # Check for accall in the same directory as the executable.
#     accall_path = os.path.join(naccess_dir, "accall")
#     if not os.path.exists(accall_path):
#         print(f"Warning: 'accall' not found at {accall_path}. NACCESS may not work.")
#     else:
#         print(f"Found accall at: {accall_path}")
#     # Build the command line.
#     cmd = f"{naccess_path} {pdb_id}.pdb"
#     if vdw_path:
#         cmd += f" -r {vdw_path}"
#     if std_path:
#         cmd += f" -s {std_path}"
#     env = os.environ.copy()
#     env["PATH"] = f"{naccess_dir}:{env.get('PATH', '')}"
#     try:
#         print(f"Running NACCESS command: {cmd}")
#         result = subprocess.run(cmd, env=env, shell=True, check=True, capture_output=True, text=True)
#         if result.returncode != 0:
#             print(f"Warning: NACCESS failed with error: {result.stderr}")
#             return False
#         print("NACCESS completed successfully")
#         return True
#     except Exception as e:
#         print(f"Warning: NACCESS failed with error: {e}")
#         return False

# def compute_asa(pdb_id):
#     """Extract ASA data from NACCESS output files (.asa or .rsa)"""
#     asa_data = {}
#     for ext in ['.asa', '.rsa']:
#         asa_file_path = f"{pdb_id}{ext}"
#         if os.path.exists(asa_file_path):
#             try:
#                 with open(asa_file_path) as asa_file:
#                     for line in asa_file:
#                         if line.startswith("RES") or line.startswith("CHAIN"):
#                             parts = line.split()
#                             if len(parts) >= 5:
#                                 chain_id = parts[2]
#                                 if parts[0] == "CHAIN":
#                                     asa = float(parts[4])
#                                     asa_data[chain_id] = asa
#                                 elif chain_id not in asa_data:
#                                     asa = float(parts[4])
#                                     asa_data[chain_id] = asa
#                                 else:
#                                     asa_data[chain_id] += float(parts[4])
#             except Exception as e:
#                 print(f"Warning: Error reading ASA file: {e}")
#             break
#     return asa_data

# def compute_molecular_weight(sequence):
#     """Calculate molecular weight of a protein sequence"""
#     weight = 0.0
#     for aa in sequence:
#         if aa in AA_WEIGHTS:
#             weight += AA_WEIGHTS[aa]
#         else:
#             print(f"Warning: Non-standard amino acid '{aa}' encountered. Using average weight.")
#             weight += sum(AA_WEIGHTS.values()) / len(AA_WEIGHTS)
#     weight += 18.02  # weight for terminal groups
#     return weight

# def parse_arguments():
#     """Parse command line arguments"""
#     parser = argparse.ArgumentParser(description='Analyze protein structures from PDB.')
#     parser.add_argument('pdb_id', nargs='?', type=str, help='PDB ID (e.g., 1ABC)')
#     parser.add_argument('--roll', type=str, help='Your roll number')
#     return parser.parse_args()

# def main():
#     default_pdb_id = "2UV0"
#     default_roll_no = "22CS10020"  # Replace with your actual roll number
#     args = parse_arguments()
#     pdb_id = args.pdb_id if args.pdb_id else default_pdb_id
#     roll_no = args.roll if args.roll else default_roll_no
#     print(f"Using PDB ID: {pdb_id}, Roll No: {roll_no}")
#     downloaded_pdb_id = download_pdb(pdb_id)
#     if not downloaded_pdb_id:
#         print("Error: Could not download required files. Exiting.")
#         return
#     pdb_sequences, fasta_sequences = extract_sequences(pdb_id)
#     chain_breaks = identify_chain_breaks(pdb_sequences, fasta_sequences)
#     naccess_success = run_naccess(pdb_id)
#     asa_data = {}
#     if naccess_success:
#         asa_data = compute_asa(pdb_id)
#     output_filename = f"A1_{pdb_id}.txt"
#     with open(output_filename, "w") as output_file:
#         output_file.write(f"Analysis of PDB ID: {pdb_id}\n")
#         output_file.write("=" * 50 + "\n\n")
#         output_file.write(f"Number of chains identified: {len(pdb_sequences)}\n\n")
#         output_file.write("Chain Information:\n")
#         output_file.write("-" * 50 + "\n")
#         for chain_id, sequence in pdb_sequences.items():
#             num_aa = len(sequence)
#             mol_weight = compute_molecular_weight(sequence)
#             asa = asa_data.get(chain_id, 0.0)
#             output_file.write(f"Chain {chain_id}: {num_aa} amino acids, Molecular Weight: {mol_weight:.2f} Da, ASA: {asa:.2f} Å²\n")
#         output_file.write("\n")
#         if chain_breaks:
#             output_file.write("Chain breaks detected:\n")
#             output_file.write("-" * 50 + "\n")
#             for chain_id, breaks in chain_breaks.items():
#                 output_file.write(f"Chain {chain_id}: {len(breaks)} break(s)\n")
#                 for pos, desc in breaks:
#                     output_file.write(f"  Position {pos}: {desc}\n")
#         else:
#             output_file.write("No chain breaks detected.\n")
#     print(f"Analysis completed. Results written to {output_filename}")
#     dirname = f"A1{roll_no}"
#     os.makedirs(dirname, exist_ok=True)
#     for filename in [output_filename, f"{pdb_id}.pdb", f"{pdb_id}.fasta"]:
#         if os.path.exists(filename):
#             try:
#                 subprocess.run(["cp", filename, f"{dirname}/"], check=True)
#             except Exception as e:
#                 print(f"Error copying {filename}: {e}")
#                 os.system(f"cp {filename} {dirname}/")
#     readme_content = (
#         "PDB Structure Analysis Program\n"
#         "==============================\n\n"
#         "This program analyzes protein structures downloaded from the RCSB Protein Data Bank.\n\n"
#         "Usage:\n"
#         "    python3 A1_22CS10020.py <pdb_id> --roll <your_roll_number>\n\n"
#         "Example:\n"
#         "    python3 A1_22CS10020.py 2UV0 --roll 22CS10020\n\n"
#         "Requirements:\n"
#         "    - Python 3.6+\n"
#         "    - Biopython (pip install biopython)\n"
#         "    - Requests (pip install requests)\n"
#         "    - NACCESS (optional, for ASA calculations)\n"
#     )
#     with open(f"{dirname}/README.txt", "w") as readme_file:
#         readme_file.write(readme_content)
#     try:
#         script_path = os.path.abspath(__file__)
#         if os.path.exists(script_path):
#             subprocess.run(["cp", script_path, f"{dirname}/A1_{roll_no}.py"], check=True)
#             subprocess.run(["cp", script_path, f"{dirname}/main.py"], check=True)
#         else:
#             subprocess.run(["cp", "A1_22CS10020.py", f"{dirname}/A1_{roll_no}.py"], check=True)
#             subprocess.run(["cp", "A1_22CS10020.py", f"{dirname}/main.py"], check=True)
#     except Exception as e:
#         print(f"Error copying script: {e}")
#         os.system(f"cp A1_22CS10020.py {dirname}/A1_{roll_no}.py")
#         os.system(f"cp A1_22CS10020.py {dirname}/main.py")
#     try:
#         subprocess.run(["zip", "-r", f"A1{roll_no}.zip", dirname], check=True)
#     except Exception as e:
#         print(f"Error creating zip file: {e}")
#         os.system(f"zip -r A1{roll_no}.zip {dirname}")
#     print(f"All files organized and zipped as A1{roll_no}.zip")

# if __name__ == "__main__":
#     main()


import os
import sys
import requests
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

def run_naccess(pdb_id):
    """Run NACCESS to calculate solvent accessible surface area"""
    print("Running NACCESS for surface area calculations...")
    candidates = [
        "naccess",
        "/usr/local/bin/naccess",
        "/usr/bin/naccess",
        "./naccess",
        "../naccess",
        os.path.expanduser("~/Downloads/naccess/naccess/naccess"),
        os.path.expanduser("~/Downloads/naccess/naccess")
    ]
    naccess_path = None
    for loc in candidates:
        candidate = loc
        if os.path.isdir(candidate):
            candidate_exe = os.path.join(candidate, "naccess")
            if os.path.exists(candidate_exe) and os.access(candidate_exe, os.X_OK):
                candidate = candidate_exe
        if os.path.exists(candidate) and os.access(candidate, os.X_OK):
            naccess_path = candidate
            break
        else:
            try:
                which_output = subprocess.check_output(["which", loc], universal_newlines=True).strip()
                if which_output:
                    naccess_path = which_output
                    break
            except Exception:
                continue
    if not naccess_path:
        print("Warning: NACCESS executable not found. ASA calculations will be skipped.")
        return False
    print(f"Found NACCESS at: {naccess_path}")
    naccess_dir = os.path.dirname(naccess_path)
    vdw_locations = [
        "vdw.radii",
        "/usr/local/share/naccess/vdw.radii",
        "/usr/share/naccess/vdw.radii",
        "./vdw.radii",
        "../vdw.radii",
        os.path.expanduser("~/Downloads/naccess/naccess/vdw.radii")
    ]
    std_locations = [
        "standard.data",
        "/usr/local/share/naccess/standard.data",
        "/usr/share/naccess/standard.data",
        "./standard.data",
        "../standard.data",
        os.path.expanduser("~/Downloads/naccess/naccess/standard.data")
    ]
    vdw_locations.append(os.path.join(naccess_dir, "vdw.radii"))
    std_locations.append(os.path.join(naccess_dir, "standard.data"))
    vdw_path = None
    for loc in vdw_locations:
        if os.path.exists(loc):
            vdw_path = loc
            break
    std_path = None
    for loc in std_locations:
        if os.path.exists(loc):
            std_path = loc
            break
    if not vdw_path:
        print("Warning: vdw.radii file not found. NACCESS may not work correctly.")
    else:
        print(f"Found vdw.radii at: {vdw_path}")
    if not std_path:
        print("Warning: standard.data file not found. NACCESS may not work correctly.")
    else:
        print(f"Found standard.data at: {std_path}")
    accall_path = os.path.join(naccess_dir, "accall")
    if not os.path.exists(accall_path):
        print(f"Warning: 'accall' not found at {accall_path}. NACCESS may not work.")
    else:
        print(f"Found accall at: {accall_path}")
    cmd = f"{naccess_path} {pdb_id}.pdb"
    if vdw_path:
        cmd += f" -r {vdw_path}"
    if std_path:
        cmd += f" -s {std_path}"
    env = os.environ.copy()
    env["PATH"] = f"{naccess_dir}:{env.get('PATH', '')}"
    try:
        print(f"Running NACCESS command: {cmd}")
        result = subprocess.run(cmd, env=env, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Warning: NACCESS failed with exit status {result.returncode}")
            print("Error output from NACCESS:")
            print(result.stderr)
            return False
        print("NACCESS completed successfully")
        return True
    except Exception as e:
        print(f"Warning: NACCESS failed with error: {e}")
        return False

def compute_asa(pdb_id):
    """Extract ASA data from NACCESS output files (.asa or .rsa)"""
    asa_data = {}
    for ext in ['.asa', '.rsa']:
        asa_file_path = f"{pdb_id}{ext}"
        if os.path.exists(asa_file_path):
            try:
                with open(asa_file_path) as asa_file:
                    for line in asa_file:
                        if line.startswith("RES") or line.startswith("CHAIN"):
                            parts = line.split()
                            if len(parts) >= 5:
                                chain_id = parts[2]
                                if parts[0] == "CHAIN":
                                    asa = float(parts[4])
                                    asa_data[chain_id] = asa
                                elif chain_id not in asa_data:
                                    asa = float(parts[4])
                                    asa_data[chain_id] = asa
                                else:
                                    asa_data[chain_id] += float(parts[4])
            except Exception as e:
                print(f"Warning: Error reading ASA file: {e}")
            break
    return asa_data

def compute_molecular_weight(sequence):
    """Calculate molecular weight of a protein sequence"""
    weight = 0.0
    for aa in sequence:
        if aa in AA_WEIGHTS:
            weight += AA_WEIGHTS[aa]
        else:
            print(f"Warning: Non-standard amino acid '{aa}' encountered. Using average weight.")
            weight += sum(AA_WEIGHTS.values()) / len(AA_WEIGHTS)
    weight += 18.02  # Add terminal groups
    return weight

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Analyze protein structures from PDB.')
    parser.add_argument('pdb_id', nargs='?', type=str, help='PDB ID (e.g., 1ABC)')
    parser.add_argument('--roll', type=str, help='Your roll number')
    return parser.parse_args()

def main():
    default_pdb_id = "2UV0"
    default_roll_no = "22CS10020"
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
    dirname = f"A1{roll_no}"
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
        "    python3 A1_22CS10020.py <pdb_id> --roll <your_roll_number>\n\n"
        "Example:\n"
        "    python3 A1_22CS10020.py 2UV0 --roll 22CS10020\n\n"
        "Requirements:\n"
        "    - Python 3.6+\n"
        "    - Biopython (pip install biopython)\n"
        "    - Requests (pip install requests)\n"
        "    - NACCESS (optional, for ASA calculations)\n"
    )
    with open(f"{dirname}/README.txt", "w") as readme_file:
        readme_file.write(readme_content)
    try:
        script_path = os.path.abspath(__file__)
        if os.path.exists(script_path):
            subprocess.run(["cp", script_path, f"{dirname}/A1_{roll_no}.py"], check=True)
            subprocess.run(["cp", script_path, f"{dirname}/main.py"], check=True)
        else:
            subprocess.run(["cp", "A1_22CS10020.py", f"{dirname}/A1_{roll_no}.py"], check=True)
            subprocess.run(["cp", "A1_22CS10020.py", f"{dirname}/main.py"], check=True)
    except Exception as e:
        print(f"Error copying script: {e}")
        os.system(f"cp A1_22CS10020.py {dirname}/A1_{roll_no}.py")
        os.system(f"cp A1_22CS10020.py {dirname}/main.py")
    try:
        subprocess.run(["zip", "-r", f"A1{roll_no}.zip", dirname], check=True)
    except Exception as e:
        print(f"Error creating zip file: {e}")
        os.system(f"zip -r A1{roll_no}.zip {dirname}")
    print(f"All files organized and zipped as A1{roll_no}.zip")

if __name__ == "__main__":
    main()