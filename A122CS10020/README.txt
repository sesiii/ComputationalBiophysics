PDB Structure Analysis Program
==============================

This program analyzes protein structures downloaded from the RCSB Protein Data Bank.

Usage:
    python3 A1_22CS10020.py <pdb_id> --roll <your_roll_number>

Example:
    python3 A1_22CS10020.py 2UV0 --roll 22CS10020

Requirements:
    - Python 3.6+
    - Biopython (pip install biopython)
    - Requests (pip install requests)
    - NACCESS (optional, for ASA calculations)
