PDB Structure Analysis Program
==============================

This program analyzes protein structures downloaded from the RCSB Protein Data Bank.

Usage:
    python3 A1_22CS10020.py <pdb_id> --roll <roll_number>

Example:
    python3 A1_22CS10020.py 2UV0 --roll 22CS10020

Requirements:
    - Python 3.6+
    - Biopython (pip install biopython)
    - Requests (pip install requests)
    - NACCESS (optional, for ASA calculations)

Note:
    Ensure that the NACCESS executable is located in one of the following paths:
    - naccess
    - /mnt/c/ComputationalBiophysics/for_student
    - or any other directory specified in the candidates list in the script.
