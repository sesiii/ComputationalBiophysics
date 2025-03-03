# Computational Biophysics Assignment

This repository contains a Python program that performs various tasks related to protein structures from the RCSB Protein Data Bank (PDB). The program takes a PDB ID as input and performs the following tasks:

1. Downloads the corresponding PDB structure and the FASTA sequence file from RCSB Protein Data Bank (rcsb.org).
2. Extracts the protein sequence from the PDB file and compares it with the FASTA sequence file.
3. Identifies the number of chains in the protein structure.
4. Checks for chain breaks in the structure and reports the chain ID and residue positions where breaks occur, if any.
5. Runs the NACCESS program to compute the solvent accessible surface area (ASA) of each chain.
6. Reports the chain ID, number of amino acids present, molecular weight of each chain, and ASA.

## File Structure

- `A1_<RollNo>.py`: The main Python program.
- `A1_<PDBID>.txt`: Example output file.
- `README.md`: This file, containing instructions on how to run the program.

## How to Run the Program

1. Ensure you have Python installed on your system.
2. Install the required Python packages using the following command:
    ```bash
    pip install requests biopython
    ```
3. Download and install the NACCESS program from [here](http://www.bioinf.manchester.ac.uk/naccess/).
4. Place the NACCESS executable in your system's PATH.
5. Run the Python program with the following command:
    ```bash
    python A1_<RollNo>.py <PDB_ID>
    ```
   Replace `<RollNo>` with your roll number and `<PDB_ID>` with the PDB ID of the protein structure you want to analyze.

## Example

To run the program for the PDB ID `1ABC`:
```bash
python A1_<RollNo>.py 1ABC
```

This will generate an output file named `A1_1ABC.txt` containing the results of the analysis.

## Notes

- Ensure you have an active internet connection to download the PDB and FASTA files.
- The program assumes the NACCESS executable is available in the system's PATH.

## References

- [RCSB Protein Data Bank](https://www.rcsb.org/)
- [Thermo Fischer Scientific - Proteins and Amino Acids](https://www.thermofisher.com/in/en/home/references/ambion-tech-support/rna-tools-and-calculators/proteins-and-amino-acids.html)
