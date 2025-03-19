This script processes protein structure files to calculate various metrics such as interface area, solvation energy, and docking metrics (LRMSD, IRMSD, and Fnat). The results are written to a file named Score.txt

## Usage
```bash
python3 file.py
```
### Requirements
Python 3.12 \
Biopython library

### Usage
Place the target protein structure file as `{target.pdb}` in the same directory as the script `file.py`.

Place the decoy protein structure files in a directory named `Decoys` in the same directory as the script. The decoy files should have the `.pdb` extension.

### Run the script:

The script will process the target and decoy files, calculate the metrics, and write the results to `Score.txt` file.

### Output
The output file `Score.txt` will contain the following columns:

- Filename
- Interface Area
- Solvation Energy
- LRMSD
- IRMSD
- Fnat score

### Functions
- `calculate_solvation_energy(structure)`: Calculates the solvation energy of a given structure.
- `calculate_interface_area(complex_struct, target_struct)`: Calculates the interface area between the complex and target structures.
- `calculate_rmsd(coords1, coords2)`: Calculates the Root Mean Square Deviation (RMSD) between two sets of coordinates.
- `calculate_docking_metrics(target_struct, decoy_struct)`: Calculates the docking metrics (LRMSD, IRMSD, and Fnat) between the target and decoy structures.
- `main()`: Main function that loads the target and decoy structures, calculates the metrics, and writes the results to `Score.txt`.
