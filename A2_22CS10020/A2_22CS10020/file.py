
from Bio.PDB import PDBParser, ShrakeRupley, Superimposer
import numpy as np
import os

# Solvation parameters from SI.pdf Table S1 (cal Å⁻² mol⁻¹)
SOLVATION_PARAMS = {
    'C': 16,   # Carbon
    'N': -6,   # Neutral nitrogen
    'O': -6,# Neutral oxygen
    'N+': -50, # Charged nitrogen of ARG, LYS
    'O-': -24, # Charged oxygen of ASP, GLU
    'S': -21 # Sulfur
}

def calculate_solvation_energy(structure):
    try:
        sr = ShrakeRupley()
        sr.compute(structure, level="A")
        solvation_energy = 0
        for atom in structure.get_atoms():
            residue = atom.get_parent()
            residue_name = residue.resname
            atom_name = atom.name
            if "C" in atom_name:
                atom_type = 'C'
            elif "N" in atom_name:
                atom_type = 'N+' if residue_name in ["ARG", "LYS"] else 'N'
            elif "O" in atom_name:
                atom_type = 'O-' if residue_name in ["ASP", "GLU"] else 'O'
            elif "S" in atom_name:
                atom_type = 'S'
            else:
                continue
            if hasattr(atom, 'sasa') and atom.sasa is not None and atom_type in SOLVATION_PARAMS:
                solvation_energy += SOLVATION_PARAMS[atom_type] * atom.sasa
        return solvation_energy
    except Exception as e:
        print(f"Error in solvation energy: {e}")
        return 0

def calculate_interface_area(complex_struct, target_struct):
    try:
        sr = ShrakeRupley()
        sr.compute(complex_struct, level="A")
        asa_complex = sum(atom.sasa for atom in complex_struct.get_atoms() if hasattr(atom, 'sasa'))
        chain_A = target_struct[0]["A"]
        sr.compute(chain_A, level="A")
        asa_A = sum(atom.sasa for atom in chain_A.get_atoms() if hasattr(atom, 'sasa'))
        chain_B = target_struct[0]["B"]
        sr.compute(chain_B, level="A")
        asa_B = sum(atom.sasa for atom in chain_B.get_atoms() if hasattr(atom, 'sasa'))
        interface_area = (asa_A + asa_B - asa_complex) / 2
        return max(interface_area, 0)
    except Exception as e:
        print(f"Error in interface area: {e}")
        return 0

def calculate_rmsd(coords1, coords2):
    try:
        diff = np.array(coords1) - np.array(coords2)
        return np.sqrt(np.mean(np.sum(diff * diff, axis=1)))
    except Exception as e:
        print(f"Error in RMSD calculation: {e}")
        return None

def calculate_docking_metrics(target_struct, decoy_struct):
    try:
        target_chain_A = target_struct[0]["A"]
        target_chain_B = target_struct[0]["B"]
        decoy_chain_A = decoy_struct[0]["A"]
        decoy_chain_B = decoy_struct[0]["B"]
    except KeyError as e:
        print(f"Error: Missing chain {e}")
        return None, None, None

    try:
        sup = Superimposer()
        target_A_atoms = [atom for atom in target_chain_A.get_atoms() if atom.name == "CA"]
        decoy_A_atoms = [atom for atom in decoy_chain_A.get_atoms() if atom.name == "CA"]
        if len(target_A_atoms) != len(decoy_A_atoms):
            print(f"Warning: Mismatched CA atoms in chain A (target: {len(target_A_atoms)}, decoy: {len(decoy_A_atoms)})")
        sup.set_atoms(target_A_atoms, decoy_A_atoms)
        sup.apply(decoy_struct[0].get_atoms())
    except Exception as e:
        print(f"Error in superposition: {e}")
        return None, None, None
    
    try:
        target_B_coords = [atom.coord for atom in target_chain_B.get_atoms() if atom.name == "CA"]
        decoy_B_coords = [atom.coord for atom in decoy_chain_B.get_atoms() if atom.name == "CA"]
        if len(target_B_coords) != len(decoy_B_coords):
            print(f"Warning: Mismatched CA atoms in chain B (target: {len(target_B_coords)}, decoy: {len(decoy_B_coords)})")
        lrmsd = calculate_rmsd(target_B_coords, decoy_B_coords)
    except Exception as e:
        print(f"Error in LRMSD: {e}")
        lrmsd = None
    
    # IRMSD: RMSD of interface residues (within 10 Å in target)
    try:
        interface_dist = 10.0
        target_interface_res = set()
        for a in target_chain_A.get_atoms():
            for b in target_chain_B.get_atoms():
                if a - b < interface_dist and a.name == "CA":
                    target_interface_res.add((a.get_parent().get_id()[1], 'A'))
                if a - b < interface_dist and b.name == "CA":
                    target_interface_res.add((b.get_parent().get_id()[1], 'B'))
        
        target_coords = []
        decoy_coords = []
        for res_id, chain_id in target_interface_res:
            try:
                if chain_id == 'A':
                    target_res = target_chain_A[res_id]
                    decoy_res = decoy_chain_A[res_id]
                else:
                    target_res = target_chain_B[res_id]
                    decoy_res = decoy_chain_B[res_id]
                target_ca = [atom for atom in target_res if atom.name == "CA"][0]
                decoy_ca = [atom for atom in decoy_res if atom.name == "CA"][0]
                target_coords.append(target_ca.coord)
                decoy_coords.append(decoy_ca.coord)
            except (KeyError, IndexError):
                continue  # Skiping th missing residues( CA atoms)
        irmsd = calculate_rmsd(target_coords, decoy_coords) if target_coords else 0
    except Exception as e:
        print(f"Error in IRMSD: {e}")
        irmsd = None
    
    # Fnat: Fraction of native contacts (<5 Å)
    try:
        contact_dist = 5.0
        target_contacts = set()
        decoy_contacts = set()
        for a in target_chain_A.get_atoms():
            for b in target_chain_B.get_atoms():
                if a - b < contact_dist and a.name == "CA" and b.name == "CA":
                    target_contacts.add((a.get_parent().id[1], b.get_parent().id[1]))
        for a in decoy_chain_A.get_atoms():
            for b in decoy_chain_B.get_atoms():
                if a - b < contact_dist and a.name == "CA" and b.name == "CA":
                    decoy_contacts.add((a.get_parent().id[1], b.get_parent().id[1]))
        fnat = len(target_contacts & decoy_contacts) / len(target_contacts) if target_contacts else 0
    except Exception as e:
        print(f"Error in Fnat: {e}")
        fnat = None
    
    return lrmsd, irmsd, fnat

def main():
    parser = PDBParser(QUIET=True)
    try:
        target_struct = parser.get_structure("target", "target.pdb")
        print("Loaded target.pdb successfully")
    except Exception as e:
        print(f"Error loading target.pdb: {e}")
        return
    
    decoy_dir = "Decoys"
    if not os.path.exists(decoy_dir):
        print(f"Error: Directory {decoy_dir} not found")
        return
    
    results = []
    decoy_files = [f for f in os.listdir(decoy_dir) if f.endswith(".pdb")]
    if not decoy_files:
        print("Error: No .pdb files found in Decoys directory")
        return
    
    for decoy_file in sorted(decoy_files):
        decoy_path = os.path.join(decoy_dir, decoy_file)
        try:
            decoy_struct = parser.get_structure("decoy", decoy_path)
            print(f"Loaded {decoy_file}")
            
            interface_area = calculate_interface_area(decoy_struct, target_struct) #interface area
            solvation_energy = calculate_solvation_energy(decoy_struct) # solvation energy
            lrmsd, irmsd, fnat = calculate_docking_metrics(target_struct, decoy_struct) # lrmsd, irmsd, fnat
            
            normalized_name = decoy_file.replace("complex.", "complex")
            results.append((normalized_name, interface_area, solvation_energy, lrmsd, irmsd, fnat))
            print(f"Processed {decoy_file}: IA={interface_area:.2f}, SE={solvation_energy:.2f}, LRMSD={lrmsd}, IRMSD={irmsd}, Fnat={fnat}")
        except Exception as e:
            print(f"Error processing {decoy_file}: {e}")
            normalized_name = decoy_file.replace("complex.", "complex")
            results.append((normalized_name, 0, 0, None, None, None))
    with open("Score.txt", "w") as f:
        f.write("Filename           | Interface Area | Solvation Energy | LRMSD | IRMSD | Fnat score\n")
        f.write("-" * 90 + "\n")
        for result in results:
            filename, ia, se, lrmsd, irmsd, fnat = result
            ia_str = f"{ia:.2f}" if ia != 0 else "0.00"
            se_str = f"{se:.2f}" if se != 0 else "0.00"
            lrmsd_str = f"{lrmsd:.2f}" if lrmsd is not None else "N/A"
            irmsd_str = f"{irmsd:.2f}" if irmsd is not None else "N/A"
            fnat_str = f"{fnat:.4f}" if fnat is not None else "N/A"
            f.write(f"{filename:<18} | {ia_str:>14} | {se_str:>16} | {lrmsd_str:>5} | {irmsd_str:>5} | {fnat_str:>10}\n")

if __name__ == "__main__":
    main()  
