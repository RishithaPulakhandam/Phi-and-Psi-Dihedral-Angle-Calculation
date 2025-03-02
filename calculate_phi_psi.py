import numpy as np

def parse_pdb(pdb_file):
    coordinates = {}
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith("ATOM"):
                atom_type = line[13:15].strip()
                if atom_type in ["N", "CA", "C"]:
                    residue_index = int(line[22:26].strip())
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    if residue_index not in coordinates:
                        coordinates[residue_index] = {}
                    coordinates[residue_index][atom_type] = np.array([x, y, z])
    return coordinates

def calculate_dihedral(p1, p2, p3, p4):
    b1 = p1 - p2
    b2 = p2 - p3
    b3 = p3 - p4
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    n1 /= np.linalg.norm(n1)
    n2 /= np.linalg.norm(n2)
    x = np.dot(n1, n2)
    m = np.cross(n1, n2)
    angle = -np.arccos(x) if np.dot(b2, m) > 0 else np.arccos(x)
    return np.degrees(angle)

def calculate_phi_psi(pdb_file):
    coordinates = parse_pdb(pdb_file)
    phi_angles = {}
    psi_angles = {}
    for res_index in sorted(coordinates.keys()):
        if res_index - 1 in coordinates and res_index + 1 in coordinates:
            try:
                c_prev = coordinates[res_index - 1]["C"]
                n = coordinates[res_index]["N"]
                ca = coordinates[res_index]["CA"]
                c = coordinates[res_index]["C"]
                n_next = coordinates[res_index + 1]["N"]
                phi_angles[res_index] = calculate_dihedral(c_prev, n, ca, c)
                psi_angles[res_index] = calculate_dihedral(n, ca, c, n_next)
            except KeyError:
                continue
    return phi_angles, psi_angles

def classify_secondary_structure(phi_angles, psi_angles):
    classification = {}
    for res_index in phi_angles.keys():
        phi = phi_angles[res_index]
        psi = psi_angles.get(res_index, 0)
        if -70 <= phi <= -50 and -60 <= psi <= -30:
            classification[res_index] = 'alpha'
        elif -140 <= phi <= -120 and 120 <= psi <= 150:
            classification[res_index] = 'beta'
        else:
            classification[res_index] = 'other'
    return classification

def write_angles_to_file(phi_angles, psi_angles, classification, filename="dihedral_angles.txt", classification_filename="secondary_structures.txt"):
    with open(filename, 'w') as out_file:
        for res_index in sorted(phi_angles.keys()):
            out_file.write(f"Residue {res_index} - Phi: {phi_angles[res_index]:.2f}, Psi: {psi_angles[res_index]:.2f}\n")
    with open(classification_filename, 'w') as class_file:
        for res_index in sorted(classification.keys()):
            class_file.write(f"Residue {res_index}: {classification[res_index]}\n")

if __name__ == "__main__":
    pdb_file_path = input("Enter the path to the PDB file: ")
    phi_angles, psi_angles = calculate_phi_psi(pdb_file_path)
    classification = classify_secondary_structure(phi_angles, psi_angles)
    write_angles_to_file(phi_angles, psi_angles, classification)
    print("Phi and Psi angles along with secondary structures have been calculated and written to files.")
