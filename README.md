# Phi-and-Psi-Dihedral-Angle-Calculation-
Phi and Psi Dihedral Angle Calculation from PDB Files

## **Overview**

This project provides a Python tool for computing phi (φ) and psi (ψ) dihedral angles in protein structures using atomic coordinates extracted from Protein Data Bank (PDB) files. It classifies secondary structures based on calculated dihedral angles and validates them against known structures.

## **Methodology**

1. **Parsing the PDB File:** Extracts atomic coordinates for backbone atoms (N, Cα, C) for each residue.
2. **Dihedral Angle Calculation:**
   - Defines backbone vectors between consecutive atoms.
   - Computes normal vectors to the planes formed by these atoms.
   - Uses the dot product and arccosine function to determine dihedral angles.
   - Determines sign using vector cross-product.
3. **Secondary Structure Classification:**
   - Classifies residues as alpha-helix, beta-sheet, or other based on angle thresholds.
4. **Validation:**
   - Compares computed angles with DSSP database annotations.
   - Visualizes results using Chimera.

## **Installation and Requirements**

- Python 3.7+
- NumPy
- Chimera (for visualization, optional)

## **Usage**

1. Place your PDB file in the working directory. (eg; 3GHG.pdb)
2. Run the script:
   ```bash
   python calculate_phi_psi.py
   ```
3. Enter the PDB file path when prompted.
4. Outputs:
   - `dihedral_angles.txt`: Contains computed phi and psi angles per residue.
   - `secondary_structures.txt`: Contains classified secondary structures.

## **Validation**

- Results were verified by comparison with DSSP annotations.
- Ramachandran plot visualization in Chimera confirmed expected angle distributions.

## **Example Output**

```
Residue 10 - Phi: -60.23, Psi: -40.12
Residue 11 - Phi: -58.90, Psi: -45.34
...
Residue 10: Alpha-helix
Residue 11: Alpha-helix
...
```
  <img width="274" alt="image" src="https://github.com/user-attachments/assets/bb8230f1-cba6-463a-8390-1dd0f6379239" />
  
I validated the result against the DSSP file and visualized the protein structure with calculated dihedral
angles using chimera.

   <img width="643" alt="image" src="https://github.com/user-attachments/assets/e7e122e6-b227-4a68-bda5-22027cc51ce9" />
   <img width="799" alt="image" src="https://github.com/user-attachments/assets/b192d154-ef11-4a1a-a5ef-5a3c4561fa4e" />
   <img width="614" alt="image" src="https://github.com/user-attachments/assets/9097497f-62d8-4601-b464-8a0d14914ebe" />

## **Future Improvements**

- Implementing a more refined machine-learning-based classification.
- Supporting additional file formats (e.g., mmCIF).
- Integrating visualization directly in the script.

## **Authors**

Developed by Rishitha Pulakhandam. 

