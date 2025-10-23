# 🧬 Protein-Ligand Interaction Analyzer

This package analyzes atomic-level interactions between proteins and ligands to identify optimal drug candidates. It calculates bond types between amino acid atoms and ligand atoms, then utilize the **AutoDock Vina scoring function** to evaluate binding affinity and rank potential protein inhibitor.

---

## 🧠 Project Overview

- **Goal**: Identify the most compatible ligand for a given protein binding site.
- **Method**: Analyze atomic interactions and apply Vina scoring to evaluate binding strength.
- **Applications**: Drug discovery, molecular docking, structure-based screening.

---

## 🔬 Core Features

### 1. Amino Acid Atom Extraction
- Parses protein structure files (e.g., PDB format).
- Extracts atomic coordinates of amino acids in the binding groove.

### 2. Ligand Atom Parsing
- Loads ligand molecules and extracts atomic data.
- Supports common formats like MOL2, SDF, or PDBQT.

### 3. Bond Type Calculation
- Computes distances between atoms to infer:
  - Hydrogen bonds
  - Ionic interactions
  - Hydrophobic contacts
  - Van der Waals forces

### 4. Vina Scoring Integration
- Applies the **AutoDock Vina scoring function** to quantify binding affinity.
- Combines geometric fit and interaction energy to rank ligands.

