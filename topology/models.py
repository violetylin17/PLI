from ..guess import guess_atom_name

all_atom_model_atom_selection = {
    'ALA': '(name N HN CA HA CB HB1 HB2 HB3 C O)',
    'ARG': '(name N HN CA HA CB HB1 HB2 CG HG1 HG2 CD HD1 HD2 NE HE CZ NH1 HH11 HH12 NH2 HH21 HH22 C O)',
    'ASP': '(name N HN CA HA CB HB1 HB2 CG OD1 OD2 C O)',
    'ASN': '(name N HN CA HA CB HB1 HB2 CG OD1 ND2 HD21 HD22 C O)',
    'CYS': '(name N HN CA HA CB HB1 HB2 SG HG1 C O)',
    'GLN': '(name N HN CA HA CB HB1 HB2 CG HG1 HG2 CD OE1 NE2 HE21 HE22 C O)',
    'GLU': '(name N HN CA HA CB HB1 HB2 CG HG1 HG2 CD OE1 OE2 C O)',
    'GLY': '(name N HN CA HA1 HA2 C O)',
    'HS2': '(name CE1 HE1 ND1 CG CB HB1 HB2 NE2 HE2 CD2 HD2)',
    'HIS': '(name N HN CA HA CB HB1 HB2 CG ND1 HD1 CD2 NE2 HE2 CE1 C O)',
    'HSD': '(name N HN CA HA CB HB1 HB2 ND1 HD1 CG CE1 HE1 NE2 CD2 HD2 C O)',
    'HSE': '(name N HN CA HA CB HB1 HB2 ND1 CG CE1 HE1 NE2 HE2 CD2 HD2 C O)',
    'HSP': '(name N HN CA HA CB HB1 HB2 CD2 HD2 CG NE2 HE2 ND1 HD1 CE1 HE1 C O )',
    'ILE': '(name N HN CA HA CB HB CG2 HG21 HG22 HG23 CG1 HG11 HG12 CD HD1 HD2 HD3 C O)',
    'LEU': '(name N HN CA HA CB HB1 HB2 CG HG CD1 HD11 HD12 HD13 CD2 HD21 HD22 HD23 C O)',
    'LYS': '(name N HN CA HA CB HB1 HB2 CG HG1 HG2 CD HD1 HD2 CE HE1 HE2 NZ HZ1 HZ2 HZ3 C O)',
    'MET': '(name N HN CA HA CB HB1 HB2 CG HG1 HG2 SD CE HE1 HE2 HE3 C O)',
    'PHE': '(name N HN CA HA CB HB1 HB2 CG CD1 HD1 CE1 HE1 CZ HZ CD2 HD2 CE2 HE2 C O)',
    'PRO': '(name N CD HD1 HD2 CA HA CB HB1 HB2 CG HG1 HG2 C O)',
    'SER': '(name N HN CA HA CB HB1 HB2 OG HG1 C O)',
    'THR': '(name N HN CA HA CB HB OG1 HG1 CG2 HG21 HG22 HG23 C O OXT)',
    'TRP': '(name N HN CA HA CB HB1 HB2 CG CD1 HD1 NE1 HE1 CE2 CD2 CE3 HE3 CZ3 HZ3 CZ2 HZ2 CH2 HH2 C O)',
    'TYR': '(name N HN CA HA CB HB1 HB2 CG CD1 HD1 CE1 HE1 CZ OH HH CD2 HD2 CE2 HE2 C O)',
    'VAL': '(name N HN CA HA CB HB CG1 HG11 HG12 HG13 CG2 HG21 HG22 HG23 C O)',
    'NCTER': '(name HT1 HT2 HT3 OT1 OT2)',
    'NTER': '(name N HT1 HT2 HT3 CA HA)',
    'CTER': '(name C OT1 OT2)'
}

# CHARMM c43b1 topall36 force field            
all_atom_model_info = {                         #(atomtype, hbond, hydro, charge)
     'ALA': {
        "N": ("NH1", "D", "", -0.47),       # (atom_type, hbond_type, charge)
        "HN": ("H", "HD", "", 0.31),
        "CA": ("CT1", "none", "H", 0.07),
        "HA": ("CT1", "none", "H", 0.09),
        "CB": ("CT3", "none", "H", -0.27),
        "HB1": ("HA3", "", "", 0.09),
        "HB2": ("HA3", "", "", 0.09),
        "HB3": ("HA3", "", "", 0.09),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'ARG': {
        "N": ("NH1", "D", "", -0.47),
        "HN": ("H", "HD", "", 0.31),
        "CA": ("CT1", "none", "H", 0.07),
        "HA": ("HB1", "", "", 0.09),
        "CB": ("CT2", "none", "H", -0.18),
        "HB1": ("HA2", "", "", 0.09),
        "HB2": ("HA2", "", "", 0.09),
        "CG": ("CT2", "none", "H", -0.18),
        "HG1": ("HA2", "none", "", 0.09),
        "HG2": ("HA2", "none", "", 0.09),
        "CD": ("CT2", "none", "P", 0.20),
        "HD1": ("HA2", "", "", 0.09),
        "HD2": ("HA2", "", "", 0.09),
        "NE": ("NC2", "D", "", -0.70),
        "HE": ("HC", "HD", "", 0.44),
        "CZ": ("C", "none", "P", 0.64),
        "NH1": ("NC2", "D", "", -0.80),
        "HH11": ("HC", "HD", "", 0.46),
        "HH12": ("HC", "HD", "", 0.46),
        "NH2": ("NC2", "D", "", -0.80),
        "HH21": ("HC", "HD", "", 0.46),
        "HH22": ("HC", "HD", "", 0.46),
        "C": ("C", "none", "H", 0.51),
        "O": ("O", "A", "", -0.51)},
    'ASP': {
        "N": ("NH1", "D", "", -0.47),
        "HN": ("H", "HD", "", 0.31),
        "CA": ("CH1E", "none", "H", 0.07),
        "HA": ("HB1", "", "", 0.09),
        "CB": ("CH2E", "none", "H", -0.28),
        "HB1": ("HA2", "", "", 0.09),
        "HB2": ("HA2", "", "", 0.09),
        "CG": ("C", "none", "P", 0.62),
        "OD1": ("OC", "A", "", -0.76),
        "OD2": ("OC", "A", "", -0.76),
        "C": ("C", "none", "H", 0.51),
        "O": ("O", "A", "", -0.51)},
    'ASN': {
        "N": ("NH1", "D", "", -0.47),
        "HN": ("H", "HD", "", 0.31),
        "CA": ("CT1", "none", "H", 0.07),
        "HA": ("HB1", "", "", 0.09),
        "CB": ("CH2E", "none", "H", -0.18),
        "HB1": ("HA2", "", "", 0.09),
        "HB2": ("HA2", "", "", 0.09),
        "CG": ("C", "none", "P", 0.55),
        "OD1": ("O", "A", "", -0.55),
        "ND2": ("NH2", "D", "", -0.62),
        "HD21": ("H", "HD", "", 0.32),
        "HD22": ("H", "HD", "", 0.30),
        "C": ("C", "none", "H", 0.51),
        "O": ("O", "A", "", -0.51)},
    'CYS': {
        "N": ("NH1", "D", "", -0.47),
        "HN": ("H", "HD", "", 0.31),
        "CA": ("CT1", "none", "H", 0.10),
        "HA": ("HB1", "", "", 0.09),
        "CB": ("CH2E", "none", "P", 0.19),
        "HB1": ("HA2", "", "", 0.09),
        "HB2": ("HA2", "", "", 0.09),
        "SG": ("S", "D", "P", -0.19),
        "HG1": ("HS", "HD", "", 0.16),
        "C": ("C", "none", "H", 0.51),
        "O": ("O", "A", "", -0.51)},
    'GLN': {
        "N": ("NH1", "N", "", -0.47),
        "HN": ("H", "HN", "", 0.31),
        "CA": ("CT1", "none", "H", 0.07),
        "HA": ("HB1", "", "", 0.09),
        "CB": ("CT2", "none", "H", -0.18),
        "HB1": ("HA2", "", "", 0.09),
        "HB2": ("HA2", "", "", 0.09),
        "CG": ("CT2", "none", "H", -0.18),
        "HG1": ("HA2", "HD", "", 0.09),
        "HG2": ("HA2", "HD", "", 0.09),
        "CD": ("CC", "none", "P", 0.55),
        "OE1": ("O", "A", "", -0.55),
        "NE2": ("NH2", "D", "", -0.62),
        "HE21": ("H", "HD", "", 0.32),
        "HE22": ("H", "HD", "", 0.30),
        "C": ("C", "none", "H", 0.51),
        "O": ("O", "A", "", -0.51)},
    'GLU': {
        "N": ("NH1", "D", "", -0.47),
        "HN": ("H", "HD", "", 0.31),
        "CA": ("CT1", "none", "H", 0.07),
        "HA": ("HB1", "", "", 0.09),
        "CB": ("CT2A", "none", "H", -0.18),
        "HB1": ("HA2", "", "", 0.09),
        "HB2": ("HA2", "", "", 0.09),
        "CG": ("CT2", "none", "H", -0.28),
        "HG1": ("HA2", "HD", "", 0.09),
        "HG2": ("HA2", "HD", "", 0.09),
        "CD": ("CC", "none", "P", 0.62),
        "OE1": ("OC", "A", "", -0.76),
        "OE2": ("OC", "A", "", -0.76),
        "C": ("C", "none", "H", 0.51),
        "O": ("O", "A", "", -0.51)},
    'GLY': {
        "N": ("NH1", "D", "", -0.47),
        "HN": ("H", "HD", "", 0.31),
        "CA": ("CT2", "none", "H", -0.02),
        "HA1": ("HB2", "", "", 0.09),
        "HA2": ("HB2", "", "", 0.09),
        "C": ("C", "none", "H", 0.51),
        "O": ("O", "A", "", -0.51)},
    'HIS': {
        "N": ("NH1", "D", "", -0.35),
        "HN": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "HA": ("HB1", "", "", 0.09),
        "CB": ("CH2E", "none", "H", 0.00),
        "HB1": ("HA2", "", "", 0.09),
        "HB2": ("HA2", "", "", 0.09),
        "CG": ("C", "none", "P", 0.10),
        "ND1": ("NR", "B", "", -0.40),
        "HD1": ("H", "HD", "", 0.30),
        "CD2": ("CR1E", "none", "P", 0.30),      # Assumption that there is no hydrogen on ND1
        "NE2": ("NR", "B", "", -0.40),
        "HE2": ("H", "HD", "", 0.30),
        "CE1": ("CR1E", "none", "P", 0.30),      # Assumption that there is no hydrogen on NE2
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'HS2': { # CE1 HE1 ND1 CG CB HB1 HB2 NE2 HE2 CD2 HD2)',
        "CE1": ("", "", "", 0.25),
        "HE1": ("", "", "", 0.13),
        "ND1": ("NR", "B", "", -0.70),
        "CG": ("C", "none", "P", 0.22),
        "CB": ("CH2E", "none", "H", -0.08),
        "HB1": ("HA2", "", "", 0.09),
        "HB2": ("HA2", "", "", 0.09),
        "NE2": ("NR", "B", "", -0.36),
        "HE2": ("H", "HD", "", 0.32),
        "CD2": ("CR1E", "none", "P", -0.05),      # Assumption that there is no hydrogen on ND1
        "HD2": ("H", "HD", "", 0.09)},        
    'HSD': {    # HIS in force field
        "N": ("NH1", "D", "", -0.47),
        "HN": ("H", "HD", "", 0.31),
        "CA": ("CT1", "none", "H", 0.07),
        "HA": ("HB1", "", "", 0.09),
        "CB": ("CT2", "none", "H", -0.09),
        "HB1": ("HA2", "", "", 0.09),
        "HB2": ("HA2", "", "", 0.09),      
        "ND1": ("NR1", "D", "", -0.36),
        "HD1": ("H", "HD", "", 0.32),
        "CG": ("CPH1", "none", "P", -0.05),
        "CE1": ("CPH2", "none", "P", 0.25),
        "HE1": ("HR1", "HD", "", 0.13),
        "NE2": ("NR2", "A", "", -0.70),
        "CD2": ("CPH1", "none", "P", 0.22),
        "HD2": ("HR3", "HD", "", 0.10),
        "C": ("C", "none", "H", 0.51),
        "O": ("O", "A", "", -0.51)},
    'HSE': {    # HSD in force field
        "N": ("NH1", "D", "", -0.47),
        "HN": ("H", "HD", "", 0.31),
        "CA": ("CT1", "none", "H", 0.07),
        "HA": ("HB1", "", "", 0.09),
        "CB": ("CT2", "none", "H", -0.08),
        "HB1": ("HA2", "", "", 0.09),
        "HB2": ("HA2", "", "", 0.09),
        "CG": ("CPH1", "none", "P", 0.22),
        "ND1": ("NR2", "A", "", -0.70),
        "CD2": ("CPH1", "none", "P", -0.05),
        "HD2": ("HR3", "", "", 0.09),
        "NE2": ("NR1", "D", "", -0.36),
        "HE2": ("H", "HD", "", 0.32),
        "CE1": ("CPH2", "none", "P", 0.25),
        "HE1": ("HR1", "", "", 0.13),
        "C": ("C", "none", "H", 0.51),
        "O": ("O", "A", "", -0.51)},
    'HSP': {    # HSC in force field
        "N": ("NH1", "D", "", -0.47),
        "HN": ("H", "HD", "", 0.31),
        "CA": ("CT1", "none", "H", 0.07),
        "HA": ("HB1", "", "", 0.09),
        "CB": ("CT2A", "none", "H", -0.05),
        "HB1": ("HA2", "", "", 0.09),
        "HB2": ("HA2", "", "", 0.09),
        "CG": ("CPH1", "none", "P", 0.19),
        "ND1": ("NR3", "D", "", -0.51),
        "HD1": ("H", "HD", "", 0.44),
        "CD2": ("CPH1", "none", "P", 0.19),
        "HD2": ("HR1", "", "", 0.13),
        "NE2": ("NR3", "D", "", -0.51),
        "HE2": ("H", "HD", "", 0.44),
        "CE1": ("CPH2", "none", "P", 0.32),
        "HE1": ("HR2", "HD", "", 0.18),
        "C": ("C", "none", "H", 0.51),
        "O": ("O", "A", "", -0.51)},
    'ILE': {
        "N": ("NH1", "D", "", -0.47),
        "HN": ("H", "HD", "", 0.31),
        "CA": ("CT1", "none", "H", 0.07),
        "HA": ("HB1", "", "", 0.09),
        "CB": ("CT1", "none", "H", -0.09),
        "HB": ("HA1", "", "", 0.09),
        "CG1": ("CT2", "none", "H", -0.18),
        "HG11": ("HA2", "", "", 0.09),
        "HG12": ("HA2", "", "", 0.09),
        "CG2": ("CT3", "none", "H", -0.27),
        "HG21": ("HA3", "", "", 0.09),
        "HG22": ("HA3", "", "", 0.09),
        "HG23": ("HA3", "", "", 0.09),
        "CD": ("CT3", "none", "H", -0.27),
        "HD1": ("HA3", "", "", 0.09),
        "HD2": ("HA3", "", "", 0.09),
        "HD3": ("HA3", "", "", 0.09),
        "C": ("C", "none", "H", 0.51),
        "O": ("O", "A", "", -0.51)},
    'LEU': {
        "N": ("NH1", "D", "", -0.47),
        "HN": ("H", "HD", "", 0.31),
        "CA": ("CT1", "none", "H", 0.07),
        "HA": ("HB1", "", "", 0.09),
        "CB": ("CT2", "none", "H", -0.18),
        "HB1": ("HA2", "", "", 0.09),
        "HB2": ("HA2", "", "", 0.09),
        "CG": ("CT1", "none", "H", -0.09),
        "HG": ("HA1", "", "", 0.09),
        "CD1": ("CT3", "none", "H", -0.27),
        "HD11": ("HA3", "", "", 0.09),
        "HD12": ("HA3", "", "", 0.09),
        "HD13": ("HA3", "", "", 0.09),
        "CD2": ("CT3", "none", "H", -0.27),
        "HD21": ("HA3", "", "", 0.09),
        "HD22": ("HA3", "", "", 0.09),
        "HD23": ("HA3", "", "", 0.09),
        "C": ("C", "none", "H", 0.51),
        "O": ("O", "A", "", -0.51)},
    'LYS': {
        "N": ("NH1", "D", "", -0.47),
        "HN": ("H", "HD", "", 0.31),
        "CA": ("CT1", "none", "H", 0.07),
        "HA": ("HB1", "", "", 0.09),
        "CB": ("CT2", "none", "H", -0.18),
        "HB1": ("HA2", "", "", 0.09),
        "HB2": ("HA2", "", "", 0.09),
        "CG": ("CT2", "none", "H", -0.18),
        "HG1": ("HA2", "", "", 0.09),
        "HG2": ("HA2", "", "", 0.09),
        "CD": ("CT2", "none", "H", -0.18),
        "HD1": ("HA2", "", "", 0.09),
        "HD2": ("HA2", "", "", 0.09),
        "CE": ("CT2", "none", "P", 0.21),
        "HE1": ("HA2", "", "", 0.05),
        "HE2": ("HA2", "", "", 0.05),
        "NZ": ("NH3", "D", "", -0.30),
        "HZ1": ("HC", "HD", "", 0.33),
        "HZ2": ("HC", "HD", "", 0.33),
        "HZ3": ("HC", "HD", "", 0.33),
        "C": ("C", "none", "H", 0.51),
        "O": ("O", "A", "", -0.51)},
    'MET': {
        "N": ("NH1", "D", "", -0.47),
        "HN": ("H", "HD", "", 0.31),
        "CA": ("CT1", "none", "H", 0.07),
        "HA": ("HB1", "", "", 0.09),
        "CB": ("CT2", "none", "H", -0.18),
        "HB1": ("HA2", "", "", 0.09),
        "HB2": ("HA2", "", "", 0.09),
        "CG": ("CT2", "none", "P", -0.14),
        "HG1": ("HA2", "", "", 0.09),
        "HG2": ("HA2", "", "", 0.09),
        "SD": ("S", "A", "P", -0.09),
        "CE": ("CT3", "none", "P", -0.22),
        "HE1": ("HA3", "", "", 0.09),
        "HE2": ("HA3", "", "", 0.09),
        "HE3": ("HA3", "", "", 0.09),
        "C": ("C", "none", "H", 0.51),
        "O": ("O", "A", "", -0.51)},
    'PHE': {
        "N": ("NH1", "D", "", -0.47),
        "HN": ("H", "HD", "", 0.31),
        "CA": ("CT1", "none", "H", 0.07),
        "HA":("HB1", "", "", 0.09),
        "CB": ("CT2", "none", "H", -0.18),
        "HB1":("HA2", "", "",0.09),
        "HB2":("HA2", "", "",0.09),
        "CG": ("CA", "none", "H", 0.00),
        "CD1": ("CA", "none", "H", -0.115),
        "HD1":("HP", "", "",0.115),
        "CD2": ("CA", "none", "H", -0.115),
        "HD2":("HP", "", "",0.115),
        "CE1": ("CA", "none", "H", -0.115),
        "HE1":("HP", "", "",0.115),
        "CE2": ("CA", "none", "H", -0.115),
        "HE2":("HP", "", "",0.115),
        "CZ": ("CA", "none", "H", -0.115),
        "HZ":("HP", "", "",0.115),
        "C": ("C", "none", "H", 0.51),
        "O": ("O", "A", "", -0.51)},
    'PRO': {
        "N": ("N", "D", "", -0.29),
        "CA": ("CP1", "none", "H", 0.02),
        "HA":("HB1", "", "",0.09),
        "CB": ("CP2", "none", "H", -0.18),
        "HB1":("HA2", "", "",0.09),
        "HB2":("HA2", "", "",0.09),
        "CG": ("CP2", "none", "H", -0.18),
        "HG1":("HA2", "", "",0.09),
        "HG2":("HA2", "", "",0.09),
        "CD": ("CP3", "none", "H", 0.00),
        "HD1":("HA2", "", "",0.09),
        "HD2":("HA2", "", "",0.09),
        "C": ("C", "none", "H", 0.51),
        "O": ("O", "A", "", -0.51)},
    'SER': {
        "N": ("NH1", "D", "", -0.47),
        "HN": ("H", "HD", "", 0.31),
        "CA": ("CT1", "none", "H", 0.07),
        "HA":("HB1", "", "",0.09),
        "CB": ("CB2", "none", "P", 0.05),
        "HB1":("HA2", "", "",0.09),
        "HB2":("HA2", "", "",0.09),
        "OG": ("OH1", "B", "", -0.66),
        "HG1": ("H", "HD", "", 0.43),
        "C": ("C", "none", "H", 0.51),
        "O": ("O", "A", "", -0.51)},
    'THR': {
        "N": ("NH1", "D", "", -0.47),
        "HN": ("H", "HD", "", 0.31),
        "CA": ("CT1", "none", "H", 0.07),
        "HA": ("HB1", "", "",0.09),
        "CB": ("CT1", "none", "P", 0.14),
        "HB": ("HA1", "", "",0.09),
        "OG1": ("OH1", "B", "", -0.66),
        "OXT": (),                                #need to confirm what it is...
        "HG1": ("H", "HD", "", 0.43),
        "CG2": ("CT3", "none", "H", -0.27),
        "HG21": ("HA3", "", "",0.09),
        "HG22": ("HA3", "", "",0.09),
        "HG23": ("HA3", "", "",0.09),
        "C": ("C", "none", "H", 0.51),
        "O": ("O", "A", "", -0.51)},
    'TRP': {
        "N": ("NH1", "D", "", -0.47),
        "HN": ("H", "HD", "", 0.31),
        "CA": ("CT1", "none", "H", 0.07),
        "HA": ("HB1", "", "",0.09),
        "CB": ("CT2", "none", "H", -0.18),
        "HB1": ("HA2", "", "",0.09),
        "HB2": ("HA2", "", "",0.09),
        "CG": ("CY", "none", "H", -0.03),
        "CD2": ("CPT", "none", "H", 0.11),
        "CE2": ("CPT", "none", "P", 0.24),
        "CE3": ("CAI", "none", "H", -0.25),
        "HE3": ("HP", "HD", "", 0.17),
        "CD1": ("CA", "none", "P", -0.15),
        "HD1": ("HP", "", "",0.22),
        "NE1": ("NY", "D", "", -0.51),
        "HE1": ("H", "HD", "", 0.37),
        "CZ2": ("CAI", "none", "H", -0.27),
        "HZ2": ("HP", "HD", "", 0.16),
        "CZ3": ("CA", "none", "H", -0.20),
        "HZ3": ("HP", "HD", "", 0.14),
        "CH2": ("CA", "none", "H", -0.14),
        "HH2": ("HP", "", "", 0.14),
        "C": ("C", "none", "H", 0.51),
        "O": ("O", "A", "", -0.51)},
    'TYR': {
        "N": ("NH1", "D", "", -0.47),
        "HN": ("H", "HD", "", 0.31),
        "CA": ("CT1", "none", "H", 0.07),
        "HA": ("HB1", "", "",0.09),
        "CB": ("CT2", "none", "H", -0.18),
        "HB1": ("HA2", "", "",0.09),
        "HB2": ("HA2", "", "",0.09),
        "CG": ("CA", "none", "H", 0.00),
        "CD1": ("CA", "none", "H", -0.115),
        "HD1": ("HP", "", "", 0.115),
        "CD2": ("CA", "none", "H", -0.115),
        "HD2": ("HP", "HD", "", 0.115),
        "CE1": ("CA", "none", "H", -0.115),
        "HE1": ("HP", "HD", "", 0.115),
        "CE2": ("CA", "none", "H", -0.115),
        "HE2": ("HP", "HD", "", 0.115),
        "CZ": ("CA", "none", "P", 0.11),
        "OH": ("OH1", "B", "", -0.54),
        "HH": ("H", "HD", "", 0.43),
        "C": ("C", "none", "H", 0.51),
        "O": ("O", "A", "", -0.51)},
    'VAL': {
        "N": ("NH1", "D", "", -0.47),
        "HN": ("H", "HD", "", 0.31),
        "CA": ("CT1", "none", "H", 0.07),
        "HA": ("HB1", "", "",0.09),
        "CB": ("CT1", "none", "H", -0.09),
        "HB": ("HA1", "", "",0.09),
        "CG1": ("CT3", "none", "H", -0.27),
        "HG11": ("HA3", "", "", 0.09),
        "HG12": ("HA3", "", "", 0.09),
        "HG13": ("HA3", "", "", 0.09),
        "CG2": ("CT3", "none", "H", -0.27),
        "HG21": ("HA3", "", "", 0.09),
        "HG22": ("HA3", "", "", 0.09),
        "HG23": ("HA3", "", "", 0.09),
        "C": ("C", "none", "H", 0.51),
        "O": ("O", "A", "", -0.51)},
    'NTER': {
        "N": ("NH3", "D", "", -0.30),
        "HT1": ("HC", "HD", "", 0.33),
        "HT2": ("HC", "HD", "", 0.33),
        "HT3": ("HC", "HD", "", 0.33),
        "CA": ("CT1", "HD", "H", 0.21),
        "HA": ("HB1", "", "", 0.10)},
    'CTER': {
        "C": ("CC", "none", "H", 0.34),
        "OT1": ("OC", "A", "", -0.67),
        "OT2": ("OC", "A", "", -0.67)}
}


non_amino_info = {
    'MG': {"MG": ("Met", "D", "", -2.0)},
    'ZN': {"ZN": ("Met", "D", "", -2.0)}
}

def get_residue_ring_atoms_unity_model(resname):
    ring = []
    if resname == "PHE":
        ring = ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]
        return ring
    elif resname == "TYR":
        ring = ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]
        return ring
    elif resname == "TRP":
        ring1 = ["CG", "CD1", "CD2", "NE1", "CE2"]
        ring2 = ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"]
        return ring1, ring2
    elif resname == "HIS" or resname == "HSD" or resname == "HSE" or resname == "HSP":
        ring = ["CG", "CD2", "ND1", "CE1", "NE2"]
        return ring
    else:
        return ring


def get_residue_main_and_side_chain_unity_model(resname):
    mainchain_atoms = ["N", "HN", "C", "CA", "O", "OT1", "OT2", "HT1", "HT2", "HT3"]
    mainchain_atoms += ["ZN", "MG", "Ca"]
    sidechain_atoms = []
    if resname == "ALA":
        sidechain_atoms = ["CB"]
    if resname == "ARG":
        sidechain_atoms = ["CB", "CG", "CD", "NE", "HE", "CZ", "NH1", "HH11", "HH12", "NH2", "HH21", "HH22"]
    if resname == "ASP":
        sidechain_atoms = ["CB", "CG", "OD1", "OD2"]
    if resname == "ASN":
        sidechain_atoms = ["CB", "CG", "OD1", "ND2", "HD21", "HD22"]
    if resname == "CYS":
        sidechain_atoms = ["CB", "SG"]
    if resname == "GLN":
        sidechain_atoms = ["CB", "CG", "CD", "OE1", "NE2", "HE21", "HE22"]
    if resname == "GLU":
        sidechain_atoms = ["CB", "CG", "CD", "OE1", "OE2"]
    if resname == "GLY":
        sidechain_atoms = []
    if resname == "HIS" or resname == "HSD" or resname == "HSE" or resname == "HSP":
        sidechain_atoms = ["CB", "CG", "ND1", "HD1", "CD2", "NE2", "HE2", "CE1"]
    if resname == "ILE":
        sidechain_atoms = ["CB", "CG1", "CG2", "CD"]
    if resname == "LEU":
        sidechain_atoms = ["CB", "CG1", "CD1", "CD2"]
    if resname == "LYS":
        sidechain_atoms = ["CB", "CG", "CD", "CE", "NZ", "HZ1", "HZ2", "HZ3"]
    if resname == "MET":
        sidechain_atoms = ["CB", "CG", "SD", "CE"]
    if resname == "PHE":
        sidechain_atoms = ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"]
    if resname == "PRO":
        sidechain_atoms = ["CB", "CG", "CD"]
    if resname == "SER":
        sidechain_atoms = ["CB", "OG", "HG1"]
    if resname == "THR":
        sidechain_atoms = ["CB", "OG1", "HG1", "CG2"]
    if resname == "TRP":
        sidechain_atoms = ["CB", "CG", "CD2", "CE2", "CE3", "CD1", "NE1", "HE1", "CZ2", "CZ3", "CH2"]
    if resname == "TYR":
        sidechain_atoms = ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH", "HH"]
    if resname == "VAL":
        sidechain_atoms = ["CB", "CG1", "CG2"]

    return mainchain_atoms, sidechain_atoms


def vdw_radius(atomname):
    atomname = guess_atom_name(atomname)

    if atomname == "H": return 1.0
    if atomname == "C": return 1.9
    if atomname == "N": return 1.8
    if atomname == "O": return 1.7
    if atomname == "F": return 1.5
    if atomname == "P": return 2.1
    if atomname == "S": return 2.0
    if atomname == "CL": return 1.8
    if atomname == "BR": return 2.0
    if atomname == "I": return 2.2
    if atomname in ["MG", "ZN", "CA"]: return 1.2