from ..guess import guess_atom_name


unity_model_atom_selection_old = {
    'ALA': '(name N or name CB or name O)',
    'ARG': '(name N or name NE or name NH1 or name NH2 or name O)',
    'ASP': '(name N or name OD1 or name OD2 or name O)',
    'ASN': '(name N or name OD1 or name ND2 or name O)',
    'CYS': '(name N or name SG or name O)',
    'GLN': '(name N or name OE1 or name NE2 or name O)',
    'GLU': '(name N or name OE1 or name OE2 or name O)',
    'GLY': '(name N or name O)',
    'HIS': '(name N or name CG or name ND1 or name CE1 or name NE2 or name CD2 or name O)',
    'HSD': '(name N or name CG or name ND1 or name CE1 or name NE2 or name CD2 or name O)',
    'HSE': '(name N or name CG or name ND1 or name CE1 or name NE2 or name CD2 or name O)',
    'HSP': '(name N or name CG or name ND1 or name CE1 or name NE2 or name CD2 or name O)',
    'ILE': '(name N or name CG2 or name CD or name O)',
    'LEU': '(name N or name CG or name CD1 or name CD2 or name O)',
    'LYS': '(name N or name NZ or name O)',
    'MET': '(name N or name SD or name CE or name O)',
    'PHE': '(name N or name CD1 or name CE1 or name CZ or name CE2 or name CD2 or name CG or name O)',
    'PRO': '(name N or name CA or name CB or name CD or name CG or name O)',
    'SER': '(name N or name OG or name O)',
    'THR': '(name N or name OG1 or name OG2 or name O)',
    'TRP': '(name N or name CG or name CD1 or name NE1 or name CE2 or name CD2 or name CE3 or name CZ3 or name CH2 or'
           ' name CZ2 or name O)',
    'TYR': '(name N or name CG or name CD1 or name CE1 or name CZ or name CE2 or name CD2 or name OH or name O)',
    'VAL': '(name N or name CG1 or name CG2  or name O)'
}

unity_model_atom_selection = {
    'ALA': '(name N or name HN or name CA or name C or name O '      # Backbone
           'or name CB)',
    'ARG': '(name N or name HN or name CA or name C or name O '      # Backbone
           'or name CB or name CG or name CD or name NE or name HE or name CZ or name NH1 or name HNH11 or name HNH12 '
           'or name NH2 or name HNH21 or name HNH22 or name O)',
    'ASP': '(name N or name HN or name CA or name C or name O '      # Backbone
           'or name CB or name CG or name OD1 or name OD2)',
    'ASN': '(name N or name HN or name CA or name C or name O '      # Backbone
           'or name CB or name CG or name OD1 or name ND2 or name HND21 or name HND22)',
    'CYS': '(name N or name HN or name CA or name C or name O '      # Backbone
           'or name CB or name SG)',
    'GLN': '(name N or name HN or name CA or name C or name O '      # Backbone
           'or name CB or name CG or name CD or name OE1 or name NE2 or name HNE21 or name HNE22)',
    'GLU': '(name N or name HN or name CA or name C or name O '      # Backbone
           'or name CB or name CG or name CD or name OE1 or name OE2)',
    'GLY': '(name N or name HN or name CA or name C or name O)',     # Backbone
    'HIS': '(name N or name HN or name CA or name C or name O '      # Backbone    # HIS is equal HSP
           'or name CB or name CG or name ND1 or name HND1 or name CD2 or name NE2 or name HNE2 or name CE1)',
    'HSD': '(name N or name HN or name CA or name C or name O '      # Backbone
           'or name CB or name CG or name ND1 or name HND1 or name CD2 or name NE2 or name CE1)',
    'HSE': '(name N or name HN or name CA or name C or name O '      # Backbone
           'or name CB or name CG or name ND1 or name CD2 or name NE2 or name HNE2 or name CE1)',
    'HSP': '(name N or name HN or name CA or name C or name O '      # Backbone
           'or name CB or name CG or name ND1 or name HND1 or name CD2 or name NE2 or name HNE2 or name CE1)',
    'ILE': '(name N or name HN or name CA or name C or name O '      # Backbone
           'or name CB or name CG2 or name CG1 or name CD)',
    'LEU': '(name N or name HN or name CA or name C or name O '      # Backbone
           'or name CB or name CG or name CD1 or name CD2)',
    'LYS': '(name N or name HN or name CA or name C or name O '      # Backbone
           'or name CB or name CG or name CD or name CE or name NZ or name HNZ1 or name HNZ2 or name HNZ3)',
    'MET': '(name N or name HN or name CA or name C or name O '      # Backbone
           'or name CB or name CG or name SD or name CE)',
    'PHE': '(name N or name HN or name CA or name C or name O '      # Backbone
           'or name CB or name CG or name CD1 or name CD2 or name CE1 or name CE2 or name CZ)',
    'PRO': '(name N or name CA or name C or name O '                # Backbone
           'or name CB or name CG or name CD)',
    'SER': '(name N or name HN or name CA or name C or name O '      # Backbone
           'or name CB or name OG or name HG1)',
    'THR': '(name N or name HN or name CA or name C or name O '      # Backbone
           'or name CB or name OG1 or name HNG1 or name CG2)',
    'TRP': '(name N or name HN or name CA or name C or name O '      # Backbone
           'or name CB or name CG or name CD2 or name CE2 or name CE3 or name CD1 or name NE1 or name HNE1 '
           'or name CZ2 or name CZ3 or name CH2)',
    'TYR': '(name N or name HN or name CA or name C or name O '      # Backbone
           'or name CB or name CG or name CD1 or name CE1 or name CD2 or name CE2 or name CZ or name OH or name HNH)',
    'VAL': '(name N or name HN or name CA or name C or name O '      # Backbone
           'or name CB or name CG1 or name CG2)'
}

# CHARMM c41b1 toph19 force field
unity_model_info_old = {
    'ALA': {
        "N": ("NH1", "D", "", -0.35),       # (atom_type, hbond_type, charge)
        "H": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH3E", "none", "H", 0.00),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'ARG': {
        "N": ("NH1", "D", "", -0.35),
        "H": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", 0.00),
        "CG": ("CH2E", "none", "H", 0.00),
        "CD": ("CH2E", "none", "P", 0.10),
        "NE": ("NH1", "D", "", -0.40),
        "HE": ("H", "HD", "", 0.30),
        "CZ": ("C", "none", "P", 0.50),
        "NH1": ("NC2", "D", "", -0.45+0.70),
        "HH11": ("HC", "HD", "", 0.0),
        "HH12": ("HC", "HD", "", 0.0),
        "NH2": ("NC2", "D", "", -0.45+0.70),
        "HH21": ("HC", "HD", "", 0.0),
        "HH22": ("HC", "HD", "", 0.0),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'ASP': {
        "N": ("NH1", "D", "", -0.35),
        "H": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", -0.16),
        "CG": ("C", "none", "P", 0.36),
        "OD1": ("OC", "A", "", -0.60),
        "OD2": ("OC", "A", "", -0.60),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'ASN': {
        "N": ("NH1", "D", "", -0.35),
        "H": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", 0.00),
        "CG": ("C", "none", "P", 0.55),
        "OD1": ("O", "A", "", -0.55),
        "ND2": ("NH2", "D", "", -0.60),
        "HD21": ("H", "HD", "", 0.30),
        "HD22": ("H", "HD", "", 0.30),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'CYS': {
        "N": ("NH1", "D", "", -0.35),
        "H": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "P", 0.19),
        "SG": ("SH1E", "D", "P", -0.19),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'GLN': {
        "N": ("NH1", "D", "", -0.35),
        "H": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", 0.00),
        "CG": ("CH2E", "none", "H", 0.00),
        "CD": ("C", "none", "P", 0.55),
        "OE1": ("O", "A", "", -0.55),
        "NE2": ("NH2", "D", "", -0.60),
        "HE21": ("H", "HD", "", 0.30),
        "HE22": ("H", "HD", "", 0.30),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'GLU': {
        "N": ("NH1", "D", "", -0.35),
        "H": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", 0.00),
        "CG": ("CH2E", "none", "H", -0.16),
        "CD": ("C", "none", "P", 0.36),
        "OE1": ("OC", "A", "", -0.60),
        "OE2": ("OC", "A", "", -0.60),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'GLY': {
        "N": ("NH1", "D", "", -0.35),
        "H": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'HIS': {
        "N": ("NH1", "D", "", -0.35),
        "H": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", 0.00),
        "CG": ("C", "none", "P", 0.10),
        "ND1": ("NR", "B", "", -0.40),
        "HD1": ("H", "HD", "", 0.30),
        "CD2": ("CR1E", "none", "P", 0.30),      # Assumption that there is no hydrogen on ND1
        "NE2": ("NR", "B", "", -0.40),
        "HE2": ("H", "HD", "", 0.30),
        "CE1": ("CR1E", "none", "P", 0.30),      # Assumption that there is no hydrogen on NE2
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'HSD': {    # HIS in force field
        "N": ("NH1", "D", "", -0.35),
        "H": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", 0.00),
        "CG": ("C", "none", "P", 0.10),
        "ND1": ("NH1", "D", "", -0.40),
        "HD1": ("H", "HD", "", 0.30),
        "CD2": ("CR1E", "none", "P", 0.10),
        "NE2": ("NR", "A", "", -0.40),
        "CE1": ("CR1E", "none", "P", 0.30),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'HSE': {    # HSD in force field
        "N": ("NH1", "D", "", -0.35),
        "H": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", 0.00),
        "CG": ("C", "none", "P", 0.10),
        "ND1": ("NR", "A", "", -0.40),
        "CD2": ("CR1E", "none", "P", 0.10),
        "NE2": ("NH1", "D", "", -0.40),
        "HE2": ("H", "HD", "", 0.30),
        "CE1": ("CR1E", "none", "P", 0.30),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'HSP': {    # HSC in force field
        "N": ("NH1", "D", "", -0.35),
        "H": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", 0.00),
        "CG": ("C", "none", "P", 0.15),
        "ND1": ("NH1", "D", "", -0.30+0.35),
        "HD1": ("H", "HD", "", 0.35),
        "CD2": ("CR1E", "none", "P", 0.20),
        "NE2": ("NH1", "D", "", -0.30+0.35),
        "HE2": ("H", "HD", "", 0.35),
        "CE1": ("CR1E", "none", "P", 0.45),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'ILE': {
        "N": ("NH1", "D", "", -0.35),
        "H": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH1E", "none", "H", 0.00),
        "CG1": ("CH3E", "none", "H", 0.00),
        "CG2": ("CH3E", "none", "H", 0.00),
        "CD": ("CH3E", "none", "H", 0.00),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'LEU': {
        "N": ("NH1", "D", "", -0.35),
        "H": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", 0.00),
        "CG": ("CH1E", "none", "H", 0.00),
        "CD1": ("CH3E", "none", "H", 0.00),
        "CD2": ("CH3E", "none", "H", 0.00),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'LYS': {
        "N": ("NH1", "D", "", -0.35),
        "H": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", 0.00),
        "CG": ("CH2E", "none", "H", 0.00),
        "CD": ("CH2E", "none", "H", 0.00),
        "CE": ("CH2E", "none", "P", 0.25),
        "NZ": ("NH3", "D", "", -0.30+0.70),
        "HZ1": ("HC", "HD", "", 0.35),
        "HZ2": ("HC", "HD", "", 0.35),
        "HZ3": ("HC", "HD", "", 0.35),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'MET': {
        "N": ("NH1", "D", "", -0.35),
        "H": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", 0.00),
        "CG": ("CH2E", "none", "P", 0.06),
        "SD": ("S", "A", "P", -0.12),
        "CE": ("CH3E", "none", "P", 0.06),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'PHE': {
        "N": ("NH1", "D", "", -0.35),
        "H": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", 0.00),
        "CG": ("C", "none", "H", 0.00),
        "CD1": ("CR1E", "none", "H", 0.00),
        "CD2": ("CR1E", "none", "H", 0.00),
        "CE1": ("CR1E", "none", "H", 0.00),
        "CE2": ("CR1E", "none", "H", 0.00),
        "CZ": ("CR1E", "none", "H", 0.00),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'PRO': {
        "N": ("N", "D", "", -0.20),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", 0.00),
        "CG": ("CH2E", "none", "H", 0.00),
        "CD": ("CH2E", "none", "H", 0.10),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'SER': {
        "N": ("NH1", "D", "", -0.35),
        "H": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "P", 0.25),
        "OG": ("OH1", "B", "", -0.65),
        "HG": ("H", "HD", "", 0.40),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'THR': {
        "N": ("NH1", "D", "", -0.35),
        "H": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH1E", "none", "P", 0.25),
        "OG1": ("OH1", "B", "", -0.65),
        "HG1": ("H", "HD", "", 0.40),
        "CG2": ("CH3E", "none", "H", 0.00),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'TRP': {
        "N": ("NH1", "D", "", -0.35),
        "H": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", 0.00),
        "CG": ("C", "none", "H", -0.03),
        "CD2": ("C", "none", "H", 0.10),
        "CE2": ("C", "none", "P", -0.04),
        "CE3": ("CR1E", "none", "H", -0.03),
        "CD1": ("CR1E", "none", "P", 0.06),
        "NE1": ("NH1", "D", "", -0.36),
        "HE1": ("H", "HD", "", 0.30),
        "CZ2": ("CR1E", "none", "H", -0.00),
        "CZ3": ("CR1E", "none", "H", -0.00),
        "CH2": ("CR1E", "none", "H", -0.00),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'TYR': {
        "N": ("NH1", "D", "", -0.35),
        "H": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", 0.00),
        "CG": ("C", "none", "H", 0.00),
        "CD1": ("CR1E", "none", "H", 0.00),
        "CD2": ("CR1E", "none", "H", 0.00),
        "CE1": ("CR1E", "none", "H", 0.00),
        "CE2": ("CR1E", "none", "H", 0.00),
        "CZ": ("CR1E", "none", "P", 0.25),
        "OH": ("OH1", "B", "", -0.65),
        "HH": ("H", "HD", "", 0.40),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'VAL': {
        "N": ("NH1", "D", "", -0.35),
        "H": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH1E", "none", "H", 0.00),
        "CG1": ("CH3E", "none", "H", 0.00),
        "CG2": ("CH3E", "none", "H", 0.00),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'NTER': {
        "N": ("NH3", "D", "", -0.30),
        "HT1": ("HC", "HD", "", 0.35),
        "HT2": ("HC", "HD", "", 0.35),
        "HT3": ("HC", "HD", "", 0.35),
        "CA": ("CH1E", "HD", "H", 0.25)},
    'CTER': {
        "C": ("C", "none", "H", 0.14),
        "OT1": ("OC", "B", "", -0.57),
        "OT2": ("OC", "B", "", -0.57)
    }
}

# CHARMM c41b1 toph19 force field
unity_model_info = {
    'ALA': {
        "N": ("NH1", "D", "", -0.35),       # (atom_type, hbond_type, charge)
        "HN": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH3E", "none", "H", 0.00),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'ARG': {
        "N": ("NH1", "D", "", -0.35),
        "HN": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", 0.00),
        "CG": ("CH2E", "none", "H", 0.00),
        "CD": ("CH2E", "none", "P", 0.10),
        "NE": ("NH1", "D", "", -0.40),
        "HE": ("H", "HD", "", 0.30),
        "CZ": ("C", "none", "P", 0.50),
        "NH1": ("NC2", "D", "", -0.45+0.70),
        "HH11": ("HC", "HD", "", 0.0),
        "HH12": ("HC", "HD", "", 0.0),
        "NH2": ("NC2", "D", "", -0.45+0.70),
        "HH21": ("HC", "HD", "", 0.0),
        "HH22": ("HC", "HD", "", 0.0),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'ASP': {
        "N": ("NH1", "D", "", -0.35),
        "HN": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", -0.16),
        "CG": ("C", "none", "P", 0.36),
        "OD1": ("OC", "A", "", -0.60),
        "OD2": ("OC", "A", "", -0.60),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'ASN': {
        "N": ("NH1", "D", "", -0.35),
        "HN": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", 0.00),
        "CG": ("C", "none", "P", 0.55),
        "OD1": ("O", "A", "", -0.55),
        "ND2": ("NH2", "D", "", -0.60),
        "HD21": ("H", "HD", "", 0.30),
        "HD22": ("H", "HD", "", 0.30),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'CYS': {
        "N": ("NH1", "D", "", -0.35),
        "HN": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "P", 0.19),
        "SG": ("SH1E", "D", "P", -0.19),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'GLN': {
        "N": ("NH1", "D", "", -0.35),
        "HN": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", 0.00),
        "CG": ("CH2E", "none", "H", 0.00),
        "CD": ("C", "none", "P", 0.55),
        "OE1": ("O", "A", "", -0.55),
        "NE2": ("NH2", "D", "", -0.60),
        "HE21": ("H", "HD", "", 0.30),
        "HE22": ("H", "HD", "", 0.30),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'GLU': {
        "N": ("NH1", "D", "", -0.35),
        "HN": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", 0.00),
        "CG": ("CH2E", "none", "H", -0.16),
        "CD": ("C", "none", "P", 0.36),
        "OE1": ("OC", "A", "", -0.60),
        "OE2": ("OC", "A", "", -0.60),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'GLY': {
        "N": ("NH1", "D", "", -0.35),
        "HN": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'HIS': {
        "N": ("NH1", "D", "", -0.35),
        "HN": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", 0.00),
        "CG": ("C", "none", "P", 0.10),
        "ND1": ("NR", "B", "", -0.40),
        "HD1": ("H", "HD", "", 0.30),
        "CD2": ("CR1E", "none", "P", 0.30),      # Assumption that there is no hydrogen on ND1
        "NE2": ("NR", "B", "", -0.40),
        "HE2": ("H", "HD", "", 0.30),
        "CE1": ("CR1E", "none", "P", 0.30),      # Assumption that there is no hydrogen on NE2
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'HSD': {    # HIS in force field
        "N": ("NH1", "D", "", -0.35),
        "HN": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", 0.00),
        "CG": ("C", "none", "P", 0.10),
        "ND1": ("NH1", "D", "", -0.40),
        "HD1": ("H", "HD", "", 0.30),
        "CD2": ("CR1E", "none", "P", 0.10),
        "NE2": ("NR", "A", "", -0.40),
        "CE1": ("CR1E", "none", "P", 0.30),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'HSE': {    # HSD in force field
        "N": ("NH1", "D", "", -0.35),
        "HN": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", 0.00),
        "CG": ("C", "none", "P", 0.10),
        "ND1": ("NR", "A", "", -0.40),
        "CD2": ("CR1E", "none", "P", 0.10),
        "NE2": ("NH1", "D", "", -0.40),
        "HE2": ("H", "HD", "", 0.30),
        "CE1": ("CR1E", "none", "P", 0.30),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'HSP': {    # HSC in force field
        "N": ("NH1", "D", "", -0.35),
        "HN": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", 0.00),
        "CG": ("C", "none", "P", 0.15),
        "ND1": ("NH1", "D", "", -0.30+0.35),
        "HD1": ("H", "HD", "", 0.35),
        "CD2": ("CR1E", "none", "P", 0.20),
        "NE2": ("NH1", "D", "", -0.30+0.35),
        "HE2": ("H", "HD", "", 0.35),
        "CE1": ("CR1E", "none", "P", 0.45),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'ILE': {
        "N": ("NH1", "D", "", -0.35),
        "HN": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH1E", "none", "H", 0.00),
        "CG1": ("CH3E", "none", "H", 0.00),
        "CG2": ("CH3E", "none", "H", 0.00),
        "CD": ("CH3E", "none", "H", 0.00),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'LEU': {
        "N": ("NH1", "D", "", -0.35),
        "HN": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", 0.00),
        "CG": ("CH1E", "none", "H", 0.00),
        "CD1": ("CH3E", "none", "H", 0.00),
        "CD2": ("CH3E", "none", "H", 0.00),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'LYS': {
        "N": ("NH1", "D", "", -0.35),
        "HN": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", 0.00),
        "CG": ("CH2E", "none", "H", 0.00),
        "CD": ("CH2E", "none", "H", 0.00),
        "CE": ("CH2E", "none", "P", 0.25),
        "NZ": ("NH3", "D", "", -0.30+0.70),
        "HZ1": ("HC", "HD", "", 0.35),
        "HZ2": ("HC", "HD", "", 0.35),
        "HZ3": ("HC", "HD", "", 0.35),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'MET': {
        "N": ("NH1", "D", "", -0.35),
        "HN": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", 0.00),
        "CG": ("CH2E", "none", "P", 0.06),
        "SD": ("S", "A", "P", -0.12),
        "CE": ("CH3E", "none", "P", 0.06),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'PHE': {
        "N": ("NH1", "D", "", -0.35),
        "HN": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", 0.00),
        "CG": ("C", "none", "H", 0.00),
        "CD1": ("CR1E", "none", "H", 0.00),
        "CD2": ("CR1E", "none", "H", 0.00),
        "CE1": ("CR1E", "none", "H", 0.00),
        "CE2": ("CR1E", "none", "H", 0.00),
        "CZ": ("CR1E", "none", "H", 0.00),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'PRO': {
        "N": ("N", "D", "", -0.20),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", 0.00),
        "CG": ("CH2E", "none", "H", 0.00),
        "CD": ("CH2E", "none", "H", 0.10),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'SER': {
        "N": ("NH1", "D", "", -0.35),
        "HN": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "P", 0.25),
        "OG": ("OH1", "B", "", -0.65),
        "HG1": ("H", "HD", "", 0.40),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'THR': {
        "N": ("NH1", "D", "", -0.35),
        "HN": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH1E", "none", "P", 0.25),
        "OG1": ("OH1", "B", "", -0.65),
        "HG1": ("H", "HD", "", 0.40),
        "CG2": ("CH3E", "none", "H", 0.00),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'TRP': {
        "N": ("NH1", "D", "", -0.35),
        "HN": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", 0.00),
        "CG": ("C", "none", "H", -0.03),
        "CD2": ("C", "none", "H", 0.10),
        "CE2": ("C", "none", "P", -0.04),
        "CE3": ("CR1E", "none", "H", -0.03),
        "CD1": ("CR1E", "none", "P", 0.06),
        "NE1": ("NH1", "D", "", -0.36),
        "HE1": ("H", "HD", "", 0.30),
        "CZ2": ("CR1E", "none", "H", -0.00),
        "CZ3": ("CR1E", "none", "H", -0.00),
        "CH2": ("CR1E", "none", "H", -0.00),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'TYR': {
        "N": ("NH1", "D", "", -0.35),
        "HN": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH2E", "none", "H", 0.00),
        "CG": ("C", "none", "H", 0.00),
        "CD1": ("CR1E", "none", "H", 0.00),
        "CD2": ("CR1E", "none", "H", 0.00),
        "CE1": ("CR1E", "none", "H", 0.00),
        "CE2": ("CR1E", "none", "H", 0.00),
        "CZ": ("CR1E", "none", "P", 0.25),
        "OH": ("OH1", "B", "", -0.65),
        "HH": ("H", "HD", "", 0.40),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'VAL': {
        "N": ("NH1", "D", "", -0.35),
        "HN": ("H", "HD", "", 0.25),
        "CA": ("CH1E", "none", "H", 0.10),
        "CB": ("CH1E", "none", "H", 0.00),
        "CG1": ("CH3E", "none", "H", 0.00),
        "CG2": ("CH3E", "none", "H", 0.00),
        "C": ("C", "none", "H", 0.55),
        "O": ("O", "A", "", -0.55)},
    'NTER': {
        "N": ("NH3", "D", "", -0.30),
        "HT1": ("HC", "HD", "", 0.35),
        "HT2": ("HC", "HD", "", 0.35),
        "HT3": ("HC", "HD", "", 0.35),
        "CA": ("CH1E", "HD", "H", 0.25)},
    'CTER': {
        "C": ("C", "none", "H", 0.14),
        "OT1": ("OC", "B", "", -0.57),
        "OT2": ("OC", "B", "", -0.57)
    }
}

cg_model_atom_selection = {
    'ALA': '(name N or name CB or name O)',
    'ARG': '(name N or name NE or name NH1 or name NH2 or name O)',
    'ASP': '(name N or name OD1 or name OD2 or name O)',
    'ASN': '(name N or name OD1 or name ND2 or name O)',
    'CYS': '(name N or name SG or name O)',
    'GLN': '(name N or name OE1 or name NE2 or name O)',
    'GLU': '(name N or name OE1 or name OE2 or name O)',
    'GLY': '(name N or name O)',
    'HIS': '(name N or name CG or name ND1 or name CE1 or name NE2 or name CD2 or name O)',
    'HSD': '(name N or name CG or name ND1 or name CE1 or name NE2 or name CD2 or name O)',
    'HSE': '(name N or name CG or name ND1 or name CE1 or name NE2 or name CD2 or name O)',
    'HSP': '(name N or name CG or name ND1 or name CE1 or name NE2 or name CD2 or name O)',
    'ILE': '(name N or name CG2 or name CD or name O)',
    'LEU': '(name N or name CG or name CD1 or name CD2 or name O)',
    'LYS': '(name N or name NZ or name O)',
    'MET': '(name N or name SD or name CE or name O)',
    'PHE': '(name N or name CD1 or name CE1 or name CZ or name CE2 or name CD2 or name CG or name O)',
    'PRO': '(name N or name CA or name CB or name CD or name CG or name O)',
    'SER': '(name N or name OG or name O)',
    'THR': '(name N or name OG1 or name OG2 or name O)',
    'TRP': '(name N or name CG or name CD1 or name NE1 or name CE2 or name CD2 or name CE3 or name CZ3 or name CH2 or'
           ' name CZ2 or name O)',
    'TYR': '(name N or name CG or name CD1 or name CE1 or name CZ or name CE2 or name CD2 or name OH or name O)',
    'VAL': '(name N or name CG1 or name CG2  or name O)'
}


# CHARMM c41b1 top_all36_prot force field
cg_model_info = {
    'ALA': [("NH1", "D", -0.47), ("CT3", "none", -0.27), ("O", "A", -0.51)],
    'ARG': [("NH1", "D", -0.47), ("NC2", "D", -0.70), ("NC2", "D", -0.80), ("NC2", "D", -0.80),
            ("O", "A", -0.51)],
    'ASP': [("NH1", "D", -0.47), ("OC", "A", -0.76), ("OC", "A", -0.76), ("O", "A", -0.51)],
    'ASN': [("NH1", "D", -0.47), ("O", "A", -0.55), ("NH2", "D", -0.62), ("O", "A", -0.51)],
    'CYS': [("NH1", "D", -0.47), ("S", "D", -0.23), ("O", "A", -0.51)],
    'GLN': [("NH1", "D", -0.47), ("O", "A", -0.55), ("NH2", "D", -0.62), ("O", "A", -0.51)],
    'GLU': [("NH1", "D", -0.47), ("OC", "A", -0.76), ("OC", "A", -0.76), ("O", "A", -0.51)],
    'GLY': [("NH1", "D", -0.47), ("O", "A", -0.51)],
    'HIS': [("NH1", "D", -0.47), ("NR1", "B", -0.36), ("CPH1", "none", -0.05), ("CPH2", "none", -0.25),
            ("NR2", "B", -0.70), ("CPH1", "none", 0.22), ("O", "A", -0.51)],
    'HSD': [("NH1", "D", -0.47), ("NR1", "D", -0.36), ("CPH1", "none", -0.05), ("CPH2", "none", -0.25),
            ("NR2", "A", -0.70), ("CPH1", "none", 0.22), ("O", "A", -0.51)],
    'HSE': [("NH1", "D", -0.47), ("NR2", "A", -0.70), ("CPH1", "none", 0.22), ("CPH2", "none", -0.25),
            ("NR1", "D", -0.36), ("CPH1", "none", -0.05), ("O", "A", -0.51)],
    'HSP': [("NH1", "D", -0.47), ("NR3", "D", -0.51), ("CPH1", "none", 0.19), ("CPH2", "none", 0.32),
            ("NR3", "D", -0.51), ("CPH1", "none", 0.19), ("O", "A", -0.51)],
    'ILE': [("NH1", "D", -0.47), ("CT3", "none", -0.27), ("CT3", "none", -0.27), ("O", "A", -0.51)],
    'LEU': [("NH1", "D", -0.47), ("CT1", "none", -0.09), ("CT3", "none", -0.27), ("CT3", "none", -0.27),
            ("O", "A", -0.51)],
    'LYS': [("NH1", "D", -0.47), ("NH3", "D", -0.30), ("O", "A", -0.51)],
    'MET': [("NH1", "D", -0.47), ("S", "A", -0.09), ("CT3", "none", -0.22), ("O", "A", -0.51)],
    'PHE': [("NH1", "D", -0.47), ("CA", "none", -0.115), ("CA", "none", -0.115), ("CA", "none", -0.115),
            ("CA", "none", -0.115), ("CA", "none", -0.115), ("CA", "none", -0.00), ("O", "A", -0.51)],
    'PRO': [("N", "D", -0.29), ("CP3", "none", 0.00), ("CP1", "none", 0.02), ("CP2", "none", -0.18),
            ("CP2", "none", -0.18), ("O", "A", -0.51)],
    'SER': [("NH1", "D", -0.47), ("OH1", "B", -0.66), ("O", "A", -0.51)],
    'THR': [("NH1", "D", -0.47), ("OH1", "B", -0.66), ("CT3", "none", -0.27), ("O", "A", -0.51)],
    'TRP': [("NH1", "D", -0.47), ("CY", "none", -0.03), ("CA", "none", -0.15), ("NY", "D", -0.51),
            ("CPT", "none", 0.24), ("CPT", "none", 0.11), ("CAI", "none", -0.25),
            ("CA", "none", -0.2), ("CAI", "none", -0.27), ("CA", "none", -0.14), ("O", "A", -0.51)],
    'TYR': [("NH1", "D", -0.47), ("CA", "none", 0.00), ("CA", "none", -0.115), ("CA", "none", -0.115),
            ("CA", "none", 0.11), ("OH1", "B", -0.54), ("CA", "none", -0.115), ("CA", "none", -0.115),
            ("O", "A", -0.51)],
    'VAL': [("NH1", "D", -0.47), ("CT3", "none", -0.27), ("CT3", "none", -0.27), ("O", "A", -0.51)]
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
    mainchain_atoms = ["N", "C", "CA", "O", "OT1", "OT2"]
    sidechain_atoms = []
    if resname == "ALA":
        sidechain_atoms = ["CB"]
    if resname == "ARG":
        sidechain_atoms = ["CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"]
    if resname == "ASP":
        sidechain_atoms = ["CB", "CG", "OD1", "OD2"]
    if resname == "ASN":
        sidechain_atoms = ["CB", "CG", "OD1", "ND2"]
    if resname == "CYS":
        sidechain_atoms = ["CB", "SG"]
    if resname == "GLN":
        sidechain_atoms = ["CB", "CG", "CD", "OE1", "NE2"]
    if resname == "GLU":
        sidechain_atoms = ["CB", "CG", "CD", "OE1", "OE2"]
    if resname == "GLY":
        sidechain_atoms = []
    if resname == "HIS" or resname == "HSD" or resname == "HSE" or resname == "HSP":
        sidechain_atoms = ["CB", "CG", "ND1", "CD2", "NE2", "CE1"]
    if resname == "ILE":
        sidechain_atoms = ["CB", "CG1", "CG2", "CD"]
    if resname == "LEU":
        sidechain_atoms = ["CB", "CG1", "CD1", "CD2"]
    if resname == "LYS":
        sidechain_atoms = ["CB", "CG", "CD", "CE", "NZ"]
    if resname == "MET":
        sidechain_atoms = ["CB", "CG", "SD", "CE"]
    if resname == "PHE":
        sidechain_atoms = ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"]
    if resname == "PRO":
        sidechain_atoms = ["CB", "CG", "CD"]
    if resname == "SER":
        sidechain_atoms = ["CB", "OG"]
    if resname == "THR":
        sidechain_atoms = ["CB", "OG1", "CG2"]
    if resname == "TRP":
        sidechain_atoms = ["CB", "CG", "CD2", "CE2", "CE3", "CD1", "NE1", "CZ2", "CZ3", "CH2", "NZ"]
    if resname == "TYR":
        sidechain_atoms = ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"]
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



