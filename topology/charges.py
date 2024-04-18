def atom_charge(atom):

    charge = 0.0
    if atom == "C": charge = 0.0
    if atom == "N": charge = -1.0
    if atom == "O": charge = -1.0
    if atom == "P": charge = -1.0
    if atom == "S": charge = -1.0
    if atom == "F": charge = -1.0
    if atom == "CL": charge = -1.0
    if atom == "BR": charge = -1.0
    if atom == "I": charge = -1.0
    if atom == "MG": charge = 2.0

    return charge
