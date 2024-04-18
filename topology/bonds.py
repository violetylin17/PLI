def bond_length(atom1, atom2):

    """
    Bond lengths taken from Handbook of Chemistry and Physics. The information provided there was very specific,
    so I tried to pick representative examples and used the bond lengths from those. Situations could arise where these
    lengths would be incorrect, probably slight errors (<0.06) in the hundreds.
    """

    distance = 0.0
    if atom1 == "C" and atom2 == "C": distance = 1.53
    if atom1 == "N" and atom2 == "N": distance = 1.425
    if atom1 == "O" and atom2 == "O": distance = 1.469
    if atom1 == "S" and atom2 == "S": distance = 2.048
    if (atom1 == "C" and atom2 == "H") or (atom1 == "H" and atom2 == "C"): distance = 1.059
    if (atom1 == "C" and atom2 == "N") or (atom1 == "N" and atom2 == "C"): distance = 1.469
    if (atom1 == "C" and atom2 == "O") or (atom1 == "O" and atom2 == "C"): distance = 1.413
    if (atom1 == "C" and atom2 == "S") or (atom1 == "S" and atom2 == "C"): distance = 1.819
    if (atom1 == "N" and atom2 == "H") or (atom1 == "H" and atom2 == "N"): distance = 1.009
    if (atom1 == "N" and atom2 == "O") or (atom1 == "O" and atom2 == "N"): distance = 1.463
    if (atom1 == "O" and atom2 == "S") or (atom1 == "S" and atom2 == "O"): distance = 1.577
    if (atom1 == "O" and atom2 == "H") or (atom1 == "H" and atom2 == "O"): distance = 0.967
    if (atom1 == "S" and atom2 == "H") or (atom1 == "H" and atom2 == "S"): distance = 2.025 / 1.5
    # This one not from source sited above. Not sure where it's from, but it wouldn't ever be used in the
    # current context ("AutoGrow")
    if (atom1 == "S" and atom2 == "N") or (atom1 == "H" and atom2 == "N"): distance = 1.633

    if (atom1 == "C" and atom2 == "F") or (atom1 == "F" and atom2 == "C"): distance = 1.399
    if (atom1 == "C" and atom2 == "CL") or (atom1 == "CL" and atom2 == "C"): distance = 1.790
    if (atom1 == "C" and atom2 == "BR") or (atom1 == "BR" and atom2 == "C"): distance = 1.910
    if (atom1 == "C" and atom2 == "I") or (atom1 == "I" and atom2 == "C"): distance = 2.162

    if (atom1 == "S" and atom2 == "BR") or (atom1 == "BR" and atom2 == "S"): distance = 2.321
    if (atom1 == "S" and atom2 == "CL") or (atom1 == "CL" and atom2 == "S"): distance = 2.283
    if (atom1 == "S" and atom2 == "F") or (atom1 == "F" and atom2 == "S"): distance = 1.640
    if (atom1 == "S" and atom2 == "I") or (atom1 == "I" and atom2 == "S"): distance = 2.687

    if (atom1 == "P" and atom2 == "BR") or (atom1 == "BR" and atom2 == "P"): distance = 2.366
    if (atom1 == "P" and atom2 == "CL") or (atom1 == "CL" and atom2 == "P"): distance = 2.008
    if (atom1 == "P" and atom2 == "F") or (atom1 == "F" and atom2 == "P"): distance = 1.495
    if (atom1 == "P" and atom2 == "I") or (atom1 == "I" and atom2 == "P"): distance = 2.490
    if (atom1 == "P" and atom2 == "O") or (atom1 == "O" and atom2 == "P"): distance = 1.6
    # estimate based on eye balling Handbook of Chemistry and Physics

    if (atom1 == "N" and atom2 == "BR") or (atom1 == "BR" and atom2 == "N"): distance = 1.843
    if (atom1 == "N" and atom2 == "CL") or (atom1 == "CL" and atom2 == "N"): distance = 1.743
    if (atom1 == "N" and atom2 == "F") or (atom1 == "F" and atom2 == "N"): distance = 1.406
    if (atom1 == "N" and atom2 == "I") or (atom1 == "I" and atom2 == "N"): distance = 2.2

    if (atom1 == "SI" and atom2 == "BR") or (atom1 == "BR" and atom2 == "SI"): distance = 2.284
    if (atom1 == "SI" and atom2 == "CL") or (atom1 == "CL" and atom2 == "SI"): distance = 2.072
    if (atom1 == "SI" and atom2 == "F") or (atom1 == "F" and atom2 == "SI"): distance = 1.636
    if (atom1 == "SI" and atom2 == "P") or (atom1 == "P" and atom2 == "SI"): distance = 2.264
    if (atom1 == "SI" and atom2 == "S") or (atom1 == "S" and atom2 == "SI"): distance = 2.145
    if (atom1 == "SI" and atom2 == "SI") or (atom1 == "SI" and atom2 == "SI"): distance = 2.359
    if (atom1 == "SI" and atom2 == "C") or (atom1 == "C" and atom2 == "SI"): distance = 1.888
    if (atom1 == "SI" and atom2 == "N") or (atom1 == "N" and atom2 == "SI"): distance = 1.743
    if (atom1 == "SI" and atom2 == "O") or (atom1 == "O" and atom2 == "SI"): distance = 1.631

    return distance
