import os
import numpy as np
from numpy import exp
from numpy.linalg import norm

from . import amino_acid
from . import guess


def list_of_interaction():
    list_of_type = {
        "hbond": ["H-bond(M)", "H-bond(S)"],
        "steric": ["Steric(M)", "Steric(S)"],
        "hydrophobic": ["Hydrophobic(M)", "Hydrophobic(S)"],
        "electrostatic": ["Electrostatic(M)", "Electrostatic(S)"],
        "pipi": ["Pi-pi stack", "Pi-t stack", "Pi-cation stack"]}

    return list_of_type


def near_active_site(lig, active_site_center=np.zeros(3), dist_criteria=5.0):
    """
    Description :
        Calculate the distance between active-site and ligand center to select target ligand

    :param lig: structure information of ligand
    :param active_site_center: center of target protein
    :param dist_criteria: Selection criteria of the distance between active-site and ligand center
    :return:
    """

    # Calculate center of ligand
    lig_coor = np.zeros(3)
    for atom in lig.atoms:
        lig_coor += atom.coordinate
    lig_center = lig_coor/len(lig.atoms)

    dist = norm(lig_center - active_site_center)
    if dist < dist_criteria:
        return True, dist
    else:
        return False, dist


def interaction_type(atomi, atomj, dist=float('inf')):
    """
    :param atomi:
    :param atomj:
    :param dist:
    :return:
    """

    hbond_criteria = 4.0        # D-H ... A  ==> maximum distance between A and D atoms
    non_polar_criteria = 4.0      # maximum interaction distance between two atoms

    atomname1 = guess.guess_atom_name(atomi.name)
    atomname2 = guess.guess_atom_name(atomj.name)

    if dist < non_polar_criteria:
        if atomname1 == 'N' or atomname1 == 'O' or atomname1 == 'F':
            if atomname2 == 'N' or atomname2 == 'O' or atomname2 == 'F':
                if dist < hbond_criteria:
                    return 1     # H-bond interaction
        elif atomname1 == 'C' or atomname1 == 'S':
            if atomname2 == 'C' or atomname2 == 'S':
                return 2            # non-polar interaction
    return 0


def old_interaction_type(atomi, atomj, dist=float('inf')):
    """
    :param atomi:
    :param atomj:
    :param dist:
    :return:
    """

    hbond_criteria = 4.0        # D-H ... A  ==> maximum distance between A and D atoms
    non_polar_criteria = 6.0      # maximum interaction distance between two atoms

    if dist < non_polar_criteria:
        if atomi.name == 'N' or atomi.name == 'O' or atomi.name == 'F':
            if atomj.name == 'N' or atomj.name == 'O' or atomj.name == 'F':
                if dist < hbond_criteria:
                    return 1     # H-bond interaction
                else:
                    return 2    # non-polar interaction
            else:
                return 2        # non-polar interaction
        elif atomi.type == 'C':
            return 2            # non-polar interaction
        elif atomi.type == 'S':
            return 2            # non-polar interaction
    else:
        return 0


def nearest_interaction(res, lig):
    """
    :param res:
    :param lig:
    :return:
    """

    interaction_list = []
    for atomi in res.atoms:
        min_dist = float('inf')
        atom1 = amino_acid.Atom()
        atom2 = amino_acid.Atom()
        for atomj in lig.atoms:
            dist = abs(norm(atomi.coordinate - atomj.coordinate))
            if dist < min_dist:
                atom1 = atomi
                atom2 = atomj
                min_dist = dist
        interaction_list.append((interaction_type(atom1, atom2, min_dist), min_dist))

    return interaction_list


def nearest_atom(res, lig, min_dist=float('-inf'), max_dist=float('inf')):

    main_nearest_atoms = []
    side_nearest_atoms = []
    for index, atomi in enumerate(res.atoms):
        flag = 0
        mini_dist = float('inf')
        for atomj in lig.atoms:
            dist = abs(norm(atomi.coordinate - atomj.coordinate))
            if dist < mini_dist and min_dist <= dist <= max_dist:
                flag = 1
                mini_dist = dist
                tmp = atomi, atomj, mini_dist

        if flag:
            if index == 0 or index == len(res.atoms)-1:
                main_nearest_atoms.append(tmp)
            else:
                side_nearest_atoms.append(tmp)

    return main_nearest_atoms, side_nearest_atoms


def contact_atoms(res, lig, min_dist=float('-inf'), max_dist=float('inf')):

    contact = []
    for atomi in res:
        for atomj in lig:
            dist = abs(norm(atomi.coordinate - atomj.coordinate))
            if min_dist <= dist <= max_dist:
                tmp = (atomi, atomj, dist)
                contact.append(tmp)

    return contact

"""
def electrostatic_scoring_function(dist, max_dist=4.0, min_dist=2.2):
    if dist < max_dist:
        score = (max_dist - dist) * (1.0 / (max_dist - min_dist))
    else:
        score = 0

    return score


def electrostatic_interaction(prot, lig, log=False):

    max_distance = 4.0
    min_distance = 2.2

    output_info1 = "\n============= Electrostatic interaction =============\n"
    output_info1 += "\nSelectoin criteria:\n"
    output_info1 += "Minimum distance > {0:2.2f} A\n".format(min_distance)
    output_info1 += "Maximum distance < {0:2.2f} A\n".format(max_distance)
    output_info1 += "\nRaw data:\n"

    # Main chain information
    main_info = "\n[ MAIN CHAIN ]\n\n"
    main_info += "                 protein   ::   ligand                      |   distance   |    score    \n"
    main_info += "----------------------------------------------------------- | ------------ | ----------- \n"

    # Side chain information
    side_info = "\n[ SIDE CHAIN ]\n\n"
    side_info += "                 protein   ::   ligand                      |   distance   |    score    \n"
    side_info += "----------------------------------------------------------- | ------------ | ----------- \n"

    main_interaction = []
    side_interaction = []
    for res in prot:
        # Main chain interaction
        main_contacts = contact_atoms(res.main_chain(), lig.atoms, min_dist=min_distance, max_dist=max_distance)

        flag = 0
        main_chain = []
        for contact in main_contacts:
            atomname1 = guess.guess_atom_name(contact[0].name)
            atomname2 = guess.guess_atom_name(contact[1].name)
            if atomname1 in ["C"] or atomname2 in ["C"]:
                continue

            force = contact[0].charge * contact[1].charge
            if force < 0 and contact[2] < max_distance:
                flag = 1
                info = contact[0].name, contact[1].name + str(contact[1].index), electrostatic_scoring_function(contact[2], min_dist=min_distance, max_dist=max_distance)
                main_chain.append(info)

                if info[2]:
                    o_info1 = " {0:4s} {1:4d}{2:4s} {3:4s} {4:6s} {5:6d}".format(res.segid, res.resid, res.type,
                                                                                 contact[0].name, contact[0].type,
                                                                                 contact[0].index)
                    o_info1 += "  ::   LIG   {0:4s} {1:2s} {2:3d} ".format(contact[1].name, contact[1].type,
                                                                           contact[1].index)
                    o_info1 += "  |   {0:^8.2f}   |   {1:^8.2f}\n".format(contact[2], info[2])
                    main_info += o_info1

        if flag:
            main_interaction.append((res.segid, res.resid, res.type, main_chain))

        # Side chain interaction
        side_contacts = contact_atoms(res.side_chain(), lig.atoms, min_dist=min_distance, max_dist=max_distance)

        flag = 0
        side_chain = []
        for contact in side_contacts:
            atomname1 = guess.guess_atom_name(contact[0].name)
            atomname2 = guess.guess_atom_name(contact[1].name)
            if atomname1 in ["C"] or atomname2 in ["C"]:
                continue

            force = contact[0].charge * contact[1].charge
            if force < 0 and contact[2] < max_distance:
                flag = 1
                info = contact[0].name, contact[1].name + str(contact[1].index), electrostatic_scoring_function(contact[2], min_dist=min_distance, max_dist=max_distance)
                side_chain.append(info)

                if info[2]:
                    o_info1 = " {0:4s} {1:4d}{2:4s} {3:4s} {4:6s} {5:6d}".format(res.segid, res.resid, res.type,
                                                                                 contact[0].name, contact[0].type,
                                                                                 contact[0].index)
                    o_info1 += "  ::   LIG   {0:4s} {1:2s} {2:3d} ".format(contact[1].name, contact[1].type,
                                                                           contact[1].index)
                    o_info1 += "  |   {0:^8.2f}   |   {1:^8.2f}\n".format(contact[2], info[2])
                    side_info += o_info1

        if flag:
            side_interaction.append((res.segid, res.resid, res.type, side_chain))

    output_info2 = "\n\n"

    if log:
        message = output_info1 + main_info + side_info + output_info2
        return main_interaction, side_interaction, message
    else:
        return main_interaction, side_interaction


def h_bond_scoring_function(dist, max_dist=4.0, min_dist=2.2):
    if dist < max_dist:
        score = (max_dist - dist) * (1.0 / (max_dist - min_dist))
    else:
        score = 0
    return score


def h_bond_interaction(prot, lig, log=True):

    max_distance = 4.0
    min_distance = 2.2

    output_info1 = "\n============= Hydrogen bond interaction =============\n"
    output_info1 += "\nSelectoin criteria:\n"
    output_info1 += "Minimum distance > {0:2.2f} A\n".format(min_distance)
    output_info1 += "Maximum distance < {0:2.2f} A\n".format(max_distance)
    output_info1 += "\nRaw data:\n"
    output_info1 += "\nA: Acceptor  D: Donor  B: Both\n"

    # Main chain information
    main_info = "\n[ MAIN CHAIN ]\n\n"
    main_info += "                     protein   ::   ligand                          |   distance   |    score    \n"
    main_info += "------------------------------------------------------------------- | ------------ | ----------- \n"

    # Side chain information
    side_info = "\n[ SIDE CHAIN ]\n\n"
    side_info += "                     protein   ::   ligand                          |   distance   |    score    \n"
    side_info += "------------------------------------------------------------------- | ------------ | ----------- \n"

    main_interaction = []
    side_interaction = []
    for res in prot:
        # Main chain interaction
        main_contacts = contact_atoms(res.main_chain(), lig.atoms, min_dist=min_distance, max_dist=max_distance)

        flag = 0
        main_chain = []
        for contact in main_contacts:
            if contact[0].hbond == "none" or contact[1].hbond == "none" \
                    or (contact[0].hbond == "D" and contact[1].hbond == "D") \
                    or (contact[0].hbond == "A" and contact[1].hbond == "A"):
                continue

            if interaction_type(contact[0], contact[1], contact[2]) == 1:       # If exist any polar interaction
                flag = 1
                info = contact[0].name, contact[1].name + str(contact[1].index), h_bond_scoring_function(contact[2], min_dist=min_distance, max_dist=max_distance)
                main_chain.append(info)

                if info[2]:
                    o_info1 = " {0:4s} {1:4d}{2:4s} {3:4s} {4:6s} {5:6d} ({6})".format(res.segid, res.resid, res.type,
                                                                                       contact[0].name,
                                                                                       contact[0].type,
                                                                                       contact[0].index,
                                                                                       contact[0].hbond)
                    o_info1 += "  ::  ({0}) LIG   {1:4s} {2:3s} {3:3d} ".format(contact[1].hbond, contact[1].name,
                                                                                contact[1].type, contact[1].index)
                    o_info1 += "  |   {0:^8.2f}   |   {1:^8.2f}\n".format(contact[2], info[2])
                    main_info += o_info1

        if flag:
            main_interaction.append((res.segid, res.resid, res.type, main_chain))

        # Side chain interaction
        side_contacts = contact_atoms(res.side_chain(), lig.atoms, min_dist=min_distance, max_dist=max_distance)

        flag = 0
        side_chain = []
        for contact in side_contacts:
            if contact[0].hbond == "none" or contact[1].hbond == "none" \
                    or (contact[0].hbond == "D" and contact[1].hbond == "D") \
                    or (contact[0].hbond == "A" and contact[1].hbond == "A"):
                continue

            if interaction_type(contact[0], contact[1], contact[2]) == 1:       # If exist any polar interaction
                flag = 1
                info = contact[0].name, contact[1].name + str(contact[1].index), h_bond_scoring_function(contact[2], min_dist=min_distance, max_dist=max_distance)
                side_chain.append(info)

                if info[2]:
                    o_info1 = " {0:4s} {1:4d}{2:4s} {3:4s} {4:6s} {5:6d} ({6})".format(res.segid, res.resid, res.type,
                                                                                       contact[0].name,
                                                                                       contact[0].type,
                                                                                       contact[0].index,
                                                                                       contact[0].hbond)
                    o_info1 += "  ::  ({0}) LIG   {1:4s} {2:3s} {3:3d} ".format(contact[1].hbond, contact[1].name,
                                                                                contact[1].type, contact[1].index)
                    o_info1 += "  |   {0:^8.2f}   |   {1:^8.2f}\n".format(contact[2], info[2])
                    side_info += o_info1

        if flag:
            side_interaction.append((res.segid, res.resid, res.type, side_chain))

    output_info2 = "\n\n"

    if log:
        message = output_info1 + main_info + side_info + output_info2
        return main_interaction, side_interaction, message
    else:
        return main_interaction, side_interaction


def hydrophobic_scoring_function(dist, max_dist=4.0, min_dist=3.2):
    if dist < max_dist:
        score = (max_dist - dist) * (1.0 / (max_dist - min_dist))
    else:
        score = 0
    return score
"""


def contact_pair(res, lig, cutoff=8.0):
    """
    Description:
        Calculate all possible contact pairs between protein and ligand in cutoff criteria

    :param res: <Residue class object>  target residue
    :param lig: <AtomGroup class object>  target ligand
    :param cutoff:  <Numerical> Maximum distance between interaction pair
    :return: <List> list of contact pair (atomI, atomJ,surface d, dist)
    """

    contact = []
    for atomi in res.atoms:
        for atomj in lig.atoms:
            dist = norm(atomi.coordinate - atomj.coordinate)
            if dist < cutoff:
                d = dist - atomi.vdw_r - atomj.vdw_r
                contact.append((atomi, atomj, d, dist))

    return contact


def steric_scoring_function(d, repul_max=0.0, g1=-0.035579, g2=-0.005156, repul=0.840245):
    """
    Description:
        Steric scoring function reference to AutoDock Vina

    :param d: <Numerical> Distance between interaction pair
    :param repul_max: <Numerical> Maximum distance criteria of repulsion score
    :param g1: <Numerical> Weight coefficient of gauss1
    :param g2: <Numerical> Weight coefficient of gauss2
    :param repul: <Numerical> Weight coefficient of replusion
    :return: <Numerical> Predict score
    """

    def gauss(dist, offset=0, width=1):
        return exp(-((dist-offset) / width)**2)

    def repulsion(dist, vmax):
        if dist < vmax:
            return dist ** 2
        else:
            return 0

    return g1*gauss(d, 0, 0.5) + g2*gauss(d, 3, 2.0) + repul*repulsion(d, repul_max)


def steric_interaction(prot, lig, cutoff=8.0, log=False):
    """
    Description:
        Calculate steric type interaction between target protein and ligand

    :param prot: <Protein class object>  target protein
    :param lig: <AtomGroup class object>  target ligand
    :param cutoff: <Numerical> Maximum distance between interaction pair
    :param log: <Bool> Enable/Disable output a log file
    :return: <Tuple> (main_chain, side_chain, log)
    """

    output_info1 = "\n============= Steric interaction =============\n"
    output_info1 += "\nSelection criteria:\n"
    output_info1 += "Cutoff distance > {0:2.2f} A\n".format(cutoff)
    output_info1 += "\nRaw data:\n"

    # Main chain information
    main_info = "\n[ MAIN CHAIN ]\n\n"
    main_info += "                 Protein   ::   Ligand                      |  Distance  |    Energy    |\n"
    main_info += " SEGID  RESID  NAME TYPE    INDEX  ::   LIG   I    J  INDEX |      A     |   kcal/mole  |\n"
    main_info += "----------------------------------------------------------- | ---------- | ------------ |\n"

    # Side chain information
    side_info = "\n[ SIDE CHAIN ]\n\n"
    side_info += "                 Protein   ::   Ligand                      |  Distance  |    Energy    |\n"
    side_info += " SEGID  RESID  NAME TYPE    INDEX  ::   LIG   I    J  INDEX |      A     |   kcal/mole  |\n"
    side_info += "----------------------------------------------------------- | ---------- | ------------ |\n"

    main_chain = []
    side_chain = []
    for res in prot.residues:
        # Main chain interaction
        main_contacts = contact_pair(res.main_chain(), lig, cutoff=cutoff)

        for contact in main_contacts:
            if contact[0].type in ["H", "HC"] or contact[0].type in ["H", "HC"]:
                continue

            # contact ==> atomI, atomJ, dist
            # contact ==> atomI, atomJ, dist, surface d by weulun
            #info = res.segid, res.resid, res.resn, contact[0].name, contact[1].name + str(contact[1].index), \
            #       steric_scoring_function(contact[2])
            info = res.segid, res.resid, res.resn, contact[0].name, contact[1].name + str(contact[1].id), \
                   steric_scoring_function(contact[2]), contact[2]  ## 2022/0817 amended by Bruce
            
            main_chain.append(info)

            # Log file message of main chain interaction
            o_info1 = " {0:4s} {1:4d}{2:4s} {3:4s} {4:6s} {5:6d}".format(res.segid, res.resid, res.resn, contact[0].name,
                                                                         contact[0].type, contact[0].id)  ## 2022/0817 amended by Bruce
            o_info1 += "  ::   LIG   {0:4s} {1:2s} {2:3d} ".format(contact[1].name, contact[1].type, contact[1].id)  ## 2022/0817 amended by Bruce
            o_info1 += "  |  {0:^8.2f}  |   {1:^8.5f}   |\n".format(contact[3], info[5])
            main_info += o_info1

        # Side chain interaction
        side_contacts = contact_pair(res.side_chain(), lig, cutoff=cutoff)

        for contact in side_contacts:
            if contact[0].type in ["H", "HC"] or contact[0].type in ["H", "HC"]:
                continue

            # contact ==> atomI, atomJ, dist
            # contact ==> atomI, atomJ, dist, d
            #info = res.segid, res.resid, res.resn, contact[0].name, contact[1].name + str(contact[1].index), \
            #       steric_scoring_function(contact[2])
            info = res.segid, res.resid, res.resn, contact[0].name, contact[1].name + str(contact[1].id), \
                   steric_scoring_function(contact[2]), contact[2]  ## 2022/0817 amended by Bruce
            side_chain.append(info)

            # Log file message of side chain interaction
            o_info1 = " {0:4s} {1:4d}{2:4s} {3:4s} {4:6s} {5:6d}".format(res.segid, res.resid, res.resn, contact[0].name,
                                                                         contact[0].type, contact[0].id) ## 2022/0817 amended by Bruce
            o_info1 += "  ::   LIG   {0:4s} {1:2s} {2:3d} ".format(contact[1].name, contact[1].type, contact[1].id)  ## 2022/0817 amended by Bruce
            o_info1 += "  |  {0:^8.2f}  |   {1:^8.5f}   |\n".format(contact[3], info[5])
            side_info += o_info1

    output_info2 = "\n\n"

    if log:
        message = output_info1 + main_info + side_info + output_info2
        return main_chain, side_chain, message
    else:
        return main_chain, side_chain


def hydrophobic_scoring_function(d, min_dist=0.5, max_dist=1.5, weight=-0.035069):
    """
    Description:
        Hydrophobic scoring function reference to AutoDock Vina

    :param d: <Numerical> Distance between interaction pair
    :param min_dist: <Numerical> Minimum distance criteria of scoring function
    :param max_dist: <Numerical> Maximum distance criteria of scoring function
    :param weight: <Numerical> Weight coefficient
    :return: <Numerical> Predict score
    """

    if d < min_dist:
        return weight * 1
    elif d > max_dist:
        return 0
    else:
        return weight * (max_dist - d) * (1.0 / (max_dist - min_dist))


def hydrophobic_interaction(prot, lig, min_dist=0.5, max_dist=1.5, cutoff=8.0, log=False):
    """
    Description:
            Calculate hydrophobic type interaction between target protein and ligand

    :param prot: <Protein class object>  target protein
    :param lig: <AtomGroup class object>  target ligand
    :param min_dist: <Numerical> Maximum distance criteria of scoring function
    :param max_dist: <Numerical> Minimum distance criteria of scoring function
    :param cutoff: <Numerical> Maximum distance between interaction pair
    :param log: <Bool> Enable/Disable output a log file
    :return: <Tuple> (main_chain, side_chain, log)
    """

    output_info1 = "\n============= Hydrophobic interaction =============\n"
    output_info1 += "\nSelection criteria:\n"
    output_info1 += "Minimum distance > {0:2.2f} A\n".format(min_dist)
    output_info1 += "Maximum distance < {0:2.2f} A\n".format(max_dist)
    output_info1 += "Cutoff distance > {0:2.2f} A\n".format(cutoff)

    output_info1 += "\nRaw data:\n"

    # Main chain information
    main_info = "\n[ MAIN CHAIN ]\n\n"
    main_info += "                 Protein   ::   Ligand                      |  Distance  |    Energy    |\n"
    main_info += " SEGID  RESID  NAME TYPE    INDEX  ::   LIG   I    J  INDEX |      A     |   kcal/mole  |\n"
    main_info += "----------------------------------------------------------- | ---------- | ------------ |\n"

    # Side chain information
    side_info = "\n[ SIDE CHAIN ]\n\n"
    side_info += "                 Protein   ::   Ligand                      |  Distance  |    Energy    |\n"
    side_info += " SEGID  RESID  NAME TYPE    INDEX  ::   LIG   I    J  INDEX |      A     |   kcal/mole  |\n"
    side_info += "----------------------------------------------------------- | ---------- | ------------ |\n"

    main_chain = []
    side_chain = []
    for res in prot.residues:
        # Main chain interaction
        main_contacts = contact_pair(res.main_chain(), lig, cutoff=cutoff)

        for contact in main_contacts:
            if contact[2] >= max_dist:
                continue

            # Hydrophobic selection criteria
            atomi_name = guess.guess_atom_name(contact[0].name)
            atomj_name = guess.guess_atom_name(contact[1].name)
            if atomi_name not in ["C", "F", "CL", "BR", "I"] or atomj_name not in ["C", "F", "CL", "BR", "I"]:
                continue

            if (contact[0].hydro, contact[1].hydro) != ("H", "H"):
                continue

            # contact ==> atomI, atomJ, dist
            # contact ==> atomI, atomJ, dist, d
            #info = res.segid, res.resid, res.resn, contact[0].name, contact[1].name + str(contact[1].index), \
            #    hydrophobic_scoring_function(contact[2], min_dist, max_dist)
            info = res.segid, res.resid, res.resn, contact[0].name, contact[1].name + str(contact[1].id), \
                hydrophobic_scoring_function(contact[2], min_dist, max_dist), contact[2]  ## 2022/0817 amended by Bruce
            main_chain.append(info)

            # Log file message of main chain interaction
            o_info1 = " {0:4s} {1:4d}{2:4s} {3:4s} {4:6s} {5:6d}".format(res.segid, res.resid, res.resn, contact[0].name,
                                                                         contact[0].type, contact[0].id)  ## 2022/0817 amended by Bruce
            o_info1 += "  ::   LIG   {0:4s} {1:2s} {2:3d} ".format(contact[1].name, contact[1].type, contact[1].id)  ## 2022/0817 amended by Bruce
            o_info1 += "  |  {0:^8.2f}  |   {1:^8.5f}   |\n".format(contact[3], info[5])
            main_info += o_info1

        # Side chain interaction
        side_contacts = contact_pair(res.side_chain(), lig, cutoff=cutoff)

        for contact in side_contacts:
            if contact[2] >= max_dist:
                continue

            # Hydrophobic selection criteria
            atomi_name = guess.guess_atom_name(contact[0].name)
            atomj_name = guess.guess_atom_name(contact[1].name)
            if atomi_name not in ["C", "F", "CL", "BR", "I"] or atomj_name not in ["C", "F", "CL", "BR", "I"]:
                continue

            if (contact[0].hydro, contact[1].hydro) != ("H", "H"):
                continue

            # contact ==> atomI, atomJ, dist
            # contact ==> atomI, atomJ, dist, d
            #info = res.segid, res.resid, res.resn, contact[0].name, contact[1].name + str(contact[1].index), \
            #    hydrophobic_scoring_function(contact[2], min_dist, max_dist)
            info = res.segid, res.resid, res.resn, contact[0].name, contact[1].name + str(contact[1].id), \
                hydrophobic_scoring_function(contact[2], min_dist, max_dist) ,contact[2]  ## 2022/0817 amended by Bruce 
            side_chain.append(info)

            # Log file message of side chain interaction
            o_info1 = " {0:4s} {1:4d}{2:4s} {3:4s} {4:6s} {5:6d}".format(res.segid, res.resid, res.resn, contact[0].name,
                                                                         contact[0].type, contact[0].id)  ## 2022/0817 amended by Bruce
            o_info1 += "  ::   LIG   {0:4s} {1:2s} {2:3d} ".format(contact[1].name, contact[1].type, contact[1].id)  ## 2022/0817 amended by Bruce
            o_info1 += "  |  {0:^8.2f}  |   {1:^8.5f}   |\n".format(contact[3], info[5])
            side_info += o_info1

    output_info2 = "\n\n"

    if log:
        message = output_info1 + main_info + side_info + output_info2
        return main_chain, side_chain, message
    else:
        return main_chain, side_chain


def hbond_scoring_function(d, min_dist=-0.7, max_dist=0.0, weight=-0.587439):
    """
    Description:
        Hbond scoring function reference to AutoDock Vina

    :param d: <Numerical> Distance between interaction pair
    :param min_dist: <Numerical> Minimum distance criteria of scoring function
    :param max_dist: <Numerical> Maximum distance criteria of scoring function
    :param weight: <Numerical> Weight coefficient
    :return: <Numerical> Predict score
    """

    if d < min_dist:
        return weight * 1
    elif d > max_dist:
        return 0
    else:
        return weight * (max_dist - d) * (1.0 / (max_dist - min_dist))


def hbond_interaction(prot, lig, min_dist=-0.7, max_dist=0.0, cutoff=8.0, log=False):
    """
    Description:
        Calculate hbond type interaction between target protein and ligand

    :param prot: <Protein class object>  target protein
    :param lig: <AtomGroup class object>  target ligand
    :param min_dist: <Numerical> Maximum distance criteria of scoring function
    :param max_dist: <Numerical> Minimum distance criteria of scoring function
    :param cutoff: <Numerical> Maximum distance between interaction pair
    :param log: <Bool> Enable/Disable output a log file
    :return: <Tuple> (main_chain, side_chain, log)
    """

    output_info1 = "\n============= Hydrogen bond interaction =============\n"
    output_info1 += "\nSelection criteria:\n"
    output_info1 += "Minimum distance > {0:2.2f} A\n".format(min_dist)
    output_info1 += "Maximum distance < {0:2.2f} A\n".format(max_dist)
    output_info1 += "Cutoff distance > {0:2.2f} A\n".format(cutoff)

    output_info1 += "\nRaw data:\n"

    # Main chain information
    main_info = "\n[ MAIN CHAIN ]\n\n"
    main_info += "                 Protein   ::   Ligand                      |  Distance  |    Energy    |\n"
    main_info += " SEGID  RESID  NAME TYPE    INDEX  ::   LIG   I    J  INDEX |      A     |   kcal/mole  |\n"
    main_info += "----------------------------------------------------------- | ---------- | ------------ |\n"

    # Side chain information
    side_info = "\n[ SIDE CHAIN ]\n\n"
    side_info += "                 Protein   ::   Ligand                      |  Distance  |    Energy    |\n"
    side_info += " SEGID  RESID  NAME TYPE    INDEX  ::   LIG   I    J  INDEX |      A     |   kcal/mole  |\n"
    side_info += "----------------------------------------------------------- | ---------- | ------------ |\n"

    main_chain = []
    side_chain = []
    for res in prot.residues:
        # Main chain interaction
        main_contacts = contact_pair(res.main_chain(), lig, cutoff=cutoff)

        for contact in main_contacts:
            if contact[2] >= max_dist:
                continue

            # Hbond selection criteria
            if contact[0].hbond == "none" or contact[1].hbond == "none" \
                    or contact[0].hbond == "HD" or contact[1].hbond == "HD" \
                    or (contact[0].hbond, contact[1].hbond) in [("D", "D"), ("A", "A")]:
                continue

            # contact ==> atomI, atomJ, dist
            # contact ==> atomI, atomJ, dist, d
            #info = res.segid, res.resid, res.resn, contact[0].name, contact[1].name + str(contact[1].index), \
            #    hbond_scoring_function(contact[2], min_dist, max_dist)
            #info = res.segid, res.resid, res.resn, contact[0].name, contact[1].name + str(contact[1].index), \
            info = res.segid, res.resid, res.resn, contact[0].name, contact[1].name + str(contact[1].id), \
                hbond_scoring_function(contact[2], min_dist, max_dist), contact[2]  ## 2022/08/16 amended by Bruce
            main_chain.append(info)

            # Log file message of main chain interaction
            o_info1 = " {0:4s} {1:4d}{2:4s} {3:4s} {4:6s} {5:6d}".format(res.segid, res.resid, res.resn, contact[0].name,
                                                                         contact[0].type, contact[0].id)  ## 2022/08/16 amended by Bruce
            #o_info1 += "  ::   LIG   {0:4s} {1:2s} {2:3d} ".format(contact[1].name, contact[1].type, contact[1].index)
            o_info1 += "  ::   LIG   {0:4s} {1:2s} {2:3d} ".format(contact[1].name, contact[1].type, contact[1].id)  ## 2022/08/16 amended by Bruce
            o_info1 += "  |  {0:^8.2f}  |   {1:^8.5f}   |\n".format(contact[3], info[5])
            main_info += o_info1

        # Side chain interaction
        side_contacts = contact_pair(res.side_chain(), lig, cutoff=cutoff)

        for contact in side_contacts:
            if contact[2] >= max_dist:
                continue

            # Hbond selection criteria
            if contact[0].hbond == "none" or contact[1].hbond == "none" \
                    or contact[0].hbond == "HD" or contact[1].hbond == "HD" \
                    or (contact[0].hbond, contact[1].hbond) in [("D", "D"), ("A", "A")]:
                continue

            # contact ==> atomI, atomJ, dist
            # contact ==> atomI, atomJ, dist, d
            #info = res.segid, res.resid, res.resn, contact[0].name, contact[1].name + str(contact[1].index), \
            #    hbond_scoring_function(contact[2], min_dist, max_dist)
            #info = res.segid, res.resid, res.resn, contact[0].name, contact[1].name + str(contact[1].index), \
            info = res.segid, res.resid, res.resn, contact[0].name, contact[1].name + str(contact[1].id), \
                hbond_scoring_function(contact[2], min_dist, max_dist), contact[2]  ## 2022/08/16 amended by Bruce
            side_chain.append(info)

            # Log file message of side chain interaction
            o_info1 = " {0:4s} {1:4d}{2:4s} {3:4s} {4:6s} {5:6d}".format(res.segid, res.resid, res.resn, contact[0].name,
                                                                         contact[0].type, contact[0].id)  ## 2022/08/16 amended by Bruce
            #o_info1 += "  ::   LIG   {0:4s} {1:2s} {2:3d} ".format(contact[1].name, contact[1].type, contact[1].index)
            o_info1 += "  ::   LIG   {0:4s} {1:2s} {2:3d} ".format(contact[1].name, contact[1].type, contact[1].id)  ## 2022/08/16 amended by Bruce
            o_info1 += "  |  {0:^8.2f}  |   {1:^8.5f}   |\n".format(contact[3], info[5])
            side_info += o_info1

    output_info2 = "\n\n"

    if log:
        message = output_info1 + main_info + side_info + output_info2
        return main_chain, side_chain, message
    else:
        return main_chain, side_chain

"""
def hydrophobic_interaction(prot, lig, log=False):

    max_distance = 5.0
    min_distance = 3.2

    output_info1 = "\n============= Hydrophobic interaction =============\n"
    output_info1 += "\nSelectoin criteria:\n"
    output_info1 += "Minimum distance > {0:2.2f} A\n".format(min_distance)
    output_info1 += "Maximum distance < {0:2.2f} A\n".format(max_distance)
    output_info1 += "\nRaw data:\n"

    # Main chain information
    main_info = "\n[ MAIN CHAIN ]\n\n"
    main_info += "                 protein   ::   ligand                      |   distance   |    score    \n"
    main_info += "----------------------------------------------------------- | ------------ | ----------- \n"

    # Side chain information
    side_info = "\n[ SIDE CHAIN ]\n\n"
    side_info += "                 protein   ::   ligand                      |   distance   |    score    \n"
    side_info += "----------------------------------------------------------- | ------------ | ----------- \n"

    main_interaction = []
    side_interaction = []
    for res in prot:
        # Main chain interaction
        main_contacts = contact_atoms(res.main_chain(), lig.atoms, min_dist=min_distance, max_dist=max_distance)

        flag = 0
        main_chain = []
        for contact in main_contacts:
            if interaction_type(contact[0], contact[1], contact[2]) == 2:       # If exist any non-polar interaction
                flag = 1
                #dist = contact[2] - vdw_radius(contact[0]) - vdw_radius(contact[1])
                #info = contact[0].name, contact[1].name, hydrophobic_scoring_function(dist)
                info = contact[0].name, contact[1].name + str(contact[1].index), hydrophobic_scoring_function(contact[2], min_dist=min_distance, max_dist=max_distance)
                main_chain.append(info)

                if info[2]:
                    o_info1 = " {0:4s} {1:4d}{2:4s} {3:4s} {4:6s} {5:6d}".format(res.segid, res.resid, res.type,
                                                                                 contact[0].name, contact[0].type,
                                                                                 contact[0].index)
                    o_info1 += "  ::   LIG   {0:4s} {1:2s} {2:3d} ".format(contact[1].name, contact[1].type,
                                                                           contact[1].index)
                    o_info1 += "  |   {0:^8.2f}   |   {1:^8.2f}\n".format(contact[2], info[2])
                    main_info += o_info1

        if flag:
            main_interaction.append((res.segid, res.resid, res.type, main_chain))

        # Side chain interaction
        side_contacts = contact_atoms(res.side_chain(), lig.atoms, min_dist=min_distance, max_dist=max_distance)

        flag = 0
        side_chain = []
        for contact in side_contacts:
            if interaction_type(contact[0], contact[1], contact[2]) == 2:      # If exist any non-polar interaction
                flag = 1
                #dist = contact[2] - vdw_radius(contact[0]) - vdw_radius(contact[1])
                #info = contact[0].name, contact[1].name, hydrophobic_scoring_function(dist)
                info = contact[0].name, contact[1].name + str(contact[1].index), hydrophobic_scoring_function(contact[2], min_dist=min_distance, max_dist=max_distance)
                side_chain.append(info)

                if info[2]:
                    o_info1 = " {0:4s} {1:4d}{2:4s} {3:4s} {4:6s} {5:6d}".format(res.segid, res.resid, res.type,
                                                                                 contact[0].name, contact[0].type,
                                                                                 contact[0].index)
                    o_info1 += "  ::   LIG   {0:4s} {1:2s} {2:3d} ".format(contact[1].name, contact[1].type,
                                                                           contact[1].index)
                    o_info1 += "  |   {0:^8.2f}   |   {1:^8.2f}\n".format(contact[2], info[2])
                    side_info += o_info1

        if flag:
            side_interaction.append((res.segid, res.resid, res.type, side_chain))

    output_info2 = "\n\n"

    if log:
        message = output_info1 + main_info + side_info + output_info2
        return main_interaction, side_interaction, message
    else:
        return main_interaction, side_interaction
"""

def project_point_onto_plane(atom_coor, ring):
    # The plane_coefficients are [a,b,c,d], where the plane is ax + by + cz = d
    # First, define a plane using coefficients a, b, c, d such that ax + by + cz = d
    a, b, c, d = ring.plane_coeff

    # Second, define a point in space (s,u,v)
    s, u, v = atom_coor

    # the formula of a line perpendicular to the plan passing through (s,u,v) is:
    # x = s + at
    # y = u + bt
    # z = v + ct
    t = (d - a * s - b * u - c * v) / (a * a + b * b + c * c)

    # here's the point closest on the plane
    x = s + a * t
    y = u + b * t
    z = v + c * t
    project_point = np.array([x, y, z])

    return project_point


def pi_pi_stacking_interaction(prot, lig, log=False):
    # Interaction criteria
    pi_stack_cutoff = 7.5
    pi_stack_angle_tolerance = 30
    pi_padding_dist = 0.75
    t_stack_angle_tolerance = 30
    t_stack_closest_cutoff = 5
    cation_stack_closest_cutoff = 4

    output_info1 = "\n============= Aromatic interaction =============\n"
    output_info1 += "\nSelectoin criteria:\n"
    output_info1 += "Pi-pi stack cutoff > {0:2.2f} A\n".format(pi_stack_cutoff)
    output_info1 += "Pi-pi stack angle tolerance < {0:2.2f} A\n".format(pi_stack_angle_tolerance)
    output_info1 += "Pi-T stack cutoff > {0:2.2f} A\n".format(t_stack_closest_cutoff)
    output_info1 += "Pi-T stack angle tolerance < {0:2.2f} A\n".format(t_stack_angle_tolerance)

    pi_pi_info = "\n[ PI-PI STACKING ]\n"
    pi_t_info = "\n[ PI-T STACKING ]\n"
    pi_cation_info = "\n[ PI-CATION STACKING ]\n"

    output_info1 += "\nLigand Ring Structure Information:\n"
    for index, ring in enumerate(lig.rings):
        output_info1 += "\nIndex: {}\n".format(index)
        output_info1 += "{} atoms ring\n".format(len(ring))
        for atom in ring.atoms:
            output_info1 += str(atom) + "\n"

    output_info1 += "\nRaw data:\n"

    # Calculate pi-pi and pi-T interaction

    pi_pi_interaction = []
    pi_t_interaction = []
    for res in prot:
        pi_pi = []
        pi_t = []
        for ring1 in res.rings:
            for ring2 in lig.rings:
                # Check distance between two rings' center < 7.5 angstrom
                ring_center_dist = abs(norm(ring1.center - ring2.center))
                if ring_center_dist <= pi_stack_cutoff:
                    # Identify pi-pi interaction
                    # Check two rings are roughly parallel (angle is smaller than 30 degree)
                    ring1_normal_vector = np.array(ring1.plane_coeff[0:3]) / norm(ring1.plane_coeff[0:3])
                    ring2_normal_vector = np.array(ring2.plane_coeff[0:3]) / norm(ring2.plane_coeff[0:3])
                    angle_between_planes = np.arccos(np.dot(ring1_normal_vector, ring2_normal_vector)) * 180 / np.pi
                    if abs(angle_between_planes - 0) <= pi_stack_angle_tolerance or \
                            abs(angle_between_planes - 180) <= pi_stack_angle_tolerance:
                        flag = False
                        for atom in ring1.atoms:
                            # Project the receptor atom onto the plane of the ligand ring
                            point_on_plane = project_point_onto_plane(atom.coordinate, ring2)
                            center_to_project_point = abs(norm(point_on_plane - ring2.center))

                            if center_to_project_point <= ring2.radius + pi_padding_dist:
                                flag = True
                                break

                        # if you've already determined it's a pi-pi stacking interaction, no need to keep trying
                        if flag is False:
                            for atom in ring2.atoms:
                                # Project the receptor atom onto the plane of the ligand ring
                                point_on_plane = project_point_onto_plane(atom.coordinate, ring1)
                                center_to_project_point = abs(norm(point_on_plane - ring1.center))
                                if center_to_project_point <= ring1.radius + pi_padding_dist:
                                    flag = True
                                    break

                        if flag is True:
                            pi_pi.append(("Ring " + str(ring1.index), "Ring " + str(ring2.index), 1))
                            pi_pi_info += "\n# {}{} Ring({}) :: LIG Ring({})\n".format(res.resid, res.type, ring1.index,
                                                                                       ring2.index)
                            pi_pi_info += "\n{}{} Ring({})\n".format(res.resid, res.type, ring1.index)
                            for atom in ring1.atoms:
                                pi_pi_info += "  " + str(atom) + "\n"

                            pi_pi_info += "\nLIG Ring({})\n".format(ring2.index)
                            for atom in ring2.atoms:
                                pi_pi_info += "  " + str(atom) + "\n"

                    elif abs(angle_between_planes - 90) <= t_stack_angle_tolerance or \
                            abs(angle_between_planes - 270) <= t_stack_angle_tolerance:
                        # Identify T-pi interaction
                        flag = False
                        for atom1 in ring1.atoms:
                            for atom2 in ring2.atoms:
                                dist = abs(norm(atom1.coordinate - atom2.coordinate))
                                if dist < t_stack_closest_cutoff:
                                    flag = True
                                    break

                        # The two rings come within 5 angstrom of each other
                        if flag:
                            point_on_receptor_plane = project_point_onto_plane(ring2.center, ring1)
                            point_on_ligand_plane = project_point_onto_plane(ring1.center, ring2)

                            # If it's a true pi-T interaction, this projected point should fall within the ring whose
                            # plane it's been projected into.
                            if abs(norm(point_on_receptor_plane - ring1.center)) <= ring1.radius + pi_padding_dist or \
                                    abs(norm(point_on_ligand_plane - ring2.center)) <= ring2.radius + pi_padding_dist:
                                pi_t.append(("Ring " + str(ring1.index), "Ring " + str(ring2.index), 1))
                                pi_t_info += "\n# {}{} Ring({}) :: LIG Ring({})\n".format(res.resid, res.type,
                                                                                          ring1.index, ring2.index)
                                pi_t_info += "\n{}{} Ring({})\n".format(res.resid, res.type, ring1.index)
                                for atom in ring1.atoms:
                                    pi_t_info += "  " + str(atom) + "\n"
                                pi_t_info += "\nLIG Ring({})\n".format(ring2.index)
                                for atom in ring2.atoms:
                                    pi_t_info += "  " + str(atom) + "\n"

        if pi_pi:
            pi_pi_interaction.append((res.segid, res.resid, res.type, pi_pi))
        if pi_t:
            pi_t_interaction.append((res.segid, res.resid, res.type, pi_t))

    # Calculate pi-cation interaction
    pi_cation_interaction = []
    for res in prot:
        pi_cation = []
        # Analyze the ring on protein
        if res.type in ["PHE", "TYR", "TRP", "HIS", "HSD", "HSE", "HSP"]:
            flag = 0
            for ring in res.rings:
                for atom in lig.atoms:
                    atomname = guess.guess_atom_name(atom.name)
                    if atom.charge > 0 and atomname in ["NA", "MG", "CA"]:
                        if abs(norm(ring.center - atom.coordinate)) <= cation_stack_closest_cutoff:
                            flag = 1
                            pi_cation.append(("Ring " + str(ring.index), atom, 1))
                            pi_cation_info += "\n# {} {}{} Ring({}) :: LIG Atom({})\n".format(res.segid, res.resid, res.type, ring.index, atom.name)
                            pi_cation_info += "\n{} {}{} Ring({})\n".format(res.segid, res.resid, res.type, ring.index)
                            for ring_atom in ring.atoms:
                                pi_cation_info += "  " + str(ring_atom) + "\n"
                            pi_cation_info += "\nLIG Atom({})\n".format(atom.name)
                            pi_cation_info += "  " + str(atom) + "\n"
            if flag:
                pi_cation_interaction.append((res.segid, res.resid, res.type, pi_cation))
        # Analyze the ring on ligand
        elif res.type in ["LYS", "ARG"]:
            flag = 0
            for ring in lig.rings:
                for atom in res.atoms:
                    if atom.name in ["NH1", "NH2", "NZ"]:
                        if abs(norm(ring.center - atom.coordinate)) <= cation_stack_closest_cutoff:
                            project_point = project_point_onto_plane(atom.coordinate, ring)
                            if abs(norm(project_point - ring.center)) <= ring.radius + pi_padding_dist:
                                flag = 1
                                pi_cation.append((atom, "Ring " + str(ring.index), 1))
                                pi_cation_info += "\n# {} {}{} Atom({}) :: LIG Ring({})\n".format(res.segid, res.resid, res.type, atom.index, ring.index)
                                pi_cation_info += "\n{} {}{} Atom({})\n".format(res.segid, res.resid, res.type,
                                                                                atom.index)
                                pi_cation_info += "  " + str() + "\n"
                                pi_cation_info += "\nLIG Ring({})\n".format(ring.index)
                                for ring_atom in ring.atoms:
                                    pi_cation_info += "  " + str(ring_atom) + "\n"
            if flag:
                pi_cation_interaction.append((res.segid, res.resid, res.type, pi_cation))

    if log:
        message = output_info1 + pi_pi_info + pi_t_info + pi_cation_info
        return pi_pi_interaction, pi_t_interaction, pi_cation_interaction, message
    else:
        return pi_pi_interaction, pi_t_interaction, pi_cation_interaction


def geometry_filter(lig):
    """
    :param lig:
    :return:
    """

    lig_coor = np.zeros(3)
    for atom in lig.atoms:
        lig_coor += atom.position
    lig_center = lig_coor/len(lig.atoms)

    n = 0   # total number of atoms
    total = 0

    for atom in lig.atoms:
        dist = abs(norm(lig_center - atom.position))
        n += 1
        total += dist
    variety = total / n

    f = open('ligand_geometry.txt', 'a')
    f.write(os.path.basename(lig.filename) + "  " + str(variety) + "\n")
    f.close()

    return os.path.basename(lig.filename), variety

