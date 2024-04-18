import os
import re
import numpy as np
import MDAnalysis as mda
from . import guess
from .topology.models import *
from .topology.bonds import *
from .topology.masses import get_atom_mass


class Atom:
    def __init__(self, segid='X', resid='N/A', resn='N/A', atomname='C', atomtype='N/A',
                 coordinate=np.array([0., 0., 0.]), mass=12.0, charge=0.0, hbond="none", hydro="", index=0, atom_id=0): ## 2022/08/16 bruce amended 
        self.__segid = segid
        self.__resid = resid
        self.__resn = resn
        self.__type = atomtype
        self.__name = atomname
        self.__charge = charge
        self.__mass = mass
        self.__vdw_r = vdw_radius(self.__name)
        self.__index = index
        self.__id = atom_id  ## 2022/08/16 bruce amended 
        self.__list_of_connecting_atoms = []
        if isinstance(coordinate, np.ndarray):
            self.__coordinate = coordinate
        else:
            self.__coordinate = np.array([0., 0., 0.])
        self.hbond = hbond
        self.hydro = hydro

    def __repr__(self):
        return "<Atom> segid " + str(self.__segid) + ", resid " + str(self.__resid) + ", resn " + str(self.__resn) + \
               ", index " + str(self.__index) + ", id " + str(self.__id) + ", name " + self.__name + ", type " + self.__type + \
               ", coordinate " + \
               "[{0:2.3f}, {1:2.3f}, {2:2.3f}]".format(self.coordinate[0], self.coordinate[1], self.coordinate[2]) + \
               ", charge " + str(self.__charge)  ## 2022/08/16 amended by Bruce

    def __str__(self):
        return "<Atom> segid " + str(self.__segid) + ", resid " + str(self.__resid) + ", resn " + str(self.__resn) + \
               ", index " + str(self.__index) + ", id " + str(self.__id) + ", name " + self.__name + ", type " + self.__type + \
               ", coordinate " + \
               "[{0:2.3f}, {1:2.3f}, {2:2.3f}]".format(self.coordinate[0], self.coordinate[1], self.coordinate[2]) + \
               ", charge " + str(self.__charge)  ## 2022/08/16 amended by Bruce

    @property
    def type(self):
        return self.__type

    @type.setter
    def type(self, value):
        self.__type = value

    @property
    def name(self):
        return self.__name

    @name.setter
    def name(self, value):
        self.__name = value

    @property
    def charge(self):
        return self.__charge

    @charge.setter
    def charge(self, value):
        self.__charge = value

    @property
    def coordinate(self):
        return self.__coordinate

    @coordinate.setter
    def coordinate(self, new_coor):
        if isinstance(new_coor, np.ndarray):
            self.__coordinate = new_coor
        else:
            print("Please input a np.ndarray data")

    @property
    def mass(self):
        return self.__mass

    @property
    def vdw_r(self):
        return self.__vdw_r

    @property
    def index(self):
        return self.__index
    
    @property  ## 2022/08/16 bruce amended 
    def id(self):
        return self.__id 
    
    @property
    def center_of_mass(self):
        return self.__coordinate

    @property
    def list_of_connecting_atoms(self):
        return self.__list_of_connecting_atoms

    def add_neighbor_atom(self, neighbor_atom):
        if not (neighbor_atom in self.__list_of_connecting_atoms):
            self.__list_of_connecting_atoms.append(neighbor_atom)



class AromaticRing:
    def __init__(self, ring=[], index=0):
        self.atoms = ring
        self.index = index
        self.center = self.ring_center()
        self.radius = self.radius()
        self.plane_coeff = self.plane_coeff()

    def __len__(self):
        return len(self.atoms)

    def __repr__(self):
        return "<AromaticRing> " + str(len(self.atoms)) + " atoms ring"

    def __str__(self):
        return "<AromaticRing> " + str(len(self.atoms)) + " atoms ring"

    def ring_center(self):
        count = 0
        coor = np.array([0., 0., 0.])
        for atom in self.atoms:
            coor += atom.coordinate
            count += 1
        return coor / count

    def radius(self):
        radius = float('-inf')
        for atom in self.atoms:
            dist = abs(np.linalg.norm(atom.coordinate - self.center))
            if dist > radius:
                radius = dist
        return radius

    def plane_coeff(self):
        if len(self.atoms) > 4:
            a = self.atoms[0].coordinate
            b = self.atoms[2].coordinate
            c = self.atoms[4].coordinate
        else:
            return np.zeros(4)

        normal_vector = np.cross(b-a, c-a)
        d = np.dot(a, normal_vector)
        coeff = list(normal_vector)
        coeff.append(d)
        return coeff


class Residue:
    def __init__(self, resn='N/A', resid='N/A', segid='N/A', atoms=[]):
        self.resn = resn
        self.resid = resid
        self.segid = segid
        self.atoms = atoms
        self.__rings = self.add_aromatic_ring()

    def __len__(self):
        return len(self.atoms)

    def __repr__(self):
        return "<Residue> segid " + self.segid + ", resid " + str(self.resid) + ", resn " + str(self.resn) + " with " \
               + str(len(self.atoms)) + " atoms"

    def __str__(self):
        return "<Residue> segid " + self.segid + ", resid " + str(self.resid) + ", resn " + str(self.resn) + " with " \
               + str(len(self.atoms)) + " atoms"

    @property
    def coordinates(self):
        coor = []
        for atom in self.atoms:
            coor.append(atom.coordinate)

        return coor

    @property
    def rings(self):
        return self.__rings

    @property
    def center_of_mass(self):
        coor = np.zeros(3)
        total_mass = 0
        for atom in self.atoms:
            total_mass += atom.mass
            coor[0] += atom.mass * atom.coordinate[0]
            coor[1] += atom.mass * atom.coordinate[1]
            coor[2] += atom.mass * atom.coordinate[2]
        return coor/total_mass

    def main_chain(self):
        main_sele, _ = get_residue_main_and_side_chain_unity_model(self.resn)

        main_chain = []
        for atom in self.atoms:
            if atom.name in main_sele:
                main_chain.append(atom)

        return Residue(resn=self.resn + "(S)", atoms=main_chain)

    def side_chain(self):
        _, side_sele = get_residue_main_and_side_chain_unity_model(self.resn)

        side_chain = []
        for atom in self.atoms:
            if atom.name in side_sele:
                side_chain.append(atom)

        return Residue(resn=self.resn + "(S)", atoms=side_chain)

    def assign_aromatic_rings(self):
        rings = []
        select = get_residue_ring_atoms_unity_model(self.resn)

        if select:
            if isinstance(select[0], str):      # Determine residue has one ring or more
                select = select,

            for ring_atom in select:
                ring = []
                for atom in self.atoms:
                    if atom.name in ring_atom:
                        ring.append(atom)
                rings.append(ring)

        return rings

    def add_aromatic_ring(self):
        rings = self.assign_aromatic_rings()

        aromatics = []
        for index, ring in enumerate(rings):
            aromatics.append(AromaticRing(ring, index))

        return aromatics

    def add_hydrophobic_label(self):
        for atom in self.atoms:
            if guess_atom_name(atom.name) in ["C", "F", "CL", "BR"]:
                self.__parse_hydrophobic_label(atom)

    def __parse_hydrophobic_label(self, atom):
        name = guess_atom_name(atom.name)
        if name == "C":
            atom.hydro = "H"
            for atomi in atom.list_of_connecting_atoms:
                if guess_atom_name(atomi.name) != "C":
                    atom.hydro = "P"
                    break
        elif name in ["F", "CL", "BR"]:
            atom.hydro = "H"
        elif name in ["S", "P"]:
            atom.hydro = "P"


class Protein:
    def __init__(self, prot_name="None", res=[]):
        self.name = prot_name
        self.__residues = res

    def __len__(self):
        return len(self.residues)

    def __repr__(self):
        return "<Protein> Name " + str(self.name) + " with " + str(len(self.__residues)) + " residues and " + \
               str(len(self.atoms)) + " atoms"

    def __str__(self):
        return "<Protein> Name " + str(self.name) + " with " + str(len(self.__residues)) + " residues and " + \
               str(len(self.atoms)) + " atoms"

    @property
    def residues(self):
        return self.__residues

    @property
    def atoms(self):
        atoms = []
        for res in self.__residues:
            atoms += res.atoms

        return atoms

    @property
    def coordinates(self):
        coor = []
        for res in self.__residues:
            coor += res.coordinates

        return coor


class AtomGroup:
    def __init__(self, group_name='None', atoms=[]):
        self.__name = os.path.basename(group_name)
        self.__atoms = atoms
        self.__rings = self.add_aromatic_ring()
        self.add_hbond_label()
        self.add_hydrophobic_label()

    def __len__(self):
        return len(self.__atoms)

    def __repr__(self):
        return "<AtomGroup> Compound " + os.path.basename(self.__name) + " with " + str(len(self)) + " atoms"

    def __str__(self):
        return "<AtomGroup> Compound " + os.path.basename(self.__name) + " with " + str(len(self)) + " atoms"

    @property
    def atoms(self):
        return self.__atoms

    @property
    def name(self):
        return os.path.basename(self.__name)

    @property
    def rings(self):
        return self.__rings

    @property
    def center_of_mass(self):
        coor = np.zeros(3)
        total_mass = 0
        for atom in self.__atoms:
            total_mass += atom.mass
            coor[0] += atom.mass * atom.coordinate[0]
            coor[1] += atom.mass * atom.coordinate[1]
            coor[2] += atom.mass * atom.coordinate[2]
        return coor/total_mass

    @property
    def coordinates(self):
        coor = []
        for atom in self.atoms:
            coor.append(atom.coordinate)

        return coor

    def read_cgenff_property(self, file_path=None):
        if file_path is None:
            file_path, _ = self.name.replace(".pdb", ".str")

        try:
            file = open(file_path, 'r')
        except IOError:
            print("Warning: There is no cgenff information of " + self.name)
        else:
            m = re.compile(r'^ATOM')
            all_info = []
            for line in file:
                if m.match(line) is not None:
                    info = []
                    for s in line.split(' '):
                        if s != '' and s != '!':
                            info.append(s)
                    all_info.append(info)

            for index, atom in enumerate(self.__atoms):
                _, _, atom.type, atom.charge, _ = all_info[index]

            file.close()

    def ring_recursive(self, atom, origin_atom, ring=[], all_rings=[], index=1):
        # Discard ring which contains more than six atoms
        if index > 6:
            return

        for connected_atom in atom.list_of_connecting_atoms:
            temp = ring.copy()
            if connected_atom not in ring:
                temp.append(connected_atom)
                self.ring_recursive(connected_atom, origin_atom, temp, all_rings, index+1)
            if connected_atom is origin_atom and index > 4:     # Consider the ring which contains more than four atoms
                all_rings.append(ring)

    def assign_aromatic_rings(self):

        # Get all the rings containing each atom in ligand
        all_rings = []
        for atom in self.atoms:
            self.ring_recursive(atom, atom, ring=[atom], all_rings=all_rings)

        # Delete redundant rings
        for i, ring1 in enumerate(all_rings):
            if len(ring1) != 0:
                for j, ring2 in enumerate(all_rings[i+1::]):
                    if len(ring2) != 0:
                        set1 = set(ring1)
                        set2 = set(ring2)
                        if set2.issubset(set1):
                            all_rings[i+j+1] = []

        while [] in all_rings:
            all_rings.remove([])

        # Figure out which of these rings are planar
        for i, ring in enumerate(all_rings):
            is_flat = True

            if len(ring) > 4:
                times = len(ring) - 1
            else:
                break

            for j in range(times):

                # first, let's see if the atom is a carbon connected to four atoms.
                # That would be a quick way of telling this is not an aromatic ring
                if ring[j].name == "C" and len(ring[j].list_of_connecting_atoms) > 4:
                    is_flat = False
                    break

                # now check the dihedral between the ring atoms to see if it's flat
                #  A
                #   \
                #     B - C
                #           \
                #            D
                pt1 = ring[j-3]
                pt2 = ring[j-2]
                pt3 = ring[j-1]
                pt4 = ring[j]
                dih = dihedral_angle(pt1.coordinate, pt2.coordinate, pt3.coordinate, pt4.coordinate)

                if (dih > -165 and dih < -15) or (dih > 15 and dih < 165):
                    is_flat = False
                    break

                if len(ring[j].list_of_connecting_atoms) == 3:
                    pt1 = ring[j].list_of_connecting_atoms[0]
                    pt2 = ring[j].list_of_connecting_atoms[1]
                    pt3 = ring[j].list_of_connecting_atoms[2]
                    dih = dihedral_angle(pt1.coordinate, pt4.coordinate, pt3.coordinate, pt2.coordinate)

                    if (dih > -165 and dih < -15) or (dih > 15 and dih < 165):
                        is_flat = False
                        break

                if not is_flat:
                    break

            if not is_flat:
                all_rings[i] = []

        # Delete non-plane rings
        while [] in all_rings:
            all_rings.remove([])

        return all_rings

    def add_aromatic_ring(self):
        rings = self.assign_aromatic_rings()

        aromatics = []
        for index, ring in enumerate(rings):
            aromatics.append(AromaticRing(ring, index=index))

        return aromatics

    def add_hbond_label(self):
        for atom in self.atoms:
            if atom.name in ["N", "O", "H", "F", "CL", "BR"]:
                self.__parse_hbond_label(atom)

    def __parse_hbond_label(self, atom):
        # Assign hydrogen bond label
        if atom.type in ["N", "O"]:
            atom.hbond = "none"
            for atomi in atom.list_of_connecting_atoms:
                if atomi.type == "HD":
                    atom.hbond = "D"
                    break
        elif atom.type in ["NA", "OA"]:
            atom.hbond = "A"
            for atomi in atom.list_of_connecting_atoms:
                if atomi.type == "HD":
                    atom.hbond = "B"
                    break
        elif atom.type == "HD":
            atom.hbond = "D"
        elif atom.type == "F" or atom.type == "CL" or atom.type == "BR":
            atom.hbond = "A"
        else:
            atom.hbond = "none"

    def add_hydrophobic_label(self):
        for atom in self.atoms:
            if atom.name in ["C", "F", "CL", "BR"]:
                self.__parse_hydrophobic_label(atom)

    def __parse_hydrophobic_label(self, atom):
        if atom.name == "C":
            atom.hydro = "H"
            for atomi in atom.list_of_connecting_atoms:
                if atomi.name != "C":
                    atom.hydro = "P"
                    break
        elif atom.name in ["F", "CL", "BR"]:
            atom.hydro = "H"
        elif atom.name in ["S", "P"]:
            atom.hydro = "P"


def three_point_angle(pt1, pt2, pt3):
    v12 = pt1 - pt2
    u12 = v12 / np.linalg.norm(v12)
    v32 = pt3 - pt2
    u32 = v32 / np.linalg.norm(v32)
    angle = abs(np.arccos(np.dot(u12, u32))) * 180 / np.pi

    return angle


def dihedral_angle(pt1, pt2, pt3, pt4):
    v12 = pt1 - pt2
    v32 = pt3 - pt2
    v23 = pt2 - pt3
    v43 = pt4 - pt3
    v12Xv32 = np.cross(v12, v32)
    u1 = v12Xv32 / np.linalg.norm(v12Xv32)
    v23Xv43 = np.cross(v23, v43)
    u2 = v23Xv43 / np.linalg.norm(v23Xv43)
    cos = np.dot(u1, u2)
    if cos > 1 or cos < -1:
        if abs(cos) - 1 > 0.005:
            print("Warning: Dihedral angle calculation have something wrong!!")
            dihedral = 0
        else:
            if cos > 1:
                cos = 1
            else:
                cos = -1
            dihedral = abs(np.arccos(cos)) * 180 / np.pi
    else:
        dihedral = abs(np.arccos(cos)) * 180 / np.pi

    return dihedral


def geometry_center(u, residues):
    """
    Description:
        Calculate the geometry center of selected amino acids

    :param u: universe of protein
    :param residues: list of residue number
    :return: coordinate of geometry center
    """
    count = 0
    coor = np.array([0., 0., 0.])
    for resid in residues:
        atoms = u.select_atoms('resid ' + str(resid))
        for atom in atoms:
            coor += atom.position
            count += 1
    return coor / count


def read_protein_unity_model(prot_path, select_segid=[], select_active_site=[]):
    u = mda.Universe(prot_path, prot_path)

    if select_active_site:
        active_site_sele = " ".join(str(i) for i in select_active_site)
        cave = u.select_atoms('byres sphzone 20 resid ' + active_site_sele)
    else:
        cave = u

    if select_segid:
        segid_sele = " ".join(i for i in select_segid)
        cave = cave.select_atoms('segid ' + segid_sele)

    seq = []
    for segment in cave.segments:
        segid = segment.segid
        seg = cave.select_atoms('segid ' + segid)
        for i, res in enumerate(seg.residues):
            resn = res.resname
            resi = res.resid

            if resn in ["F", "CL", "BR", "I", "MG", "ZN", "CA"]:
                atoms = seg.select_atoms('resid ' + str(resi))
                res_info = non_amino_info[resn].copy()
            else:
                res_info = all_atom_model_info[res.resname].copy()
                if i == 0:
                    for key in all_atom_model_info["NTER"].keys():
                        res_info[key] = all_atom_model_info["NTER"][key]
                elif i == len(seg.residues) - 1:
                    for key in all_atom_model_info["CTER"].keys():
                        res_info[key] = all_atom_model_info["CTER"][key]

                atoms = seg.select_atoms('resid ' + str(resi) + ' and (' + all_atom_model_atom_selection[resn] + \
                                         ' or ' + all_atom_model_atom_selection["NCTER"] + ")")

            info = []
            for atom in atoms:
                atomindex = atom.index  ## 2022/08/16 amended by Bruce 
                atom_id = atom.id  ## 2022/08/16 amended by Bruce
                atomname = atom.name
                atomtype, hbond, hydro, charge = res_info[atomname]
                coor = atom.position.astype('float64')
                mass = atom.mass

                a = Atom(segid=segid, resid=resi, resn=resn, atomname=atomname, atomtype=atomtype, coordinate=coor, \
                          mass=mass, charge=charge, hbond=hbond, hydro=hydro, index=atomindex, atom_id=atom_id)  ## 2022/08/16 amended
                info.append(a)

            seq.append(Residue(resn, resi, segid, info))

    if select_active_site:
        return Protein(res=seq), geometry_center(u, select_active_site)
    else:
        return Protein(res=seq)

def read_protein_all_atom_model(prot_path, select_segid=[], select_active_site=[]):
    u = mda.Universe(prot_path, prot_path)

    if select_active_site:
        active_site_sele = " ".join(str(i) for i in select_active_site)
        cave = u.select_atoms('byres sphzone 20 resid ' + active_site_sele)
    else:
        cave = u

    if select_segid:
        segid_sele = " ".join(i for i in select_segid)
        cave = cave.select_atoms('segid ' + segid_sele)

    seq = []
    for segment in cave.segments:
        segid = segment.segid
        seg = cave.select_atoms('segid ' + segid)
        for i, res in enumerate(seg.residues):
            resn = res.resname
            resi = res.resid

            if resn in ["F", "CL", "BR", "I", "MG", "ZN", "CA"]:
                atoms = seg.select_atoms('resid ' + str(resi))
                res_info = non_amino_info[resn].copy()
            else:
                # res_info = unity_model_info[res.resname].copy()
                res_info = all_atom_model_info[resn].copy()
                if i == 0:
                    #for key in unity_model_info["NTER"].keys():
                    #    res_info[key] = unity_model_info["NTER"][key]
                    for key in all_atom_model_info["NTER"].keys():
                        res_info[key] = all_atom_model_info["NTER"][key]
                
                elif i == len(seg.residues) - 1:
                    #or key in unity_model_info["CTER"].keys():
                    #    res_info[key] = unity_model_info["CTER"][key]
                    for key in all_atom_model_info["CTER"].keys():
                        res_info[key] = all_atom_model_info["CTER"][key]
                
                #atoms = seg.select_atoms('resid ' + str(resi) + ' and (' + unity_model_atom_selection[resn] +
                #                         ' or ' + unity_model_atom_selection["NCTER"] + ")")
                atoms = seg.select_atoms('resid ' + str(resi) + ' and (' + all_atom_model_atom_selection[resn] + ")") 
                                         #+' or ' + all_atom_model_atom_selection["NCTER"] + ")")

            info = []
            for atom in atoms:
                #atomindex = atom.index
                atomindex = atom.index  ## 2022/08/16 amended by Bruce 
                atom_id = atom.id  ## 2022/08/16 amended by Bruce
                atomname = atom.name
                atomtype, hbond, hydro, charge = res_info[atomname]
                coor = atom.position.astype('float64')
                mass = atom.mass

                a = Atom(segid=segid, resid=resi, resn=resn, atomname=atomname, atomtype=atomtype, coordinate=coor,
                         mass=mass, charge=charge, hbond=hbond, hydro=hydro, index=atomindex, atom_id=atom_id)  ## 2022/08/16 amended by Bruce
                info.append(a)

            seq.append(Residue(resn, resi, segid, info))

    if select_active_site:
        return Protein(res=seq), geometry_center(u, select_active_site)
    else:
        return Protein(res=seq)

def read_protein_cg(prot_path, active_site=[]):
    u = mda.Universe(prot_path, prot_path)

    if active_site:
        active_site_sele = " ".join(str(i) for i in active_site)
        cave = u.select_atoms('byres sphzone 20 resid ' + active_site_sele)
    else:
        cave = u

    seq = []
    for res in cave.residues:
        resn = res.resname
        resi = res.resid

        atoms = cave.select_atoms('resid ' + str(resi) + ' and ' + cg_model_atom_selection[resn])
        res_info = cg_model_info[res.resname]

        info = []
        for index, atom in enumerate(atoms):
            a = Atom(atomname=atom.type, atomtype=res_info[index][0], coordinate=atom.position, mass=atom.mass,
                     charge=res_info[index][3], hbond=res_info[index][1], index=index)
            info.append(a)

        seq.append(Residue(resn, resi, info))

    if active_site:
        return seq, geometry_center(u, active_site)
    else:
        return seq


def read_atom_group(file_path):
    u = mda.Universe(file_path, file_path)

    name = u.filename
    #atoms = u.select_atoms('not name H*')

    info = []
    for index, atom in enumerate(u.atoms):
        atomname = guess_atom_name(atom.name) ## 2022/09/21 amended by Bruce
        atomtype = atom.type
        charge = atom.charge
        coor = atom.position.astype('float64')
        atom_index = atom.index  ## 2022/08/16 amended by Bruce
        atom_id = atom.id  ## 2022/08/16 amended by Bruce
        if atom.mass == 0.00:
            mass = get_atom_mass(atomname)
        else:
            mass = atom.mass

        a = Atom(atomname=atomname, atomtype=atomtype, coordinate=coor, mass=mass,
                 charge=charge, index=atom_index, atom_id=atom_id)
        info.append(a)

    add_bond_information(info)
    group = AtomGroup(name, info)

    return group


read_ligand = read_atom_group   # alias of function


def add_bond_information(atoms):
    for atom1 in atoms:
        for atom2 in atoms:
            if atom1 is not atom2:
                diff = atom1.coordinate - atom2.coordinate
                dist = np.linalg.norm(diff.astype("float64"))

                atomname1 = guess.guess_atom_name(atom1.name)
                atomname2 = guess.guess_atom_name(atom2.name)

                if dist <= bond_length(atomname1, atomname2) * 1.2:
                    atoms[atom1.index].add_neighbor_atom(atoms[atom2.index])
                    atoms[atom2.index].add_neighbor_atom(atoms[atom1.index])


def add_hbond_information(info):
    for atom in info:
        type_info = "none"
        # None
        if atom.name == "C": type_info = "none"

        # Donor
        if atom.name == "H":
            atom2 = atom.list_of_connecting_atoms[0]
            if atom2 == "O" or atom2 == "N":
                type_info = "D"
            else:
                type_info = "none"

        # Acceptor
        if atom.name == "N":
            if len(atom.list_of_connecting_atoms) <= 2:
                flag = 0
                for atom2 in atom.list_of_connecting_atoms:
                    dist = abs(np.linalg.norm(atom2.coordinate - atom.coordinate))
                    if dist < 1.30:     # N=C
                        flag = 1
                        break

                if flag:
                    type_info = "A"
                else:
                    type_info = "D"
            else:
                type_info = "A"
        if atom.name == "P": type_info = "A"
        if atom.name == "S": type_info = "A"
        if atom.name == "CL": type_info = "A"
        if atom.name == "BR": type_info = "A"
        if atom.name == "O":
            if len(atom.list_of_connecting_atoms) == 1:
                num_of_connected_carbon = 0
                num_of_connected_oxygen = 0
                num_of_others = 0
                for atom2 in atom.list_of_connecting_atoms[0].list_of_connecting_atoms:
                    if atom2.name == "C":
                        num_of_connected_carbon += 1
                    elif atom2.name == "O":
                        num_of_connected_oxygen += 1
                    else:
                        num_of_others += 1

                if num_of_connected_carbon + num_of_others == 2:
                    dist = abs(np.linalg.norm(atom.coordinate - atom.list_of_connecting_atoms[0].coordinate))
                    if dist < 1.30:
                        type_info = "A"  # Keton  RC(=O)R'
                    else:
                        type_info = "D"
                elif num_of_connected_carbon == 1:
                    if num_of_connected_oxygen == 2:
                        type_info = "B"  # Carboxyl group C(=O)OH
                    else:
                        dist = abs(np.linalg.norm(atom.coordinate - atom.list_of_connecting_atoms[0].coordinate))
                        if dist < 1.30:
                            type_info = "A"  # C=O
                        else:
                            type_info = "B"  # C-OH
            else:
                type_info = "A"     # Ether R-O-R

        atom.hbond = type_info
