# -*- coding: utf-8 -*-


"""Molecular constitutional and topological indices."""

from rdkit import Chem
from rdkit.Chem import Lipinski as LPK


class Constitution:
    """Constitutional descriptors."""

    @staticmethod
    def calculate_heavy_mol_weight(mol: Chem.Mol) -> float:
        """calculate_ molecular weight of heavy atoms."""
        MolWeight = 0
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 1:
                MolWeight += atom.GetMass()
        return MolWeight

    @staticmethod
    def calculate_average_mol_weight(mol: Chem.Mol) -> float:
        """calculate_ average molecular weight of heavy atoms."""
        MolWeight = 0
        for atom in mol.GetAtoms():
            MolWeight += atom.GetMass()
        return MolWeight / mol.GetNumAtoms()

    @staticmethod
    def calculate_hydrogen_number(mol: Chem.Mol) -> float:
        """calculate_ number of Hydrogens."""
        i = 0
        Hmol = Chem.AddHs(mol)
        for atom in Hmol.GetAtoms():
            if atom.GetAtomicNum() == 1:
                i += 1
        return i

    @staticmethod
    def calculate_halogen_number(mol: Chem.Mol) -> float:
        """calculate_ number of Halogens."""
        i = 0
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() in [9, 17, 35, 53]:
                i += 1
        return i

    @staticmethod
    def calculate_hetero_number(mol: Chem.Mol) -> float:
        """calculate_ number of Heteroatoms."""
        i = 0
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() not in [1, 6]:
                i += 1
        return mol.GetNumAtoms() - i

    @staticmethod
    def calculate_heavy_atom_number(mol: Chem.Mol) -> float:
        """calculate_ number of Heavy atoms."""
        return mol.GetNumHeavyAtoms()

    @staticmethod
    def _calculate_element_number(mol: Chem.Mol, Atomic_number=6) -> float:
        """calculate_ number of atoms with specified element."""
        i = 0
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == Atomic_number:
                i += 1
        return i

    @staticmethod
    def calculate_fluorine_number(mol: Chem.Mol) -> float:
        """calculate_ number of Fluorine atoms."""
        return Constitution._calculate_element_number(mol, Atomic_number=9)

    @staticmethod
    def calculate_chlorine_number(mol: Chem.Mol) -> float:
        """calculate_ number of Fluorine atoms."""
        return Constitution._calculate_element_number(mol, Atomic_number=17)


    def calculate_bromine_number(mol: Chem.Mol) -> float:
        """calculate_ number of Bromine atoms."""
        return Constitution._calculate_element_number(mol, Atomic_number=35)

    @staticmethod
    def calculate_iodine_number(mol: Chem.Mol) -> float:
        """calculate_ number of Iodine atoms."""
        return Constitution._calculate_element_number(mol, Atomic_number=53)

    @staticmethod
    def calculate_carbon_number(mol: Chem.Mol) -> float:
        """calculate_ number of Carbon atoms."""
        return Constitution._calculate_element_number(mol, Atomic_number=6)

    @staticmethod
    def calculate_phosphor_number(mol: Chem.Mol) -> float:
        """Calcualtion number of Phosphor atoms."""
        return Constitution._calculate_element_number(mol, Atomic_number=15)

    @staticmethod
    def calculate_sulfur_number(mol: Chem.Mol) -> float:
        """calculate_ number of Sulfur atoms."""
        return Constitution._calculate_element_number(mol, Atomic_number=16)

    @staticmethod
    def calculate_oxygen_number(mol: Chem.Mol) -> float:
        """calculate_ number of Oxygen atoms."""
        return Constitution._calculate_element_number(mol, Atomic_number=8)

    @staticmethod
    def calculate_nitrogen_number(mol: Chem.Mol) -> float:
        """calculate_ number of Nitrogen atoms."""
        return Constitution._calculate_element_number(mol, Atomic_number=7)

    @staticmethod
    def calculate_ring_number(mol: Chem.Mol) -> float:
        """calculate_ number of rings."""
        return len(Chem.GetSSSR(mol))

    @staticmethod
    def calculate_rot_bond_number(mol: Chem.Mol) -> float:
        """calculate_ number of rotatable bonds."""
        return LPK.NumRotatableBonds(mol)

    @staticmethod
    def calculate_hdonor_number(mol: Chem.Mol) -> float:
        """calculate_ number of Hydrongen bond donors."""
        return LPK.NumHDonors(mol)

    @staticmethod
    def calculate_hacceptor_number(mol: Chem.Mol) -> float:
        """calculate_ number of Hydrogen bond acceptors."""
        return LPK.NumHAcceptors(mol)

    @staticmethod
    def calculate_singlebond_number(mol: Chem.Mol) -> float:
        """calculate_ number of single bonds."""
        i = 0
        for bond in mol.GetBonds():
            if bond.GetBondType().name == 'SINGLE':
                i += 1
        return i

    @staticmethod
    def calculate_doublebond_number(mol: Chem.Mol) -> float:
        """calculate_ number of double bonds."""
        i = 0
        for bond in mol.GetBonds():
            if bond.GetBondType().name == 'DOUBLE':
                i += 1
        return i

    @staticmethod
    def calculate_triplebond_number(mol: Chem.Mol) -> float:
        """calculate_ number of triple bonds."""
        i = 0
        for bond in mol.GetBonds():
            if bond.GetBondType().name == 'TRIPLE':
                i += 1
        return i

    @staticmethod
    def calculate_aromaticbond_number(mol: Chem.Mol) -> float:
        """calculate_ number of aromatic bonds."""
        i = 0
        for bond in mol.GetBonds():
            if bond.GetBondType().name == 'AROMATIC':
                i += 1
        return i

    @staticmethod
    def calculate_allatom_number(mol: Chem.Mol) -> float:
        """calculate_ number of all atoms."""
        return Chem.AddHs(mol).GetNumAtoms()

    @staticmethod
    def _calculate_path_n(mol: Chem.Mol, path_Length=2) -> float:
        """calculate_ number of path of length N."""
        return len(Chem.FindAllPathsOfLengthN(mol, path_Length, useBonds=1))

    @staticmethod
    def calculate_path_1(mol: Chem.Mol) -> float:
        """calculate_ number of path length of 1."""
        return Constitution._calculate_path_n(mol, 1)

    @staticmethod
    def calculate_path_2(mol: Chem.Mol) -> float:
        """calculate_ number of path length of 2."""
        return Constitution._calculate_path_n(mol, 2)

    @staticmethod
    def calculate_path_3(mol: Chem.Mol) -> float:
        """calculate_ number of path length of 3."""
        return Constitution._calculate_path_n(mol, 3)

    @staticmethod
    def calculate_path_4(mol: Chem.Mol) -> float:
        """calculate_ number of path length of 4."""
        return Constitution._calculate_path_n(mol, 4)

    @staticmethod
    def calculate_path_5(mol: Chem.Mol) -> float:
        """calculate_ number of path length of 5."""
        return Constitution._calculate_path_n(mol, 5)

    @staticmethod
    def calculate_path_6(mol: Chem.Mol) -> float:
        """calculate_ number of path length of 6."""
        return Constitution._calculate_path_n(mol, 6)

    @staticmethod
    def get_all(mol: Chem.Mol) -> dict:
        """Get all (30) constitutional descriptors."""
        result = {}
        for DesLabel in _constitutional.keys():
            result[DesLabel] = _constitutional[DesLabel](mol)
        return result


_constitutional = {'Weight': Constitution.calculate_heavy_mol_weight,
                   'nH': Constitution.calculate_hydrogen_number,
                   'nHal': Constitution.calculate_halogen_number,
                   'nHet': Constitution.calculate_hetero_number,
                   'nHA': Constitution.calculate_heavy_atom_number,
                   'nF': Constitution.calculate_fluorine_number,
                   'nCl': Constitution.calculate_chlorine_number,
                   'nBr': Constitution.calculate_bromine_number,
                   'nI': Constitution.calculate_iodine_number,
                   'nC': Constitution.calculate_carbon_number,
                   'nP': Constitution.calculate_phosphor_number,
                   'nS': Constitution.calculate_oxygen_number,
                   'nO': Constitution.calculate_oxygen_number,
                   'nN': Constitution.calculate_nitrogen_number,
                   'nRing': Constitution.calculate_ring_number,
                   'nRotB': Constitution.calculate_rot_bond_number,
                   'nHBD': Constitution.calculate_hdonor_number,
                   'nHBA': Constitution.calculate_hacceptor_number,
                   'nSBond': Constitution.calculate_singlebond_number,
                   'nDBond': Constitution.calculate_doublebond_number,
                   'nTBond': Constitution.calculate_triplebond_number,
                   'nAroBond': Constitution.calculate_aromaticbond_number,
                   'nAtom': Constitution.calculate_allatom_number,
                   'AWeight': Constitution.calculate_average_mol_weight,
                   'path_L1': Constitution.calculate_path_1,
                   'path_L2': Constitution.calculate_path_2,
                   'path_L3': Constitution.calculate_path_3,
                   'path_L4': Constitution.calculate_path_4,
                   'path_L5': Constitution.calculate_path_5,
                   'path_L6': Constitution.calculate_path_6,
                   }
