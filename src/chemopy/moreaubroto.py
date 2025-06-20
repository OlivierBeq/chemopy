# -*- coding: utf-8 -*-


"""Moreau-Broto autocorrelation descriptors."""


import numpy as np
from rdkit import Chem

from .atom_property import get_relative_atomic_property


class MoreauBroto:
    """Moreau-Broto autocorrelation descriptors."""

    @staticmethod
    def _calculate_moreau_broto_autocorrelation(mol: Chem.Mol, lag: int = 1, propertylabel: str = 'm') -> float:
        """Calculate weighted Moreau-Broto autocorrelation descriptors.

        :param lag: topological distance between atom i and atom j.
        :param propertylabel: type of weighted property
        """
        Natom = mol.GetNumAtoms()
        distance_matrix = Chem.GetDistanceMatrix(mol)
        res = 0.0
        for i in range(Natom):
            for j in range(Natom):
                if distance_matrix[i, j] == lag:
                    atom1 = mol.GetAtomWithIdx(i)
                    atom2 = mol.GetAtomWithIdx(j)
                    temp1 = get_relative_atomic_property(element=atom1.GetSymbol(), propertyname=propertylabel)
                    temp2 = get_relative_atomic_property(element=atom2.GetSymbol(), propertyname=propertylabel)
                    res = res + temp1 * temp2
                else:
                    res = res + 0.0
        return np.log(res / 2 + 1)

    @staticmethod
    def calculate_moreau_broto_auto_mass(mol: Chem.Mol) -> dict:
        """Calculate Moreau-Broto autocorrelation with carbon-scaled atomic mass."""
        res = {}
        for i in range(8):
            res[f'ATSm{i + 1}'] = MoreauBroto._calculate_moreau_broto_autocorrelation(mol, lag=i + 1, propertylabel='m')
        return res

    @staticmethod
    def calculate_moreau_broto_auto_volume(mol: Chem.Mol) -> dict:
        """Calculate Moreau-Broto autocorrelation with carbon-scaled atomic van der Waals volume."""
        res = {}
        for i in range(8):
            res[f'ATSv{i + 1}'] = MoreauBroto._calculate_moreau_broto_autocorrelation(mol, lag=i + 1, propertylabel='V')
        return res

    @staticmethod
    def calculate_moreau_broto_auto_electronegativity(mol: Chem.Mol) -> dict:
        """Calculate Moreau-Broto autocorrelation with carbon-scaled atomic Sanderson electronegativity."""
        res = {}
        for i in range(8):
            res[f'ATSe{i + 1}'] = MoreauBroto._calculate_moreau_broto_autocorrelation(mol, lag=i + 1, propertylabel='En')
        return res

    @staticmethod
    def calculate_moreau_broto_auto_polarizability(mol: Chem.Mol) -> dict:
        """Calculate Moreau-Broto autocorrelation with carbon-scaled atomic polarizability."""
        res = {}
        for i in range(8):
            res[f'ATSp{i + 1}'] = MoreauBroto._calculate_moreau_broto_autocorrelation(mol, lag=i + 1, propertylabel='alapha')
        return res

    @staticmethod
    def get_all(mol: Chem.Mol) -> dict:
        """Calculate all (32) Moreau-Broto autocorrelation descriptors."""
        res = {}
        res.update(MoreauBroto.calculate_moreau_broto_auto_mass(mol))
        res.update(MoreauBroto.calculate_moreau_broto_auto_volume(mol))
        res.update(MoreauBroto.calculate_moreau_broto_auto_electronegativity(mol))
        res.update(MoreauBroto.calculate_moreau_broto_auto_polarizability(mol))
        return res
