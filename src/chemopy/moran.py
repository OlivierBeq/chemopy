# -*- coding: utf-8 -*-


"""Moran autocorrelation descriptors."""

import numpy as np
from rdkit import Chem

from .atom_property import get_relative_atomic_property


class Moran:
    """Moran autocorrelation descriptors."""

    @staticmethod
    def _calculate_moran_autocorrelation(mol: Chem.Mol, lag: int = 1, propertylabel: str = 'm') -> float:
        """Calculate weighted Moran autocorrelation descriptors.

        :param lag: topological distance between atom i and atom j.
        :param propertylabel: type of weighted property
        """
        Natom = mol.GetNumAtoms()
        prolist = []
        for i in mol.GetAtoms():
            temp = get_relative_atomic_property(i.GetSymbol(), propertyname=propertylabel)
            prolist.append(temp)
        aveweight = sum(prolist) / Natom
        tempp = [np.square(x - aveweight) for x in prolist]
        GetDistanceMatrix = Chem.GetDistanceMatrix(mol)
        res = 0.0
        index = 0
        for i in range(Natom):
            for j in range(Natom):
                if GetDistanceMatrix[i, j] == lag:
                    atom1 = mol.GetAtomWithIdx(i)
                    atom2 = mol.GetAtomWithIdx(j)
                    temp1 = get_relative_atomic_property(element=atom1.GetSymbol(), propertyname=propertylabel)
                    temp2 = get_relative_atomic_property(element=atom2.GetSymbol(), propertyname=propertylabel)
                    res = res + (temp1 - aveweight) * (temp2 - aveweight)
                    index += 1
                else:
                    res = res + 0.0
        if sum(tempp) == 0 or index == 0:
            result = 0
        else:
            result = (res / index) / (sum(tempp) / Natom)
        return result

    @staticmethod
    def calculate_moran_auto_mass(mol: Chem.Mol) -> dict:
        """Calculate Moran autocorrelation with carbon-scaled atomic mass."""
        res = {}
        for i in range(8):
            res[f'MATSm{i + 1}'] = Moran._calculate_moran_autocorrelation(mol, lag=i + 1, propertylabel='m')
        return res

    @staticmethod
    def calculate_moran_auto_volume(mol: Chem.Mol) -> dict:
        """Calculate Moran autocorrelation with carbon-scaled atomic van der Waals volume."""
        res = {}
        for i in range(8):
            res[f'MATSv{i + 1}'] = Moran._calculate_moran_autocorrelation(mol, lag=i + 1, propertylabel='V')
        return res

    @staticmethod
    def calculate_moran_auto_electronegativity(mol: Chem.Mol) -> dict:
        """Calculate Moran autocorrelation with carbon-scaled atomic Sanderson electronegativity."""
        res = {}
        for i in range(8):
            res[f'MATSe{i + 1}'] = Moran._calculate_moran_autocorrelation(mol, lag=i + 1, propertylabel='En')
        return res

    @staticmethod
    def calculate_moran_auto_polarizability(mol: Chem.Mol) -> dict:
        """Calculate Moran autocorrelation with carbon-scaled atomic polarizability."""
        res = {}
        for i in range(8):
            res[f'MATSp{i + 1}'] = Moran._calculate_moran_autocorrelation(mol, lag=i + 1, propertylabel='alapha')
        return res

    @staticmethod
    def get_all(mol: Chem.Mol) -> dict:
        """Calcualate all (32) Moran autocorrelation descriptors."""
        res = {}
        res.update(Moran.calculate_moran_auto_mass(mol))
        res.update(Moran.calculate_moran_auto_volume(mol))
        res.update(Moran.calculate_moran_auto_electronegativity(mol))
        res.update(Moran.calculate_moran_auto_polarizability(mol))
        return res
