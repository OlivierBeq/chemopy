# -*- coding: utf-8 -*-


"""Geary autocorrelation indices."""

import numpy as np
from rdkit import Chem

from .atom_property import get_relative_atomic_property


class Geary:
    """Geary autocorrelation indices."""

    @staticmethod
    def _calculate_geary_autocorrelation(mol: Chem.Mol, lag: int = 1, propertylabel: str = "m") -> float:
        """Calculate weighted Geary autocorrelation descriptors.

        :param lag: topological distance between atom i and atom j.
        :param propertylabel: weighted property.
        """
        Natom = mol.GetNumAtoms()
        prop = np.array(
            [get_relative_atomic_property(atom.GetSymbol(), propertyname=propertylabel) for atom in mol.GetAtoms()]
        )
        aveweight = prop.sum() / Natom
        tempp_sum = np.square(prop - aveweight).sum()
        distance_matrix = Chem.GetDistanceMatrix(mol)
        mask = distance_matrix == lag
        index = int(mask.sum())
        res = np.square(prop[:, None] - prop[None, :])[mask].sum()
        if tempp_sum == 0 or index == 0:
            result = 0
        else:
            result = (res / index / 2) / (tempp_sum / (Natom - 1))
        return result

    @staticmethod
    def calculate_geary_auto_mass(mol: Chem.Mol) -> dict:
        """Calculate Geary autocorrelation descriptors from carbon-scaled atomic mass."""
        res = {}
        for i in range(8):
            res[f"GATSm{i + 1}"] = Geary._calculate_geary_autocorrelation(mol, lag=i + 1, propertylabel="m")
        return res

    @staticmethod
    def calculate_geary_auto_volume(mol: Chem.Mol) -> dict:
        """Calculate Geary autocorrelation descriptors from carbon-scaled atomic van der Waals volume."""
        res = {}
        for i in range(8):
            res[f"GATSv{i + 1}"] = Geary._calculate_geary_autocorrelation(mol, lag=i + 1, propertylabel="V")
        return res

    @staticmethod
    def calculate_geary_auto_electronegativity(mol: Chem.Mol) -> dict:
        """Calculate Geary autocorrelation descriptors from carbon-scaled atomic Sanderson electronegativity."""
        res = {}
        for i in range(8):
            res[f"GATSe{i + 1}"] = Geary._calculate_geary_autocorrelation(mol, lag=i + 1, propertylabel="En")
        return res

    @staticmethod
    def calculate_geary_auto_polarizability(mol: Chem.Mol) -> dict:
        """Calculate Geary autocorrelation descriptors from carbon-scaled atomic polarizability."""
        res = {}
        for i in range(8):
            res[f"GATSp{i + 1}"] = Geary._calculate_geary_autocorrelation(mol, lag=i + 1, propertylabel="alapha")
        return res

    @staticmethod
    def get_all(mol: Chem.Mol) -> dict:
        """Calcualate all (32) Geary autocorrelation descriptors.

        Carbon-scaled weigthed schemes: atomic mass, atomic van der Waals volume,
                                        Sanderson electronegativity, atomic polarizability
        """
        res = {}
        res.update(Geary.calculate_geary_auto_mass(mol))
        res.update(Geary.calculate_geary_auto_volume(mol))
        res.update(Geary.calculate_geary_auto_electronegativity(mol))
        res.update(Geary.calculate_geary_auto_polarizability(mol))
        return res
