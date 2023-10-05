# -*- coding: utf-8 -*-


"""Molecular physical/chemical properties."""


import math
import os
import subprocess  # noqa: S404
import tempfile
from platform import architecture
from sys import platform
from typing import List, Union

from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit.Chem import MolSurf

from .geo_opt import clean_mopac_files


class MolecularProperties:
    """Molecular physical/chemical properties."""

    @staticmethod
    def calculate_mollogp(mol: Chem.Mol) -> float:
        """Cacluate of MolLogP.

        From Wildman and G. M. Crippen JCICS _39_ 868-873 (1999).
        """
        return Crippen._pyMolLogP(mol)

    @staticmethod
    def calculate_mollogp2(mol: Chem.Mol) -> float:
        """Cacluate MolLogP^2.

        From Wildman and G. M. Crippen JCICS _39_ 868-873 (1999).
        """
        res = Crippen._pyMolLogP(mol)
        return res * res

    @staticmethod
    def calculate_molmr(mol: Chem.Mol) -> float:
        """Cacluate molecular refraction.

        From Wildman and G. M. Crippen JCICS _39_ 868-873 (1999).
        """
        return Crippen._pyMolMR(mol)

    @staticmethod
    def calculate_tpsa(mol: Chem.Mol) -> float:
        """Calculate the topological polar surface area.

        From Ertl P. et al., J.Med.Chem. (2000), 43,3714-3717.
        """
        return MolSurf.TPSA(mol)

    @staticmethod
    def _calculate_bond_number(mol: Chem.Mol, bondtype: str = 'SINGLE') -> float:
        """Calculate number of bond of specified type.

        :param bondtype: can be SINGLE, DOUBLE, TRIPLE or AROMATIC.
        """
        i = 0
        for bond in mol.GetBonds():
            if bond.GetBondType().name == bondtype:
                i += 1
        return i

    @staticmethod
    def calculate_unsaturation_index(mol: Chem.Mol) -> float:
        """Calculate unsaturation index."""
        nd = MolecularProperties._calculate_bond_number(mol, bondtype='DOUBLE')
        nt = MolecularProperties._calculate_bond_number(mol, bondtype='TRIPLE')
        na = MolecularProperties._calculate_bond_number(mol, bondtype='AROMATIC')
        res = math.log((1 + nd + nt + na), 2)
        return res

    @staticmethod
    def calculate_hydrophilicity_factor(mol: Chem.Mol) -> float:
        """Calculate hydrophilicity factor.

        From Todeschini R. et al., SAR QSAR Environ Res (1997), 7,173-193.
        """
        nheavy = mol.GetNumHeavyAtoms()
        nc = 0
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:
                nc += 1
        nhy = 0
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 7 or atom.GetAtomicNum() == 8 or atom.GetAtomicNum() == 16:
                atomn = atom.GetNeighbors()
                for i in atomn:
                    if i.GetAtomicNum() == 1:
                        nhy += 1
        res = ((1 + nhy) * math.log((1 + nhy), 2) + nc * (1.0 / nheavy * math.log(1.0 / nheavy, 2)) + math.sqrt(
            (nhy + 0.0) / (nheavy * nheavy))) / math.log(1.0 + nheavy)
        return res

    @staticmethod
    def get_all(mol: Chem.Mol) -> dict:
        """Get all (6) constitutional descriptors."""
        result = {}
        for des_label in molecular_property.keys():
            result[des_label] = molecular_property[des_label](mol)
        return result


molecular_property = {'MR': MolecularProperties.calculate_molmr,
                      'LogP': MolecularProperties.calculate_mollogp,
                      'LogP2': MolecularProperties.calculate_mollogp2,
                      'TPSA': MolecularProperties.calculate_tpsa,
                      'UI': MolecularProperties.calculate_unsaturation_index,
                      'Hy': MolecularProperties.calculate_hydrophilicity_factor,
                      }

