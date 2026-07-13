# -*- coding: utf-8 -*-


"""Burden eigvenvalue descriptors."""

import numpy as np
import numpy.linalg
from rdkit import Chem
from natsort import natsorted

from .atom_property import get_relative_atomic_property


class BCUT:
    """BCUT descriptors."""

    _BOND_ORDER_SQRT = {"SINGLE": 1.0, "DOUBLE": np.sqrt(2), "TRIPLE": np.sqrt(3), "AROMATIC": np.sqrt(1.5)}

    @staticmethod
    def _get_burden_base_matrix(mol: Chem.Mol) -> np.array:
        """Build the propertylabel-independent part of the Burden matrix (off-diagonal + fill).

        Only the diagonal (`_get_burden_matrix`) depends on the chosen atomic property, so this
        (bond-order-derived off-diagonal, 0.001-filled elsewhere) is computed once per molecule
        and reused for each of the 4 property variants instead of being rebuilt identically 4 times.
        """
        mol = Chem.AddHs(mol)
        AdMatrix = Chem.GetAdjacencyMatrix(mol)
        base = np.full(AdMatrix.shape, 0.001, dtype=np.float32)
        np.fill_diagonal(base, 0.0)
        for i0, i1 in np.argwhere(AdMatrix):
            bond = mol.GetBondBetweenAtoms(int(i0), int(i1))
            base[i0, i1] = BCUT._BOND_ORDER_SQRT.get(bond.GetBondType().name, base[i0, i1])
        return base

    @staticmethod
    def _get_burden_matrix(mol: Chem.Mol, propertylabel: str = "m", base_matrix: np.array = None) -> np.array:
        """Calculate weighted Burden matrix and eigenvalues."""
        mol = Chem.AddHs(mol)
        if base_matrix is None:
            base_matrix = BCUT._get_burden_base_matrix(mol)
        AdMatrix1 = base_matrix.copy()
        # The diagonal elements of B, Bii, are either given by
        # the carbon normalized atomic mass,
        # van der Waals volume, Sanderson electronegativity,
        # and polarizability of atom i.
        diag = np.array(
            [
                get_relative_atomic_property(element=atom.GetSymbol(), propertyname=propertylabel)
                for atom in mol.GetAtoms()
            ]
        )
        np.fill_diagonal(AdMatrix1, diag)
        return np.real(np.linalg.eigvals(AdMatrix1))

    @staticmethod
    def calculate_burden_mass(mol: Chem.Mol, base_matrix: np.array = None) -> dict:
        """Calculate (16) Burden descriptors from atomic mass."""
        temp = BCUT._get_burden_matrix(mol, propertylabel="m", base_matrix=base_matrix)
        temp1 = np.sort(temp[temp >= 0])
        temp2 = np.sort(np.abs(temp[temp < 0]))
        if len(temp1) < 8:
            temp1 = np.concatenate((np.zeros(8), temp1))
        if len(temp2) < 8:
            temp2 = np.concatenate((np.zeros(8), temp2))
        bcut = [
            "bcutm16",
            "bcutm15",
            "bcutm14",
            "bcutm13",
            "bcutm12",
            "bcutm11",
            "bcutm10",
            "bcutm9",
            "bcutm8",
            "bcutm7",
            "bcutm6",
            "bcutm5",
            "bcutm4",
            "bcutm3",
            "bcutm2",
            "bcutm1",
        ]
        bcutvalue = np.concatenate((temp2[-8:], temp1[-8:]))
        bcutvalue = [i for i in bcutvalue]
        res = dict(natsorted(dict(zip(bcut, bcutvalue)).items()))
        return res

    @staticmethod
    def calculate_burden_vdw(mol: Chem.Mol, base_matrix: np.array = None) -> dict:
        """Calculate (16) Burden descriptors from atomic volumes."""
        temp = BCUT._get_burden_matrix(mol, propertylabel="V", base_matrix=base_matrix)
        temp1 = np.sort(temp[temp >= 0])
        temp2 = np.sort(np.abs(temp[temp < 0]))
        if len(temp1) < 8:
            temp1 = np.concatenate((np.zeros(8), temp1))
        if len(temp2) < 8:
            temp2 = np.concatenate((np.zeros(8), temp2))
        bcut = [
            "bcutv16",
            "bcutv15",
            "bcutv14",
            "bcutv13",
            "bcutv12",
            "bcutv11",
            "bcutv10",
            "bcutv9",
            "bcutv8",
            "bcutv7",
            "bcutv6",
            "bcutv5",
            "bcutv4",
            "bcutv3",
            "bcutv2",
            "bcutv1",
        ]
        bcutvalue = np.concatenate((temp2[-8:], temp1[-8:]))
        bcutvalue = [i for i in bcutvalue]
        res = dict(natsorted(dict(zip(bcut, bcutvalue)).items()))
        return res

    @staticmethod
    def calculate_burden_electronegativity(mol: Chem.Mol, base_matrix: np.array = None) -> dict:
        """Calculate (16) Burden descriptors from atomic electronegativity."""
        temp = BCUT._get_burden_matrix(mol, propertylabel="En", base_matrix=base_matrix)
        temp1 = np.sort(temp[temp >= 0])
        temp2 = np.sort(np.abs(temp[temp < 0]))
        if len(temp1) < 8:
            temp1 = np.concatenate((np.zeros(8), temp1))
        if len(temp2) < 8:
            temp2 = np.concatenate((np.zeros(8), temp2))
        bcut = [
            "bcute16",
            "bcute15",
            "bcute14",
            "bcute13",
            "bcute12",
            "bcute11",
            "bcute10",
            "bcute9",
            "bcute8",
            "bcute7",
            "bcute6",
            "bcute5",
            "bcute4",
            "bcute3",
            "bcute2",
            "bcute1",
        ]
        bcutvalue = np.concatenate((temp2[-8:], temp1[-8:]))
        bcutvalue = [i for i in bcutvalue]
        res = dict(natsorted(dict(zip(bcut, bcutvalue)).items()))
        return res

    @staticmethod
    def calculate_burden_polarizability(mol: Chem.Mol, base_matrix: np.array = None) -> dict:
        """Calculate (16) Burden descriptors from atomic polarizability."""
        temp = BCUT._get_burden_matrix(mol, propertylabel="alapha", base_matrix=base_matrix)
        temp1 = np.sort(temp[temp >= 0])
        temp2 = np.sort(np.abs(temp[temp < 0]))
        if len(temp1) < 8:
            temp1 = np.concatenate((np.zeros(8), temp1))
        if len(temp2) < 8:
            temp2 = np.concatenate((np.zeros(8), temp2))
        bcut = [
            "bcutp16",
            "bcutp15",
            "bcutp14",
            "bcutp13",
            "bcutp12",
            "bcutp11",
            "bcutp10",
            "bcutp9",
            "bcutp8",
            "bcutp7",
            "bcutp6",
            "bcutp5",
            "bcutp4",
            "bcutp3",
            "bcutp2",
            "bcutp1",
        ]
        bcutvalue = np.concatenate((temp2[-8:], temp1[-8:]))
        bcutvalue = [i for i in bcutvalue]
        res = dict(natsorted(dict(zip(bcut, bcutvalue)).items()))
        return res

    @staticmethod
    def get_all(mol: Chem.Mol) -> dict:
        """Calculate all (64) Burden descriptors."""
        base_matrix = BCUT._get_burden_base_matrix(mol)
        bcut = {}
        bcut.update(BCUT.calculate_burden_mass(mol, base_matrix=base_matrix))
        bcut.update(BCUT.calculate_burden_vdw(mol, base_matrix=base_matrix))
        bcut.update(BCUT.calculate_burden_electronegativity(mol, base_matrix=base_matrix))
        bcut.update(BCUT.calculate_burden_polarizability(mol, base_matrix=base_matrix))
        return bcut
