# -*- coding: utf-8 -*-


"""Burden eigvenvalue descriptors."""

import numpy as np
import numpy.linalg
from rdkit import Chem
from natsort import natsorted

from .atomproperty import get_relative_atomic_property


class BCUT:
    """BCUT descriptors."""
    
    @staticmethod
    def _get_burden_matrix(mol: Chem.Mol, propertylabel: str = 'm') -> np.array:
        """Calculate weighted Burden matrix and eigenvalues."""
        mol = Chem.AddHs(mol)
        Natom = mol.GetNumAtoms()
        AdMatrix = Chem.GetAdjacencyMatrix(mol)
        bondindex = np.argwhere(AdMatrix)
        AdMatrix1 = np.array(AdMatrix, dtype=np.float32)
        # The diagonal elements of B, Bii, are either given by
        # the carbon normalized atomic mass,
        # van der Waals volume, Sanderson electronegativity,
        # and polarizability of atom i.
        for i in range(Natom):
            atom = mol.GetAtomWithIdx(i)
            temp = get_relative_atomic_property(element=atom.GetSymbol(), propertyname=propertylabel)
            AdMatrix1[i, i] = temp
        # The element of B connecting atoms i and j, Bij,
        # is equal to the square root of the bond
        # order between atoms i and j.
        for i in bondindex:
            bond = mol.GetBondBetweenAtoms(int(i[0]), int(i[1]))
            if bond.GetBondType().name == 'SINGLE':
                AdMatrix1[i[0], i[1]] = np.sqrt(1)
            if bond.GetBondType().name == "DOUBLE":
                AdMatrix1[i[0], i[1]] = np.sqrt(2)
            if bond.GetBondType().name == "TRIPLE":
                AdMatrix1[i[0], i[1]] = np.sqrt(3)
            if bond.GetBondType().name == "AROMATIC":
                AdMatrix1[i[0], i[1]] = np.sqrt(1.5)
        # All other elements of B (corresponding non-bonded
        # atom pairs) are set to 0.001
        bondnonindex = np.argwhere(AdMatrix == 0)
        for i in bondnonindex:
            if i[0] != i[1]:
                AdMatrix1[i[0], i[1]] = 0.001
        return np.real(np.linalg.eigvals(AdMatrix1))

    @staticmethod
    def calculate_burden_mass(mol: Chem.Mol) -> dict:
        """Calculate (16) Burden descriptors from atomic mass."""
        temp =  BCUT._get_burden_matrix(mol, propertylabel='m')
        temp1 = np.sort(temp[temp >= 0])
        temp2 = np.sort(np.abs(temp[temp < 0]))
        if len(temp1) < 8:
            temp1 = np.concatenate((np.zeros(8), temp1))
        if len(temp2) < 8:
            temp2 = np.concatenate((np.zeros(8), temp2))
        bcut = ["bcutm16", "bcutm15", "bcutm14", "bcutm13", "bcutm12", "bcutm11", "bcutm10",
                "bcutm9", "bcutm8", "bcutm7", "bcutm6", "bcutm5", "bcutm4", "bcutm3",
                "bcutm2", "bcutm1"]
        bcutvalue = np.concatenate((temp2[-8:], temp1[-8:]))
        bcutvalue = [i for i in bcutvalue]
        res = dict(natsorted(dict(zip(bcut, bcutvalue)).items()))
        return res

    @staticmethod
    def calculate_burden_vdw(mol: Chem.Mol) -> dict:
        """Calculate (16) Burden descriptors from atomic volumes."""
        temp = BCUT._get_burden_matrix(mol, propertylabel='V')
        temp1 = np.sort(temp[temp >= 0])
        temp2 = np.sort(np.abs(temp[temp < 0]))
        if len(temp1) < 8:
            temp1 = np.concatenate((np.zeros(8), temp1))
        if len(temp2) < 8:
            temp2 = np.concatenate((np.zeros(8), temp2))
        bcut = ["bcutv16", "bcutv15", "bcutv14", "bcutv13", "bcutv12", "bcutv11", "bcutv10",
                "bcutv9", "bcutv8", "bcutv7", "bcutv6", "bcutv5", "bcutv4", "bcutv3",
                "bcutv2", "bcutv1"]
        bcutvalue = np.concatenate((temp2[-8:], temp1[-8:]))
        bcutvalue = [i for i in bcutvalue]
        res = dict(natsorted(dict(zip(bcut, bcutvalue)).items()))
        return res

    @staticmethod
    def calculate_burden_electronegativity(mol: Chem.Mol) -> dict:
        """Calculate (16) Burden descriptors from atomic electronegativity."""
        temp = BCUT._get_burden_matrix(mol, propertylabel='En')
        temp1 = np.sort(temp[temp >= 0])
        temp2 = np.sort(np.abs(temp[temp < 0]))
        if len(temp1) < 8:
            temp1 = np.concatenate((np.zeros(8), temp1))
        if len(temp2) < 8:
            temp2 = np.concatenate((np.zeros(8), temp2))
        bcut = ["bcute16", "bcute15", "bcute14", "bcute13", "bcute12", "bcute11", "bcute10",
                "bcute9", "bcute8", "bcute7", "bcute6", "bcute5", "bcute4", "bcute3",
                "bcute2", "bcute1"]
        bcutvalue = np.concatenate((temp2[-8:], temp1[-8:]))
        bcutvalue = [i for i in bcutvalue]
        res = dict(natsorted(dict(zip(bcut, bcutvalue)).items()))
        return res

    @staticmethod
    def calculate_burden_polarizability(mol: Chem.Mol) -> dict:
        """Calculate (16) Burden descriptors from atomic polarizability."""
        temp = BCUT._get_burden_matrix(mol, propertylabel='alapha')
        temp1 = np.sort(temp[temp >= 0])
        temp2 = np.sort(np.abs(temp[temp < 0]))
        if len(temp1) < 8:
            temp1 = np.concatenate((np.zeros(8), temp1))
        if len(temp2) < 8:
            temp2 = np.concatenate((np.zeros(8), temp2))
        bcut = ["bcutp16", "bcutp15", "bcutp14", "bcutp13", "bcutp12", "bcutp11", "bcutp10",
                "bcutp9", "bcutp8", "bcutp7", "bcutp6", "bcutp5", "bcutp4", "bcutp3",
                "bcutp2", "bcutp1"]
        bcutvalue = np.concatenate((temp2[-8:], temp1[-8:]))
        bcutvalue = [i for i in bcutvalue]
        res = dict(natsorted(dict(zip(bcut, bcutvalue)).items()))
        return res

    @staticmethod
    def get_all(mol: Chem.Mol) -> dict:
        """Calculate all (64) Burden descriptors."""
        bcut = {}
        bcut.update(BCUT.calculate_burden_mass(mol))
        bcut.update(BCUT.calculate_burden_vdw(mol))
        bcut.update(BCUT.calculate_burden_electronegativity(mol))
        bcut.update(BCUT.calculate_burden_polarizability(mol))
        return bcut
