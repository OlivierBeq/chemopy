# -*- coding: utf-8 -*-


"""Burden eigvenvalue descriptors."""

import numpy
import numpy.linalg
from rdkit import Chem

from chemopy.AtomProperty import GetRelativeAtomicProperty


def _GetBurdenMatrix(mol: Chem.Mol, propertylabel: str = 'm') -> numpy.matrix:
    """Calculate weighted Burden matrix and eigenvalues."""
    mol = Chem.AddHs(mol)
    Natom = mol.GetNumAtoms()
    AdMatrix = Chem.GetAdjacencyMatrix(mol)
    bondindex = numpy.argwhere(AdMatrix)
    AdMatrix1 = numpy.array(AdMatrix, dtype=numpy.float32)
    # The diagonal elements of B, Bii, are either given by
    # the carbon normalized atomic mass,
    # van der Waals volume, Sanderson electronegativity,
    # and polarizability of atom i.
    for i in range(Natom):
        atom = mol.GetAtomWithIdx(i)
        temp = GetRelativeAtomicProperty(element=atom.GetSymbol(), propertyname=propertylabel)
        AdMatrix1[i, i] = round(temp, 3)
    # The element of B connecting atoms i and j, Bij,
    # is equal to the square root of the bond
    # order between atoms i and j.
    for i in bondindex:
        bond = mol.GetBondBetweenAtoms(int(i[0]), int(i[1]))
        if bond.GetBondType().name == 'SINGLE':
            AdMatrix1[i[0], i[1]] = round(numpy.sqrt(1), 3)
        if bond.GetBondType().name == "DOUBLE":
            AdMatrix1[i[0], i[1]] = round(numpy.sqrt(2), 3)
        if bond.GetBondType().name == "TRIPLE":
            AdMatrix1[i[0], i[1]] = round(numpy.sqrt(3), 3)
        if bond.GetBondType().name == "AROMATIC":
            AdMatrix1[i[0], i[1]] = round(numpy.sqrt(1.5), 3)
    # All other elements of B (corresponding non bonded
    # atom pairs) are set to 0.001
    bondnonindex = numpy.argwhere(AdMatrix == 0)
    for i in bondnonindex:
        if i[0] != i[1]:
            AdMatrix1[i[0], i[1]] = 0.001
    return numpy.real(numpy.linalg.eigvals(AdMatrix1))


def CalculateBurdenMass(mol: Chem.Mol) -> dict:
    """Calculate (16) Burden descriptors from atomic mass."""
    temp = _GetBurdenMatrix(mol, propertylabel='m')
    temp1 = numpy.sort(temp[temp >= 0])
    temp2 = numpy.sort(numpy.abs(temp[temp < 0]))
    if len(temp1) < 8:
        temp1 = numpy.concatenate((numpy.zeros(8), temp1))
    if len(temp2) < 8:
        temp2 = numpy.concatenate((numpy.zeros(8), temp2))
    bcut = ["bcutm16", "bcutm15", "bcutm14", "bcutm13", "bcutm12", "bcutm11", "bcutm10",
            "bcutm9", "bcutm8", "bcutm7", "bcutm6", "bcutm5", "bcutm4", "bcutm3",
            "bcutm2", "bcutm1"]
    bcutvalue = numpy.concatenate((temp2[-8:], temp1[-8:]))
    bcutvalue = [round(i, 3) for i in bcutvalue]
    res = dict(zip(bcut, bcutvalue))
    return res


def CalculateBurdenVDW(mol: Chem.Mol) -> dict:
    """Calculate (16) Burden descriptors from atomic volumes."""
    temp = _GetBurdenMatrix(mol, propertylabel='V')
    temp1 = numpy.sort(temp[temp >= 0])
    temp2 = numpy.sort(numpy.abs(temp[temp < 0]))
    if len(temp1) < 8:
        temp1 = numpy.concatenate((numpy.zeros(8), temp1))
    if len(temp2) < 8:
        temp2 = numpy.concatenate((numpy.zeros(8), temp2))
    bcut = ["bcutv16", "bcutv15", "bcutv14", "bcutv13", "bcutv12", "bcutv11", "bcutv10",
            "bcutv9", "bcutv8", "bcutv7", "bcutv6", "bcutv5", "bcutv4", "bcutv3",
            "bcutv2", "bcutv1"]
    bcutvalue = numpy.concatenate((temp2[-8:], temp1[-8:]))
    bcutvalue = [round(i, 3) for i in bcutvalue]
    res = dict(zip(bcut, bcutvalue))
    return res


def CalculateBurdenElectronegativity(mol: Chem.Mol) -> dict:
    """Calculate (16) Burden descriptors from atomic electronegativity."""
    temp = _GetBurdenMatrix(mol, propertylabel='En')
    temp1 = numpy.sort(temp[temp >= 0])
    temp2 = numpy.sort(numpy.abs(temp[temp < 0]))
    if len(temp1) < 8:
        temp1 = numpy.concatenate((numpy.zeros(8), temp1))
    if len(temp2) < 8:
        temp2 = numpy.concatenate((numpy.zeros(8), temp2))
    bcut = ["bcute16", "bcute15", "bcute14", "bcute13", "bcute12", "bcute11", "bcute10",
            "bcute9", "bcute8", "bcute7", "bcute6", "bcute5", "bcute4", "bcute3",
            "bcute2", "bcute1"]
    bcutvalue = numpy.concatenate((temp2[-8:], temp1[-8:]))
    bcutvalue = [round(i, 3) for i in bcutvalue]
    res = dict(zip(bcut, bcutvalue))
    return res


def CalculateBurdenPolarizability(mol: Chem.Mol) -> dict:
    """Calculate (16) Burden descriptors from atomic polarizability."""
    temp = _GetBurdenMatrix(mol, propertylabel='alapha')
    temp1 = numpy.sort(temp[temp >= 0])
    temp2 = numpy.sort(numpy.abs(temp[temp < 0]))
    if len(temp1) < 8:
        temp1 = numpy.concatenate((numpy.zeros(8), temp1))
    if len(temp2) < 8:
        temp2 = numpy.concatenate((numpy.zeros(8), temp2))
    bcut = ["bcutp16", "bcutp15", "bcutp14", "bcutp13", "bcutp12", "bcutp11", "bcutp10",
            "bcutp9", "bcutp8", "bcutp7", "bcutp6", "bcutp5", "bcutp4", "bcutp3",
            "bcutp2", "bcutp1"]
    bcutvalue = numpy.concatenate((temp2[-8:], temp1[-8:]))
    bcutvalue = [round(i, 3) for i in bcutvalue]
    res = dict(zip(bcut, bcutvalue))
    return res


def GetBurden(mol: Chem.Mol) -> dict:
    """Calculate all (64) Burden descriptors."""
    bcut = {}
    bcut.update(CalculateBurdenMass(mol))
    bcut.update(CalculateBurdenVDW(mol))
    bcut.update(CalculateBurdenElectronegativity(mol))
    bcut.update(CalculateBurdenPolarizability(mol))
    return bcut


########################################################################
# if __name__ =='__main__':
#     smi5=['CCOCCC','CCC(C)CC','CC(C)CCC','CC(C)C(C)C','CCCCCN','c1ccccc1N','C']
#     for index, smi in enumerate(smi5):
#         m = Chem.MolFromSmiles(smi)
#         print(index+1)
#         print(smi, '\n')
#         print(GetBurden(m))
#         print(len(GetBurden(m)))
