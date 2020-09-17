# -*- coding: utf-8 -*-


"""Geary autocorrelation indices."""

import numpy
from rdkit import Chem

from pychem.AtomProperty import GetRelativeAtomicProperty


def _CalculateGearyAutocorrelation(mol: Chem.Mol, lag: int = 1, propertylabel: str = 'm') -> float:
    """Calculate weighted Geary autocorrelation descriptors.

    :param lag: topological distance between atom i and atom j.
    :param propertylabel: weighted property.
    """
    Natom = mol.GetNumAtoms()
    prolist = []
    for i in mol.GetAtoms():
        temp = GetRelativeAtomicProperty(i.GetSymbol(), propertyname=propertylabel)
        prolist.append(temp)
    aveweight = sum(prolist) / Natom
    tempp = [numpy.square(x - aveweight) for x in prolist]
    GetDistanceMatrix = Chem.GetDistanceMatrix(mol)
    res = 0.0
    index = 0
    for i in range(Natom):
        for j in range(Natom):
            if GetDistanceMatrix[i, j] == lag:
                atom1 = mol.GetAtomWithIdx(i)
                atom2 = mol.GetAtomWithIdx(j)
                temp1 = GetRelativeAtomicProperty(element=atom1.GetSymbol(), propertyname=propertylabel)
                temp2 = GetRelativeAtomicProperty(element=atom2.GetSymbol(), propertyname=propertylabel)
                res += numpy.square(temp1 - temp2)
                index += 1
            else:
                res += 0.0
    if sum(tempp) == 0 or index == 0:
        result = 0
    else:
        result = (res / index / 2) / (sum(tempp) / (Natom - 1))
    return round(result, 3)


def CalculateGearyAutoMass(mol: Chem.Mol) -> dict:
    """Calculate Geary autocorrelation descriptors from carbon-scaled atomic mass."""
    res = {}
    for i in range(8):
        res[f'GATSm{i + 1}'] = _CalculateGearyAutocorrelation(mol, lag=i + 1, propertylabel='m')
    return res


def CalculateGearyAutoVolume(mol: Chem.Mol) -> dict:
    """Calculate Geary autocorrelation descriptors from carbon-scaled atomic van der Waals volume."""
    res = {}
    for i in range(8):
        res[f'GATSv{i + 1}'] = _CalculateGearyAutocorrelation(mol, lag=i + 1, propertylabel='V')
    return res


def CalculateGearyAutoElectronegativity(mol: Chem.Mol) -> dict:
    """Calculate Geary autocorrelation descriptors from carbon-scaled atomic Sanderson electronegativity."""
    res = {}
    for i in range(8):
        res[f'GATSe{i + 1}'] = _CalculateGearyAutocorrelation(mol, lag=i + 1, propertylabel='En')
    return res


def CalculateGearyAutoPolarizability(mol: Chem.Mol) -> dict:
    """Calculate Geary autocorrelation descriptors from carbon-scaled atomic polarizability."""
    res = {}
    for i in range(8):
        res[f'GATSp{i + 1}'] = _CalculateGearyAutocorrelation(mol, lag=i + 1, propertylabel='alapha')
    return res


def GetGearyAuto(mol: Chem.Mol) -> dict:
    """Calcualate all (32) Geary autocorrelation descriptors.

    Carbon-scaled weigthed schemes: atomic mass, atomic van der Waals volume,
                                    Sanderson electronegativity, atomic polarizability
    """
    res = {}
    res.update(CalculateGearyAutoMass(mol))
    res.update(CalculateGearyAutoVolume(mol))
    res.update(CalculateGearyAutoElectronegativity(mol))
    res.update(CalculateGearyAutoPolarizability(mol))
    return res


###########################################################################
# if __name__=='__main__':
#     smi5=['COCCCC','CCC(C)CC','CC(C)CCC','CC(C)C(C)C','CCOCCN','c1ccccc1N']
#     smis = ['CCCC','CCCCC','CCCCCC','CC(N)C(=O)O','CC(N)C(=O)[O-].[Na+]']
#     for index, smi in enumerate(smi5):
#         m = Chem.MolFromSmiles(smi)
#         print(index+1)
#         print(smi)
# ##        print('\t',CalculateEstateFingerprint(m))
# ##        print('\t',CalculateEstateValue(m))
# ##        print('\t',CalculateMaxAtomTypeEState(m))
# ##        print('\t', CalculateMinAtomTypeEState(m))
#         print(GetGearyAuto(m))
