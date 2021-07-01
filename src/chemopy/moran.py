# -*- coding: utf-8 -*-


"""Moran autocorrelation descriptors."""

import numpy
from rdkit import Chem

from chemopy.AtomProperty import GetRelativeAtomicProperty


def _CalculateMoranAutocorrelation(mol: Chem.Mol, lag: int = 1, propertylabel: str = 'm') -> float:
    """Calculate weighted Moran autocorrelation descriptors.

    :param lag: topological distance between atom i and atom j.
    :param propertylabel: type of weighted property
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
                res = res + (temp1 - aveweight) * (temp2 - aveweight)
                index += 1
            else:
                res = res + 0.0
    if sum(tempp) == 0 or index == 0:
        result = 0
    else:
        result = (res / index) / (sum(tempp) / Natom)
    return round(result, 3)


def CalculateMoranAutoMass(mol: Chem.Mol) -> dict:
    """Calculate Moran autocorrelation with carbon-scaled atomic mass."""
    res = {}
    for i in range(8):
        res[f'MATSm{i}'] = _CalculateMoranAutocorrelation(mol, lag=i + 1, propertylabel='m')
    return res


def CalculateMoranAutoVolume(mol: Chem.Mol) -> dict:
    """Calculate Moran autocorrelation with carbon-scaled atomic van der Waals volume."""
    res = {}
    for i in range(8):
        res[f'MATSv{i}'] = _CalculateMoranAutocorrelation(mol, lag=i + 1, propertylabel='V')
    return res


def CalculateMoranAutoElectronegativity(mol: Chem.Mol) -> dict:
    """Calculate Moran autocorrelation with carbon-scaled atomic Sanderson electronegativity."""
    res = {}
    for i in range(8):
        res[f'MATSe{i}'] = _CalculateMoranAutocorrelation(mol, lag=i + 1, propertylabel='En')
    return res


def CalculateMoranAutoPolarizability(mol: Chem.Mol) -> dict:
    """Calculate Moran autocorrelation with carbon-scaled atomic polarizability."""
    res = {}
    for i in range(8):
        res[f'MATSp{i}'] = _CalculateMoranAutocorrelation(mol, lag=i + 1, propertylabel='alapha')
    return res


def GetMoranAuto(mol: Chem.Mol) -> dict:
    """Calcualate all (32) Moran autocorrelation descriptors."""
    res = {}
    res.update(CalculateMoranAutoMass(mol))
    res.update(CalculateMoranAutoVolume(mol))
    res.update(CalculateMoranAutoElectronegativity(mol))
    res.update(CalculateMoranAutoPolarizability(mol))
    return res


# ###########################################################################
# if __name__=='__main__':
#     smi5=['COCCCC','CCC(C)CC','CC(C)CCC','CC(C)C(C)C','CCOCCN','c1ccccc1N']
#     smis = ['CCCC','CCCCC','CCCCCC','CC(N)C(=O)O','CC(N)C(=O)[O-].[Na+]']
#     for index, smi in enumerate(smis):
#         m = Chem.MolFromSmiles(smi)
#         print(index+1)
#         print(smi)
# ##        print('\t',CalculateEstateFingerprint(m))
# ##        print('\t',CalculateEstateValue(m))
# ##        print('\t',CalculateMaxAtomTypeEState(m))
# ##        print('\t', CalculateMinAtomTypeEState(m))
#         print(GetMoranAuto(m))
