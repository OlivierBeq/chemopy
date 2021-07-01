# -*- coding: utf-8 -*-


"""Moreau-Broto autocorrelation descriptors."""


import numpy
from rdkit import Chem

from chemopy.AtomProperty import GetRelativeAtomicProperty


def _CalculateMoreauBrotoAutocorrelation(mol: Chem.Mol, lag: int = 1, propertylabel: str = 'm') -> float:
    """Calculate weighted Moreau-Broto autocorrelation descriptors.

    :param lag: topological distance between atom i and atom j.
    :param propertylabel: type of weighted property
    """
    Natom = mol.GetNumAtoms()
    GetDistanceMatrix = Chem.GetDistanceMatrix(mol)
    res = 0.0
    for i in range(Natom):
        for j in range(Natom):
            if GetDistanceMatrix[i, j] == lag:
                atom1 = mol.GetAtomWithIdx(i)
                atom2 = mol.GetAtomWithIdx(j)
                temp1 = GetRelativeAtomicProperty(element=atom1.GetSymbol(), propertyname=propertylabel)
                temp2 = GetRelativeAtomicProperty(element=atom2.GetSymbol(), propertyname=propertylabel)
                res = res + temp1 * temp2
            else:
                res = res + 0.0
    return round(numpy.log(res / 2 + 1), 3)


def CalculateMoreauBrotoAutoMass(mol: Chem.Mol) -> dict:
    """Calculate Moreau-Broto autocorrelation with carbon-scaled atomic mass."""
    res = {}
    for i in range(8):
        res[f'ATSm{i}'] = _CalculateMoreauBrotoAutocorrelation(mol, lag=i + 1, propertylabel='m')
    return res


def CalculateMoreauBrotoAutoVolume(mol: Chem.Mol) -> dict:
    """Calculate Moreau-Broto autocorrelation with carbon-scaled atomic van der Waals volume."""
    res = {}
    for i in range(8):
        res[f'ATSv{i}'] = _CalculateMoreauBrotoAutocorrelation(mol, lag=i + 1, propertylabel='V')
    return res


def CalculateMoreauBrotoAutoElectronegativity(mol: Chem.Mol) -> dict:
    """Calculate Moreau-Broto autocorrelation with carbon-scaled atomic Sanderson electronegativity."""
    res = {}
    for i in range(8):
        res[f'ATSe{i}'] = _CalculateMoreauBrotoAutocorrelation(mol, lag=i + 1, propertylabel='En')
    return res


def CalculateMoreauBrotoAutoPolarizability(mol: Chem.Mol) -> dict:
    """Calculate Moreau-Broto autocorrelation with carbon-scaled atomic polarizability."""
    res = {}
    for i in range(8):
        res[f'ATSp{i}'] = _CalculateMoreauBrotoAutocorrelation(mol, lag=i + 1, propertylabel='alapha')
    return res


def GetMoreauBrotoAuto(mol: Chem.Mol) -> dict:
    """Calculate all (32) Moreau-Broto autocorrelation descriptors."""
    res = {}
    res.update(CalculateMoreauBrotoAutoMass(mol))
    res.update(CalculateMoreauBrotoAutoVolume(mol))
    res.update(CalculateMoreauBrotoAutoElectronegativity(mol))
    res.update(CalculateMoreauBrotoAutoPolarizability(mol))
    return res


###########################################################################
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
#         print(len(GetMoreauBrotoAuto(m)))
