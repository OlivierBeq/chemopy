# -*- coding: utf-8 -*-


"""Estate fingerprints and values.

Based on Hall, Money, Kier, J. Chem. Inf. Comput. Sci. (1991)
doi: 10.1021/ci00001a012
"""


import numpy
from rdkit import Chem
from rdkit.Chem.EState import Fingerprinter as ESFP

import chemopy.AtomTypes as ATEstate


def _CalculateEState(mol: Chem.Mol, skipH: bool = True) -> float:
    """Get the EState value of each atom in the molecule."""
    mol = Chem.AddHs(mol)
    if skipH:
        mol = Chem.RemoveHs(mol)
    tb1 = Chem.GetPeriodicTable()
    nAtoms = mol.GetNumAtoms()
    Is = numpy.zeros(nAtoms, numpy.float)
    for i in range(nAtoms):
        at = mol.GetAtomWithIdx(i)
        atNum = at.GetAtomicNum()
        d = at.GetDegree()
        if d > 0:
            h = at.GetTotalNumHs()
            dv = tb1.GetNOuterElecs(atNum) - h
            # dv=numpy.array(_AtomHKDeltas(at),'d')
            N = _GetPrincipleQuantumNumber(atNum)
            Is[i] = (4.0 / (N * N) * dv + 1) / d
    dists = Chem.GetDistanceMatrix(mol, useBO=0, useAtomWts=0)
    dists += 1
    accum = numpy.zeros(nAtoms, numpy.float)
    for i in range(nAtoms):
        for j in range(i + 1, nAtoms):
            p = dists[i, j]
            if p < 1e6:
                temp = (Is[i] - Is[j]) / (p * p)
                accum[i] += temp
                accum[j] -= temp
    res = accum + Is
    return res


def _GetPrincipleQuantumNumber(atNum: int) -> int:
    """Get the principle quantum number of atom from atomic number."""
    if atNum <= 2:
        return 1
    elif atNum <= 10:
        return 2
    elif atNum <= 18:
        return 3
    elif atNum <= 36:
        return 4
    elif atNum <= 54:
        return 5
    elif atNum <= 86:
        return 6
    else:
        return 7


def CalculateHeavyAtomEState(mol: Chem.Mol) -> float:
    """Calculate sum of the EState indices over all heavy atoms."""
    return round(sum(_CalculateEState(mol)), 3)


def _CalculateAtomEState(mol: Chem.Mol, AtomicNum=6) -> float:
    """Calculate the sum of the EState indices over all atoms with specified atomic number."""
    nAtoms = mol.GetNumAtoms()
    Is = numpy.zeros(nAtoms, numpy.float)
    Estate = _CalculateEState(mol)
    for i in range(nAtoms):
        at = mol.GetAtomWithIdx(i)
        atNum = at.GetAtomicNum()
        if atNum == AtomicNum:
            Is[i] = Estate[i]
    res = sum(Is)
    return res


def CalculateCAtomEState(mol: Chem.Mol) -> float:
    """Calculate the sum of the EState indices over all C atoms."""
    return _CalculateAtomEState(mol, AtomicNum=6)


def CalculateHalogenEState(mol: Chem.Mol) -> float:
    """Calculate the sum of the EState indices over all Halogen atoms."""
    Nf = _CalculateAtomEState(mol, AtomicNum=9)
    Ncl = _CalculateAtomEState(mol, AtomicNum=17)
    Nbr = _CalculateAtomEState(mol, AtomicNum=35)
    Ni = _CalculateAtomEState(mol, AtomicNum=53)
    return round(Nf + Ncl + Nbr + Ni, 3)


def CalculateHeteroEState(mol: Chem.Mol) -> float:
    """Calculate the sum of the EState indices over all hetero atoms."""
    Ntotal = sum(_CalculateEState(mol))
    NC = _CalculateAtomEState(mol, AtomicNum=6)
    NH = _CalculateAtomEState(mol, AtomicNum=1)
    return round(Ntotal - NC - NH, 3)


def CalculateAverageEState(mol: Chem.Mol) -> float:
    """Calculate the ratio of the sum of the EState indices over heavy atoms and the number of non-hydrogen atoms."""
    N = mol.GetNumAtoms()
    return round(sum(_CalculateEState(mol)) / N, 3)


def CalculateMaxEState(mol: Chem.Mol) -> float:
    """Calculate the maximal Estate value in all atoms."""
    return round(max(_CalculateEState(mol)), 3)


def CalculateMinEState(mol: Chem.Mol) -> float:
    """Calculate the minimal Estate value in all atoms."""
    return round(min(_CalculateEState(mol)), 3)


def CalculateDiffMaxMinEState(mol: Chem.Mol) -> float:
    """Calculate the difference between Smax and Smin."""
    return round(max(_CalculateEState(mol)) - min(_CalculateEState(mol)), 3)


def CalculateEstateFingerprint(mol: Chem.Mol, implementation='rdkit') -> dict:
    """Calculate the sum of EState values for each EState atom type.

    :param implementation: either rdkit or chemopy. chemopy rounds
                           to the third decimal place but not rdkit.
    """
    if implementation not in ['rdkit', 'chemopy']:
        raise ValueError('Implementation of AtomTypeEState must be either rdkit or chemopy.')
    if implementation == 'chemopy':
        AT = ATEstate.GetAtomLabel(mol)
        Estate = _CalculateEState(mol)
        res = []
        for i in AT:
            if i == []:
                res.append(0)
            else:
                res.append(sum(Estate[k] for k in i))
        ESresult = {}
        for n, es in enumerate(res):
            ESresult[f'S{n+1}'] = round(es, 3)
        return ESresult
    else:  # RDKit with more decimals than chemopy
        temp = ESFP.FingerprintMol(mol)
        res = {}
        for i, j in enumerate(temp[1]):
            res[f'S{i + 1}'] = round(j, 3)
        return res


def CalculateAtomTypeEstateFingerprint(mol: Chem.Mol) -> dict:
    """Calculate EState Fingerprints.

    This is the counts of each EState atom type in the molecule.
    """
    temp = ESFP.FingerprintMol(mol)
    res = {}
    for i, j in enumerate(temp[0]):
        res[f'Sfinger{i + 1}'] = j
    return res


def CalculateMaxAtomTypeEState(mol: Chem.Mol) -> dict:
    """Calculate the maximum of EState value."""
    AT = ATEstate.GetAtomLabel(mol)
    Estate = _CalculateEState(mol)
    res = []
    for i in AT:
        if i == []:
            res.append(0)
        else:
            res.append(max(Estate[k] for k in i))
    ESresult = {}
    for n, es in enumerate(res):
        ESresult[f'Smax{n + 1}'] = round(es, 3)
    return ESresult


def CalculateMinAtomTypeEState(mol: Chem.Mol) -> dict:
    """Calculate the minimum of EState value."""
    AT = ATEstate.GetAtomLabel(mol)
    Estate = _CalculateEState(mol)
    res = []
    for i in AT:
        if i == []:
            res.append(0)
        else:
            res.append(min(Estate[k] for k in i))
    ESresult = {}
    for n, es in enumerate(res):
        ESresult[f'Smin{n + 1}'] = round(es, 3)
    return ESresult


def GetEstateDescriptors(mol: Chem.Mol) -> dict:
    """Calculate all (8) EState descriptors."""
    result = {}
    result.update({'Shev': CalculateHeavyAtomEState(mol)})
    result.update({'Scar': CalculateCAtomEState(mol)})
    result.update({'Shal': CalculateHalogenEState(mol)})
    result.update({'Shet': CalculateHeteroEState(mol)})
    result.update({'Save': CalculateAverageEState(mol)})
    result.update({'Smax': CalculateMaxEState(mol)})
    result.update({'Smin': CalculateMinEState(mol)})
    result.update({'DS': CalculateDiffMaxMinEState(mol)})
    return result


def GetEstateFingerprints(mol: Chem.Mol) -> dict:
    """Calculate all (316) EState descriptors."""
    result = {}
    result.update(CalculateEstateFingerprint(mol))
    result.update(CalculateAtomTypeEstateFingerprint(mol))
    result.update(CalculateMaxAtomTypeEState(mol))
    result.update(CalculateMinAtomTypeEState(mol))
    return result


################################################################
# if __name__=='__main__':
#     smi5=['COCCCC','CCC(C)CC','CC(C)CCC','CC(C)C(C)C','CCOCCN','c1ccccc1N']
#     smis = ['CCCC','CCCCC','CCCCCC','CC(N)C(=O)O','CC(N)C(=O)[O-].[Na+]']
#     for index, smi in enumerate(smis):
#         m = Chem.MolFromSmiles(smi)
#         print(index+1)
#         print(smi )
# ##        print('\t',CalculateEstateFingerprint(m))
# ##        print('\t',CalculateEstateValue(m))
# ##        print('\t',CalculateMaxAtomTypeEState(m))
# ##        print('\t', CalculateMinAtomTypeEState(m))
#         print(GetEstate(m))
#         print(len(GetEstate(m)))
