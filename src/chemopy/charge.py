# -*- coding: utf-8 -*-


"""Charge descriptors based on Gasteiger/Marseli partial charges."""

import numpy
from rdkit import Chem
from rdkit.Chem import rdPartialCharges as GMCharge

iter_step = 12


def _CalculateElementMaxPCharge(mol: Chem.Mol, AtomicNum: int = 6):
    """Get the most positive charge of atom with specified atomic number."""
    Hmol = Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
    res = []
    for atom in Hmol.GetAtoms():
        if atom.GetAtomicNum() == AtomicNum:
            res.append(float(atom.GetProp('_GasteigerCharge')))
    if res == []:
        return 0
    else:
        return round(max(res), 3)


def _CalculateElementMaxNCharge(mol: Chem.Mol, AtomicNum: int = 6) -> float:
    """Get the most negative charge of atom with specified atomic number."""
    Hmol = Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
    res = []
    for atom in Hmol.GetAtoms():
        if atom.GetAtomicNum() == AtomicNum:
            res.append(float(atom.GetProp('_GasteigerCharge')))
    if res == []:
        return 0
    else:
        return round(min(res), 3)


def CalculateHMaxPCharge(mol: Chem.Mol) -> float:
    """Get most positive charge of all hydrogen atoms."""
    return _CalculateElementMaxPCharge(mol, AtomicNum=1)


def CalculateCMaxPCharge(mol: Chem.Mol) -> float:
    """Get most positive charge of all carbon atoms."""
    return _CalculateElementMaxPCharge(mol, AtomicNum=6)


def CalculateNMaxPCharge(mol: Chem.Mol) -> float:
    """Get most positive charge of all nitrogen atoms."""
    return _CalculateElementMaxPCharge(mol, AtomicNum=7)


def CalculateOMaxPCharge(mol: Chem.Mol) -> float:
    """Get most positive charge of all oxygen atoms."""
    return _CalculateElementMaxPCharge(mol, AtomicNum=8)


def CalculateHMaxNCharge(mol) -> float:
    """Get most negative charge of all hydrogen atoms."""
    return _CalculateElementMaxNCharge(mol, AtomicNum=1)


def CalculateCMaxNCharge(mol: Chem.Mol) -> float:
    """Get most negative charge of all carbon atoms."""
    return _CalculateElementMaxNCharge(mol, AtomicNum=6)


def CalculateNMaxNCharge(mol: Chem.Mol) -> float:
    """Get most negative charge of all nitrogen atoms."""
    return _CalculateElementMaxNCharge(mol, AtomicNum=7)


def CalculateOMaxNCharge(mol: Chem.Mol) -> float:
    """Get most negative charge of all oxygen atoms."""
    return _CalculateElementMaxNCharge(mol, AtomicNum=8)


def CalculateAllMaxPCharge(mol: Chem.Mol) -> float:
    """Get most positive charge of all atoms."""
    Hmol = Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
    res = []
    for atom in Hmol.GetAtoms():
        res.append(float(atom.GetProp('_GasteigerCharge')))
    if res == []:
        return 0
    else:
        return round(max(res), 3)


def CalculateAllMaxNCharge(mol: Chem.Mol) -> float:
    """Get most negative charge of all atoms."""
    Hmol = Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
    res = []
    for atom in Hmol.GetAtoms():
        res.append(float(atom.GetProp('_GasteigerCharge')))
    if res == []:
        return 0
    else:
        return round(min(res), 3)


def _CalculateElementSumSquareCharge(mol: Chem.Mol, AtomicNum: int = 6) -> float:
    """Get the sum of square charges of atoms with specified atomic number."""
    Hmol = Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
    res = []
    for atom in Hmol.GetAtoms():
        if atom.GetAtomicNum() == AtomicNum:
            res.append(float(atom.GetProp('_GasteigerCharge')))
    if res == []:
        return 0
    else:
        return round(sum(numpy.square(res)), 3)


def CalculateHSumSquareCharge(mol: Chem.Mol) -> float:
    """Get the sum of square charges of hydrogen atoms."""
    return _CalculateElementSumSquareCharge(mol, AtomicNum=1)


def CalculateCSumSquareCharge(mol: Chem.Mol) -> float:
    """Get the sum of square charges of carbon atoms."""
    return _CalculateElementSumSquareCharge(mol, AtomicNum=6)


def CalculateNSumSquareCharge(mol: Chem.Mol) -> float:
    """Get the sum of square charges of nitrogen atoms."""
    return _CalculateElementSumSquareCharge(mol, AtomicNum=7)


def CalculateOSumSquareCharge(mol: Chem.Mol):
    """Get the sum of square charges of oxygen atoms."""
    return _CalculateElementSumSquareCharge(mol, AtomicNum=8)


def CalculateAllSumSquareCharge(mol: Chem.Mol):
    """Get the sum of square charges of all atoms."""
    Hmol = Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
    res = []
    for atom in Hmol.GetAtoms():
        res.append(float(atom.GetProp('_GasteigerCharge')))
    if res == []:
        return 0
    else:
        return round(sum(numpy.square(res)), 3)


def CalculateTotalPCharge(mol: Chem.Mol) -> float:
    """Get the total positive charge."""
    Hmol = Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
    res = []
    for atom in Hmol.GetAtoms():
        res.append(float(atom.GetProp('_GasteigerCharge')))
    if res == []:
        return 0
    else:
        cc = numpy.array(res, 'd')
        return round(sum(cc[cc > 0]), 3)


def CalculateMeanPCharge(mol: Chem.Mol) -> float:
    """Get the average positive charge."""
    Hmol = Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
    res = []
    for atom in Hmol.GetAtoms():
        res.append(float(atom.GetProp('_GasteigerCharge')))
    if res == []:
        return 0
    else:
        cc = numpy.array(res, 'd')
        return round(numpy.mean(cc[cc > 0]), 3)


def CalculateTotalNCharge(mol: Chem.Mol) -> float:
    """Ge the total negative charge."""
    Hmol = Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
    res = []
    for atom in Hmol.GetAtoms():
        res.append(float(atom.GetProp('_GasteigerCharge')))
    if res == []:
        return 0
    else:
        cc = numpy.array(res, 'd')
        return round(sum(cc[cc < 0]), 3)


def CalculateMeanNCharge(mol: Chem.Mol) -> float:
    """Get the average negative charge."""
    Hmol = Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
    res = []
    for atom in Hmol.GetAtoms():
        res.append(float(atom.GetProp('_GasteigerCharge')))
    if res == []:
        return 0
    else:
        cc = numpy.array(res, 'd')
        return round(numpy.mean(cc[cc < 0]), 3)


def CalculateTotalAbsoulteCharge(mol: Chem.Mol):
    """Get the total absolute charge."""
    Hmol = Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
    res = []
    for atom in Hmol.GetAtoms():
        res.append(float(atom.GetProp('_GasteigerCharge')))
    if res == []:
        return 0
    else:
        cc = numpy.array(res, 'd')
        return round(sum(numpy.absolute(cc)), 3)


def CalculateMeanAbsoulteCharge(mol: Chem.Mol) -> float:
    """Get the average absolute charge."""
    Hmol = Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
    res = []
    for atom in Hmol.GetAtoms():
        res.append(float(atom.GetProp('_GasteigerCharge')))
    if res == []:
        return 0
    else:
        cc = numpy.array(res, 'd')
        return round(numpy.mean(numpy.absolute(cc)), 3)


def CalculateRelativePCharge(mol: Chem.Mol) -> float:
    """Get the ratio between the most positive partial charge and the total positive charge."""
    Hmol = Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
    res = []
    for atom in Hmol.GetAtoms():
        res.append(float(atom.GetProp('_GasteigerCharge')))
    if res == []:
        return 0
    else:
        cc = numpy.array(res, 'd')
        if sum(cc[cc > 0]) == 0:
            return 0
        else:
            return round(max(res) / sum(cc[cc > 0]), 3)


def CalculateRelativeNCharge(mol: Chem.Mol) -> float:
    """Get the ratio between the most negative partial charge and the total negative charge."""
    Hmol = Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
    res = []
    for atom in Hmol.GetAtoms():
        res.append(float(atom.GetProp('_GasteigerCharge')))
    if res == []:
        return 0
    else:
        cc = numpy.array(res, 'd')
        if sum(cc[cc < 0]) == 0:
            return 0
        else:
            return round(min(res) / sum(cc[cc < 0]), 3)


def CalculateLocalDipoleIndex(mol: Chem.Mol) -> float:
    """Calculate the local dipole index (D)."""
    GMCharge.ComputeGasteigerCharges(mol, iter_step)
    res = []
    for atom in mol.GetAtoms():
        res.append(float(atom.GetProp('_GasteigerCharge')))
    cc = [numpy.absolute(res[x.GetBeginAtom().GetIdx()] - res[x.GetEndAtom().GetIdx()]) for x in mol.GetBonds()]
    B = len(mol.GetBonds())
    return 0 if len(cc) == 0.0 else round(sum(cc) / B, 3)


def CalculateSubmolPolarityPara(mol: Chem.Mol) -> float:
    """Calculate the submolecular polarity parameter (SPP)."""
    return round(CalculateAllMaxPCharge(mol) - CalculateAllMaxNCharge(mol), 3)


_Charge = {'SPP': CalculateSubmolPolarityPara,
           'LDI': CalculateLocalDipoleIndex,
           'Rnc': CalculateRelativeNCharge,
           'Rpc': CalculateRelativePCharge,
           'Mac': CalculateMeanAbsoulteCharge,
           'Tac': CalculateTotalAbsoulteCharge,
           'Mnc': CalculateMeanNCharge,
           'Tnc': CalculateTotalNCharge,
           'Mpc': CalculateMeanPCharge,
           'Tpc': CalculateTotalPCharge,
           'Qass': CalculateAllSumSquareCharge,
           'QOss': CalculateOSumSquareCharge,
           'QNss': CalculateNSumSquareCharge,
           'QCss': CalculateCSumSquareCharge,
           'QHss': CalculateHSumSquareCharge,
           'Qmin': CalculateAllMaxNCharge,
           'Qmax': CalculateAllMaxPCharge,
           'QOmin': CalculateOMaxNCharge,
           'QNmin': CalculateNMaxNCharge,
           'QCmin': CalculateCMaxNCharge,
           'QHmin': CalculateHMaxNCharge,
           'QOmax': CalculateOMaxPCharge,
           'QNmax': CalculateNMaxPCharge,
           'QCmax': CalculateCMaxPCharge,
           'QHmax': CalculateHMaxPCharge}


def GetCharge(mol: Chem.Mol) -> dict:
    """Get all (25) constitutional descriptors."""
    result = {}
    for DesLabel in _Charge.keys():
        result[DesLabel] = _Charge[DesLabel](mol)
    return result


##############################################################################
# if __name__ =='__main__':
#     smis = ['CCCC','CCCCC','CCCCCC','CC(N)C(=O)O','CC(N)C(=O)[O-].[Na+]']
#     smi5=['CCCCCC','CCC(C)CC','CC(C)CCC','CC(C)C(C)C','CCCCCN','c1ccccc1N']
#     for index, smi in enumerate(smis):
#         m = Chem.MolFromSmiles(smi)
#         print(index+1)
#         print(smi)
#         print('\t',GetCharge(m))
#         print(len(GetCharge(m)))
