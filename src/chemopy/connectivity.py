# -*- coding: utf-8 -*-

"""Molecular connectivity topological indices."""

from typing import List

import numpy
from rdkit import Chem
from rdkit.Chem import rdchem

from chemopy.topology import _HallKierDeltas

periodicTable = rdchem.GetPeriodicTable()


def CalculateChi0(mol: Chem.Mol) -> float:
    """Calculate molecular connectivity chi index for path order 0."""
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    while 0 in deltas:
        deltas.remove(0)
    deltas = numpy.array(deltas, 'd')
    res = sum(numpy.sqrt(1. / deltas))
    return res


def CalculateChi1(mol: Chem.Mol):
    """Calculate molecular connectivity chi index for path order 1."""
    cc = [x.GetBeginAtom().GetDegree() * x.GetEndAtom().GetDegree() for x in mol.GetBonds()]
    if len(cc) == 0:
        return 0.0
    while 0 in cc:
        cc.remove(0)
    cc = numpy.array(cc, 'd')
    res = sum(numpy.sqrt(1. / cc))
    return res


def CalculateMeanRandic(mol: Chem.Mol) -> float:
    """Calculate mean chi1 (Randic) connectivity index."""
    cc = [x.GetBeginAtom().GetDegree() * x.GetEndAtom().GetDegree() for x in mol.GetBonds()]
    if len(cc) == 0:
        return 0.0
    while 0 in cc:
        cc.remove(0)
    cc = numpy.array(cc, 'd')
    res = numpy.mean(numpy.sqrt(1. / cc))
    return res


def _CalculateChinp(mol: Chem.Mol, NumPath: int = 2) -> float:
    """Calculate molecular connectivity chi index for path order 2."""
    accum = 0.0
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    for path in Chem.FindAllPathsOfLengthN(mol, NumPath + 1, useBonds=0):
        cAccum = 1.0
        for idx in path:
            cAccum *= deltas[idx]
        if cAccum:
            accum += 1. / numpy.sqrt(cAccum)
    return accum


def CalculateChi2(mol: Chem.Mol) -> float:
    """Calculate molecular connectivity chi index for path order 2."""
    return _CalculateChinp(mol, NumPath=2)


def CalculateChi3p(mol: Chem.Mol) -> float:
    """Calculate molecular connectivity chi index for path order 3."""
    return _CalculateChinp(mol, NumPath=3)


def CalculateChi4p(mol: Chem.Mol) -> float:
    """Calculate molecular connectivity chi index for path order 4."""
    return _CalculateChinp(mol, NumPath=4)


def CalculateChi5p(mol: Chem.Mol) -> float:
    """Calculate molecular connectivity chi index for path order 5."""
    return _CalculateChinp(mol, NumPath=5)


def CalculateChi6p(mol: Chem.Mol) -> float:
    """Calculate molecular connectivity chi index for path order 6."""
    return _CalculateChinp(mol, NumPath=6)


def CalculateChi7p(mol: Chem.Mol) -> float:
    """Calculate molecular connectivity chi index for path order 7."""
    return _CalculateChinp(mol, NumPath=7)


def CalculateChi8p(mol: Chem.Mol) -> float:
    """Calculate molecular connectivity chi index for path order 8."""
    return _CalculateChinp(mol, NumPath=8)


def CalculateChi9p(mol: Chem.Mol) -> float:
    """Calculate molecular connectivity chi index for path order 9."""
    return _CalculateChinp(mol, NumPath=9)


def CalculateChi10p(mol: Chem.Mol) -> float:
    """Calculate molecular connectivity chi index for path order 10."""
    return _CalculateChinp(mol, NumPath=10)


def CalculateChi3c(mol: Chem.Mol) -> float:
    """Calculate molecular connectivity chi index for cluster."""
    accum = 0.0
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    patt = Chem.MolFromSmarts('*~*(~*)~*')
    HPatt = mol.GetSubstructMatches(patt)
    for cluster in HPatt:
        deltas = [mol.GetAtomWithIdx(x).GetDegree() for x in cluster]
        while 0 in deltas:
            deltas.remove(0)
        if deltas != []:
            deltas1 = numpy.array(deltas, numpy.float)
            accum = accum + 1. / numpy.sqrt(deltas1.prod())
    return accum


def CalculateChi4c(mol: Chem.Mol) -> float:
    """Calculate molecular connectivity chi index for cluster."""
    accum = 0.0
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    patt = Chem.MolFromSmarts('*~*(~*)(~*)~*')
    HPatt = mol.GetSubstructMatches(patt)
    for cluster in HPatt:
        deltas = [mol.GetAtomWithIdx(x).GetDegree() for x in cluster]
        while 0 in deltas:
            deltas.remove(0)
        if deltas != []:
            deltas1 = numpy.array(deltas, numpy.float)
            accum = accum + 1. / numpy.sqrt(deltas1.prod())
    return accum


def CalculateChi4pc(mol: Chem.Mol) -> float:
    """Calculate molecular connectivity chi index for path/cluster."""
    accum = 0.0
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    patt = Chem.MolFromSmarts('*~*(~*)~*~*')
    HPatt = mol.GetSubstructMatches(patt)
    for cluster in HPatt:
        deltas = [mol.GetAtomWithIdx(x).GetDegree() for x in cluster]
        while 0 in deltas:
            deltas.remove(0)
        if deltas != []:
            deltas1 = numpy.array(deltas, numpy.float)
            accum = accum + 1. / numpy.sqrt(deltas1.prod())
    return accum


def CalculateDeltaChi3c4pc(mol: Chem.Mol) -> float:
    """Calculate the difference between chi3c and chi4pc."""
    return abs(CalculateChi3c(mol) - CalculateChi4pc(mol))


def _CalculateChinch(mol: Chem.Mol, NumCycle=3):
    """Calculate molecular connectivity chi index for cycles of n."""
    accum = 0.0
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    for tup in mol.GetRingInfo().AtomRings():
        cAccum = 1.0
        if len(tup) == NumCycle:
            for idx in tup:
                cAccum *= deltas[idx]
            if cAccum:
                accum += 1. / numpy.sqrt(cAccum)
    return accum


def CalculateChi3ch(mol: Chem.Mol) -> float:
    """Calculate molecular connectivity chi index for cycles of 3."""
    return _CalculateChinch(mol, NumCycle=3)


def CalculateChi4ch(mol: Chem.Mol) -> float:
    """Calculate molecular connectivity chi index for cycles of 4."""
    return _CalculateChinch(mol, NumCycle=4)


def CalculateChi5ch(mol: Chem.Mol) -> float:
    """Calculate molecular connectivity chi index for cycles of 5."""
    return _CalculateChinch(mol, NumCycle=5)


def CalculateChi6ch(mol: Chem.Mol) -> float:
    """Calculate molecular connectivity chi index for cycles of 6."""
    return _CalculateChinch(mol, NumCycle=6)


def CalculateChiv0(mol: Chem.Mol) -> float:
    """Calculate valence molecular connectivity chi index for path order 0."""
    deltas = _HallKierDeltas(mol, skipHs=0)
    while 0 in deltas:
        deltas.remove(0)
    deltas = numpy.array(deltas, 'd')
    res = sum(numpy.sqrt(1. / deltas))
    return res


def _CalculateChivnp(mol: Chem.Mol, NumPath: int = 1) -> float:
    """Calculate valence molecular connectivity chi index for path order 1."""
    accum = 0.0
    deltas = _HallKierDeltas(mol, skipHs=0)
    for path in Chem.FindAllPathsOfLengthN(mol, NumPath + 1, useBonds=0):
        cAccum = 1.0
        for idx in path:
            cAccum *= deltas[idx]
        if cAccum:
            accum += 1. / numpy.sqrt(cAccum)
    return accum


def CalculateChiv1(mol: Chem.Mol) -> float:
    """Calculate valence molecular connectivity chi index for path order 1."""
    return _CalculateChivnp(mol, NumPath=1)


def CalculateChiv2(mol: Chem.Mol) -> float:
    """Calculate valence molecular connectivity chi index for path order 2."""
    return _CalculateChivnp(mol, NumPath=2)


def CalculateChiv3p(mol: Chem.Mol) -> float:
    """Calculate valence molecular connectivity chi index for path order 3."""
    return _CalculateChivnp(mol, NumPath=3)


def CalculateChiv4p(mol: Chem.Mol) -> float:
    """Calculate valence molecular connectivity chi index for path order 4."""
    return _CalculateChivnp(mol, NumPath=4)


def CalculateChiv5p(mol: Chem.Mol) -> float:
    """Calculate valence molecular connectivity chi index for path order 5."""
    return _CalculateChivnp(mol, NumPath=5)


def CalculateChiv6p(mol: Chem.Mol) -> float:
    """Calculate valence molecular connectivity chi index for path order 6."""
    return _CalculateChivnp(mol, NumPath=6)


def CalculateChiv7p(mol: Chem.Mol) -> float:
    """Calculate valence molecular connectivity chi index for path order 7."""
    return _CalculateChivnp(mol, NumPath=7)


def CalculateChiv8p(mol: Chem.Mol) -> float:
    """Calculate valence molecular connectivity chi index for path order 8."""
    return _CalculateChivnp(mol, NumPath=8)


def CalculateChiv9p(mol: Chem.Mol) -> float:
    """Calculate valence molecular connectivity chi index for path order 9."""
    return _CalculateChivnp(mol, NumPath=9)


def CalculateChiv10p(mol: Chem.Mol) -> float:
    """Calculate valence molecular connectivity chi index for path order 10."""
    return _CalculateChivnp(mol, NumPath=10)


def CalculateDeltaChi0(mol: Chem.Mol) -> float:
    """Calculate the difference between chi0v and chi0."""
    return abs(CalculateChiv0(mol) - CalculateChi0(mol))


def CalculateDeltaChi1(mol: Chem.Mol) -> float:
    """Calculate the difference between chi1v and chi1."""
    return abs(CalculateChiv1(mol) - CalculateChi1(mol))


def CalculateDeltaChi2(mol: Chem.Mol) -> float:
    """Calculate the difference between chi2v and chi2."""
    return abs(_CalculateChivnp(mol, NumPath=2) - _CalculateChinp(mol, NumPath=2))


def CalculateDeltaChi3(mol: Chem.Mol) -> float:
    """Calculate the difference between chi3v and chi3."""
    return abs(_CalculateChivnp(mol, NumPath=3) - _CalculateChinp(mol, NumPath=3))


def CalculateDeltaChi4(mol: Chem.Mol) -> float:
    """Calculate the difference between chi4v and chi4."""
    return abs(_CalculateChivnp(mol, NumPath=4) - _CalculateChinp(mol, NumPath=4))


def _AtomHallKierDeltas(atom: Chem.Atom, skipHs: bool = False) -> List[float]:
    """Calculate Kier & Hall atomic valence delta-values for molecular connectivity.

    From Kier L. and Hall L., J. Pharm. Sci. (1983), 72(10),1170-1173.
    """
    global periodicTable
    res = []
    n = atom.GetAtomicNum()
    if n > 1:
        nV = periodicTable.GetNOuterElecs(n)
        nHs = atom.GetTotalNumHs()
        if n < 10:
            res.append(float(nV - nHs))
        else:
            res.append(float(nV - nHs) / float(n - nV - 1))
    elif not skipHs:
        res.append(0.0)
    return res


def CalculateChiv3c(mol: Chem.Mol) -> float:
    """Calculate valence molecular connectivity chi index for cluster."""
    accum = 0.0
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    patt = Chem.MolFromSmarts('*~*(~*)~*')
    HPatt = mol.GetSubstructMatches(patt)
    for cluster in HPatt:
        deltas = [_AtomHallKierDeltas(mol.GetAtomWithIdx(x)) for x in cluster]
        while 0 in deltas:
            deltas.remove(0)
        if deltas != []:
            deltas1 = numpy.array(deltas, numpy.float)
            accum = accum + 1. / numpy.sqrt(deltas1.prod())
    return accum


def CalculateChiv4c(mol: Chem.Mol) -> float:
    """Calculate valence molecular connectivity chi index for cluster."""
    accum = 0.0
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    patt = Chem.MolFromSmarts('*~*(~*)(~*)~*')
    HPatt = mol.GetSubstructMatches(patt)
    for cluster in HPatt:
        deltas = [_AtomHallKierDeltas(mol.GetAtomWithIdx(x)) for x in cluster]
        while 0 in deltas:
            deltas.remove(0)
        if deltas != []:
            deltas1 = numpy.array(deltas, numpy.float)
            accum = accum + 1. / numpy.sqrt(deltas1.prod())
    return accum


def CalculateChiv4pc(mol: Chem.Mol) -> float:
    """Calculate valence molecular connectivity chi index for path/cluster."""
    accum = 0.0
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    patt = Chem.MolFromSmarts('*~*(~*)~*~*')
    HPatt = mol.GetSubstructMatches(patt)
    for cluster in HPatt:
        deltas = [_AtomHallKierDeltas(mol.GetAtomWithIdx(x)) for x in cluster]
        while 0 in deltas:
            deltas.remove(0)
        if deltas != []:
            deltas1 = numpy.array(deltas, numpy.float)
            accum = accum + 1. / numpy.sqrt(deltas1.prod())
    return accum


def CalculateDeltaChiv3c4pc(mol: Chem.Mol) -> float:
    """Calculate the difference between chiv3c and chiv4pc."""
    return abs(CalculateChiv3c(mol) - CalculateChiv4pc(mol))


def _CalculateChivnch(mol: Chem.Mol, NumCyc=3):
    """Calculate valence molecular connectivity chi index for cycles of n."""
    accum = 0.0
    deltas = _HallKierDeltas(mol, skipHs=0)
    for tup in mol.GetRingInfo().AtomRings():
        cAccum = 1.0
        if len(tup) == NumCyc:
            for idx in tup:
                cAccum *= deltas[idx]
            if cAccum:
                accum += 1. / numpy.sqrt(cAccum)
    return accum


def CalculateChiv3ch(mol: Chem.Mol) -> float:
    """Calculate valence molecular connectivity chi index for cycles of 3."""
    return _CalculateChivnch(mol, NumCyc=3)


def CalculateChiv4ch(mol: Chem.Mol) -> float:
    """Calculate valence molecular connectivity chi index for cycles of 4."""
    return _CalculateChivnch(mol, NumCyc=4)


def CalculateChiv5ch(mol: Chem.Mol) -> float:
    """Calculate valence molecular connectivity chi index for cycles of 5."""
    return _CalculateChivnch(mol, NumCyc=5)


def CalculateChiv6ch(mol: Chem.Mol) -> float:
    """Calculate valence molecular connectivity chi index for cycles of 6."""
    return _CalculateChivnch(mol, NumCyc=6)


_connectivity = {'Chi0': CalculateChi0,
                 'Chi1': CalculateChi1,
                 'mChi1': CalculateMeanRandic,
                 'Chi2': CalculateChi2,
                 'Chi3': CalculateChi3p,
                 'Chi4': CalculateChi4p,
                 'Chi5': CalculateChi5p,
                 'Chi6': CalculateChi6p,
                 'Chi7': CalculateChi7p,
                 'Chi8': CalculateChi8p,
                 'Chi9': CalculateChi9p,
                 'Chi10': CalculateChi10p,
                 'Chi3c': CalculateChi3c,
                 'Chi4c': CalculateChi4c,
                 'Chi4pc': CalculateChi4pc,
                 'Chi3ch': CalculateChi3ch,
                 'Chi4ch': CalculateChi4ch,
                 'Chi5ch': CalculateChi5ch,
                 'Chi6ch': CalculateChi6ch,
                 'knotp': CalculateDeltaChi3c4pc,
                 'Chiv0': CalculateChiv0,
                 'Chiv1': CalculateChiv1,
                 'Chiv2': CalculateChiv2,
                 'Chiv3': CalculateChiv3p,
                 'Chiv4': CalculateChiv4p,
                 'Chiv5': CalculateChiv5p,
                 'Chiv6': CalculateChiv6p,
                 'Chiv7': CalculateChiv7p,
                 'Chiv8': CalculateChiv8p,
                 'Chiv9': CalculateChiv9p,
                 'Chiv10': CalculateChiv10p,
                 'dchi0': CalculateDeltaChi0,
                 'dchi1': CalculateDeltaChi1,
                 'dchi2': CalculateDeltaChi2,
                 'dchi3': CalculateDeltaChi3,
                 'dchi4': CalculateDeltaChi4,
                 'Chiv3c': CalculateChiv3c,
                 'Chiv4c': CalculateChiv4c,
                 'Chiv4pc': CalculateChiv4pc,
                 'Chiv3ch': CalculateChiv3ch,
                 'Chiv4ch': CalculateChiv4ch,
                 'Chiv5ch': CalculateChiv5ch,
                 'Chiv6ch': CalculateChiv6ch,
                 'knotpv': CalculateDeltaChiv3c4pc}


def GetConnectivity(mol: Chem.Mol) -> dict:
    """Get all (44) connectivity descriptors."""
    result = {}
    for DesLabel in _connectivity.keys():
        result[DesLabel] = round(_connectivity[DesLabel](mol), 3)
    return result


###############################################################################
# if __name__ =='__main__':
#     smis = ['CCCC','CCCCC','CCCCCC','CC(N)C(=O)O','CC(N)C(=O)[O-].[Na+]']
#     smi5=['CCCCCC','CCC(C)CC','CC(C)CCC','CC(C)C(C)C','CCCCCN','c1ccccc1N']
#     for index, smi in enumerate(smis):
#         m = Chem.MolFromSmiles(smi)
#         print(index+1)
#         print(smi)
#         print('\t',GetConnectivity(m))
#         print('\t',len(GetConnectivity(m)))
