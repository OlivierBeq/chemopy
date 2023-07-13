# -*- coding: utf-8 -*-

"""Molecular connectivity topological indices."""

from typing import List

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdchem

from .topology import Topology

periodicTable = rdchem.GetPeriodicTable()


class Connectivity:
    """Topological indices"""

    @staticmethod
    def calculate_chi0(mol: Chem.Mol) -> float:
        """Calculate molecular connectivity chi index for path order 0."""
        deltas = [x.GetDegree() for x in mol.GetAtoms()]
        while [0] in deltas:
            deltas.remove([0])
        deltas = np.array(deltas, 'd')
        res = sum(np.sqrt(1. / deltas))
        return res

    @staticmethod
    def calculate_chi1(mol: Chem.Mol) -> float:
        """Calculate molecular connectivity chi index for path order 1."""
        cc = [x.GetBeginAtom().GetDegree() * x.GetEndAtom().GetDegree() for x in mol.GetBonds()]
        if len(cc) == 0:
            return 0.0
        while 0 in cc:
            cc.remove(0)
        cc = np.array(cc, 'd')
        res = sum(np.sqrt(1. / cc))
        return res

    @staticmethod
    def calculate_mean_randic(mol: Chem.Mol) -> float:
        """Calculate mean chi1 (Randic) connectivity index."""
        cc = [x.GetBeginAtom().GetDegree() * x.GetEndAtom().GetDegree() for x in mol.GetBonds()]
        if len(cc) == 0:
            return 0.0
        while 0 in cc:
            cc.remove(0)
        cc = np.array(cc, 'd')
        res = np.mean(np.sqrt(1. / cc))
        return res

    @staticmethod
    def _calculate_chinp(mol: Chem.Mol, NumPath: int = 2) -> float:
        """Calculate molecular connectivity chi index for path order 2."""
        accum = 0.0
        deltas = [x.GetDegree() for x in mol.GetAtoms()]
        for path in Chem.FindAllPathsOfLengthN(mol, NumPath + 1, useBonds=0):
            cAccum = 1.0
            for idx in path:
                cAccum *= deltas[idx]
            if cAccum:
                accum += 1. / np.sqrt(cAccum)
        return accum

    @staticmethod
    def calculate_chi2(mol: Chem.Mol) -> float:
        """Calculate molecular connectivity chi index for path order 2."""
        return Connectivity._calculate_chinp(mol, NumPath=2)

    @staticmethod
    def calculate_chi3p(mol: Chem.Mol) -> float:
        """Calculate molecular connectivity chi index for path order 3."""
        return Connectivity._calculate_chinp(mol, NumPath=3)

    @staticmethod
    def calculate_chi4p(mol: Chem.Mol) -> float:
        """Calculate molecular connectivity chi index for path order 4."""
        return Connectivity._calculate_chinp(mol, NumPath=4)

    @staticmethod
    def calculate_chi5p(mol: Chem.Mol) -> float:
        """Calculate molecular connectivity chi index for path order 5."""
        return Connectivity._calculate_chinp(mol, NumPath=5)

    @staticmethod
    def calculate_chi6p(mol: Chem.Mol) -> float:
        """Calculate molecular connectivity chi index for path order 6."""
        return Connectivity._calculate_chinp(mol, NumPath=6)

    @staticmethod
    def calculate_chi7p(mol: Chem.Mol) -> float:
        """Calculate molecular connectivity chi index for path order 7."""
        return Connectivity._calculate_chinp(mol, NumPath=7)

    @staticmethod
    def calculate_chi8p(mol: Chem.Mol) -> float:
        """Calculate molecular connectivity chi index for path order 8."""
        return Connectivity._calculate_chinp(mol, NumPath=8)

    @staticmethod
    def calculate_chi9p(mol: Chem.Mol) -> float:
        """Calculate molecular connectivity chi index for path order 9."""
        return Connectivity._calculate_chinp(mol, NumPath=9)

    @staticmethod
    def calculate_chi10p(mol: Chem.Mol) -> float:
        """Calculate molecular connectivity chi index for path order 10."""
        return Connectivity._calculate_chinp(mol, NumPath=10)

    @staticmethod
    def calculate_chi3c(mol: Chem.Mol) -> float:
        """Calculate molecular connectivity chi index for cluster."""
        accum = 0.0
        patt = Chem.MolFromSmarts('*~*(~*)~*')
        HPatt = mol.GetSubstructMatches(patt)
        for cluster in HPatt:
            deltas = [mol.GetAtomWithIdx(x).GetDegree() for x in cluster]
            while [0] in deltas:
                deltas.remove([0])
            if deltas != []:
                deltas1 = np.array(deltas, dtype=float)
                accum = accum + 1. / np.sqrt(deltas1.prod())
        return accum

    @staticmethod
    def calculate_chi4c(mol: Chem.Mol) -> float:
        """Calculate molecular connectivity chi index for cluster."""
        accum = 0.0
        patt = Chem.MolFromSmarts('*~*(~*)(~*)~*')
        HPatt = mol.GetSubstructMatches(patt)
        for cluster in HPatt:
            deltas = [mol.GetAtomWithIdx(x).GetDegree() for x in cluster]
            while [0] in deltas:
                deltas.remove([0])
            if deltas != []:
                deltas1 = np.array(deltas, dtype=float)
                accum = accum + 1. / np.sqrt(deltas1.prod())
        return accum

    @staticmethod
    def calculate_chi4pc(mol: Chem.Mol) -> float:
        """Calculate molecular connectivity chi index for path/cluster."""
        accum = 0.0
        patt = Chem.MolFromSmarts('*~*(~*)~*~*')
        HPatt = mol.GetSubstructMatches(patt)
        for cluster in HPatt:
            deltas = [mol.GetAtomWithIdx(x).GetDegree() for x in cluster]
            while [0] in deltas:
                deltas.remove([0])
            if deltas != []:
                deltas1 = np.array(deltas, dtype=float)
                accum = accum + 1. / np.sqrt(deltas1.prod())
        return accum

    @staticmethod
    def calculate_delta_chi3c4pc(mol: Chem.Mol) -> float:
        """Calculate the difference between chi3c and chi4pc."""
        return abs(Connectivity.calculate_chi3c(mol) - Connectivity.calculate_chi4pc(mol))

    @staticmethod
    def _calculate_chinch(mol: Chem.Mol, NumCycle=3) -> float:
        """Calculate molecular connectivity chi index for cycles of n."""
        accum = 0.0
        deltas = [x.GetDegree() for x in mol.GetAtoms()]
        for tup in mol.GetRingInfo().AtomRings():
            cAccum = 1.0
            if len(tup) == NumCycle:
                for idx in tup:
                    cAccum *= deltas[idx]
                if cAccum:
                    accum += 1. / np.sqrt(cAccum)
        return accum

    @staticmethod
    def calculate_chi3ch(mol: Chem.Mol) -> float:
        """Calculate molecular connectivity chi index for cycles of 3."""
        return Connectivity._calculate_chinch(mol, NumCycle=3)

    @staticmethod
    def calculate_chi4ch(mol: Chem.Mol) -> float:
        """Calculate molecular connectivity chi index for cycles of 4."""
        return Connectivity._calculate_chinch(mol, NumCycle=4)

    @staticmethod
    def calculate_chi5ch(mol: Chem.Mol) -> float:
        """Calculate molecular connectivity chi index for cycles of 5."""
        return Connectivity._calculate_chinch(mol, NumCycle=5)

    @staticmethod
    def calculate_chi6ch(mol: Chem.Mol) -> float:
        """Calculate molecular connectivity chi index for cycles of 6."""
        return Connectivity._calculate_chinch(mol, NumCycle=6)

    @staticmethod
    def calculate_chiv0(mol: Chem.Mol) -> float:
        """Calculate valence molecular connectivity chi index for path order 0."""
        deltas = Topology._hall_kier_deltas(mol, skipHs=0)
        while 0 in deltas:
            deltas.remove(0)
        deltas = np.array(deltas, 'd')
        res = sum(np.sqrt(1. / deltas))
        return res

    @staticmethod
    def _calculate_chivnp(mol: Chem.Mol, NumPath: int = 1) -> float:
        """Calculate valence molecular connectivity chi index for path order 1."""
        accum = 0.0
        deltas = Topology._hall_kier_deltas(mol, skipHs=False)
        for path in Chem.FindAllPathsOfLengthN(mol, NumPath + 1, useBonds=0):
            cAccum = 1.0
            for idx in path:
                cAccum *= deltas[idx]
            if cAccum:
                accum += 1. / np.sqrt(cAccum)
        return accum

    @staticmethod
    def calculate_chiv1(mol: Chem.Mol) -> float:
        """Calculate valence molecular connectivity chi index for path order 1."""
        return Connectivity._calculate_chivnp(mol, NumPath=1)

    @staticmethod
    def calculate_chiv2(mol: Chem.Mol) -> float:
        """Calculate valence molecular connectivity chi index for path order 2."""
        return Connectivity._calculate_chivnp(mol, NumPath=2)

    @staticmethod
    def calculate_chiv3p(mol: Chem.Mol) -> float:
        """Calculate valence molecular connectivity chi index for path order 3."""
        return Connectivity._calculate_chivnp(mol, NumPath=3)

    @staticmethod
    def calculate_chiv4p(mol: Chem.Mol) -> float:
        """Calculate valence molecular connectivity chi index for path order 4."""
        return Connectivity._calculate_chivnp(mol, NumPath=4)

    @staticmethod
    def calculate_chiv5p(mol: Chem.Mol) -> float:
        """Calculate valence molecular connectivity chi index for path order 5."""
        return Connectivity._calculate_chivnp(mol, NumPath=5)

    @staticmethod
    def calculate_chiv6p(mol: Chem.Mol) -> float:
        """Calculate valence molecular connectivity chi index for path order 6."""
        return Connectivity._calculate_chivnp(mol, NumPath=6)

    @staticmethod
    def calculate_chiv7p(mol: Chem.Mol) -> float:
        """Calculate valence molecular connectivity chi index for path order 7."""
        return Connectivity._calculate_chivnp(mol, NumPath=7)

    @staticmethod
    def calculate_chiv8p(mol: Chem.Mol) -> float:
        """Calculate valence molecular connectivity chi index for path order 8."""
        return Connectivity._calculate_chivnp(mol, NumPath=8)

    @staticmethod
    def calculate_chiv9p(mol: Chem.Mol) -> float:
        """Calculate valence molecular connectivity chi index for path order 9."""
        return Connectivity._calculate_chivnp(mol, NumPath=9)

    @staticmethod
    def calculate_chiv10p(mol: Chem.Mol) -> float:
        """Calculate valence molecular connectivity chi index for path order 10."""
        return Connectivity._calculate_chivnp(mol, NumPath=10)

    @staticmethod
    def calculate_delta_chi0(mol: Chem.Mol) -> float:
        """Calculate the difference between chi0v and chi0."""
        return abs(Connectivity.calculate_chiv0(mol) - Connectivity.calculate_chi0(mol))

    @staticmethod
    def calculate_delta_chi1(mol: Chem.Mol) -> float:
        """Calculate the difference between chi1v and chi1."""
        return abs(Connectivity.calculate_chiv1(mol) - Connectivity.calculate_chi1(mol))

    @staticmethod
    def calculate_delta_chi2(mol: Chem.Mol) -> float:
        """Calculate the difference between chi2v and chi2."""
        return abs(Connectivity._calculate_chivnp(mol, NumPath=2) - Connectivity._calculate_chinp(mol, NumPath=2))

    @staticmethod
    def calculate_delta_chi3(mol: Chem.Mol) -> float:
        """Calculate the difference between chi3v and chi3."""
        return abs(Connectivity._calculate_chivnp(mol, NumPath=3) - Connectivity._calculate_chinp(mol, NumPath=3))

    @staticmethod
    def calculate_delta_chi4(mol: Chem.Mol) -> float:
        """Calculate the difference between chi4v and chi4."""
        return abs(Connectivity._calculate_chivnp(mol, NumPath=4) - Connectivity._calculate_chinp(mol, NumPath=4))

    @staticmethod
    def _atom_hall_kier_deltas(atom: Chem.Atom, skipHs: bool = False) -> List[float]:
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

    @staticmethod
    def calculate_chiv3c(mol: Chem.Mol) -> float:
        """Calculate valence molecular connectivity chi index for cluster."""
        accum = 0.0
        patt = Chem.MolFromSmarts('*~*(~*)~*')
        HPatt = mol.GetSubstructMatches(patt)
        for cluster in HPatt:
            deltas = [Connectivity._atom_hall_kier_deltas(mol.GetAtomWithIdx(x)) for x in cluster]
            while [0] in deltas:
                deltas.remove([0])
            if deltas != []:
                deltas1 = np.array(deltas, dtype=float)
                accum = accum + 1. / np.sqrt(deltas1.prod())
        return accum

    @staticmethod
    def calculate_chiv4c(mol: Chem.Mol) -> float:
        """Calculate valence molecular connectivity chi index for cluster."""
        accum = 0.0
        patt = Chem.MolFromSmarts('*~*(~*)(~*)~*')
        HPatt = mol.GetSubstructMatches(patt)
        for cluster in HPatt:
            deltas = [Connectivity._atom_hall_kier_deltas(mol.GetAtomWithIdx(x)) for x in cluster]
            while [0] in deltas:
                deltas.remove([0])
            if deltas != []:
                deltas1 = np.array(deltas, dtype=float)
                accum = accum + 1. / np.sqrt(deltas1.prod())
        return accum

    @staticmethod
    def calculate_chiv4pc(mol: Chem.Mol) -> float:
        """Calculate valence molecular connectivity chi index for path/cluster."""
        accum = 0.0
        patt = Chem.MolFromSmarts('*~*(~*)~*~*')
        HPatt = mol.GetSubstructMatches(patt)
        for cluster in HPatt:
            deltas = [Connectivity._atom_hall_kier_deltas(mol.GetAtomWithIdx(x)) for x in cluster]
            while [0] in deltas:
                deltas.remove([0])
            if deltas != []:
                deltas1 = np.array(deltas, dtype=float)
                accum = accum + 1. / np.sqrt(deltas1.prod())
        return accum

    @staticmethod
    def calculate_delta_chiv3c4pc(mol: Chem.Mol) -> float:
        """Calculate the difference between chiv3c and chiv4pc."""
        return abs(Connectivity.calculate_chiv3c(mol) - Connectivity.calculate_chiv4pc(mol))

    @staticmethod
    def _calculate_chivnch(mol: Chem.Mol, NumCyc=3) -> float:
        """Calculate valence molecular connectivity chi index for cycles of n."""
        accum = 0.0
        deltas = Topology._hall_kier_deltas(mol, skipHs=0)
        for tup in mol.GetRingInfo().AtomRings():
            cAccum = 1.0
            if len(tup) == NumCyc:
                for idx in tup:
                    cAccum *= deltas[idx]
                if cAccum:
                    accum += 1. / np.sqrt(cAccum)
        return accum

    @staticmethod
    def calculate_chiv3ch(mol: Chem.Mol) -> float:
        """Calculate valence molecular connectivity chi index for cycles of 3."""
        return Connectivity._calculate_chivnch(mol, NumCyc=3)

    @staticmethod
    def calculate_chiv4ch(mol: Chem.Mol) -> float:
        """Calculate valence molecular connectivity chi index for cycles of 4."""
        return Connectivity._calculate_chivnch(mol, NumCyc=4)

    @staticmethod
    def calculate_chiv5ch(mol: Chem.Mol) -> float:
        """Calculate valence molecular connectivity chi index for cycles of 5."""
        return Connectivity._calculate_chivnch(mol, NumCyc=5)

    @staticmethod
    def calculate_chiv6ch(mol: Chem.Mol) -> float:
        """Calculate valence molecular connectivity chi index for cycles of 6."""
        return Connectivity._calculate_chivnch(mol, NumCyc=6)


    _connectivity = {'Chiv0': calculate_chiv0,
                     'Chiv1': calculate_chiv1,
                     'Chiv2': calculate_chiv2,
                     'Chiv3': calculate_chiv3p,
                     'Chiv4': calculate_chiv4p,
                     'Chiv5': calculate_chiv5p,
                     'Chiv6': calculate_chiv6p,
                     'Chiv7': calculate_chiv7p,
                     'Chiv8': calculate_chiv8p,
                     'Chiv9': calculate_chiv9p,
                     'Chiv10': calculate_chiv10p,
                     'Chi3c': calculate_chi3c,
                     'Chi4c': calculate_chi4c,
                     'Chi4pc': calculate_chi4pc,
                     'Chi3ch': calculate_chi3ch,
                     'Chi4ch': calculate_chi4ch,
                     'Chi5ch': calculate_chi5ch,
                     'Chi6ch': calculate_chi6ch,
                     'Chi0': calculate_chi0,
                     'Chi1': calculate_chi1,
                     'Chi2': calculate_chi2,
                     'Chi3': calculate_chi3p,
                     'Chi4': calculate_chi4p,
                     'Chi5': calculate_chi5p,
                     'Chi6': calculate_chi6p,
                     'Chi7': calculate_chi7p,
                     'Chi8': calculate_chi8p,
                     'Chi9': calculate_chi9p,
                     'Chi10': calculate_chi10p,
                     'Chiv3c': calculate_chiv3c,
                     'Chiv4c': calculate_chiv4c,
                     'Chiv4pc': calculate_chiv4pc,
                     'Chiv3ch': calculate_chiv3ch,
                     'Chiv4ch': calculate_chiv4ch,
                     'Chiv5ch': calculate_chiv5ch,
                     'Chiv6ch': calculate_chiv6ch,
                     'mChi1': calculate_mean_randic,
                     'knotp': calculate_delta_chi3c4pc,
                     'dchi0': calculate_delta_chi0,
                     'dchi1': calculate_delta_chi1,
                     'dchi2': calculate_delta_chi2,
                     'dchi3': calculate_delta_chi3,
                     'dchi4': calculate_delta_chi4,
                     'knotpv': calculate_delta_chiv3c4pc}

    @staticmethod
    def get_all(mol: Chem.Mol) -> dict:
        """Get all (44) connectivity descriptors."""
        result = {}
        for DesLabel in Connectivity._connectivity.keys():
            result[DesLabel] = Connectivity._connectivity[DesLabel](mol)
        return result
