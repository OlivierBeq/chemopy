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
    def calculate_chi0(mol: Chem.Mol, deltas=None) -> float:
        """Calculate molecular connectivity chi index for path order 0."""
        if deltas is None:
            deltas = np.array([x.GetDegree() for x in mol.GetAtoms()], dtype=float)
        deltas = deltas[deltas != 0]
        return float(np.sum(np.sqrt(1.0 / deltas)))

    @staticmethod
    def calculate_chi1(mol: Chem.Mol) -> float:
        """Calculate molecular connectivity chi index for path order 1."""
        cc = [x.GetBeginAtom().GetDegree() * x.GetEndAtom().GetDegree() for x in mol.GetBonds()]
        if len(cc) == 0:
            return 0.0
        while 0 in cc:
            cc.remove(0)
        cc = np.array(cc, "d")
        res = sum(np.sqrt(1.0 / cc))
        return res

    @staticmethod
    def calculate_mean_randic(mol: Chem.Mol) -> float:
        """Calculate mean chi1 (Randic) connectivity index."""
        cc = [x.GetBeginAtom().GetDegree() * x.GetEndAtom().GetDegree() for x in mol.GetBonds()]
        if len(cc) == 0:
            return 0.0
        while 0 in cc:
            cc.remove(0)
        cc = np.array(cc, "d")
        res = np.mean(np.sqrt(1.0 / cc))
        return res

    @staticmethod
    def _calculate_chinp(mol: Chem.Mol, NumPath: int = 2, deltas=None) -> float:
        """Calculate molecular connectivity chi index for path order 2."""
        if deltas is None:
            deltas = np.array([x.GetDegree() for x in mol.GetAtoms()], dtype=float)
        paths = Chem.FindAllPathsOfLengthN(mol, NumPath + 1, useBonds=0)
        if len(paths) == 0:
            return 0.0
        products = deltas[np.array(paths)].prod(axis=1)
        products = products[products != 0]
        return float(np.sum(1.0 / np.sqrt(products)))

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
        patt = Chem.MolFromSmarts("*~*(~*)~*")
        HPatt = mol.GetSubstructMatches(patt)
        for cluster in HPatt:
            deltas = [mol.GetAtomWithIdx(x).GetDegree() for x in cluster]
            while [0] in deltas:
                deltas.remove([0])
            if deltas != []:
                deltas1 = np.array(deltas, dtype=float)
                accum = accum + 1.0 / np.sqrt(deltas1.prod())
        return accum

    @staticmethod
    def calculate_chi4c(mol: Chem.Mol) -> float:
        """Calculate molecular connectivity chi index for cluster."""
        accum = 0.0
        patt = Chem.MolFromSmarts("*~*(~*)(~*)~*")
        HPatt = mol.GetSubstructMatches(patt)
        for cluster in HPatt:
            deltas = [mol.GetAtomWithIdx(x).GetDegree() for x in cluster]
            while [0] in deltas:
                deltas.remove([0])
            if deltas != []:
                deltas1 = np.array(deltas, dtype=float)
                accum = accum + 1.0 / np.sqrt(deltas1.prod())
        return accum

    @staticmethod
    def calculate_chi4pc(mol: Chem.Mol) -> float:
        """Calculate molecular connectivity chi index for path/cluster."""
        accum = 0.0
        patt = Chem.MolFromSmarts("*~*(~*)~*~*")
        HPatt = mol.GetSubstructMatches(patt)
        for cluster in HPatt:
            deltas = [mol.GetAtomWithIdx(x).GetDegree() for x in cluster]
            while [0] in deltas:
                deltas.remove([0])
            if deltas != []:
                deltas1 = np.array(deltas, dtype=float)
                accum = accum + 1.0 / np.sqrt(deltas1.prod())
        return accum

    @staticmethod
    def calculate_delta_chi3c4pc(mol: Chem.Mol) -> float:
        """Calculate the difference between chi3c and chi4pc."""
        return abs(Connectivity.calculate_chi3c(mol) - Connectivity.calculate_chi4pc(mol))

    @staticmethod
    def _calculate_chinch(mol: Chem.Mol, NumCycle=3, deltas=None) -> float:
        """Calculate molecular connectivity chi index for cycles of n."""
        if deltas is None:
            deltas = np.array([x.GetDegree() for x in mol.GetAtoms()], dtype=float)
        rings = [tup for tup in mol.GetRingInfo().AtomRings() if len(tup) == NumCycle]
        if not rings:
            return 0.0
        products = deltas[np.array(rings)].prod(axis=1)
        products = products[products != 0]
        return float(np.sum(1.0 / np.sqrt(products)))

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
    def calculate_chiv0(mol: Chem.Mol, deltas=None) -> float:
        """Calculate valence molecular connectivity chi index for path order 0."""
        if deltas is None:
            deltas = np.array(Topology._hall_kier_deltas(mol, skipHs=0), dtype=float)
        deltas = deltas[deltas != 0]
        return float(np.sum(np.sqrt(1.0 / deltas)))

    @staticmethod
    def _calculate_chivnp(mol: Chem.Mol, NumPath: int = 1, deltas=None) -> float:
        """Calculate valence molecular connectivity chi index for path order 1."""
        if deltas is None:
            deltas = np.array(Topology._hall_kier_deltas(mol, skipHs=False), dtype=float)
        paths = Chem.FindAllPathsOfLengthN(mol, NumPath + 1, useBonds=0)
        if len(paths) == 0:
            return 0.0
        products = deltas[np.array(paths)].prod(axis=1)
        products = products[products != 0]
        return float(np.sum(1.0 / np.sqrt(products)))

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
        patt = Chem.MolFromSmarts("*~*(~*)~*")
        HPatt = mol.GetSubstructMatches(patt)
        for cluster in HPatt:
            deltas = [Connectivity._atom_hall_kier_deltas(mol.GetAtomWithIdx(x)) for x in cluster]
            while [0] in deltas:
                deltas.remove([0])
            if deltas != []:
                deltas1 = np.array(deltas, dtype=float)
                accum = accum + 1.0 / np.sqrt(deltas1.prod())
        return accum

    @staticmethod
    def calculate_chiv4c(mol: Chem.Mol) -> float:
        """Calculate valence molecular connectivity chi index for cluster."""
        accum = 0.0
        patt = Chem.MolFromSmarts("*~*(~*)(~*)~*")
        HPatt = mol.GetSubstructMatches(patt)
        for cluster in HPatt:
            deltas = [Connectivity._atom_hall_kier_deltas(mol.GetAtomWithIdx(x)) for x in cluster]
            while [0] in deltas:
                deltas.remove([0])
            if deltas != []:
                deltas1 = np.array(deltas, dtype=float)
                accum = accum + 1.0 / np.sqrt(deltas1.prod())
        return accum

    @staticmethod
    def calculate_chiv4pc(mol: Chem.Mol) -> float:
        """Calculate valence molecular connectivity chi index for path/cluster."""
        accum = 0.0
        patt = Chem.MolFromSmarts("*~*(~*)~*~*")
        HPatt = mol.GetSubstructMatches(patt)
        for cluster in HPatt:
            deltas = [Connectivity._atom_hall_kier_deltas(mol.GetAtomWithIdx(x)) for x in cluster]
            while [0] in deltas:
                deltas.remove([0])
            if deltas != []:
                deltas1 = np.array(deltas, dtype=float)
                accum = accum + 1.0 / np.sqrt(deltas1.prod())
        return accum

    @staticmethod
    def calculate_delta_chiv3c4pc(mol: Chem.Mol) -> float:
        """Calculate the difference between chiv3c and chiv4pc."""
        return abs(Connectivity.calculate_chiv3c(mol) - Connectivity.calculate_chiv4pc(mol))

    @staticmethod
    def _calculate_chivnch(mol: Chem.Mol, NumCyc=3, deltas=None) -> float:
        """Calculate valence molecular connectivity chi index for cycles of n."""
        if deltas is None:
            deltas = np.array(Topology._hall_kier_deltas(mol, skipHs=0), dtype=float)
        rings = [tup for tup in mol.GetRingInfo().AtomRings() if len(tup) == NumCyc]
        if not rings:
            return 0.0
        products = deltas[np.array(rings)].prod(axis=1)
        products = products[products != 0]
        return float(np.sum(1.0 / np.sqrt(products)))

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

    @staticmethod
    def get_all(mol: Chem.Mol) -> dict:
        """Get all (44) connectivity descriptors.

        Several descriptors below (e.g. `dchi2`/`dchi3`/`dchi4`, `knotp`, `knotpv`) are, by
        definition, differences of other descriptors already computed here; and the plain/valence
        Hall-Kier deltas they all depend on are computed once and reused, instead of every single
        descriptor recomputing them (and refinding the same atom paths/rings) independently.
        """
        plain_deltas = np.array([x.GetDegree() for x in mol.GetAtoms()], dtype=float)
        valence_deltas = np.array(Topology._hall_kier_deltas(mol, skipHs=0), dtype=float)

        result = {}
        result["Chiv0"] = Connectivity.calculate_chiv0(mol, deltas=valence_deltas)
        chivnp = {n: Connectivity._calculate_chivnp(mol, NumPath=n, deltas=valence_deltas) for n in range(1, 11)}
        for n in range(1, 11):
            result[f"Chiv{n}"] = chivnp[n]
        result["Chi3c"] = Connectivity.calculate_chi3c(mol)
        result["Chi4c"] = Connectivity.calculate_chi4c(mol)
        result["Chi4pc"] = Connectivity.calculate_chi4pc(mol)
        result["Chi3ch"] = Connectivity._calculate_chinch(mol, NumCycle=3, deltas=plain_deltas)
        result["Chi4ch"] = Connectivity._calculate_chinch(mol, NumCycle=4, deltas=plain_deltas)
        result["Chi5ch"] = Connectivity._calculate_chinch(mol, NumCycle=5, deltas=plain_deltas)
        result["Chi6ch"] = Connectivity._calculate_chinch(mol, NumCycle=6, deltas=plain_deltas)
        result["Chi0"] = Connectivity.calculate_chi0(mol, deltas=plain_deltas)
        result["Chi1"] = Connectivity.calculate_chi1(mol)
        chinp = {n: Connectivity._calculate_chinp(mol, NumPath=n, deltas=plain_deltas) for n in range(2, 11)}
        for n in range(2, 11):
            result[f"Chi{n}"] = chinp[n]
        result["Chiv3c"] = Connectivity.calculate_chiv3c(mol)
        result["Chiv4c"] = Connectivity.calculate_chiv4c(mol)
        result["Chiv4pc"] = Connectivity.calculate_chiv4pc(mol)
        result["Chiv3ch"] = Connectivity._calculate_chivnch(mol, NumCyc=3, deltas=valence_deltas)
        result["Chiv4ch"] = Connectivity._calculate_chivnch(mol, NumCyc=4, deltas=valence_deltas)
        result["Chiv5ch"] = Connectivity._calculate_chivnch(mol, NumCyc=5, deltas=valence_deltas)
        result["Chiv6ch"] = Connectivity._calculate_chivnch(mol, NumCyc=6, deltas=valence_deltas)
        result["mChi1"] = Connectivity.calculate_mean_randic(mol)
        result["knotp"] = abs(result["Chi3c"] - result["Chi4pc"])
        result["dchi0"] = abs(result["Chiv0"] - result["Chi0"])
        result["dchi1"] = abs(result["Chiv1"] - result["Chi1"])
        result["dchi2"] = abs(chivnp[2] - chinp[2])
        result["dchi3"] = abs(chivnp[3] - chinp[3])
        result["dchi4"] = abs(chivnp[4] - chinp[4])
        result["knotpv"] = abs(result["Chiv3c"] - result["Chiv4pc"])
        return result
