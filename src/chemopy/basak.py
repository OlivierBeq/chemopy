# -*- coding: utf-8 -*-


"""Basak information topological indices."""

from collections import Counter

import numpy as np
from rdkit import Chem

from .topology import Topology


class Basak:
    """Basak information content descriptors."""

    @staticmethod
    def calculate_basak_ic0(mol: Chem.Mol) -> float:
        """Calculate information content of order 0."""
        Hmol = Chem.AddHs(mol)
        nAtoms = Hmol.GetNumAtoms()
        IC = []
        for i in range(nAtoms):
            at = Hmol.GetAtomWithIdx(i)
            IC.append(at.GetAtomicNum())
        Unique = np.unique(IC)
        NAtomType = len(Unique)
        NTAtomType = np.zeros(NAtomType, dtype=float)
        for i in range(NAtomType):
            NTAtomType[i] = IC.count(Unique[i])
        if nAtoms != 0:
            BasakIC = Topology._calculate_entropy(NTAtomType / nAtoms)
        else:
            BasakIC = 0.0
        return BasakIC

    @staticmethod
    def calculate_basak_sic0(mol: Chem.Mol) -> float:
        """Calculate the structural information content of order 0."""
        Hmol = Chem.AddHs(mol)
        nAtoms = Hmol.GetNumAtoms()
        IC = Basak.calculate_basak_ic0(mol)
        if nAtoms <= 1:
            BasakSIC = 0.0
        else:
            BasakSIC = IC / np.log2(nAtoms)
        return BasakSIC

    @staticmethod
    def calculate_basak_cic0(mol: Chem.Mol) -> float:
        """Calculate the complementary information content of order 0."""
        Hmol = Chem.AddHs(mol)
        nAtoms = Hmol.GetNumAtoms()
        IC = Basak.calculate_basak_ic0(mol)
        if nAtoms <= 1:
            BasakCIC = 0.0
        else:
            BasakCIC = np.log2(nAtoms) - IC
        return BasakCIC

    @staticmethod
    def _calculate_basak_ic_n(mol: Chem.Mol, NumPath=1) -> float:
        """Calculate the information content of order n."""
        Hmol = Chem.AddHs(mol)
        nAtoms = Hmol.GetNumAtoms()
        TotalPath = Chem.FindAllPathsOfLengthN(Hmol, NumPath, useBonds=0, useHs=1)
        if len(TotalPath) == 0:
            return 0.0
        # Atomic numbers are static per molecule: look them up once instead of
        # repeating GetAtomWithIdx()/GetAtomicNum() calls inside the path loop below.
        atomic_nums = [atom.GetAtomicNum() for atom in Hmol.GetAtoms()]
        # Group each path's atomic-number sequence by the atom(s) it starts/ends at,
        # in a single pass over TotalPath, instead of rescanning all paths per atom.
        per_atom_paths = [[] for _ in range(nAtoms)]
        for index in TotalPath:
            start, end = index[0], index[-1]
            per_atom_paths[start].append(tuple(atomic_nums[kk] for kk in index[1:]))
            if end == start:
                # Preserve the reference semantics: a path that starts and ends at the
                # same atom contributes twice (two independent `if` checks, not `elif`).
                per_atom_paths[end].append(tuple(atomic_nums[kk] for kk in reversed(index[:-1])))
            else:
                per_atom_paths[end].append(tuple(atomic_nums[kk] for kk in reversed(index[:-1])))
        # Build each atom's canonical (hashable) environment signature and count
        # duplicates directly, instead of an O(nAtoms^2) manual list-equality scan.
        signatures = ((atomic_nums[i], *per_atom_paths[i]) for i in range(nAtoms))
        res = list(Counter(signatures).values())
        return Topology._calculate_entropy(np.array(res, dtype=float) / sum(res))

    @staticmethod
    def calculate_basak_ic1(mol: Chem.Mol) -> float:
        """Calculate the information content of order 1."""
        return Basak._calculate_basak_ic_n(mol, NumPath=2)

    @staticmethod
    def calculate_basak_ic2(mol: Chem.Mol) -> float:
        """Calculate the information content of order 2."""
        return Basak._calculate_basak_ic_n(mol, NumPath=3)

    @staticmethod
    def calculate_basak_ic3(mol: Chem.Mol) -> float:
        """Calculate the information content of order 3."""
        return Basak._calculate_basak_ic_n(mol, NumPath=4)

    @staticmethod
    def calculate_basak_ic4(mol: Chem.Mol) -> float:
        """Calculate the information content of order 4."""
        return Basak._calculate_basak_ic_n(mol, NumPath=5)

    @staticmethod
    def calculate_basak_ic5(mol: Chem.Mol) -> float:
        """Calculate the information content of order 5."""
        return Basak._calculate_basak_ic_n(mol, NumPath=6)

    @staticmethod
    def calculate_basak_ic6(mol: Chem.Mol) -> float:
        """Calculate the information content of order 6."""
        return Basak._calculate_basak_ic_n(mol, NumPath=7)

    @staticmethod
    def calculate_basak_sic1(mol: Chem.Mol) -> float:
        """Calculate the structural information content of order 1."""
        Hmol = Chem.AddHs(mol)
        nAtoms = Hmol.GetNumAtoms()
        IC = Basak.calculate_basak_ic1(mol)
        if nAtoms <= 1:
            BasakSIC = 0.0
        else:
            BasakSIC = IC / np.log2(nAtoms)
        return BasakSIC

    @staticmethod
    def calculate_basak_sic2(mol: Chem.Mol) -> float:
        """Calculate the structural information content of order 2."""
        Hmol = Chem.AddHs(mol)
        nAtoms = Hmol.GetNumAtoms()
        IC = Basak.calculate_basak_ic2(mol)
        if nAtoms <= 1:
            BasakSIC = 0.0
        else:
            BasakSIC = IC / np.log2(nAtoms)
        return BasakSIC

    @staticmethod
    def calculate_basak_sic3(mol: Chem.Mol) -> float:
        """Calculate the structural information content of order 3."""
        Hmol = Chem.AddHs(mol)
        nAtoms = Hmol.GetNumAtoms()
        IC = Basak.calculate_basak_ic3(mol)
        if nAtoms <= 1:
            BasakSIC = 0.0
        else:
            BasakSIC = IC / np.log2(nAtoms)
        return BasakSIC

    @staticmethod
    def calculate_basak_sic4(mol: Chem.Mol) -> float:
        """Calculate the structural information content of order 4."""
        Hmol = Chem.AddHs(mol)
        nAtoms = Hmol.GetNumAtoms()
        IC = Basak.calculate_basak_ic4(mol)
        if nAtoms <= 1:
            BasakSIC = 0.0
        else:
            BasakSIC = IC / np.log2(nAtoms)
        return BasakSIC

    @staticmethod
    def calculate_basak_sic5(mol: Chem.Mol) -> float:
        """Calculate the structural information content of order 5."""
        Hmol = Chem.AddHs(mol)
        nAtoms = Hmol.GetNumAtoms()
        IC = Basak.calculate_basak_ic5(mol)
        if nAtoms <= 1:
            BasakSIC = 0.0
        else:
            BasakSIC = IC / np.log2(nAtoms)
        return BasakSIC

    @staticmethod
    def calculate_basak_sic6(mol: Chem.Mol) -> float:
        """Calculate the structural information content of order 6."""
        Hmol = Chem.AddHs(mol)
        nAtoms = Hmol.GetNumAtoms()
        IC = Basak.calculate_basak_ic6(mol)
        if nAtoms <= 1:
            BasakSIC = 0.0
        else:
            BasakSIC = IC / np.log2(nAtoms)
        return BasakSIC

    @staticmethod
    def calculate_basak_cic1(mol: Chem.Mol) -> float:
        """Calculate the complementary information content of order 1."""
        Hmol = Chem.AddHs(mol)
        nAtoms = Hmol.GetNumAtoms()
        IC = Basak.calculate_basak_ic1(mol)
        if nAtoms <= 1:
            BasakCIC = 0.0
        else:
            BasakCIC = np.log2(nAtoms) - IC
        return BasakCIC

    @staticmethod
    def calculate_basak_cic2(mol: Chem.Mol) -> float:
        """Calculate the complementary information content of order 2."""
        Hmol = Chem.AddHs(mol)
        nAtoms = Hmol.GetNumAtoms()
        IC = Basak.calculate_basak_ic2(mol)
        if nAtoms <= 1:
            BasakCIC = 0.0
        else:
            BasakCIC = np.log2(nAtoms) - IC
        return BasakCIC

    @staticmethod
    def calculate_basak_cic3(mol: Chem.Mol) -> float:
        """Calculate the complementary information content of order 3."""
        Hmol = Chem.AddHs(mol)
        nAtoms = Hmol.GetNumAtoms()
        IC = Basak.calculate_basak_ic3(mol)
        if nAtoms <= 1:
            BasakCIC = 0.0
        else:
            BasakCIC = np.log2(nAtoms) - IC
        return BasakCIC

    @staticmethod
    def calculate_basak_cic4(mol: Chem.Mol) -> float:
        """Calculate the complementary information content of order 4."""
        Hmol = Chem.AddHs(mol)
        nAtoms = Hmol.GetNumAtoms()
        IC = Basak.calculate_basak_ic4(mol)
        if nAtoms <= 1:
            BasakCIC = 0.0
        else:
            BasakCIC = np.log2(nAtoms) - IC
        return BasakCIC

    @staticmethod
    def calculate_basak_cic5(mol: Chem.Mol) -> float:
        """Calculate the complementary information content of order 5."""
        Hmol = Chem.AddHs(mol)
        nAtoms = Hmol.GetNumAtoms()
        IC = Basak.calculate_basak_ic5(mol)
        if nAtoms <= 1:
            BasakCIC = 0.0
        else:
            BasakCIC = np.log2(nAtoms) - IC
        return BasakCIC

    @staticmethod
    def calculate_basak_cic6(mol: Chem.Mol) -> float:
        """Calculate the complementary information content of order 6."""
        Hmol = Chem.AddHs(mol)
        nAtoms = Hmol.GetNumAtoms()
        IC = Basak.calculate_basak_ic6(mol)
        if nAtoms <= 1:
            BasakCIC = 0.0
        else:
            BasakCIC = np.log2(nAtoms) - IC
        return BasakCIC

    @staticmethod
    def get_all(mol: Chem.Mol) -> dict:
        """Calculate all (21) Basak descriptors."""
        # calculate_basak_{s,c}ic{1..6} each independently recompute _calculate_basak_ic_n;
        # compute each order once here and derive SIC/CIC from it instead.
        nAtoms = Chem.AddHs(mol).GetNumAtoms()
        log2n = np.log2(nAtoms) if nAtoms > 1 else None
        result = {}
        for order in range(7):
            ic = Basak.calculate_basak_ic0(mol) if order == 0 else Basak._calculate_basak_ic_n(mol, NumPath=order + 1)
            result[f"IC{order}"] = ic
            result[f"SIC{order}"] = ic / log2n if log2n is not None else 0.0
            result[f"CIC{order}"] = log2n - ic if log2n is not None else 0.0
        return result
