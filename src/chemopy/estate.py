# -*- coding: utf-8 -*-


"""Estate fingerprints and values.

Based on Hall, Money, Kier, J. Chem. Inf. Comput. Sci. (1991)
doi: 10.1021/ci00001a012
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem.EState import Fingerprinter as ESFP

from . import atom_types as ATEstate


class EState:
    """Electrotopological descriptors and fingerprints."""

    @staticmethod
    def _calculate_estate(mol: Chem.Mol, skipH: bool = True) -> np.array:
        """Get the EState value of each atom in the molecule."""
        mol = Chem.AddHs(mol)
        if skipH:
            mol = Chem.RemoveHs(mol)
        tb1 = Chem.GetPeriodicTable()
        nAtoms = mol.GetNumAtoms()
        Is = np.zeros(nAtoms, dtype=float)
        for i in range(nAtoms):
            at = mol.GetAtomWithIdx(i)
            atNum = at.GetAtomicNum()
            d = at.GetDegree()
            if d > 0:
                h = at.GetTotalNumHs()
                dv = tb1.GetNOuterElecs(atNum) - h
                N = EState._get_principle_quantum_number(atNum)
                Is[i] = (4.0 / (N * N) * dv + 1) / d
        dists = Chem.GetDistanceMatrix(mol, useBO=0, useAtomWts=0)
        dists = dists + 1
        # Vectorized equivalent of the reference i<j double loop: the contribution matrix
        # C[i, j] = (Is[i]-Is[j]) / dists[i, j]^2 is antisymmetric (C[j, i] == -C[i, j]),
        # so summing each full row reproduces the same accum[i] the pairwise loop built.
        mask = dists < 1e6
        np.fill_diagonal(mask, False)
        diff = Is[:, None] - Is[None, :]
        contributions = np.where(mask, diff / (dists * dists), 0.0)
        accum = contributions.sum(axis=1)
        res = accum + Is
        return res

    @staticmethod
    def _get_principle_quantum_number(atNum: int) -> int:
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

    @staticmethod
    def calculate_heavy_atom_estate(mol: Chem.Mol) -> float:
        """Calculate sum of the EState indices over all heavy atoms."""
        return sum(EState._calculate_estate(mol))

    @staticmethod
    def _calculate_atom_estate(mol: Chem.Mol, atomic_num=6) -> float:
        """Calculate the sum of the EState indices over all atoms with specified atomic number."""
        nAtoms = mol.GetNumAtoms()
        Is = np.zeros(nAtoms, dtype=float)
        Estate = EState._calculate_estate(mol, skipH=False)
        for i in range(nAtoms):
            at = mol.GetAtomWithIdx(i)
            atNum = at.GetAtomicNum()
            if atNum == atomic_num:
                Is[i] = Estate[i]
        res = sum(Is)
        return res

    @staticmethod
    def calculate_c_atom_estate(mol: Chem.Mol) -> float:
        """Calculate the sum of the EState indices over all C atoms."""
        return EState._calculate_atom_estate(mol, atomic_num=6)

    @staticmethod
    def calculate_halogen_estate(mol: Chem.Mol) -> float:
        """Calculate the sum of the EState indices over all Halogen atoms."""
        Nf = EState._calculate_atom_estate(mol, atomic_num=9)
        Ncl = EState._calculate_atom_estate(mol, atomic_num=17)
        Nbr = EState._calculate_atom_estate(mol, atomic_num=35)
        Ni = EState._calculate_atom_estate(mol, atomic_num=53)
        return Nf + Ncl + Nbr + Ni

    @staticmethod
    def calculate_hetero_estate(mol: Chem.Mol) -> float:
        """Calculate the sum of the EState indices over all hetero atoms."""
        Ntotal = sum(EState._calculate_estate(mol))
        NC = EState._calculate_atom_estate(mol, atomic_num=6)
        NH = EState._calculate_atom_estate(mol, atomic_num=1)
        return Ntotal - NC - NH

    @staticmethod
    def calculate_average_estate(mol: Chem.Mol) -> float:
        """Calculate the ratio of the sum of the EState indices over heavy atoms and the number of non-hydrogen atoms."""
        N = mol.GetNumAtoms()
        return sum(EState._calculate_estate(mol)) / N

    @staticmethod
    def calculate_max_estate(mol: Chem.Mol) -> float:
        """Calculate the maximal Estate value in all atoms."""
        return max(EState._calculate_estate(mol))

    @staticmethod
    def calculate_min_estate(mol: Chem.Mol) -> float:
        """Calculate the minimal Estate value in all atoms."""
        return min(EState._calculate_estate(mol))

    @staticmethod
    def calculate_diff_max_min_estate(mol: Chem.Mol) -> float:
        """Calculate the difference between Smax and Smin."""
        return max(EState._calculate_estate(mol)) - min(EState._calculate_estate(mol))

    @staticmethod
    def calculate_estate_fingerprint(mol: Chem.Mol, implementation="rdkit", binary: bool = False) -> dict:
        """Calculate the sum of EState values for each EState atom type.

        :param implementation: either rdkit or chemopy. chemopy rounds
                               to the third decimal place but not rdkit.
        :param binary: should bineray values be returned instead of continous ones
        """
        if implementation not in ["rdkit", "chemopy"]:
            raise ValueError("Implementation of AtomTypeEState must be either rdkit or chemopy.")
        if implementation == "chemopy":
            AT = ATEstate.get_atom_label(mol)
            Estate = EState._calculate_estate(mol)
            res = []
            for i in AT:
                if i == []:
                    res.append(0)
                else:
                    res.append(sum(Estate[k] for k in i))
            ESresult = {f"S{n + 1}": es for n, es in enumerate(res)}
            if binary:
                ESresult = {f"EStateFP_{i + 1}": 0 if x == 0 else 1 for i, x in enumerate(ESresult.values())}
            return ESresult
        else:
            temp = ESFP.FingerprintMol(mol)
            if binary:
                res = {f"EStateFP_{i + 1}": 0 if x == 0 else 1 for i, x in enumerate(temp[0])}
            else:
                res = {f"S{i + 1}": j for i, j in enumerate(temp[1])}
            return res

    @staticmethod
    def calculate_atom_type_estate_fingerprint(mol: Chem.Mol) -> dict:
        """Calculate EState Fingerprints.

        This is the counts of each EState atom type in the molecule.
        """
        temp = ESFP.FingerprintMol(mol)
        res = {}
        for i, j in enumerate(temp[0]):
            res[f"Sfinger{i + 1}"] = j
        return res

    @staticmethod
    def calculate_max_atom_type_estate(mol: Chem.Mol) -> dict:
        """Calculate the maximum of EState value."""
        AT = ATEstate.get_atom_label(mol)
        Estate = EState._calculate_estate(mol)
        res = []
        for i in AT:
            if i == []:
                res.append(0)
            else:
                res.append(max(Estate[k] for k in i))
        ESresult = {}
        for n, es in enumerate(res):
            ESresult[f"Smax{n + 1}"] = es
        return ESresult

    @staticmethod
    def calculate_min_atom_type_estate(mol: Chem.Mol) -> dict:
        """Calculate the minimum of EState value."""
        AT = ATEstate.get_atom_label(mol)
        Estate = EState._calculate_estate(mol)
        res = []
        for i in AT:
            if i == []:
                res.append(0)
            else:
                res.append(min(Estate[k] for k in i))
        ESresult = {}
        for n, es in enumerate(res):
            ESresult[f"Smin{n + 1}"] = es
        return ESresult

    @staticmethod
    def get_all_descriptors(mol: Chem.Mol) -> dict:
        """Calculate all (8) EState descriptors.

        The reference implementation recomputes `_calculate_estate` (an O(nAtoms^2)
        distance-weighted accumulation) independently for each of the 8 descriptors below
        (~15 calls total, since `calculate_hetero_estate`/`calculate_halogen_estate` each
        call it multiple times themselves); here it is computed at most twice (once per
        `skipH` variant) and reused, exactly reproducing each descriptor's original formula.
        """
        estate_no_h = EState._calculate_estate(mol, skipH=True)
        estate_with_h = EState._calculate_estate(mol, skipH=False)
        atomic_nums = np.array([atom.GetAtomicNum() for atom in Chem.AddHs(mol).GetAtoms()])

        def atom_type_sum(atomic_num):
            mask = atomic_nums == atomic_num
            return float(estate_with_h[mask].sum()) if mask.any() else 0.0

        Ntotal = float(estate_no_h.sum())
        NC = atom_type_sum(6)
        NH = atom_type_sum(1)
        Smax = float(estate_no_h.max())
        Smin = float(estate_no_h.min())
        result = {
            "Shev": Ntotal,
            "Scar": NC,
            "Shal": atom_type_sum(9) + atom_type_sum(17) + atom_type_sum(35) + atom_type_sum(53),
            "Shet": Ntotal - NC - NH,
            "Save": Ntotal / mol.GetNumAtoms(),
            "Smax": Smax,
            "Smin": Smin,
            "DS": Smax - Smin,
        }
        return result

    @staticmethod
    def get_all_fps(mol: Chem.Mol) -> dict:
        """Calculate all (316) EState descriptors."""
        # calculate_max_atom_type_estate/calculate_min_atom_type_estate each independently
        # recompute get_atom_label()/_calculate_estate(); compute both once here instead.
        AT = ATEstate.get_atom_label(mol)
        Estate = EState._calculate_estate(mol)
        max_result = {}
        min_result = {}
        for n, i in enumerate(AT):
            if i == []:
                max_result[f"Smax{n + 1}"] = 0
                min_result[f"Smin{n + 1}"] = 0
            else:
                values = [Estate[k] for k in i]
                max_result[f"Smax{n + 1}"] = max(values)
                min_result[f"Smin{n + 1}"] = min(values)

        result = {}
        result.update(EState.calculate_estate_fingerprint(mol))
        result.update(max_result)
        result.update(min_result)
        return result
