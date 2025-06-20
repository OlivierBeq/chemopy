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
        dists += 1
        accum = np.zeros(nAtoms, dtype=float)
        for i in range(nAtoms):
            for j in range(i + 1, nAtoms):
                p = dists[i, j]
                if p < 1e6:
                    temp = (Is[i] - Is[j]) / (p * p)
                    accum[i] += temp
                    accum[j] -= temp
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
    def calculate_estate_fingerprint(mol: Chem.Mol, implementation='rdkit', binary: bool = False) -> dict:
        """Calculate the sum of EState values for each EState atom type.

        :param implementation: either rdkit or chemopy. chemopy rounds
                               to the third decimal place but not rdkit.
        :param binary: should bineray values be returned instead of continous ones
        """
        if implementation not in ['rdkit', 'chemopy']:
            raise ValueError('Implementation of AtomTypeEState must be either rdkit or chemopy.')
        if implementation == 'chemopy':
            AT = ATEstate.get_atom_label(mol)
            Estate = EState._calculate_estate(mol)
            res = []
            for i in AT:
                if i == []:
                    res.append(0)
                else:
                    res.append(sum(Estate[k] for k in i))
            ESresult = {f'S{n + 1}': es for n, es in enumerate(res)}
            if binary:
                ESresult = {f'EStateFP_{i + 1}' : 0 if x == 0 else 1 for i, x in enumerate(ESresult.values())}
            return ESresult
        else:
            temp = ESFP.FingerprintMol(mol)
            if binary:
                res = {f'EStateFP_{i + 1}' : 0 if x == 0 else 1 for i, x in enumerate(temp[0])}
            else:
                res = {f'S{i + 1}': j for i, j in enumerate(temp[1])}
            return res

    @staticmethod
    def calculate_atom_type_estate_fingerprint(mol: Chem.Mol) -> dict:
        """Calculate EState Fingerprints.

        This is the counts of each EState atom type in the molecule.
        """
        temp = ESFP.FingerprintMol(mol)
        res = {}
        for i, j in enumerate(temp[0]):
            res[f'Sfinger{i + 1}'] = j
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
            ESresult[f'Smax{n + 1}'] = es
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
            ESresult[f'Smin{n + 1}'] = es
        return ESresult

    @staticmethod
    def get_all_descriptors(mol: Chem.Mol) -> dict:
        """Calculate all (8) EState descriptors."""
        result = {}
        result.update({'Shev': EState.calculate_heavy_atom_estate(mol)})
        result.update({'Scar': EState.calculate_c_atom_estate(mol)})
        result.update({'Shal': EState.calculate_halogen_estate(mol)})
        result.update({'Shet': EState.calculate_hetero_estate(mol)})
        result.update({'Save': EState.calculate_average_estate(mol)})
        result.update({'Smax': EState.calculate_max_estate(mol)})
        result.update({'Smin': EState.calculate_min_estate(mol)})
        result.update({'DS': EState.calculate_diff_max_min_estate(mol)})
        return result

    @staticmethod
    def get_all_fps(mol: Chem.Mol) -> dict:
        """Calculate all (316) EState descriptors."""
        result = {}
        result.update(EState.calculate_estate_fingerprint(mol))
        result.update(EState.calculate_max_atom_type_estate(mol))
        result.update(EState.calculate_min_atom_type_estate(mol))
        return result
