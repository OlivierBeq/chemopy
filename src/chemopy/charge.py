# -*- coding: utf-8 -*-


"""Charge descriptors based on Gasteiger/Marseli partial charges."""

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdPartialCharges as GMCharge


iter_step = 12


class Charge:
    """Charge descriptors."""
    
    @staticmethod
    def _calculate_element_max_pcharge(mol: Chem.Mol, AtomicNum: int = 6) -> float:
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
            return max(res)

    @staticmethod
    def _calculate_element_max_ncharge(mol: Chem.Mol, AtomicNum: int = 6) -> float:
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
            return min(res)

    @staticmethod
    def calculate_h_max_pcharge(mol: Chem.Mol) -> float:
        """Get most positive charge of all hydrogen atoms."""
        return Charge._calculate_element_max_pcharge(mol, AtomicNum=1)

    @staticmethod
    def calculate_c_max_pcharge(mol: Chem.Mol) -> float:
        """Get most positive charge of all carbon atoms."""
        return Charge._calculate_element_max_pcharge(mol, AtomicNum=6)

    @staticmethod
    def calculate_n_max_pcharge(mol: Chem.Mol) -> float:
        """Get most positive charge of all nitrogen atoms."""
        return Charge._calculate_element_max_pcharge(mol, AtomicNum=7)

    @staticmethod
    def calculate_o_max_pcharge(mol: Chem.Mol) -> float:
        """Get most positive charge of all oxygen atoms."""
        return Charge._calculate_element_max_pcharge(mol, AtomicNum=8)

    @staticmethod
    def calculate_h_max_ncharge(mol) -> float:
        """Get most negative charge of all hydrogen atoms."""
        return Charge._calculate_element_max_ncharge(mol, AtomicNum=1)

    @staticmethod
    def calculate_c_max_ncharge(mol: Chem.Mol) -> float:
        """Get most negative charge of all carbon atoms."""
        return Charge._calculate_element_max_ncharge(mol, AtomicNum=6)

    @staticmethod
    def calculate_n_max_ncharge(mol: Chem.Mol) -> float:
        """Get most negative charge of all nitrogen atoms."""
        return Charge._calculate_element_max_ncharge(mol, AtomicNum=7)

    @staticmethod
    def calculate_o_max_ncharge(mol: Chem.Mol) -> float:
        """Get most negative charge of all oxygen atoms."""
        return Charge._calculate_element_max_ncharge(mol, AtomicNum=8)

    @staticmethod
    def calculate_all_max_pcharge(mol: Chem.Mol) -> float:
        """Get most positive charge of all atoms."""
        Hmol = Chem.AddHs(mol)
        GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
        res = []
        for atom in Hmol.GetAtoms():
            res.append(float(atom.GetProp('_GasteigerCharge')))
        if res == []:
            return 0
        else:
            return max(res)

    @staticmethod
    def calculate_all_max_ncharge(mol: Chem.Mol) -> float:
        """Get most negative charge of all atoms."""
        Hmol = Chem.AddHs(mol)
        GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
        res = []
        for atom in Hmol.GetAtoms():
            res.append(float(atom.GetProp('_GasteigerCharge')))
        if res == []:
            return 0
        else:
            return min(res)

    @staticmethod
    def _calculate_element_sum_square_charge(mol: Chem.Mol, AtomicNum: int = 6) -> float:
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
            return sum(np.square(res))

    @staticmethod
    def calculate_h_sum_square_charge(mol: Chem.Mol) -> float:
        """Get the sum of square charges of hydrogen atoms."""
        return Charge._calculate_element_sum_square_charge(mol, AtomicNum=1)

    @staticmethod
    def calculate_c_sum_square_charge(mol: Chem.Mol) -> float:
        """Get the sum of square charges of carbon atoms."""
        return Charge._calculate_element_sum_square_charge(mol, AtomicNum=6)

    @staticmethod
    def calculate_n_sum_square_charge(mol: Chem.Mol) -> float:
        """Get the sum of square charges of nitrogen atoms."""
        return Charge._calculate_element_sum_square_charge(mol, AtomicNum=7)

    @staticmethod
    def calculate_o_sum_square_charge(mol: Chem.Mol) -> float:
        """Get the sum of square charges of oxygen atoms."""
        return Charge._calculate_element_sum_square_charge(mol, AtomicNum=8)

    @staticmethod
    def calculate_all_sum_square_charge(mol: Chem.Mol)-> float:
        """Get the sum of square charges of all atoms."""
        Hmol = Chem.AddHs(mol)
        GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
        res = []
        for atom in Hmol.GetAtoms():
            res.append(float(atom.GetProp('_GasteigerCharge')))
        if res == []:
            return 0
        else:
            return sum(np.square(res))

    @staticmethod
    def calculate_total_pcharge(mol: Chem.Mol) -> float:
        """Get the total positive charge."""
        Hmol = Chem.AddHs(mol)
        GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
        res = []
        for atom in Hmol.GetAtoms():
            res.append(float(atom.GetProp('_GasteigerCharge')))
        if res == []:
            return 0
        else:
            cc = np.array(res, 'd')
            return sum(cc[cc > 0])

    @staticmethod
    def calculate_mean_pcharge(mol: Chem.Mol) -> float:
        """Get the average positive charge."""
        Hmol = Chem.AddHs(mol)
        GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
        res = []
        for atom in Hmol.GetAtoms():
            res.append(float(atom.GetProp('_GasteigerCharge')))
        if res == []:
            return 0
        else:
            cc = np.array(res, 'd')
            cc = cc[cc > 0]
            return np.mean(cc) if len(cc) else np.nan

    @staticmethod
    def calculate_total_ncharge(mol: Chem.Mol) -> float:
        """Ge the total negative charge."""
        Hmol = Chem.AddHs(mol)
        GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
        res = []
        for atom in Hmol.GetAtoms():
            res.append(float(atom.GetProp('_GasteigerCharge')))
        if res == []:
            return 0
        else:
            cc = np.array(res, 'd')
            cc = cc[cc < 0]
            return sum(cc) if len(cc) else np.nan

    @staticmethod
    def calculate_mean_ncharge(mol: Chem.Mol) -> float:
        """Get the average negative charge."""
        Hmol = Chem.AddHs(mol)
        GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
        res = []
        for atom in Hmol.GetAtoms():
            res.append(float(atom.GetProp('_GasteigerCharge')))
        if res == []:
            return 0
        else:
            cc = np.array(res, 'd')
            cc = cc[cc < 0]
            return np.mean(cc) if len(cc) else np.nan

    @staticmethod
    def calculate_total_absolute_charge(mol: Chem.Mol) -> float:
        """Get the total absolute charge."""
        Hmol = Chem.AddHs(mol)
        GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
        res = []
        for atom in Hmol.GetAtoms():
            res.append(float(atom.GetProp('_GasteigerCharge')))
        if res == []:
            return 0
        else:
            cc = np.array(res, 'd')
            return sum(np.absolute(cc))

    @staticmethod
    def calculate_mean_absolute_charge(mol: Chem.Mol) -> float:
        """Get the average absolute charge."""
        Hmol = Chem.AddHs(mol)
        GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
        res = []
        for atom in Hmol.GetAtoms():
            res.append(float(atom.GetProp('_GasteigerCharge')))
        if res == []:
            return 0
        else:
            cc = np.array(res, 'd')
            return np.mean(np.absolute(cc))

    @staticmethod
    def calculate_relative_pcharge(mol: Chem.Mol) -> float:
        """Get the ratio between the most positive partial charge and the total positive charge."""
        Hmol = Chem.AddHs(mol)
        GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
        res = []
        for atom in Hmol.GetAtoms():
            res.append(float(atom.GetProp('_GasteigerCharge')))
        if res == []:
            return 0
        else:
            cc = np.array(res, 'd')
            if sum(cc[cc > 0]) == 0:
                return 0
            else:
                return max(res) / sum(cc[cc > 0])

    @staticmethod
    def calculate_relative_ncharge(mol: Chem.Mol) -> float:
        """Get the ratio between the most negative partial charge and the total negative charge."""
        Hmol = Chem.AddHs(mol)
        GMCharge.ComputeGasteigerCharges(Hmol, iter_step)
        res = []
        for atom in Hmol.GetAtoms():
            res.append(float(atom.GetProp('_GasteigerCharge')))
        if res == []:
            return 0
        else:
            cc = np.array(res, 'd')
            if sum(cc[cc < 0]) == 0:
                return 0
            else:
                return min(res) / sum(cc[cc < 0])

    @staticmethod
    def calculate_local_dipole_index(mol: Chem.Mol) -> float:
        """calculate_ the local dipole index (D)."""
        GMCharge.ComputeGasteigerCharges(mol, iter_step)
        res = []
        for atom in mol.GetAtoms():
            res.append(float(atom.GetProp('_GasteigerCharge')))
        cc = [np.absolute(res[x.GetBeginAtom().GetIdx()] - res[x.GetEndAtom().GetIdx()]) for x in mol.GetBonds()]
        B = len(mol.GetBonds())
        return 0 if len(cc) == 0.0 else sum(cc) / B

    @staticmethod
    def calculate_submol_polarity_param(mol: Chem.Mol) -> float:
        """calculate_ the submolecular polarity parameter (SPP)."""
        return Charge.calculate_all_max_pcharge(mol) - Charge.calculate_all_max_ncharge(mol)
    
    @staticmethod
    def get_all(mol: Chem.Mol) -> dict:
        """Get all (25) constitutional descriptors."""
        result = {}
        for DesLabel in _charge.keys():
            result[DesLabel] = _charge[DesLabel](mol)
        return result


_charge = {'QHmax': Charge.calculate_h_max_pcharge,
           'QCmax': Charge.calculate_c_max_pcharge,
           'QNmax': Charge.calculate_n_max_pcharge,
           'QOmax': Charge.calculate_o_max_pcharge,
           'QHmin': Charge.calculate_h_max_ncharge,
           'QCmin': Charge.calculate_c_max_ncharge,
           'QNmin': Charge.calculate_n_max_ncharge,
           'QOmin': Charge.calculate_o_max_ncharge,
           'Qmax': Charge.calculate_all_max_pcharge,
           'Qmin': Charge.calculate_all_max_ncharge,
           'QHss': Charge.calculate_h_sum_square_charge,
           'QCss': Charge.calculate_c_sum_square_charge,
           'QNss': Charge.calculate_n_sum_square_charge,
           'QOss': Charge.calculate_o_sum_square_charge,
           'Qass': Charge.calculate_all_sum_square_charge,
           'Mpc': Charge.calculate_mean_pcharge,
           'Tpc': Charge.calculate_total_pcharge,
           'Mnc': Charge.calculate_mean_ncharge,
           'Tnc': Charge.calculate_total_ncharge,
           'Mac': Charge.calculate_mean_absolute_charge,
           'Tac': Charge.calculate_total_absolute_charge,
           'Rnc': Charge.calculate_relative_ncharge,
           'Rpc': Charge.calculate_relative_pcharge,
           'SPP': Charge.calculate_submol_polarity_param,
           'LDI': Charge.calculate_local_dipole_index,
           }
