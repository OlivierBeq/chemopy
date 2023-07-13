
# -*- coding: utf-8 -*-


"""Kier and Hall's kappa indices."""


from rdkit import Chem
from rdkit.Chem import rdchem

periodicTable = rdchem.GetPeriodicTable()

class Kappa:
    """Kier and Hall's kappa indices."""

    @staticmethod
    def calculate_kappa1(mol: Chem.Mol) -> float:
        """Calculate molecular shape index for one bonded fragment."""
        P1 = mol.GetNumBonds(onlyHeavy=1)
        A = mol.GetNumHeavyAtoms()
        denom = P1 + 0.0
        if denom:
            kappa = (A) * (A - 1) ** 2 / denom ** 2
        else:
            kappa = 0.0
        return kappa

    @staticmethod
    def calculate_kappa2(mol: Chem.Mol) -> float:
        """Calculate molecular shape index for two bonded fragments."""
        P2 = len(Chem.FindAllPathsOfLengthN(mol, 2))
        A = mol.GetNumHeavyAtoms()
        denom = P2 + 0.0
        if denom:
            kappa = (A - 1) * (A - 2) ** 2 / denom ** 2
        else:
            kappa = 0.0
        return kappa

    @staticmethod
    def calculate_kappa3(mol: Chem.Mol) -> float:
        """Calculate molecular shape index for three bonded fragments."""
        P3 = len(Chem.FindAllPathsOfLengthN(mol, 3))
        A = mol.GetNumHeavyAtoms()
        denom = P3 + 0.0
        if denom:
            if A % 2 == 1:
                kappa = (A - 1) * (A - 3) ** 2 / denom ** 2
            else:
                kappa = (A - 3) * (A - 2) ** 2 / denom ** 2
        else:
            kappa = 0.0
        return kappa

    hall_kier_alphas = {'Br': [None, None, 0.48],
                      'C': [-0.22, -0.13, 0.0],
                      'Cl': [None, None, 0.29],
                      'F': [None, None, -0.07],
                      'H': [0.0, 0.0, 0.0],
                      'I': [None, None, 0.73],
                      'N': [-0.29, -0.2, -0.04],
                      'O': [None, -0.2, -0.04],
                      'P': [None, 0.3, 0.43],
                      'S': [None, 0.22, 0.35]}

    @staticmethod
    def _hall_kier_alpha(mol: Chem.Mol) -> float:
        """Calculate Hall-Kier alpha value for a molecule."""
        alphaSum = 0.0
        rC = periodicTable.GetRb0(6)
        for atom in mol.GetAtoms():
            atNum = atom.GetAtomicNum()
            if not atNum:
                continue
            symb = atom.GetSymbol()
            alphaV = Kappa.hall_kier_alphas.get(symb, None)
            if alphaV is not None:
                hyb = atom.GetHybridization() - 2
                if hyb < len(alphaV):
                    alpha = alphaV[hyb]
                    if alpha is None:
                        alpha = alphaV[-1]
                else:
                    alpha = alphaV[-1]
            else:
                rA = periodicTable.GetRb0(atNum)
                alpha = rA / rC - 1
            alphaSum += alpha
        return alphaSum

    @staticmethod
    def calculate_kappa_alapha1(mol: Chem.Mol) -> float:
        """Calculate molecular shape index for one bonded fragment."""
        P1 = mol.GetNumBonds(onlyHeavy=1)
        A = mol.GetNumHeavyAtoms()
        alpha = Kappa._hall_kier_alpha(mol)
        denom = P1 + alpha
        if denom:
            kappa = (A + alpha) * (A + alpha - 1) ** 2 / denom ** 2
        else:
            kappa = 0.0
        return kappa

    @staticmethod
    def calculate_kappa_alapha2(mol: Chem.Mol) -> float:
        """Calculate molecular shape index for two bonded fragments."""
        P2 = len(Chem.FindAllPathsOfLengthN(mol, 2))
        A = mol.GetNumHeavyAtoms()
        alpha = Kappa._hall_kier_alpha(mol)
        denom = P2 + alpha
        if denom:
            kappa = (A + alpha - 1) * (A + alpha - 2) ** 2 / denom ** 2
        else:
            kappa = 0.0
        return kappa

    @staticmethod
    def calculate_kappa_alapha3(mol: Chem.Mol) -> float:
        """Calculate molecular shape index for three bonded fragments."""
        P3 = len(Chem.FindAllPathsOfLengthN(mol, 3))
        A = mol.GetNumHeavyAtoms()
        alpha = Kappa._hall_kier_alpha(mol)
        denom = P3 + alpha
        if denom:
            if A % 2 == 1:
                kappa = (A + alpha - 1) * (A + alpha - 3) ** 2 / denom ** 2
            else:
                kappa = (A + alpha - 3) * (A + alpha - 2) ** 2 / denom ** 2
        else:
            kappa = 0.0
        return kappa

    @staticmethod
    def calculate_flexibility(mol: Chem.Mol) -> float:
        """Calculate Kier molecular flexibility index."""
        kappa1 = Kappa.calculate_kappa_alapha1(mol)
        kappa2 = Kappa.calculate_kappa_alapha2(mol)
        A = mol.GetNumHeavyAtoms()
        phi = kappa1 * kappa2 / (A + 0.0)
        return phi

    @staticmethod
    def get_all(mol: Chem.Mol) -> dict:
        """Calculate all (7) kappa values."""
        res = {}
        res['kappam1'] = Kappa.calculate_kappa_alapha1(mol)
        res['kappam2'] = Kappa.calculate_kappa_alapha2(mol)
        res['kappam3'] = Kappa.calculate_kappa_alapha3(mol)
        res['phi'] = Kappa.calculate_flexibility(mol)
        res['kappa1'] = Kappa.calculate_kappa1(mol)
        res['kappa2'] = Kappa.calculate_kappa2(mol)
        res['kappa3'] = Kappa.calculate_kappa3(mol)
        return res
