# -*- coding: utf-8 -*-


"""Basak information topological indices."""


import copy

import numpy as np
from rdkit import Chem

from .topology import Topology


class Basak:
    """Basak information content descriptors."""
    
    @staticmethod
    def calculate_basak_ic0(mol: Chem.Mol)-> float:
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
    def calculate_basak_sic0(mol: Chem.Mol)-> float:
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
    def calculate_basak_cic0(mol: Chem.Mol)-> float:
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
    def _calculate_basak_ic_n(mol: Chem.Mol, NumPath=1)-> float:
        """Calculate the information content of order n."""
        Hmol = Chem.AddHs(mol)
        nAtoms = Hmol.GetNumAtoms()
        TotalPath = Chem.FindAllPathsOfLengthN(Hmol, NumPath, useBonds=0, useHs=1)
        if len(TotalPath) == 0:
            BasakIC = 0.0
        else:
            IC = {}
            for i in range(nAtoms):
                temp = []
                at = Hmol.GetAtomWithIdx(i)
                temp.append(at.GetAtomicNum())
                for index in TotalPath:
                    if i == index[0]:
                        temp.append([Hmol.GetAtomWithIdx(kk).GetAtomicNum() for kk in index[1:]])
                    if i == index[-1]:
                        cds = list(index)
                        cds.reverse()
                        temp.append([Hmol.GetAtomWithIdx(kk).GetAtomicNum() for kk in cds[1:]])
                IC[str(i)] = temp
            cds = []
            for value in IC.values():
                cds.append(value)
            kkk = list(range(len(cds)))
            aaa = copy.deepcopy(kkk)
            res = []
            for i in aaa:
                if i in kkk:
                    jishu = 0
                    kong = []
                    temp1 = cds[i]
                    for j in aaa:
                        if cds[j] == temp1:
                            jishu += 1
                            kong.append(j)
                    for ks in kong:
                        kkk.remove(ks)
                    res.append(jishu)
            BasakIC = Topology._calculate_entropy(np.array(res, dtype=float) / sum(res))
        return BasakIC

    @staticmethod
    def calculate_basak_ic1(mol: Chem.Mol)-> float:
        """Calculate the information content of order 1."""
        return Basak._calculate_basak_ic_n(mol, NumPath=2)

    @staticmethod
    def calculate_basak_ic2(mol: Chem.Mol)-> float:
        """Calculate the information content of order 2."""
        return Basak._calculate_basak_ic_n(mol, NumPath=3)

    @staticmethod
    def calculate_basak_ic3(mol: Chem.Mol)-> float:
        """Calculate the information content of order 3."""
        return Basak._calculate_basak_ic_n(mol, NumPath=4)

    @staticmethod
    def calculate_basak_ic4(mol: Chem.Mol)-> float:
        """Calculate the information content of order 4."""
        return Basak._calculate_basak_ic_n(mol, NumPath=5)

    @staticmethod
    def calculate_basak_ic5(mol: Chem.Mol)-> float:
        """Calculate the information content of order 5."""
        return Basak._calculate_basak_ic_n(mol, NumPath=6)

    @staticmethod
    def calculate_basak_ic6(mol: Chem.Mol)-> float:
        """Calculate the information content of order 6."""
        return Basak._calculate_basak_ic_n(mol, NumPath=7)

    @staticmethod
    def calculate_basak_sic1(mol: Chem.Mol)-> float:
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
    def calculate_basak_sic2(mol: Chem.Mol)-> float:
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
    def calculate_basak_sic3(mol: Chem.Mol)-> float:
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
    def calculate_basak_sic4(mol: Chem.Mol)-> float:
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
    def calculate_basak_sic5(mol: Chem.Mol)-> float:
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
    def calculate_basak_sic6(mol: Chem.Mol)-> float:
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
    def calculate_basak_cic1(mol: Chem.Mol)-> float:
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
    def calculate_basak_cic2(mol: Chem.Mol)-> float:
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
    def calculate_basak_cic3(mol: Chem.Mol)-> float:
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
    def calculate_basak_cic4(mol: Chem.Mol)-> float:
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
    def calculate_basak_cic5(mol: Chem.Mol)-> float:
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
    def calculate_basak_cic6(mol: Chem.Mol)-> float:
        """Calculate the complementary information content of order 6."""
        Hmol = Chem.AddHs(mol)
        nAtoms = Hmol.GetNumAtoms()
        IC = Basak.calculate_basak_ic6(mol)
        if nAtoms <= 1:
            BasakCIC = 0.0
        else:
            BasakCIC = np.log2(nAtoms) - IC
        return BasakCIC
    
    
    _basak = {'IC0': calculate_basak_ic0,
              'IC1': calculate_basak_ic1,
              'IC2': calculate_basak_ic2,
              'IC3': calculate_basak_ic3,
              'IC4': calculate_basak_ic4,
              'IC5': calculate_basak_ic5,
              'IC6': calculate_basak_ic6,
              'SIC0': calculate_basak_sic0,
              'SIC1': calculate_basak_sic1,
              'SIC2': calculate_basak_sic2,
              'SIC3': calculate_basak_sic3,
              'SIC4': calculate_basak_sic4,
              'SIC5': calculate_basak_sic5,
              'SIC6': calculate_basak_sic6,
              'CIC0': calculate_basak_cic0,
              'CIC1': calculate_basak_cic1,
              'CIC2': calculate_basak_cic2,
              'CIC3': calculate_basak_cic3,
              'CIC4': calculate_basak_cic4,
              'CIC5': calculate_basak_cic5,
              'CIC6': calculate_basak_cic6,
              }

    @staticmethod
    def get_all(mol: Chem.Mol) -> dict:
        """Calculate all (21) Basak descriptors."""
        result = {}
        for DesLabel in Basak._basak.keys():
            result[DesLabel] = Basak._basak[DesLabel](mol)
        return result
