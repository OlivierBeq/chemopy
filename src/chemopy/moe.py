# -*- coding: utf-8 -*-


"""MOE-type descriptors.

Includes LabuteASA, TPSA, slogPVSA, MRVSA,
         PEOEVSA, EstateVSA and VSAEstate.
"""

from typing import List

from rdkit import Chem
from rdkit.Chem import MolSurf
from rdkit.Chem.EState import EState_VSA as EVSA


class MOE:
    """MOE-type 2D descriptors."""

    @staticmethod
    def calculate_labute_asa(mol: Chem.Mol) -> dict:
        """Calculate Labute's Approximate Surface Area (ASA from MOE)."""
        asa = MolSurf.pyLabuteASA(mol, includeHs=1)
        return asa

    @staticmethod
    def calculate_tpsa(mol: Chem.Mol) -> dict:
        """Calculate topological polar surface area based on fragments.

        Implementation based on the Daylight contrib program TPSA.
        """
        tpsa = MolSurf.TPSA(mol)
        return tpsa

    @staticmethod
    def calculate_slogpvsa(mol: Chem.Mol, bins: List[float] = None) -> dict:
        """Get MOE-type descriptors using LogP and SA contributions.

        :param bins: interval boundaries used in the P_VSA calculation.
                     The default SLOGP bins are [-0.4,-0.2,0,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6].
        """
        temp = MolSurf.SlogP_VSA_(mol, bins, force=1)
        res = {}
        for i, j in enumerate(temp):
            res[f'slogPVSA{i}'] = j
        return res

    @staticmethod
    def calculate_smrvsa(mol: Chem.Mol, bins: List[float] = None) -> dict:
        """Get MOE-type descriptors using MR and SA contributions.

        :param bins: interval boundaries used in the P_VSA calculation.
                     The default SMR bins are [1.29, 1.82, 2.24, 2.45, 2.75, 3.05, 3.63,3.8,4.0].
        """
        temp = MolSurf.SMR_VSA_(mol, bins, force=1)
        res = {}
        for i, j in enumerate(temp):
            res[f'SMRVSA{i}'] = j
        return res

    @staticmethod
    def calculate_peoevsa(mol: Chem.Mol, bins: List[float] = None) -> dict:
        """Get MOE-type descriptors using partial charges and SA contributions.

        :param bins: interval boundaries used in the P_VSA calculation.
                     The default PEOE bins are [-.3,-.25,-.20,-.15,-.10,-.05,0,.05,.10,.15,.20,.25,.30].
        """
        temp = MolSurf.PEOE_VSA_(mol, bins, force=1)
        res = {}
        for i, j in enumerate(temp):
            res[f'PEOEVSA{i}'] = j
        return res

    @staticmethod
    def calculate_estate_vsa(mol: Chem.Mol, bins: List[float] = None) -> dict:
        """Get MOE-type descriptors using Estate indices and SA contributions.

        :param bins: interval boundaries used in the P_VSA calculation.
                     The default Estate bins are [-0.390,0.290,0.717,1.165,1.540,1.807,2.05,4.69,9.17,15.0].
        """
        temp = EVSA.EState_VSA_(mol, bins, force=1)
        res = {}
        for i, j in enumerate(temp):
            res[f'EstateVSA{i}'] = j
        return res

    @staticmethod
    def calculate_vsa_estate(mol: Chem.Mol, bins: List[float] = None) -> dict:
        """Get MOE-type descriptors using SA and Estate indices contributions.

        :param bins: interval boundaries used in the P_VSA calculation.
                     The default VSA bins are [4.78,5.00,5.410,5.740,6.00,6.07,6.45,7.00,11.0].
        """
        temp = EVSA.VSA_EState_(mol, bins, force=1)
        res = {}
        for i, j in enumerate(temp):
            res[f'VSAEstate{i}'] = j
        return res

    @staticmethod
    def get_all(mol: Chem.Mol) -> dict:
        """Calculate all (59) MOE-type descriptors."""
        result = {}
        result['TPSA1'] = MOE.calculate_tpsa(mol)
        result['LabuteASA'] = MOE.calculate_labute_asa(mol)
        result.update(MOE.calculate_slogpvsa(mol, bins=None))  # 12 values
        result.update(MOE.calculate_smrvsa(mol, bins=None))  # 10 values
        result.update(MOE.calculate_peoevsa(mol, bins=None))  # 14 values
        result.update(MOE.calculate_estate_vsa(mol, bins=None))  # 11 values
        result.update(MOE.calculate_vsa_estate(mol, bins=None))  # 10 values
        return result
