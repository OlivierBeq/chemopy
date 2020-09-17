# -*- coding: utf-8 -*-


"""MOE-type descriptors.

Includes LabuteASA, TPSA, slogPVSA, MRVSA,
         PEOEVSA, EstateVSA and VSAEstate.
"""

from typing import List

from rdkit import Chem
from rdkit.Chem import MolSurf as MOE
from rdkit.Chem.EState import EState_VSA as EVSA


def CalculateLabuteASA(mol: Chem.Mol) -> dict:
    """Calculate Labute's Approximate Surface Area (ASA from MOE)."""
    res = {}
    temp = MOE.pyLabuteASA(mol, includeHs=1)
    res['LabuteASA'] = round(temp, 3)
    return res


def CalculateTPSA(mol: Chem.Mol) -> dict:
    """Calculate topological polar surface area based on fragments.

    Implementation based on the Daylight contrib program TPSA.
    """
    res = {}
    temp = MOE.TPSA(mol)
    res['TPSA1'] = round(temp, 3)
    return res


def CalculateSLOGPVSA(mol: Chem.Mol, bins: List[float] = None) -> dict:
    """Get MOE-type descriptors using LogP and SA contributions.

    :param bins: interval boundaries used in the P_VSA calculation.
                 The default SLOGP bins are [-0.4,-0.2,0,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6].
    """
    temp = MOE.SlogP_VSA_(mol, bins, force=1)
    res = {}
    for i, j in enumerate(temp):
        res[f'slogPVSA{i}'] = round(j, 3)
    return res


def CalculateSMRVSA(mol: Chem.Mol, bins: List[float] = None) -> dict:
    """Get MOE-type descriptors using MR and SA contributions.

    :param bins: interval boundaries used in the P_VSA calculation.
                 The default SMR bins are [1.29, 1.82, 2.24, 2.45, 2.75, 3.05, 3.63,3.8,4.0].
    """
    temp = MOE.SMR_VSA_(mol, bins, force=1)
    res = {}
    for i, j in enumerate(temp):
        res[f'MRVSA{i}'] = round(j, 3)
    return res


def CalculatePEOEVSA(mol: Chem.Mol, bins: List[float] = None) -> dict:
    """Get MOE-type descriptors using partial charges and SA contributions.

    :param bins: interval boundaries used in the P_VSA calculation.
                 The default PEOE bins are [-.3,-.25,-.20,-.15,-.10,-.05,0,.05,.10,.15,.20,.25,.30].
    """
    temp = MOE.PEOE_VSA_(mol, bins, force=1)
    res = {}
    for i, j in enumerate(temp):
        res[f'PEOEVSA{i}'] = round(j, 3)
    return res


def CalculateEstateVSA(mol: Chem.Mol, bins: List[float] = None) -> dict:
    """Get MOE-type descriptors using Estate indices and SA contributions.

    :param bins: interval boundaries used in the P_VSA calculation.
                 The default Estate bins are [-0.390,0.290,0.717,1.165,1.540,1.807,2.05,4.69,9.17,15.0].
    """
    temp = EVSA.EState_VSA_(mol, bins, force=1)
    res = {}
    for i, j in enumerate(temp):
        res[f'EstateVSA{i}'] = round(j, 3)
    return res


def CalculateVSAEstate(mol: Chem.Mol, bins: List[float] = None) -> dict:
    """Get MOE-type descriptors using SA and Estate indices contributions.

    :param bins: interval boundaries used in the P_VSA calculation.
                 The default VSA bins are [4.78,5.00,5.410,5.740,6.00,6.07,6.45,7.00,11.0].
    """
    temp = EVSA.VSA_EState_(mol, bins, force=1)
    res = {}
    for i, j in enumerate(temp):
        res[f'VSAEstate{i}'] = round(j, 3)
    return res


def GetMOE(mol: Chem.Mol) -> dict:
    """Calculate all (59) MOE-type descriptors."""
    result = {}
    result.update(CalculateLabuteASA(mol))
    result.update(CalculateTPSA(mol))
    result.update(CalculateSLOGPVSA(mol, bins=None))  # 12 values
    result.update(CalculateSMRVSA(mol, bins=None))  # 10 values
    result.update(CalculatePEOEVSA(mol, bins=None))  # 14 values
    result.update(CalculateEstateVSA(mol, bins=None))  # 11 values
    result.update(CalculateVSAEstate(mol, bins=None))  # 10 values
    return result


#########################################################################
# if __name__=="__main__":
#     smi5=['COCCCC','CCC(C)CC','CC(C)CCC','CC(C)C(C)C','CCOCCN','c1ccccc1N']
#     smis = ['CCCC','CCCCC','CCCCCC','CC(N)C(=O)O','CC(N)C(=O)[O-].[Na+]']
#     for index, smi in enumerate(smis):
#         m = Chem.MolFromSmiles(smi)
#         print(index+1)
#         print(smi)
#         print('\t',GetMOE(m))
#         print('\t', len(GetMOE(m)))
