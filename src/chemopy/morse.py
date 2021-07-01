# -*- coding: utf-8 -*-


"""3D morse descriptors."""


import math
from typing import List

from rdkit import Chem

from chemopy.AtomProperty import GetRelativeAtomicProperty
from chemopy.GeoOpt import _ReadCoordinates
from chemopy.utils import GetGeometricalDistanceMatrix, GetR

# Parameter for RDF equation
_beta = 100


def CalculateUnweightMoRSE(ChargeCoordinates: List[List[float]]) -> dict:
    """Calculate unweighted 3D MoRse descriptors.

    :param ChargeCoordinates: Atomic coordinates and charges as read by chemopy.GeoOpt._ReadCoordinates
    """
    R = GetR(n=30)
    temp = []
#    ChargeCoordinates=_ReadCoordinates('temp.arc')
    for i in ChargeCoordinates:
        # if i[0]!='H':
        temp.append([float(i[1]), float(i[2]), float(i[3])])
    DM = GetGeometricalDistanceMatrix(temp)
    nAT = len(temp)
    RDFresult = {}
    for kkk, Ri in enumerate(R):
        res = 0.0
        for j in range(nAT - 1):
            for k in range(j + 1, nAT):
                res = res + math.sin(Ri * DM[j, k]) / (Ri * DM[j, k])
        RDFresult[f'MoRSEU{kkk + 1}'] = round(res, 3)
    return RDFresult


def CalculateChargeMoRSE(ChargeCoordinates: List[List[float]]) -> dict:
    """Calculate 3D MoRse descriptors from atomic charge.

    :param ChargeCoordinates: Atomic coordinates and charges as read by chemopy.GeoOpt._ReadCoordinates
    """
    R = GetR(n=30)
    temp = []
    charge = []
    for i in ChargeCoordinates:
        # if i[0]!='H':
        temp.append([float(i[1]), float(i[2]), float(i[3])])
        charge.append(float(i[4]))
    DM = GetGeometricalDistanceMatrix(temp)
    nAT = len(temp)
    RDFresult = {}
    for kkk, Ri in enumerate(R):
        res = 0.0
        for j in range(nAT - 1):
            for k in range(j + 1, nAT):
                res = res + charge[j] * charge[k] * math.sin(Ri * DM[j, k]) / (Ri * DM[j, k])
        RDFresult[f'MoRSEC{kkk + 1}'] = round(res, 3)
    return RDFresult


def CalculateMassMoRSE(mol: Chem.Mol, ChargeCoordinates: List[List[float]]) -> dict:
    """Calculate 3D MoRse descriptors from atomic mass.

    :param ChargeCoordinates: Atomic coordinates and charges as read by chemopy.GeoOpt._ReadCoordinates
    """
    R = GetR(n=30)
    temp = []
    mass = [i.atomicmass for i in mol.atoms]
    for i in ChargeCoordinates:
        # if i[0]!='H':
        temp.append([float(i[1]), float(i[2]), float(i[3])])
    DM = GetGeometricalDistanceMatrix(temp)
    nAT = len(temp)
    RDFresult = {}
    for kkk, Ri in enumerate(R):
        res = 0.0
        for j in range(nAT - 1):
            for k in range(j + 1, nAT):
                res = res + mass[j] * mass[k] * math.sin(Ri * DM[j, k]) / (Ri * DM[j, k])
        RDFresult[f'MoRSEM{kkk + 1}'] = round(res / 144, 3)
    return RDFresult


def CalculateAtomicNumberMoRSE(mol: Chem.Mol, ChargeCoordinates: List[List[float]]) -> dict:
    """Calculate 3D MoRse descriptors from atomic number.

    :param ChargeCoordinates: Atomic coordinates and charges as read by chemopy.GeoOpt._ReadCoordinates
    """
    R = GetR(n=30)
    temp = []
    mass = [i.atomicnum for i in mol.atoms]
    for i in ChargeCoordinates:
        # if i[0]!='H':
        temp.append([float(i[1]), float(i[2]), float(i[3])])
    DM = GetGeometricalDistanceMatrix(temp)
    nAT = len(temp)
    RDFresult = {}
    for kkk, Ri in enumerate(R):
        res = 0.0
        for j in range(nAT - 1):
            for k in range(j + 1, nAT):
                res = res + mass[j] * mass[k] * math.sin(Ri * DM[j, k]) / (Ri * DM[j, k])
        RDFresult[f'MoRSEN{kkk + 1}'] = round(res / 144, 3)
    return RDFresult


def CalculatePolarizabilityMoRSE(ChargeCoordinates: List[List[float]]) -> dict:
    """Calculate 3D MoRse descriptors from atomic polarizablity.

    :param ChargeCoordinates: Atomic coordinates and charges as read by chemopy.GeoOpt._ReadCoordinates
    """
    R = GetR(n=30)
    temp = []
    polarizability = []
    for i in ChargeCoordinates:
        # if i[0]!='H':
        temp.append([float(i[1]), float(i[2]), float(i[3])])
        polarizability.append(GetRelativeAtomicProperty(i[0], 'alapha'))
    DM = GetGeometricalDistanceMatrix(temp)
    nAT = len(temp)
    RDFresult = {}
    for kkk, Ri in enumerate(R):
        res = 0.0
        for j in range(nAT - 1):
            for k in range(j + 1, nAT):
                res = res + polarizability[j] * polarizability[k] * math.sin(Ri * DM[j, k]) / (Ri * DM[j, k])
        RDFresult[f'MoRSEP{kkk + 1}'] = round(res, 3)
    return RDFresult


def CalculateSandersonElectronegativityMoRSE(ChargeCoordinates: List[List[float]]) -> dict:
    """Calculate 3D MoRse descriptors from Sanderson electronegativity.

    :param ChargeCoordinates: Atomic coordinates and charges as read by chemopy.GeoOpt._ReadCoordinates
    """
    R = GetR(n=30)
    temp = []
    En = []
    for i in ChargeCoordinates:
        # if i[0]!='H':
        temp.append([float(i[1]), float(i[2]), float(i[3])])
        En.append(GetRelativeAtomicProperty(i[0], 'En'))
    DM = GetGeometricalDistanceMatrix(temp)
    nAT = len(temp)
    RDFresult = {}
    for kkk, Ri in enumerate(R):
        res = 0.0
        for j in range(nAT - 1):
            for k in range(j + 1, nAT):
                res = res + En[j] * En[k] * math.sin(Ri * DM[j, k]) / (Ri * DM[j, k])
        RDFresult[f'MoRSEE{kkk + 1}'] = round(res, 3)
    return RDFresult


def CalculateVDWVolumeMoRSE(ChargeCoordinates: List[List[float]]) -> dict:
    """Calculate 3D MoRse descriptors from van der Waals volume.

    :param ChargeCoordinates: Atomic coordinates and charges as read by chemopy.GeoOpt._ReadCoordinates
    """
    R = GetR(n=30)
    temp = []
    VDW = []
    for i in ChargeCoordinates:
        # if i[0]!='H':
        temp.append([float(i[1]), float(i[2]), float(i[3])])
        VDW.append(GetRelativeAtomicProperty(i[0], 'V'))
    DM = GetGeometricalDistanceMatrix(temp)
    nAT = len(temp)
    RDFresult = {}
    for kkk, Ri in enumerate(R):
        res = 0.0
        for j in range(nAT - 1):
            for k in range(j + 1, nAT):
                res = res + VDW[j] * VDW[k] * math.sin(Ri * DM[j, k]) / (Ri * DM[j, k])
        RDFresult[f'MoRSEV{kkk + 1}'] = round(res, 3)
    return RDFresult


def GetMoRSEUnweighted(arc_file: str) -> dict:
    """Get all unweighted 3D-Morse descriptors.

    :param arc_file: Path to MOPAC .arc file
    """
    ChargeCoordinates = _ReadCoordinates(arc_file)
    result = CalculateUnweightMoRSE(ChargeCoordinates)
    return result


def GetMoRSECharge(arc_file: str) -> dict:
    """Get all 3D-Morse descriptors from charge schemes.

    :param arc_file: Path to MOPAC .arc file
    """
    ChargeCoordinates = _ReadCoordinates(arc_file)
    result = CalculateChargeMoRSE(ChargeCoordinates)
    return result


def GetMoRSEMass(mol: Chem.Mol, arc_file: str) -> dict:
    """Get all 3D-Morse descriptors from on mass schemes.

    :param arc_file: Path to MOPAC .arc file
    """
    ChargeCoordinates = _ReadCoordinates(arc_file)
    result = CalculateMassMoRSE(mol, ChargeCoordinates)
    return result


def GetMoRSEAtomicNumber(mol: Chem.Mol, arc_file: str) -> dict:
    """Get all 3D-Morse descriptors from atomic number schemes.

    :param arc_file: Path to MOPAC .arc file
    """
    ChargeCoordinates = _ReadCoordinates(arc_file)
    result = CalculateAtomicNumberMoRSE(mol, ChargeCoordinates)
    return result


def GetMoRSEPolarizability(arc_file: str) -> dict:
    """Get all 3D-Morse descriptors from polarizability schemes.

    :param arc_file: Path to MOPAC .arc file
    """
    ChargeCoordinates = _ReadCoordinates(arc_file)
    result = CalculatePolarizabilityMoRSE(ChargeCoordinates)
    return result


def GetMoRSESandersonElectronegativity(arc_file: str) -> dict:
    """Get all 3D-Morse descriptors from Sanderson Electronegativity schemes.

    :param arc_file: Path to MOPAC .arc file
    """
    ChargeCoordinates = _ReadCoordinates(arc_file)
    result = CalculateSandersonElectronegativityMoRSE(ChargeCoordinates)
    return result


def GetMoRSEVDWVolume(arc_file: str) -> dict:
    """Get all 3D-Morse descriptors from VDW Volume schemes.

    :param arc_file: Path to MOPAC .arc file
    """
    ChargeCoordinates = _ReadCoordinates(arc_file)
    result = CalculateVDWVolumeMoRSE(ChargeCoordinates)
    return result


def GetMoRSE(mol: Chem.Mol, arc_file: str) -> dict:
    """Get all (210) 3D-Morse descriptors with different (un)weighted schemes.

    :param arc_file: Path to MOPAC .arc file
    """
    result = {}
    ChargeCoordinates = _ReadCoordinates(arc_file)
    result.update(CalculateUnweightMoRSE(ChargeCoordinates))
    result.update(CalculateChargeMoRSE(ChargeCoordinates))
    result.update(CalculateMassMoRSE(mol, ChargeCoordinates))
    result.update(CalculateAtomicNumberMoRSE(mol, ChargeCoordinates))
    result.update(CalculatePolarizabilityMoRSE(ChargeCoordinates))
    result.update(CalculateSandersonElectronegativityMoRSE(ChargeCoordinates))
    result.update(CalculateVDWVolumeMoRSE(ChargeCoordinates))
    return result


# if __name__=="__main__":
#     from GeoOpt import GetARCFile
#     mol='C1C=CCCS1'
#     inputmol=pybel.readstring('smi',mol)
#     dir = GetARCFile(inputmol)
#     #filename='temp'
#     ChargeCoordinates=_ReadCoordinates()
#     print(CalculateVDWVolumeMoRSE(ChargeCoordinates))
#     print(len(GetMoRSE(inputmol, dir)))
#     shutil.rmtree(dir, ignore_errors=True)
