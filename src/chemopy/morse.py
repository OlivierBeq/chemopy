# -*- coding: utf-8 -*-


"""3D morse descriptors."""


import math
from typing import List

import scipy
from rdkit import Chem

from pychem.AtomProperty import GetRelativeAtomicProperty
from pychem.GeoOpt import _ReadCoordinates

_beta = 100

def GetR(n: int) -> List[float]:
    """Calcuate the parameters R of the RDF equation."""
    R = []
    for i in range(2, n + 2):
        R.append(float(i * 0.5))
    return R


def _GetAtomDistance(x: List[float], y: List[float]) -> float:
    """Calculate Euclidean distance between two atomic coordinates."""
    temp = [math.pow(x[0] - y[0], 2), math.pow(x[1] - y[1], 2), math.pow(x[2] - y[2], 2)]
    res = math.sqrt(sum(temp))
    return res
    

def _GetGementricalDistanceMatrix(CoordinateList: List[List[float]]) -> scipy.matrix:
    """Calculate distance matrix from coordinate list."""
    NAtom = len(CoordinateList)
    DistanceMatrix = scipy.zeros((NAtom, NAtom))
    for i in range(NAtom - 1):
        for j in range(i + 1, NAtom):
            DistanceMatrix[i, j] = _GetAtomDistance(CoordinateList[i], CoordinateList[j])
            DistanceMatrix[j, i] = DistanceMatrix[i, j]
    return DistanceMatrix

    
def CalculateUnweightMoRSE(ChargeCoordinates: List[List[float]]) -> dict:
    """Calculate unweighted 3D MoRse descriptors.

    :param ChargeCoordinates: Atomic coordinates and charges as read by pychem.GeoOpt._ReadCoordinates
    """
    R = GetR(n=30)
    temp = []
#    ChargeCoordinates=_ReadCoordinates('temp.arc')
    for i in ChargeCoordinates:
        # if i[0]!='H':
        temp.append([float(i[1]), float(i[2]), float(i[3])])
    DM = _GetGementricalDistanceMatrix(temp)
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

    :param ChargeCoordinates: Atomic coordinates and charges as read by pychem.GeoOpt._ReadCoordinates
    """
    R = GetR(n=30)
    temp = []
    charge = []
    for i in ChargeCoordinates:
        # if i[0]!='H':
        temp.append([float(i[1]), float(i[2]), float(i[3])])
        charge.append(float(i[4]))
    DM = _GetGementricalDistanceMatrix(temp)
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

    :param ChargeCoordinates: Atomic coordinates and charges as read by pychem.GeoOpt._ReadCoordinates
    """
    R = GetR(n=30)
    temp = []
    mass = [i.atomicmass for i in mol.atoms]
    for i in ChargeCoordinates:
        # if i[0]!='H':
        temp.append([float(i[1]), float(i[2]), float(i[3])])
    DM = _GetGementricalDistanceMatrix(temp)
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

    :param ChargeCoordinates: Atomic coordinates and charges as read by pychem.GeoOpt._ReadCoordinates
    """
    R = GetR(n=30)
    temp = []
    mass = [i.atomicnum for i in mol.atoms]
    for i in ChargeCoordinates:
        # if i[0]!='H':
        temp.append([float(i[1]), float(i[2]), float(i[3])])
    DM = _GetGementricalDistanceMatrix(temp)
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

    :param ChargeCoordinates: Atomic coordinates and charges as read by pychem.GeoOpt._ReadCoordinates
    """
    R = GetR(n=30)
    temp = []
    polarizability = []
    for i in ChargeCoordinates:
        # if i[0]!='H':
        temp.append([float(i[1]), float(i[2]), float(i[3])])
        polarizability.append(GetRelativeAtomicProperty(i[0], 'alapha'))
    DM = _GetGementricalDistanceMatrix(temp)
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

    :param ChargeCoordinates: Atomic coordinates and charges as read by pychem.GeoOpt._ReadCoordinates
    """
    R = GetR(n=30)
    temp = []
    En = []
    for i in ChargeCoordinates:
        # if i[0]!='H':
        temp.append([float(i[1]), float(i[2]), float(i[3])])
        En.append(GetRelativeAtomicProperty(i[0], 'En'))
    DM = _GetGementricalDistanceMatrix(temp)
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

    :param ChargeCoordinates: Atomic coordinates and charges as read by pychem.GeoOpt._ReadCoordinates
    """
    R = GetR(n=30)
    temp = []
    VDW = []
    for i in ChargeCoordinates:
        # if i[0]!='H':
        temp.append([float(i[1]), float(i[2]), float(i[3])])
        VDW.append(GetRelativeAtomicProperty(i[0], 'V'))
    DM = _GetGementricalDistanceMatrix(temp)
    nAT = len(temp)
    RDFresult = {}
    for kkk, Ri in enumerate(R):
        res = 0.0
        for j in range(nAT - 1):
            for k in range(j + 1, nAT):
                res = res + VDW[j] * VDW[k] * math.sin(Ri * DM[j, k]) / (Ri * DM[j, k])
        RDFresult[f'MoRSEV{kkk + 1}'] = round(res, 3)
    return RDFresult



def GetMoRSEUnweighted(mol: Chem.Mol) -> dict:
    """Get all unweighted 3D-Morse descriptors."""
    filename = 'temp'
    ChargeCoordinates = _ReadCoordinates(filename) 
    result = CalculateUnweightMoRSE(ChargeCoordinates)
    return result


def GetMoRSECharge(mol: Chem.Mol) -> dict:
    """Get all 3D-Morse descriptors from charge schemes."""
    filename = 'temp'
    ChargeCoordinates = _ReadCoordinates(filename) 
    result = CalculateChargeMoRSE(ChargeCoordinates)
    return result


def GetMoRSEMass(mol: Chem.Mol) -> dict:
    """Get all 3D-Morse descriptors from on mass schemes."""
    filename = 'temp'
    ChargeCoordinates = _ReadCoordinates(filename) 
    result=CalculateMassMoRSE(mol, ChargeCoordinates)
    return result


def GetMoRSEAtomicNumber(mol: Chem.Mol) -> dict:
    """Get all 3D-Morse descriptors from atomic number schemes."""
    filename = 'temp'
    ChargeCoordinates = _ReadCoordinates(filename) 
    result = CalculateAtomicNumberMoRSE(mol, ChargeCoordinates)
    return result


def GetMoRSEPolarizability(mol: Chem.Mol) -> dict:
    """Get all 3D-Morse descriptors from polarizability schemes."""
    filename = 'temp'
    ChargeCoordinates = _ReadCoordinates(filename) 
    result = CalculatePolarizabilityMoRSE(ChargeCoordinates)
    return result


def GetMoRSESandersonElectronegativity(mol: Chem.Mol) -> dict:
    """Get all 3D-Morse descriptors from Sanderson Electronegativity schemes."""
    filename = 'temp'
    ChargeCoordinates = _ReadCoordinates(filename) 
    result = CalculateSandersonElectronegativityMoRSE(ChargeCoordinates)
    return result


def GetMoRSEVDWVolume(mol: Chem.Mol) -> dict:
    """Get all 3D-Morse descriptors from VDW Volume schemes."""
    filename = 'temp'
    ChargeCoordinates = _ReadCoordinates(filename) 
    result=CalculateVDWVolumeMoRSE(ChargeCoordinates)
    return result


def GetMoRSE(mol: Chem.Mol) -> dict:
    """Get all (210) 3D-Morse descriptors with different (un)weighted schemes."""
    result = {}
    filename = 'temp'
    ChargeCoordinates = _ReadCoordinates(filename)
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
