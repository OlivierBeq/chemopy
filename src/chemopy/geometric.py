# -*- coding: utf-8 -*-


"""Geometrical molecular descriptors.

Optimizes molecular structure using MOPAC.
"""

import math
from typing import List, Tuple

import scipy
from openbabel import pybel

from pychem.GeoOpt import _ReadCoordinates


def GetAtomDistance(x: List[float], y: List[float]) -> float:
    """Calculate Euclidean distance between two atomic coordinates."""
    temp = [math.pow(x[0] - y[0], 2), math.pow(x[1] - y[1], 2), math.pow(x[2] - y[2], 2)]
    res = math.sqrt(sum(temp))
    return res


def GetGementricalDistanceMatrix(CoordinateList: List[List[float]]) -> scipy.matrix:
    """Calculate distance matrix from coordinate list."""
    NAtom = len(CoordinateList)
    DistanceMatrix = scipy.zeros((NAtom, NAtom))
    for i in range(NAtom - 1):
        for j in range(i + 1, NAtom):
            DistanceMatrix[i, j] = GetAtomDistance(CoordinateList[i], CoordinateList[j])
            DistanceMatrix[j, i] = DistanceMatrix[i, j]
    return DistanceMatrix


           
def _GetMassCenter(MassCoordinates: List[Tuple[float, Tuple[float, float, float]]]) -> Tuple[float, float, float]:
    """Get the center of mass.

    :param MassCoordinates: list of atomic masses and coordinates
                            in the format [(atommass, (x,y,z)), ...]
    """
    res1 = 0.0
    res2 = 0.0
    res3 = 0.0
    temp = []
    for i in MassCoordinates:
        res1 += i[0] * i[1][0]
        res2 += i[0] * i[1][1]
        res3 += i[0] * i[1][2]
        temp.append(i[0])
    result = (res1 / sum(temp), res2 / sum(temp), res3 / sum(temp))
    return result


def _GetGeometricalCenter(ChargeCoordinates: List[List[float]]) -> Tuple[float, float, float]:
    """Get the geometrical center.

    :param ChargeCoordinates: Atomic coordinates and charges as read by pychem.GeoOpt._ReadCoordinates
    """
    res1 = []
    res2 = []
    res3 = []
    for i in ChargeCoordinates:
        res1.append(float(i[1]))
        res2.append(float(i[2]))
        res3.append(float(i[3]))
    result = (scipy.mean(res1), scipy.mean(res2), scipy.mean(res3))
    return result


def Calculate3DWienerWithH(ChargeCoordinates: List[List[float]]) -> float:
    """Calculate 3D Wiener index from geometrical distance matrix of a MOPAC optimized molecule (including Hs).

    :param ChargeCoordinates: Atomic coordinates and charges as read by pychem.GeoOpt._ReadCoordinates
    """
    temp = []
    for i in ChargeCoordinates:
        temp.append([float(i[1]), float(i[2]), float(i[3])])
    DistanceMatrix = GetGeometricalDistanceMatrix(temp)
    return round(scipy.sum(DistanceMatrix) / 2.0, 3)


def Calculate3DWienerWithoutH(ChargeCoordinates: List[List[float]]) -> float:
    """Calculate 3D Wiener index from geometrical distance matrix of a MOPAC optimized molecule (not including Hs).

    :param ChargeCoordinates: Atomic coordinates and charges as read by pychem.GeoOpt._ReadCoordinates
    """
    temp = []
    for i in ChargeCoordinates:
        if i[0] != 'H':
            temp.append([float(i[1]), float(i[2]), float(i[3])])
    DistanceMatrix = GetGeometricalDistanceMatrix(temp)
    return round(scipy.sum(DistanceMatrix) / 2.0, 3)


def CalculatePetitjean3DIndex(ChargeCoordinates: List[List[float]]) -> float:
    """Calculate Petitjean Index from molecular gemetrical distance matrix.

    The 3D Petitjean shape index (PJI3) is calculated
    dividing the difference between geometric diameter and
    radius by the geometric radius [P.A. Bath, A.R. Poirrette,
    P. Willett, F.H. Allen, J.Chem.Inf.Comput.Sci. 1995, 35, 714-716].
    The geometric radius of a molecule is defined as the minimum
    geometric eccentricity and the diameter is defined as the
    maximum geometric eccentricity in the molecule, the atom
    geometric eccentricity being the longest geometric distance
    from the considered atom to any other atom in the molecule.

    :param ChargeCoordinates: Atomic coordinates and charges as read by pychem.GeoOpt._ReadCoordinates
    """
    temp = []
    for i in ChargeCoordinates:
        temp.append([float(i[1]), float(i[2]), float(i[3])])
    DistanceMatrix = GetGeometricalDistanceMatrix(temp)
    temp1 = scipy.amax(DistanceMatrix, axis=0)
    return round(max(temp1) / min(temp1) - 1.0, 3)


def CalculateGeometricalDiameter(ChargeCoordinates: List[List[float]]) -> float:
    """Calculate the geometrical diameter.

    :param ChargeCoordinates: Atomic coordinates and charges as read by pychem.GeoOpt._ReadCoordinates
    """
    temp = []
    for i in ChargeCoordinates:
        temp.append([float(i[1]), float(i[2]), float(i[3])])
    DistanceMatrix = GetGeometricalDistanceMatrix(temp)
    temp1 = scipy.amax(DistanceMatrix, axis=0)
    return round(max(temp1), 3)


def CalculateTopoElectronic(ChargeCoordinates: List[List[float]]) -> float:
    """Calculate Topographic electronic descriptors.

    :param ChargeCoordinates: Atomic coordinates and charges as read by pychem.GeoOpt._ReadCoordinates
    """
    pass


def CalculateGravitational3D1(mol: pybel.Molecule, ChargeCoordinates: List[List[float]]) -> float:
    """Calculate Gravitational 3D index from all atoms.

    :param mol: molecule
    :param ChargeCoordinates: Atomic coordinates and charges as read by pychem.GeoOpt._ReadCoordinates
    """
    mol.removeh()
    mol.addh()
    temp = []
    for i, j in enumerate(ChargeCoordinates):
        temp.append([mol.atoms[i].atomicmass, [float(j[1]), float(j[2]), float(j[3])]])
    nAT = len(temp)
    result = 0.0
    for i in range(nAT - 1):
        for j in range(i + 1, nAT):
            dis = GetAtomDistance(temp[i][1], temp[j][1])
            result += temp[i][0] * temp[j][0] / scipy.power(dis, p=2)
    return round(float(result) / 100, 3)
            

def CalculateGravitational3D2((mol,ChargeCoordinates)):
    """Calculate Gravitational 3D index from bonded atoms.

    Katritzky, A.R. et al., J.Phys.Chem. 1996, 100, 10400-10407]
    :param mol: molecule
    :param ChargeCoordinates: Atomic coordinates and charges as read by pychem.GeoOpt._ReadCoordinates
    """
    pass
        


def CalculateRadiusofGyration(mol: pybel.Molecule, ChargeCoordinates: List[List[float]]) -> float:
    """Calculate Radius of gyration.

    :param mol: molecule
    :param ChargeCoordinates: Atomic coordinates and charges as read by pychem.GeoOpt._ReadCoordinates
    """
    mol.addh()
    temp = []
    for i, j in enumerate(ChargeCoordinates):
        temp.append([mol.atoms[i].atomicmass, [float(j[1]), float(j[2]), float(j[3])]])
    nAT = len(temp)
    masscenter = _GetCenterOfMass(temp)
    result = 0.0
    for i in range(nAT):
        dis = GetAtomDistance(temp[i][1], masscenter)
        result += temp[i][0] * scipy.power(dis, p=2)
    return round(scipy.sqrt(float(result / mol.molwt)), 3)


def GetInertiaMatrix(mol: pybel.Molecule, ChargeCoordinates: List[List[float]]) -> scipy.matrix:
    """Get Inertia matrix based on atomic mass and optimized coordinates.

    :param mol: molecule
    :param ChargeCoordinates: Atomic coordinates and charges as read by pychem.GeoOpt._ReadCoordinates
    """
    mol.removeh()
    mol.addh()
    temp = []
    for i, j in enumerate(ChargeCoordinates):
        temp.append([mol.atoms[i].atomicmass, [float(j[1]), float(j[2]), float(j[3])]])
    nAT = len(temp)
    InertiaMatrix = scipy.zeros((3, 3))
    res11 = 0.0
    res22 = 0.0
    res33 = 0.0
    res12 = 0.0
    res13 = 0.0
    res23 = 0.0
    for i in range(nAT):
        res11 += temp[i][0] * (math.pow(temp[i][1][1], 2) + math.pow(temp[i][1][2], 2))
        res22 += temp[i][0] * (math.pow(temp[i][1][0], 2) + math.pow(temp[i][1][2], 2))
        res33 += temp[i][0] * (math.pow(temp[i][1][0], 2) + math.pow(temp[i][1][1], 2))
        res12 += temp[i][0] * (temp[i][1][0] * temp[i][1][1])
        res13 += temp[i][0] * (temp[i][1][0] * temp[i][1][2])
        res23 += temp[i][0] * (temp[i][1][1] * temp[i][1][2])
    InertiaMatrix[0, 0] = res11
    InertiaMatrix[1, 1] = res22
    InertiaMatrix[2, 2] = res33
    InertiaMatrix[0, 1] = res12
    InertiaMatrix[0, 2] = res13
    InertiaMatrix[1, 2] = res23
    InertiaMatrix[1, 0] = res12
    InertiaMatrix[2, 0] = res13
    InertiaMatrix[2, 1] = res23
    return InertiaMatrix


def CalculatePrincipalMomentofInertia(mol: pybel.Molecule, ChargeCoordinates: List[List[float]]) -> dict:
    """Calculate X, Y and Z-principal geometric moments.

    derived from ADAPT developed by Jurs.
    :param mol: molecule
    :param ChargeCoordinates: Atomic coordinates and charges as read by pychem.GeoOpt._ReadCoordinates
    """
    InertiaMatrix = GetInertiaMatrix(mol, ChargeCoordinates)
    ma = scipy.mean(InertiaMatrix, axis=1)
    ms = scipy.std(InertiaMatrix, axis=1, ddof=1)
    bb = scipy.ones((3, 1))
    InertiaMatrix = (InertiaMatrix - bb * ma.T) / (bb * ms.T)
    u, s, v = scipy.linalg.svd(InertiaMatrix)
    res = {}
    res['IA'] = round(s[2], 3)
    res['IB'] = round(s[1], 3)
    res['IC'] = round(s[0], 3)
    return res


def CalculateRatioPMI(mol: pybel.Molecule, ChargeCoordinates: List[List[float]]) -> dict:
    """Calculate the ratio of principal moment of inertia.

    derived from ADAPT developed by Jurs.
    :param mol: molecule
    :param ChargeCoordinates: Atomic coordinates and charges as read by pychem.GeoOpt._ReadCoordinates
    """
    temp = CalculatePrincipalMomentofInertia(mol, ChargeCoordinates)
    res = {}
    res['IA/B'] = round(temp['IA'] / temp['IB'], 3)
    res['IA/C'] = round(temp['IA'] / temp['IC'], 3)
    res['IB/C'] = round(temp['IB'] / temp['IC'], 3)
    return res


def CalculateHarary3D(ChargeCoordinates: List[List[float]]) -> float:
    """Calculate 3D-Harary index as the sum of all the reciprocal geometric distances.

    :param ChargeCoordinates: Atomic coordinates and charges as read by pychem.GeoOpt._ReadCoordinates
    """
    temp = []
    for i in ChargeCoordinates:
        temp.append([float(i[1]), float(i[2]), float(i[3])])
    DistanceMatrix = GetGeometricalDistanceMatrix(temp)
    nAT = len(temp)
    res = 0.0
    for i in range(nAT - 1):
        for j in range(i + 1, nAT):
            if DistanceMatrix[i, j] == 0:
                cds = 0.0
            else:
                cds = 1. / DistanceMatrix[i, j]
            res = res + cds
    return round(res, 3)


def CalculateAverageGeometricalDistanceDegree(ChargeCoordinates: List[List[float]]) -> float:
    """Calculate the average geometric distance degree (AGDD).

    This is the ratio between the sum of all geometric distance degrees and the atoms.
    :param ChargeCoordinates: Atomic coordinates and charges as read by pychem.GeoOpt._ReadCoordinates
    """
    temp = []
    for i in ChargeCoordinates:
        temp.append([float(i[1]), float(i[2]), float(i[3])])
    DistanceMatrix = GetGementricalDistanceMatrix(temp)
    nAT = len(temp)
    res = sum(sum(DistanceMatrix)) / nAT
    return round(res, 3)


def CalculateAbsEigenvalueSumOnGeometricMatrix(ChargeCoordinates: List[List[float]]) -> float:
    """Calculate the absolute eigenvalue sum on geometry matrix (SEig).

    This is the sum of the absolute eigenvalues of the geometry matrix.
    :param ChargeCoordinates: Atomic coordinates and charges as read by pychem.GeoOpt._ReadCoordinates
    """
    temp = []
    for i in ChargeCoordinates:
        temp.append([float(i[1]), float(i[2]), float(i[3])])
    DistanceMatrix = GetGementricalDistanceMatrix(temp) 
    
    u, s, vt = scipy.linalg.svd(DistanceMatrix)
    return round(sum(abs(s)), 3)


def CalculateSPANR(mol: pybel.Molecule, ChargeCoordinates: List[List[float]]) -> float:
    """Calculate the span R.

    This is defined as the radius of the smallest sphere,
    centred on the centre of mass, completely enclosing all atoms of a molecule.

    Arteca G.A. et al., Molecular Shape Descriptors in Reviews in
    Computational Chemistry - Vol. 9, K.B. Lipkowitz, D. Boyd (Eds.),
    VCH Publishers, New York (NY), pp. 191-253, 1991.
    :param mol: molecule
    :param ChargeCoordinates: Atomic coordinates and charges as read by pychem.GeoOpt._ReadCoordinates
    """
    mol.removeh()
    mol.addh()
    temp = []
    for i, j in enumerate(ChargeCoordinates):
        temp.append([mol.atoms[i].atomicmass, [float(j[1]), float(j[2]), float(j[3])]])
    masscenter = _GetMassCenter(temp)  
    res = []
    for i in temp:
        res.append(GetAtomDistance(i[1], masscenter))
    return round(float(max(res)), 3)


def CalculateAverageSPANR(mol: pybel.Molecule, ChargeCoordinates: List[List[float]]) -> float:
    """Calculate the average span R (SPAM).

    This is the root square of the ratio of SPAN over the number of atoms.
    :param mol: molecule
    :param ChargeCoordinates: Atomic coordinates and charges as read by pychem.GeoOpt._ReadCoordinates
    """
    mol.removeh()
    mol.addh()
    temp = []
    for i, j in enumerate(ChargeCoordinates):
        temp.append([mol.atoms[i].atomicmass, [float(j[1]), float(j[2]), float(j[3])]])
    nAT = len(temp)
    masscenter = _GetMassCenter(temp)      
    
    res = []
    for i in temp:
        res.append(GetAtomDistance(i[1], masscenter))
    return round(math.pow(float(max(res)) / nAT, 0.5), 3)


def CalculateMolecularEccentricity(mol: pybel.Molecule, ChargeCoordinates: List[List[float]]) -> float:
    """Calculate molecular eccentricity.

    G.A. Arteca, Molecular Shape Descriptors in Reviews
    in Computational Chemistry - Vol. 9, K.B. Lipkowitz, D. Boyd (Eds.),
    VCH Publishers, New York (NY), pp. 191-253, 1991.
    :param mol: molecule
    :param ChargeCoordinates: Atomic coordinates and charges as read by pychem.GeoOpt._ReadCoordinates
    """
    InertiaMatrix = GetInertiaMatrix(mol, ChargeCoordinates)
    u, s, v = scipy.linalg.svd(InertiaMatrix)
    res1 = s[0]
    res3 = s[2]
    res = math.pow(res1 * res1 - res3 * res3, 1. / 2) / res1
    return round(res, 3)


def GetGeometric(mol: pybel.Molecule) -> dict:
    """Get all (20) geometrical descriptors.

    :param mol: the molecule
    """
    filename = 'temp'
    ChargeCoordinates = _ReadCoordinates(filename)
    res = {}
    res['W3DH'] = Calculate3DWienerWithH(ChargeCoordinates)
    res['W3D'] = Calculate3DWienerWithoutH(ChargeCoordinates)
    res['Petitj3D'] = CalculatePetitjean3DIndex(ChargeCoordinates)
    res['GeDi'] = CalculateGemetricalDiameter(ChargeCoordinates)
    res['grav'] = CalculateGravitational3D1(mol, ChargeCoordinates)
    res['rygr'] = CalculateRadiusofGyration(mol, ChargeCoordinates)
    res['Harary3D'] = CalculateHarary3D(ChargeCoordinates)
    res['AGDD'] = CalculateAverageGeometricalDistanceDegree(ChargeCoordinates)
    res['SEig'] = CalculateAbsEigenvalueSumOnGeometricMatrix(ChargeCoordinates)
    res['SPAN'] = CalculateSPANR(mol, ChargeCoordinates)
    res['ASPAN'] = CalculateAverageSPANR(mol, ChargeCoordinates)
    res['MEcc'] = CalculateMolecularEccentricity(mol, ChargeCoordinates)
    #res.update(CalculatePrincipalMomentofInertia(mol,ChargeCoordinates))
    #res.update(CalculateRatioPMI(mol,ChargeCoordinates))
    return res


############
# if __name__=="__main__":
#     from GeoOpt import GetARCFile
#     mol='C1C=CCCS1'
#     mol='ClC(Cl)(Cl)Cl'
#     inputmol=pybel.readstring('smi',mol)
#     dir = GetARCFile(inputmol)
#     result = GetGeometric(inputmol, dir)
#     print(result)
#     print(len(result))
#     shutil.rmtree(dir, ignore_errors=True)
