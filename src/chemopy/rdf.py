# -*- coding: utf-8 -*-


"""3D Radial Distribution Function (RDF) descriptors."""

import math

from pychem.AtomProperty import GetRelativeAtomicProperty
from pychem.GeoOpt import _ReadCoordinates
from pychem.utils import GetGeometricalDistanceMatrix, GetR

_beta = 100


def CalculateUnweightRDF(ChargeCoordinates):
    """Calculate unweighted RDF descriptors."""
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
                res = res + math.exp(-_beta * math.pow(Ri - DM[j, k], 2))
        RDFresult[f'RDFU{kkk + 1}'] = round(res, 3)
    return RDFresult


def CalculateChargeRDF(ChargeCoordinates):
    """Calculate RDF descriptors with Charge schemes."""
    R = GetR(n=30)
    temp = []
    Charge = []
#    ChargeCoordinates=_ReadCoordinates('temp.arc')
    for i in ChargeCoordinates:
        # if i[0]!='H':
        temp.append([float(i[1]), float(i[2]), float(i[3])])
        Charge.append(float(i[4]))
    DM = GetGeometricalDistanceMatrix(temp)
    nAT = len(temp)
    RDFresult = {}
    for kkk, Ri in enumerate(R):
        res = 0.0
        for j in range(nAT - 1):
            for k in range(j + 1, nAT):
                res = res + Charge[j] * Charge[k] * math.exp(-_beta * math.pow(Ri - DM[j, k], 2))
        RDFresult[f'RDFC{kkk + 1}'] = round(res, 3)
    return RDFresult


def CalculateMassRDF(mol, ChargeCoordinates):
    """Calculate RDF descriptors with Mass schemes."""
    mol.addh()
    mass = [i.atomicmass for i in mol.atoms]
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
                res = res + mass[j] * mass[k] * math.exp(-_beta * math.pow(Ri - DM[j, k], 2))
        RDFresult[f'RDFM{kkk + 1}'] = round(res / 144, 3)
    return RDFresult


def CalculatePolarizabilityRDF(ChargeCoordinates):
    """Calculate RDF descriptors with Polarizability schemes."""
    R = GetR(n=30)
    temp = []
    polarizability = []
#    ChargeCoordinates=_ReadCoordinates('temp.arc')
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
                res = res + polarizability[j] * polarizability[k] * math.exp(-_beta * math.pow(Ri - DM[j, k], 2))
        RDFresult[f'RDFP{kkk + 1}'] = round(res, 3)
    return RDFresult


def CalculateSandersonElectronegativityRDF(ChargeCoordinates):
    """Calculate RDF descriptors with Sanderson Electronegativity schemes."""
    R = GetR(n=30)
    temp = []
    EN = []
#    ChargeCoordinates=_ReadCoordinates('temp.arc')
    for i in ChargeCoordinates:
        # if i[0]!='H':
        temp.append([float(i[1]), float(i[2]), float(i[3])])
        EN.append(GetRelativeAtomicProperty(i[0], 'En'))
    DM = GetGeometricalDistanceMatrix(temp)
    nAT = len(temp)
    RDFresult = {}
    for kkk, Ri in enumerate(R):
        res = 0.0
        for j in range(nAT - 1):
            for k in range(j + 1, nAT):
                res = res + EN[j] * EN[k] * math.exp(-_beta * math.pow(Ri - DM[j, k], 2))
        RDFresult[f'RDFE{kkk + 1}'] = round(res, 3)
    return RDFresult


def CalculateVDWVolumeRDF(ChargeCoordinates):
    """Calculate RDF with atomic van der Waals volume shemes."""
    R = GetR(n=30)
    temp = []
    VDW = []
#    ChargeCoordinates=_ReadCoordinates('temp.arc')
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
                res = res + VDW[j] * VDW[k] * math.exp(-_beta * math.pow(Ri - DM[j, k], 2))
        RDFresult[f'RDFV{kkk + 1}'] = round(res, 3)
    return RDFresult


def GetRDFUnweighed(mol, arc_file):
    """Obtain all Unweighed RDF descriptors."""
    ChargeCoordinates = _ReadCoordinates(arc_file)
    result = CalculateUnweightRDF(ChargeCoordinates)
    return result


def GetRDFCharge(mol, arc_file):
    """Obtain all RDF descriptors with Charge schemes."""
    ChargeCoordinates = _ReadCoordinates(arc_file)
    result = CalculateChargeRDF(ChargeCoordinates)
    return result


def GetRDFMass(mol, arc_file):
    """Obtain all RDF descriptors with Mass schemes."""
    ChargeCoordinates = _ReadCoordinates(arc_file)
    result = CalculateMassRDF(mol, ChargeCoordinates)
    return result


def GetRDFPolarizability(mol, arc_file):
    """Obtain all RDF descriptors with Polarizability schemes."""
    ChargeCoordinates = _ReadCoordinates(arc_file)
    result = CalculatePolarizabilityRDF(ChargeCoordinates)
    return result


def GetRDFSandersonElectronegativity(mol, arc_file):
    """Obtain all RDF descriptors with Sanderson Electronegativity schemes."""
    ChargeCoordinates = _ReadCoordinates(arc_file)
    result = CalculateSandersonElectronegativityRDF(ChargeCoordinates)
    return result


def GetRDFVDWVolume(mol, arc_file):
    """Obtain all RDF descriptors with VDW Volume schemes."""
    ChargeCoordinates = _ReadCoordinates(arc_file)
    result = CalculateVDWVolumeRDF(ChargeCoordinates)
    return result


def GetRDF(mol, arc_file):
    """Obtain all (180) RDF descriptors with different (un)weighted schemes."""
    result = {}
    ChargeCoordinates = _ReadCoordinates(arc_file)
    result.update(CalculateUnweightRDF(ChargeCoordinates))
    result.update(CalculateChargeRDF(ChargeCoordinates))
    result.update(CalculateMassRDF(mol, ChargeCoordinates))
    result.update(CalculatePolarizabilityRDF(ChargeCoordinates))
    result.update(CalculateSandersonElectronegativityRDF(ChargeCoordinates))
    result.update(CalculateVDWVolumeRDF(ChargeCoordinates))
    return result


# def _GetHTMLDoc():
#     """
# #     Write HTML documentation for this module.
# #     """
#     import pydoc
#     pydoc.writedoc('rdf')
############################################################################
# if __name__=="__main__":

#     import pybel
#     from GeoOpt import GetARCFile
#     import shutil
#     mol='C1C=CCCS1'
#     inputmol=pybel.readstring('smi',mol)
#     dir_ = GetARCFile(inputmol)
#     filename = os.path.join(dir_, 'temp')
#     ChargeCoordinates=_ReadCoordinates(filename)
#     res=CalculateVDWVolumeRDF(ChargeCoordinates)
#     print(res)
#     print(len(GetRDF(inputmol, dir_)))
#     shutil.rmtree(dir_, ignore_errors=True)
