# -*- coding: utf-8 -*-


"""Charged partial surface area (CPSA) descriptors."""

from chemopy import asa
from chemopy.GeoOpt import GetAtomClassList, _ReadCoordinates


def GetChargeSA(arc_file, RadiusProbe=1.5, n_sphere_point=960):
    """Get atom symbol, charge and partial solvent-accessible surface areas for all atoms.

    :param arc_file: Path to MOPAC .arc file
    :param RadiusProbe: radius of the probe used to calculate SASA
    :param n_sphere_point: number of points per atom to calculate SASA
    """
    ChargeCoordinates = _ReadCoordinates(arc_file)
    atoms = GetAtomClassList(ChargeCoordinates)
    FASA = asa.calculate_asa(atoms, RadiusProbe, n_sphere_point)
    res = []
    for i in range(len(FASA)):
        res.append([ChargeCoordinates[i][0], ChargeCoordinates[i][4], FASA[i]])
    return res


def CalculateASA(ChargeSA):
    """Calculate solvent-accessible surface area."""
    res = 0.0
    for i in ChargeSA:
        res = res + i[2]
    return res


def CalculateMSA(arc_file):
    """Calculate molecular surface areas.

    :param arc_file: Path to MOPAC .arc file
    """
    ChargeSA = GetChargeSA(arc_file, RadiusProbe=0, n_sphere_point=960)
    res = 0.0
    for i in ChargeSA:
        res = res + i[2]
    return res


def CalculatePNSA1(ChargeSA):
    """Calculate partial negative area."""
    res = 0.0
    for i in ChargeSA:
        if float(i[1]) < 0:
            res = res + i[2]
    return res


def CalculatePPSA1(ChargeSA):
    """Calculate partial positive area."""
    res = 0.0
    for i in ChargeSA:
        if float(i[1]) > 0:
            res = res + i[2]
    return res


def CalculatePNSA2(ChargeSA):
    """Calculate total charge weighted negative surface area."""
    temp1, temp2 = 0.0, 0.0
    for i in ChargeSA:
        if float(i[1]) < 0:
            temp1 += float(i[1])
            temp2 += i[2]
    res = temp1 * temp2
    return res


def CalculatePPSA2(ChargeSA):
    """Calculate total charge weighted positive surface area."""
    temp1, temp2 = 0.0, 0.0
    for i in ChargeSA:
        if float(i[1]) > 0:
            temp1 += float(i[1])
            temp2 += i[2]
    res = temp1 * temp2
    return res


def CalculatePNSA3(ChargeSA):
    """Calculate atom charge weighted negative surface area."""
    res = 0.0
    for i in ChargeSA:
        if float(i[1]) < 0:
            res += float(i[1]) * i[2]
    return res


def CalculatePPSA3(ChargeSA):
    """Calculate atom charge weighted positive surface area."""
    res = 0.0
    for i in ChargeSA:
        if float(i[1]) > 0:
            res += float(i[1]) * i[2]
    return res


def CalculateDPSA1(ChargeSA):
    """Calculate difference in charged partial surface area."""
    return CalculatePPSA1(ChargeSA) - CalculatePNSA1(ChargeSA)


def CalculateDPSA2(ChargeSA):
    """Calculate difference in total charge weighted partial surface area."""
    return CalculatePPSA2(ChargeSA) - CalculatePNSA2(ChargeSA)


def CalculateDPSA3(ChargeSA):
    """Calculate difference in atomic charge weighted surface area."""
    return CalculatePPSA3(ChargeSA) - CalculatePNSA3(ChargeSA)


def CalculateFNSA1(ChargeSA):
    """Calculate fractional charged partial negative surface area."""
    temp = 0.0
    for i in ChargeSA:
        temp += i[2]
    return CalculatePNSA1(ChargeSA) / temp


def CalculateFNSA2(ChargeSA):
    """Calculate fractional charged total negative surface area."""
    temp = 0.0
    for i in ChargeSA:
        temp += i[2]
    return CalculatePNSA2(ChargeSA) / temp


def CalculateFNSA3(ChargeSA):
    """Calculate fractional charged atom negative surface area."""
    temp = 0.0
    for i in ChargeSA:
        temp += i[2]
    return CalculatePNSA3(ChargeSA) / temp


def CalculateFPSA1(ChargeSA):
    """Calculate fractional charged partial positive surface area."""
    temp = 0.0
    for i in ChargeSA:
        temp += i[2]
    return CalculatePPSA1(ChargeSA) / temp


def CalculateFPSA2(ChargeSA):
    """Calculate fractional charged total positive surface area."""
    temp = 0.0
    for i in ChargeSA:
        temp += i[2]
    return CalculatePPSA2(ChargeSA) / temp


def CalculateFPSA3(ChargeSA):
    """Calculate fractional charged atom positive surface area."""
    temp = 0.0
    for i in ChargeSA:
        temp += i[2]
    return CalculatePPSA3(ChargeSA) / temp


def CalculateWNSA1(ChargeSA):
    """Calculate surface weighted charged partial negative surface area."""
    temp = 0.0
    for i in ChargeSA:
        temp += i[2]
    return CalculatePNSA1(ChargeSA) * temp / 1000


def CalculateWNSA2(ChargeSA):
    """Calculate surface weighted charged total negative surface area."""
    temp = 0.0
    for i in ChargeSA:
        temp += i[2]
    return CalculatePNSA2(ChargeSA) * temp / 1000


def CalculateWNSA3(ChargeSA):
    """Calculate surface weighted charged atom negative surface area."""
    temp = 0.0
    for i in ChargeSA:
        temp += i[2]
    return CalculatePNSA3(ChargeSA) * temp / 1000


def CalculateWPSA1(ChargeSA):
    """Calculate surface weighted charged partial positive surface area."""
    temp = 0.0
    for i in ChargeSA:
        temp += i[2]
    return CalculatePPSA1(ChargeSA) * temp / 1000


def CalculateWPSA2(ChargeSA):
    """Calculate surface weighted charged total positive surface area."""
    temp = 0.0
    for i in ChargeSA:
        temp += i[2]
    return CalculatePPSA2(ChargeSA) * temp / 1000


def CalculateWPSA3(ChargeSA):
    """Calculate surface weighted charged atom positive surface area."""
    temp = 0.0
    for i in ChargeSA:
        temp += i[2]
    return CalculatePPSA3(ChargeSA) * temp / 1000


def CalculateTASA(ChargeSA):
    """Calculate total apolar (hydrophobic) surface area."""
    res = 0.0
    for i in ChargeSA:
        if abs(float(i[1])) < 0.2:
            res += i[2]
    return res


def CalculateTPSA(ChargeSA):
    """Calculate total polar surface area."""
    res = 0.0
    for i in ChargeSA:
        if abs(float(i[1])) >= 0.2:
            res += i[2]
    return res


def CalculateRatioTATP(ChargeSA):
    """Calculate ratio between TASA and TPSA (FrTATP)."""
    res = 0.0
    if CalculateTPSA(ChargeSA) == 0:
        return res
    else:
        return CalculateTASA(ChargeSA) / CalculateTPSA(ChargeSA)


def CalculateRASA(ChargeSA):
    """Calculate relative hydrophobic surface area."""
    temp = 0.0
    for i in ChargeSA:
        temp += i[2]
    return CalculateTASA(ChargeSA) / temp


def CalculateRPSA(ChargeSA):
    """Calculate relative polar surface area."""
    temp = 0.0
    for i in ChargeSA:
        temp += i[2]
    return CalculateTPSA(ChargeSA) / temp


def CalculateRNCS(ChargeSA):
    """Calculate relative negative charge surface area."""
    charge = []
    for i in ChargeSA:
        charge.append(float(i[1]))
    temp = []
    for i in ChargeSA:
        temp.append(i[2])
    RNCG = min(charge) / sum(i for i in charge if i < 0)
    return temp[charge.index(min(charge))] / RNCG


def CalculateRPCS(ChargeSA):
    """Calculate relative positive charge surface area."""
    charge = []
    for i in ChargeSA:
        charge.append(float(i[1]))
    temp = []
    for i in ChargeSA:
        temp.append(i[2])
    RPCG = max(charge) / sum(i for i in charge if i > 0)
    return temp[charge.index(min(charge))] / RPCG


def GetCPSA(arc_file):
    """Get all CPSA descriptors.

    :param arc_file: Path to MOPAC .arc file
    """
    res = {}
    ChargeSA = GetChargeSA(arc_file, RadiusProbe=1.5, n_sphere_point=5000)
    res['SASA'] = round(CalculateASA(ChargeSA), 3)
    res['MSA'] = round(CalculateMSA(arc_file), 3)
    res['PNSA1'] = round(CalculatePNSA1(ChargeSA), 3)
    res['PNSA2'] = round(CalculatePNSA2(ChargeSA), 3)
    res['PNSA3'] = round(CalculatePNSA3(ChargeSA), 3)
    res['PPSA1'] = round(CalculatePPSA1(ChargeSA), 3)
    res['PPSA2'] = round(CalculatePPSA2(ChargeSA), 3)
    res['PPSA3'] = round(CalculatePPSA3(ChargeSA), 3)
    res['DPSA1'] = round(CalculateDPSA1(ChargeSA), 3)
    res['DPSA2'] = round(CalculateDPSA2(ChargeSA), 3)
    res['DPSA3'] = round(CalculateDPSA3(ChargeSA), 3)
    res['FNSA1'] = round(CalculateFNSA1(ChargeSA), 3)
    res['FNSA2'] = round(CalculateFNSA2(ChargeSA), 3)
    res['FNSA3'] = round(CalculateFNSA3(ChargeSA), 3)
    res['FPSA1'] = round(CalculateFPSA1(ChargeSA), 3)
    res['FPSA2'] = round(CalculateFPSA2(ChargeSA), 3)
    res['FPSA3'] = round(CalculateFPSA3(ChargeSA), 3)
    res['WNSA1'] = round(CalculateWNSA1(ChargeSA), 3)
    res['WNSA2'] = round(CalculateWNSA2(ChargeSA), 3)
    res['WNSA3'] = round(CalculateWNSA3(ChargeSA), 3)
    res['WPSA1'] = round(CalculateWPSA1(ChargeSA), 3)
    res['WPSA2'] = round(CalculateWPSA2(ChargeSA), 3)
    res['WPSA3'] = round(CalculateWPSA3(ChargeSA), 3)
    res['TASA'] = round(CalculateTASA(ChargeSA), 3)
    res['TPSA'] = round(CalculateTPSA(ChargeSA), 3)
    res['RASA'] = round(CalculateRASA(ChargeSA), 3)
    res['RPSA'] = round(CalculateRPSA(ChargeSA), 3)
    res['RNCS'] = round(CalculateRNCS(ChargeSA), 3)
    res['RPCS'] = round(CalculateRPCS(ChargeSA), 3)
    res['FrTATP'] = round(CalculateRatioTATP(ChargeSA), 3)
    return res


#################################################
# if __name__=="__main__":
#     from GeoOpt import GetARCFile
#     mol='C1C=CCS1'
#     inputmol=pybel.readstring('smi',mol)
#     dir = GetARCFile(inputmol)
#     result=GetCPSA(dir)
#     print(result)
#     print(len(result))
#     shutil.rmtree(dir, ignore_errors=True)
