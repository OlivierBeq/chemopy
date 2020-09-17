# -*- coding: utf-8 -*-


"""Quantum chemistry descriptors."""

from typing import List, Tuple

import numpy


def _GetMax(x: float) -> float:
    """Get the maximum of x.

    If x is empty return 0.0.
    """
    if x == []:
        return 0.0
    else:
        return max(x)


def _GetMin(x: float) -> float:
    """Get the minimum of x.

    Ff x is empty return 0.0.
    """
    if x == []:
        return 0.0
    else:
        return min(x)


def ReadFile(filename: str) -> dict:
    """Read basic quantum chemistry descriptors from the .arc file."""
    inputdict = {}
    f = open(filename, 'r')
    for line in f.readlines():
        value = line[10:34].strip()
        if value == "HEAT OF FORMATION":
            # 1ev = 96.4853 kj/mol
            inputdict['Hf'] = float(line.strip().upper().split('=')[-1].upper().strip('KJ/MOL')) / 96.4853
        if value == "TOTAL ENERGY":
            inputdict['ET'] = float(line.strip().upper().split('=')[1].upper().strip('EV'))
        if value == "DIPOLE":
            inputdict['mu'] = float(line.strip().upper().split('=')[1].split('DEBYE')[0])
        if value.startswith("HOMO LUMO ENERGIES"):
            inputdict['EHomo'] = float(line.split('=')[1].strip().split()[0])
            inputdict['ELumo'] = float(line.split('=')[1].strip().split()[1])
        elif value.startswith("HOMO (SOMO) LUMO"):
            data = line.split('=')[1].strip().translate(str.maketrans('', '', '()')).split()
            data = list(filter(None, data))
            inputdict['EHomo'] = float(data[0])
            inputdict['ESomo'] = float(data[1])
            inputdict['ELumo'] = float(data[2])
        elif value.startswith("ALPHA SOMO LUMO"):
            data = line.split('=')[1].strip().translate(str.maketrans('', '', '()')).split()
            data = list(filter(None, data))
            inputdict['EaSomo'] = float(data[0])
            inputdict['EaLumo'] = float(data[1])
        elif value.startswith("BETA  SOMO LUMO"):
            data = line.split('=')[1].strip().translate(str.maketrans('', '', '()')).split()
            data = list(filter(None, data))
            inputdict['EbSomo'] = float(data[0])
            inputdict['EbLumo'] = float(data[1])
        if line[10:26] == "MOLECULAR WEIGHT":
            inputdict['Mw'] = float(line[-12:-1])
        elif value == "COSMO AREA":
            inputdict['CoArea'] = float(line.split('=')[1].strip().split()[0])
        elif value == "COSMO VOLUME":
            inputdict['CoVolume'] = float(line.split('=')[1].strip().split()[0])
    f.close()
    return inputdict


def _ReadCharge(arc_file: str) -> List[Tuple[str, float]]:
    """Read the charge of each atom in .arc file.

    :param arc_file: Path to MOPAC .arc file
    """
    Charge = []
    with open(arc_file, 'r') as f:
        templine = f.readlines()

    for line in range(len(templine)):
        if templine[line][-7:-1] == "CHARGE":
            k = line
            break

    for i in templine[k + 4: len(templine) - 1]:
        temp = i.split()
        Charge.append((temp[0].strip(), temp[-1].strip()))
    return Charge


def GetChargeDescriptors(arc_file: str) -> dict:
    """Calculate charge descriptors.

    :param arc_file: Path to MOPAC .arc file
    """
    res = {}
    Htemp = []
    Ctemp = []
    Ntemp = []
    Otemp = []
    temp = []
    Charge = _ReadCharge(arc_file)
    for i in Charge:
        temp.append(float(i[1]))
        if i[0] == 'H':
            Htemp.append(float(i[1]))
        if i[0] == 'C':
            Ctemp.append(float(i[1]))
        if i[0] == 'N':
            Ntemp.append(float(i[1]))
        if i[0] == 'O':
            Otemp.append(float(i[1]))
    res['QHmax'] = round(_GetMax(Htemp), 3)
    res['QCmax'] = round(_GetMax(Ctemp), 3)
    res['QNmax'] = round(_GetMax(Ntemp), 3)
    res['QOmax'] = round(_GetMax(Otemp), 3)
    res['QHmin'] = round(_GetMin(Htemp), 3)
    res['QCmin'] = round(_GetMin(Ctemp), 3)
    res['QNmin'] = round(_GetMin(Ntemp), 3)
    res['QOmin'] = round(_GetMin(Otemp), 3)
    res['Qmax'] = round(max(temp), 3)
    res['Qmin'] = round(min(temp), 3)
    res['QHss'] = round(sum(i * i for i in Htemp), 3)
    res['QCss'] = round(sum(i * i for i in Ctemp), 3)
    res['QNss'] = round(sum(i * i for i in Ntemp), 3)
    res['QOss'] = round(sum(i * i for i in Otemp), 3)
    res['Qass'] = round(sum(i * i for i in temp), 3)
    res['Mpc'] = round(numpy.mean([i for i in temp if i > 0]), 3)
    res['Tpc'] = round(sum(i for i in temp if i > 0), 3)
    res['Mnc'] = round(numpy.mean([i for i in temp if i < 0]), 3)
    res['Tnc'] = round(sum(i for i in temp if i < 0), 3)
    res['Tac'] = round(sum(numpy.abs(i) for i in temp), 3)
    res['Mac'] = round(numpy.mean([numpy.abs(i) for i in temp]), 3)
    res['Rpc'] = round(_GetMax(temp) / res['Tpc'], 3)
    res['Rnc'] = round(_GetMin(temp) / res['Tnc'], 3)
    return res


def CalculateBasicQC(inputdict: dict) -> dict:
    """Calculate between 38 and 40 quantum chemical descriptors.

    Derived from Lumo, Homo, dipole moment, enthalpy and the total energy.
    """
    if 'EHomo' in inputdict.keys():
        EHomo = inputdict['EHomo']
    else:
        EHomo = inputdict['EaSomo']
    if 'ELumo' in inputdict.keys():
        ELumo = inputdict['ELumo']
    else:
        ELumo = inputdict['EaLumo']
    dict_ = {}
    dict_.update(inputdict)
    dict_['GAP'] = ELumo - EHomo
    dict_['S'] = 2. / (ELumo - EHomo)
    dict_['eta'] = (ELumo - EHomo) / 2.0
    dict_['fHL'] = EHomo / ELumo
    dict_['IP'] = -EHomo
    dict_['EA'] = -ELumo
    dict_['xmu'] = (-ELumo - EHomo) / 2.0
    return dict_


def GetQuantumChemistry(arc_file: str) -> dict:
    """Get all quantum chemistry descriptors.

    :param arc_file: Path to MOPAC .arc file
    """
    inputdict = ReadFile(arc_file)
    res = CalculateBasicQC(inputdict)
    res.update(GetChargeDescriptors(arc_file))
    return res


# if __name__=="__main__":

#     from GeoOpt import GetARCFile
#     mol='CC(N)C(=O)O'
#     inputmol=pybel.readstring('smi',mol)
#     dir_ = GetARCFile(inputmol)
#     result=GetQuantumChemistry(dir_)
#     print(result)
#     print(len(result))
#     shutil.rmtree(dir_, ignore_errors=True)
