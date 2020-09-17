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
    f = file(filename, 'r')
    for line in f.readlines():
        if line[10:27] == "HEAT OF FORMATION":
            inputdict['Hf'] = float(line[-18:-6]) / 96.4853 ##1ev=96.4853kj/mol
        if line[10:22] == "TOTAL ENERGY":
            inputdict['ET'] = float(line[-16:-4])
        if line[10:16] == "DIPOLE":
            inputdict['mu'] = float(line[-16:-7])
        if line[10:28] == "HOMO LUMO ENERGIES":
            inputdict['EHomo'] = float(line[-7:-1])
            inputdict['ELumo'] = float(line[-19:-8])
        #if line[10:26]=="MOLECULAR WEIGHT":
        #    inputdict['Mw']=float(line[-12:-1])
        #if line[10:20]=="COSMO AREA":
        #    inputdict['CoArea']=float(line[-24:-17])
        #if line[10:22]=="COSMO VOLUME":
        #    inputdict['CoVolume']=float(line[-24:-17])
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
    ##converting the unit into suitable one
    EHomo = inputdict['EHomo']
    ELumo = inputdict['ELumo']

    
    dict = {}
    dict.update(inputdict)
    dict['GAP'] = ELumo-EHomo
    dict['S'] = 2. / (ELumo-EHomo)
    dict['eta'] = (ELumo -EHomo) / 2.0
    dict['fHL'] = EHomo / ELumo
    dict['IP'] = -EHomo
    dict['EA'] = -ELumo
    dict['xmu'] = (-ELumo - EHomo) / 2.0
    return dict


def GetQuantumChemistry():
    """Get all quantum chemistry descriptors."""
    filename = 'temp'
    inputdict = ReadFile(filename)
    # res = CalculateBasicQC(inputdict)
    res = {}
    res.update(GetChargeDescriptors(filename))
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
