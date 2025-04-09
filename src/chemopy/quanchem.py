# -*- coding: utf-8 -*-


"""Quantum chemistry descriptors."""

from typing import List, Tuple

from rdkit import Chem
import numpy as np


class QuantumChemistry:
    """Quantum chemistry descriptors."""

    @staticmethod
    def _get_max(x: float) -> float:
        """Get the maximum of x.

        If x is empty return 0.0.
        """
        if x == []:
            return 0.0
        else:
            return max(x)

    @staticmethod
    def _get_min(x: float) -> float:
        """Get the minimum of x.

        Ff x is empty return 0.0.
        """
        if x == []:
            return 0.0
        else:
            return min(x)

    @staticmethod
    def read_file(filename: str) -> dict:
        """Read basic quantum chemistry descriptors from an .arc file."""
        inputdict = {'Hf': np.nan, 'ET': np.nan, 'mu': np.nan, 'EHomo': np.nan, 'ELumo': np.nan,
                     'ESomo': np.nan, 'EaSomo': np.nan, 'EaLumo': np.nan, 'EbSomo': np.nan, 'EbLumo': np.nan,
                     'mMw': np.nan, 'CoArea': np.nan, 'CoVolume': np.nan}
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
                inputdict['mMw'] = float(line[-12:-1])
            elif value == "COSMO AREA":
                inputdict['CoArea'] = float(line.split('=')[1].strip().split()[0])
            elif value == "COSMO VOLUME":
                inputdict['CoVolume'] = float(line.split('=')[1].strip().split()[0])
        f.close()
        return inputdict

    @staticmethod
    def _read_charge(arc_file: str) -> List[Tuple[str, float]]:
        """Read the charge of each atom in .arc file.

        :param arc_file: Path to MOPAC .arc file
        """
        charge = []
        with open(arc_file, 'r') as f:
            templine = f.readlines()

        for line in range(len(templine)):
            if templine[line][-7:-1] == "CHARGE":
                k = line
                break

        for i in templine[k + 4: len(templine) - 1]:
            temp = i.split()
            charge.append((temp[0].strip(), temp[-1].strip()))
        return charge

    @staticmethod
    def get_charge_descriptors(arc_file: str) -> dict:
        """Calculate charge descriptors.

        :param arc_file: Path to MOPAC .arc file
        """
        res = {}
        Htemp = []
        Ctemp = []
        Ntemp = []
        Otemp = []
        temp = []
        Charge = QuantumChemistry._read_charge(arc_file)
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
        res['qQHmax'] = QuantumChemistry._get_max(Htemp)
        res['qQCmax'] = QuantumChemistry._get_max(Ctemp)
        res['qQNmax'] = QuantumChemistry._get_max(Ntemp)
        res['qQOmax'] = QuantumChemistry._get_max(Otemp)
        res['qQHmin'] = QuantumChemistry._get_min(Htemp)
        res['qQCmin'] = QuantumChemistry._get_min(Ctemp)
        res['qQNmin'] = QuantumChemistry._get_min(Ntemp)
        res['qQOmin'] = QuantumChemistry._get_min(Otemp)
        res['qQmax'] = max(temp)
        res['qQmin'] = min(temp)
        res['qQHss'] = sum(i * i for i in Htemp)
        res['qQCss'] = sum(i * i for i in Ctemp)
        res['qQNss'] = sum(i * i for i in Ntemp)
        res['qQOss'] = sum(i * i for i in Otemp)
        res['qQass'] = sum(i * i for i in temp)
        res['qMpc'] = np.mean([i for i in temp if i > 0])
        res['qTpc'] = sum(i for i in temp if i > 0)
        res['qMnc'] = np.mean([i for i in temp if i < 0])
        res['qTnc'] = sum(i for i in temp if i < 0)
        res['qMac'] = np.mean([np.abs(i) for i in temp])
        res['qTac'] = sum(np.abs(i) for i in temp)
        res['qRpc'] = QuantumChemistry._get_max(temp) / res['qTpc']
        res['qRnc'] = QuantumChemistry._get_min(temp) / res['qTnc']
        return res

    @staticmethod
    def calculate_basic_quantum_chemistry(inputdict: dict) -> dict:
        """Calculate between 38 and 40 quantum chemical descriptors.

        Derived from Lumo, Homo, dipole moment, enthalpy and the total energy.
        """
        if inputdict.get('EHomo') is not None:
            EHomo = inputdict['EHomo']
        else:
            EHomo = inputdict['EaSomo']
        if inputdict.get('ELumo') is not None:
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

    @staticmethod
    def get_all(mol: Chem.Mol, arc_file: str) -> dict:
        """Get all quantum chemistry descriptors.

        :param arc_file: Path to MOPAC .arc file
        """
        inputdict = QuantumChemistry.read_file(arc_file)
        res = QuantumChemistry.calculate_basic_quantum_chemistry(inputdict)
        res.update(QuantumChemistry.get_charge_descriptors(arc_file))
        return res
