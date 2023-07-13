# -*- coding: utf-8 -*-


"""3D MoRSE descriptors."""


import math
from typing import List

import numpy as np
from rdkit import Chem

from .atomproperty import get_relative_atomic_property
from .geo_opt import read_coordinates
from .utils import get_geometrical_distance_matrix, get_r

# Parameter for RDF equation
_beta = 100

class MoRSE:
    """3D MoRSE descriptors."""

    @staticmethod
    def calculate_unweight_morse(charge_coordinates: List[List[float]]) -> dict:
        """Calculate unweighted 3D MoRse descriptors.

        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        R = get_r(n=30)
        temp = []
        for i in charge_coordinates:
            # if i[0]!='H':
            temp.append([float(i[1]), float(i[2]), float(i[3])])
        DM = get_geometrical_distance_matrix(temp)
        nAT = len(temp)
        RDFresult = {}
        for kkk, Ri in enumerate(R):
            res = 0.0
            for j in range(nAT - 1):
                for k in range(j + 1, nAT):
                    res = res + math.sin(Ri * DM[j, k]) / (Ri * DM[j, k])
            RDFresult[f'MoRSEU{kkk + 1}'] = res
        return RDFresult

    @staticmethod
    def calculate_charge_morse(charge_coordinates: List[List[float]]) -> dict:
        """Calculate 3D MoRse descriptors from atomic charge.

        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        R = get_r(n=30)
        temp = []
        charge = []
        for i in charge_coordinates:
            # if i[0]!='H':
            temp.append([float(i[1]), float(i[2]), float(i[3])])
            charge.append(float(i[4]))
        DM = get_geometrical_distance_matrix(temp)
        nAT = len(temp)
        RDFresult = {}
        for kkk, Ri in enumerate(R):
            res = 0.0
            for j in range(nAT - 1):
                for k in range(j + 1, nAT):
                    res = res + charge[j] * charge[k] * math.sin(Ri * DM[j, k]) / (Ri * DM[j, k])
            RDFresult[f'MoRSEC{kkk + 1}'] = res
        return RDFresult

    @staticmethod
    def calculate_mass_morse(mol: Chem.Mol, charge_coordinates: List[List[float]]) -> dict:
        """Calculate 3D MoRse descriptors from atomic mass.

        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        R = get_r(n=30)
        temp = []
        mass = [atom.GetMass() for atom in Chem.AddHs(mol).GetAtoms()]
        for i in charge_coordinates:
            # if i[0]!='H':
            temp.append([float(i[1]), float(i[2]), float(i[3])])
        DM = get_geometrical_distance_matrix(temp)
        nAT = len(temp)
        RDFresult = {}
        for kkk, Ri in enumerate(R):
            res = 0.0
            for j in range(nAT - 1):
                for k in range(j + 1, nAT):
                    res += mass[j] * mass[k] * math.sin(Ri * DM[j, k]) / (Ri * DM[j, k])
            RDFresult[f'MoRSEM{kkk + 1}'] = res / 144
        return RDFresult

    @staticmethod
    def calculate_atomic_number_morse(mol: Chem.Mol, charge_coordinates: List[List[float]]) -> dict:
        """Calculate 3D MoRse descriptors from atomic number.

        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        R = get_r(n=30)
        temp = []
        mass = [atom.GetMass() for atom in Chem.AddHs(mol).GetAtoms()]
        for i in charge_coordinates:
            # if i[0]!='H':
            temp.append([float(i[1]), float(i[2]), float(i[3])])
        DM = get_geometrical_distance_matrix(temp)
        nAT = len(temp)
        RDFresult = {}
        for kkk, Ri in enumerate(R):
            res = 0.0
            for j in range(nAT - 1):
                for k in range(j + 1, nAT):
                    res += mass[j] * mass[k] * math.sin(Ri * DM[j, k]) / (Ri * DM[j, k])
            RDFresult[f'MoRSEN{kkk + 1}'] = res / 144
        return RDFresult

    @staticmethod
    def calculate_polarizability_morse(charge_coordinates: List[List[float]]) -> dict:
        """Calculate 3D MoRse descriptors from atomic polarizablity.

        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        R = get_r(n=30)
        temp = []
        polarizability = []
        for i in charge_coordinates:
            # if i[0]!='H':
            temp.append([float(i[1]), float(i[2]), float(i[3])])
            polarizability.append(get_relative_atomic_property(i[0], 'alapha'))
        DM = get_geometrical_distance_matrix(temp)
        nAT = len(temp)
        RDFresult = {}
        for kkk, Ri in enumerate(R):
            res = 0.0
            for j in range(nAT - 1):
                for k in range(j + 1, nAT):
                    res = res + polarizability[j] * polarizability[k] * math.sin(Ri * DM[j, k]) / (Ri * DM[j, k])
            RDFresult[f'MoRSEP{kkk + 1}'] = res
        return RDFresult

    @staticmethod
    def calculate_sanderson_electronegativity_morse(charge_coordinates: List[List[float]]) -> dict:
        """Calculate 3D MoRse descriptors from Sanderson electronegativity.

        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        R = get_r(n=30)
        temp = []
        En = []
        for i in charge_coordinates:
            # if i[0]!='H':
            temp.append([float(i[1]), float(i[2]), float(i[3])])
            En.append(get_relative_atomic_property(i[0], 'En'))
        DM = get_geometrical_distance_matrix(temp)
        nAT = len(temp)
        RDFresult = {}
        for kkk, Ri in enumerate(R):
            res = 0.0
            for j in range(nAT - 1):
                for k in range(j + 1, nAT):
                    res = res + En[j] * En[k] * math.sin(Ri * DM[j, k]) / (Ri * DM[j, k])
            RDFresult[f'MoRSEE{kkk + 1}'] = res
        return RDFresult

    @staticmethod
    def calculate_vdw_volume_morse(charge_coordinates: List[List[float]]) -> dict:
        """Calculate 3D MoRse descriptors from van der Waals volume.

        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        R = get_r(n=30)
        temp = []
        VDW = []
        for i in charge_coordinates:
            # if i[0]!='H':
            temp.append([float(i[1]), float(i[2]), float(i[3])])
            VDW.append(get_relative_atomic_property(i[0], 'V'))
        DM = get_geometrical_distance_matrix(temp)
        nAT = len(temp)
        RDFresult = {}
        for kkk, Ri in enumerate(R):
            res = 0.0
            for j in range(nAT - 1):
                for k in range(j + 1, nAT):
                    res = res + VDW[j] * VDW[k] * math.sin(Ri * DM[j, k]) / (Ri * DM[j, k])
            RDFresult[f'MoRSEV{kkk + 1}'] = res
        return RDFresult

    @staticmethod
    def get_morse_unweighted(arc_file: str) -> dict:
        """Get all unweighted 3D-Morse descriptors.

        :param arc_file: Path to MOPAC .arc file
        """
        ChargeCoordinates = read_coordinates(arc_file)
        result = MoRSE.CalculateUnweightMoRSE(ChargeCoordinates)
        return result

    @staticmethod
    def get_morse_charge(arc_file: str) -> dict:
        """Get all 3D-Morse descriptors from charge schemes.

        :param arc_file: Path to MOPAC .arc file
        """
        ChargeCoordinates = read_coordinates(arc_file)
        result = MoRSE.CalculateChargeMoRSE(ChargeCoordinates)
        return result

    @staticmethod
    def get_morse_mass(mol: Chem.Mol, arc_file: str) -> dict:
        """Get all 3D-Morse descriptors from on mass schemes.

        :param arc_file: Path to MOPAC .arc file
        """
        ChargeCoordinates = read_coordinates(arc_file)
        result = MoRSE.CalculateMassMoRSE(mol, ChargeCoordinates)
        return result

    @staticmethod
    def get_morse_atomic_number(mol: Chem.Mol, arc_file: str) -> dict:
        """Get all 3D-Morse descriptors from atomic number schemes.

        :param arc_file: Path to MOPAC .arc file
        """
        ChargeCoordinates = read_coordinates(arc_file)
        result = MoRSE.CalculateAtomicNumberMoRSE(mol, ChargeCoordinates)
        return result

    @staticmethod
    def get_morse_polarizability(arc_file: str) -> dict:
        """Get all 3D-Morse descriptors from polarizability schemes.

        :param arc_file: Path to MOPAC .arc file
        """
        ChargeCoordinates = read_coordinates(arc_file)
        result = MoRSE.CalculatePolarizabilityMoRSE(ChargeCoordinates)
        return result

    @staticmethod
    def get_morse_sanderson_electronegativity(arc_file: str) -> dict:
        """Get all 3D-Morse descriptors from Sanderson Electronegativity schemes.

        :param arc_file: Path to MOPAC .arc file
        """
        ChargeCoordinates = read_coordinates(arc_file)
        result = MoRSE.CalculateSandersonElectronegativityMoRSE(ChargeCoordinates)
        return result

    @staticmethod
    def get_morse_vdw_volume(arc_file: str) -> dict:
        """Get all 3D-Morse descriptors from VDW Volume schemes.

        :param arc_file: Path to MOPAC .arc file
        """
        ChargeCoordinates = read_coordinates(arc_file)
        result = MoRSE.CalculateVDWVolumeMoRSE(ChargeCoordinates)
        return result

    @staticmethod
    def get_all(mol: Chem.Mol, arc_file: str) -> dict:
        """Get all (210) 3D-Morse descriptors with different (un)weighted schemes.

        :param arc_file: Path to MOPAC .arc file
        """
        result = {}
        charge_coordinates = read_coordinates(arc_file)
        result.update(MoRSE.calculate_unweight_morse(charge_coordinates))
        result.update(MoRSE.calculate_charge_morse(charge_coordinates))
        result.update(MoRSE.calculate_mass_morse(mol, charge_coordinates)) if mol is not None else np.nan
        result.update(MoRSE.calculate_atomic_number_morse(mol, charge_coordinates)) if mol is not None else np.nan
        result.update(MoRSE.calculate_polarizability_morse(charge_coordinates))
        result.update(MoRSE.calculate_sanderson_electronegativity_morse(charge_coordinates))
        result.update(MoRSE.calculate_vdw_volume_morse(charge_coordinates))
        return result
