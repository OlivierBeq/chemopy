# -*- coding: utf-8 -*-


"""3D MoRSE descriptors."""

from typing import List

import numpy as np
from rdkit import Chem

from .atom_property import get_relative_atomic_property
from .geo_opt import read_coordinates
from .utils import get_geometrical_distance_matrix, get_r

# Parameter for RDF equation
_beta = 100


def _morse_sums(DM: np.array, weights: np.array, R_values: List[float]) -> np.array:
    """Vectorized equivalent of the reference triple loop.

    For each R in R_values, computes sum over atom pairs i<j of
    `weights[i] * weights[j] * sin(R * DM[i, j]) / (R * DM[i, j])`, which the reference
    implementation computed with a Python-level `for R: for j: for k>j` loop nest.
    """
    n = DM.shape[0]
    iu = np.triu_indices(n, k=1)
    d = DM[iu]
    w = weights[iu[0]] * weights[iu[1]]
    R_arr = np.asarray(R_values)[:, None]
    rd = R_arr * d[None, :]
    return np.sum(w[None, :] * np.sin(rd) / rd, axis=1)


class MoRSE:
    """3D MoRSE descriptors."""

    @staticmethod
    def calculate_unweight_morse(charge_coordinates: List[List[float]]) -> dict:
        """Calculate unweighted 3D MoRse descriptors.

        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        R = get_r(n=30)
        temp = [[float(i[1]), float(i[2]), float(i[3])] for i in charge_coordinates]
        DM = get_geometrical_distance_matrix(temp)
        weights = np.ones(len(temp))
        sums = _morse_sums(DM, weights, R)
        return {f"MoRSEU{kkk + 1}": res for kkk, res in enumerate(sums)}

    @staticmethod
    def calculate_charge_morse(charge_coordinates: List[List[float]]) -> dict:
        """Calculate 3D MoRse descriptors from atomic charge.

        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        R = get_r(n=30)
        temp = [[float(i[1]), float(i[2]), float(i[3])] for i in charge_coordinates]
        charge = np.array([float(i[4]) for i in charge_coordinates])
        DM = get_geometrical_distance_matrix(temp)
        sums = _morse_sums(DM, charge, R)
        return {f"MoRSEC{kkk + 1}": res for kkk, res in enumerate(sums)}

    @staticmethod
    def calculate_mass_morse(mol: Chem.Mol, charge_coordinates: List[List[float]]) -> dict:
        """Calculate 3D MoRse descriptors from atomic mass.

        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        R = get_r(n=30)
        temp = [[float(i[1]), float(i[2]), float(i[3])] for i in charge_coordinates]
        mass = np.array([atom.GetMass() for atom in Chem.AddHs(mol).GetAtoms()])
        DM = get_geometrical_distance_matrix(temp)
        sums = _morse_sums(DM, mass, R)
        return {f"MoRSEM{kkk + 1}": res / 144 for kkk, res in enumerate(sums)}

    @staticmethod
    def calculate_atomic_number_morse(mol: Chem.Mol, charge_coordinates: List[List[float]]) -> dict:
        """Calculate 3D MoRse descriptors from atomic number.

        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        R = get_r(n=30)
        temp = [[float(i[1]), float(i[2]), float(i[3])] for i in charge_coordinates]
        mass = np.array([atom.GetMass() for atom in Chem.AddHs(mol).GetAtoms()])
        DM = get_geometrical_distance_matrix(temp)
        sums = _morse_sums(DM, mass, R)
        return {f"MoRSEN{kkk + 1}": res / 144 for kkk, res in enumerate(sums)}

    @staticmethod
    def calculate_polarizability_morse(charge_coordinates: List[List[float]]) -> dict:
        """Calculate 3D MoRse descriptors from atomic polarizablity.

        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        R = get_r(n=30)
        temp = [[float(i[1]), float(i[2]), float(i[3])] for i in charge_coordinates]
        polarizability = np.array([get_relative_atomic_property(i[0], "alapha") for i in charge_coordinates])
        DM = get_geometrical_distance_matrix(temp)
        sums = _morse_sums(DM, polarizability, R)
        return {f"MoRSEP{kkk + 1}": res for kkk, res in enumerate(sums)}

    @staticmethod
    def calculate_sanderson_electronegativity_morse(charge_coordinates: List[List[float]]) -> dict:
        """Calculate 3D MoRse descriptors from Sanderson electronegativity.

        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        R = get_r(n=30)
        temp = [[float(i[1]), float(i[2]), float(i[3])] for i in charge_coordinates]
        En = np.array([get_relative_atomic_property(i[0], "En") for i in charge_coordinates])
        DM = get_geometrical_distance_matrix(temp)
        sums = _morse_sums(DM, En, R)
        return {f"MoRSEE{kkk + 1}": res for kkk, res in enumerate(sums)}

    @staticmethod
    def calculate_vdw_volume_morse(charge_coordinates: List[List[float]]) -> dict:
        """Calculate 3D MoRse descriptors from van der Waals volume.

        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        R = get_r(n=30)
        temp = [[float(i[1]), float(i[2]), float(i[3])] for i in charge_coordinates]
        VDW = np.array([get_relative_atomic_property(i[0], "V") for i in charge_coordinates])
        DM = get_geometrical_distance_matrix(temp)
        sums = _morse_sums(DM, VDW, R)
        return {f"MoRSEV{kkk + 1}": res for kkk, res in enumerate(sums)}

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
