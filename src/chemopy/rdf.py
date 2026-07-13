# -*- coding: utf-8 -*-


"""3D Radial Distribution Function (RDF) descriptors."""

import numpy as np
from rdkit import Chem

from .atom_property import get_relative_atomic_property
from .geo_opt import read_coordinates
from .utils import get_geometrical_distance_matrix, get_r

_beta = 100


def _rdf_sums(DM: np.array, weights: np.array, R_values) -> np.array:
    """Vectorized equivalent of the reference triple loop.

    For each R in R_values, computes sum over atom pairs i<j of
    `weights[i] * weights[j] * exp(-_beta * (R - DM[i, j])^2)`, which the reference
    implementation computed with a Python-level `for R: for j: for k>j` loop nest.
    """
    n = DM.shape[0]
    iu = np.triu_indices(n, k=1)
    d = DM[iu]
    w = weights[iu[0]] * weights[iu[1]]
    R_arr = np.asarray(R_values)[:, None]
    return np.sum(w[None, :] * np.exp(-_beta * np.power(R_arr - d[None, :], 2)), axis=1)


class RDF:
    """3D Radial Distribution Function descriptors."""

    @staticmethod
    def calculate_unweight_rdf(charge_coordinates):
        """Calculate unweighted RDF descriptors."""
        R = get_r(n=30)
        temp = [[float(i[1]), float(i[2]), float(i[3])] for i in charge_coordinates]
        DM = get_geometrical_distance_matrix(temp)
        weights = np.ones(len(temp))
        sums = _rdf_sums(DM, weights, R)
        return {f"RDFU{kkk + 1}": res for kkk, res in enumerate(sums)}

    @staticmethod
    def calculate_charge_rdf(charge_coordinates):
        """Calculate RDF descriptors with Charge schemes."""
        R = get_r(n=30)
        temp = [[float(i[1]), float(i[2]), float(i[3])] for i in charge_coordinates]
        Charge = np.array([float(i[4]) for i in charge_coordinates])
        DM = get_geometrical_distance_matrix(temp)
        sums = _rdf_sums(DM, Charge, R)
        return {f"RDFC{kkk + 1}": res for kkk, res in enumerate(sums)}

    @staticmethod
    def calculate_mass_rdf(mol: Chem.Mol, charge_coordinates):
        """Calculate RDF descriptors with Mass schemes."""
        mass = np.array([atom.GetMass() for atom in Chem.AddHs(mol).GetAtoms()])
        R = get_r(n=30)
        temp = [[float(i[1]), float(i[2]), float(i[3])] for i in charge_coordinates]
        DM = get_geometrical_distance_matrix(temp)
        sums = _rdf_sums(DM, mass, R)
        return {f"RDFM{kkk + 1}": res / 144 for kkk, res in enumerate(sums)}

    @staticmethod
    def calculate_polarizability_rdf(charge_coordinates):
        """Calculate RDF descriptors with Polarizability schemes."""
        R = get_r(n=30)
        temp = [[float(i[1]), float(i[2]), float(i[3])] for i in charge_coordinates]
        polarizability = np.array([get_relative_atomic_property(i[0], "alapha") for i in charge_coordinates])
        DM = get_geometrical_distance_matrix(temp)
        sums = _rdf_sums(DM, polarizability, R)
        return {f"RDFP{kkk + 1}": res for kkk, res in enumerate(sums)}

    @staticmethod
    def calculate_sanderson_electronegativity_rdf(charge_coordinates):
        """Calculate RDF descriptors with Sanderson Electronegativity schemes."""
        R = get_r(n=30)
        temp = [[float(i[1]), float(i[2]), float(i[3])] for i in charge_coordinates]
        EN = np.array([get_relative_atomic_property(i[0], "En") for i in charge_coordinates])
        DM = get_geometrical_distance_matrix(temp)
        sums = _rdf_sums(DM, EN, R)
        return {f"RDFE{kkk + 1}": res for kkk, res in enumerate(sums)}

    @staticmethod
    def calculate_vdw_volume_rdf(charge_coordinates):
        """Calculate RDF with atomic van der Waals volume shemes."""
        R = get_r(n=30)
        temp = [[float(i[1]), float(i[2]), float(i[3])] for i in charge_coordinates]
        VDW = np.array([get_relative_atomic_property(i[0], "V") for i in charge_coordinates])
        DM = get_geometrical_distance_matrix(temp)
        sums = _rdf_sums(DM, VDW, R)
        return {f"RDFV{kkk + 1}": res for kkk, res in enumerate(sums)}

    @staticmethod
    def get_rdf_unweighed(arc_file):
        """Obtain all Unweighed RDF descriptors."""
        ChargeCoordinates = read_coordinates(arc_file)
        result = RDF.calculate_unweight_rdf(ChargeCoordinates)
        return result

    @staticmethod
    def get_rdf_charge(arc_file):
        """Obtain all RDF descriptors with Charge schemes."""
        ChargeCoordinates = read_coordinates(arc_file)
        result = RDF.calculate_charge_rdf(ChargeCoordinates)
        return result

    @staticmethod
    def get_rdf_mass(mol: Chem.Mol, arc_file):
        """Obtain all RDF descriptors with Mass schemes."""
        ChargeCoordinates = read_coordinates(arc_file)
        result = RDF.calculate_mass_rdf(mol, ChargeCoordinates)
        return result

    @staticmethod
    def get_rdf_polarizability(arc_file):
        """Obtain all RDF descriptors with Polarizability schemes."""
        ChargeCoordinates = read_coordinates(arc_file)
        result = RDF.calculate_polarizability_rdf(ChargeCoordinates)
        return result

    @staticmethod
    def get_rdf_sanderson_electronegativity(arc_file):
        """Obtain all RDF descriptors with Sanderson Electronegativity schemes."""
        ChargeCoordinates = read_coordinates(arc_file)
        result = RDF.calculate_sanderson_electronegativity_rdf(ChargeCoordinates)
        return result

    @staticmethod
    def get_rdf_vdw_volume(arc_file):
        """Obtain all RDF descriptors with VDW Volume schemes."""
        ChargeCoordinates = read_coordinates(arc_file)
        result = RDF.calculate_vdw_volume_rdf(ChargeCoordinates)
        return result

    @staticmethod
    def get_all(mol: Chem.Mol, arc_file: str) -> dict:
        """Obtain all (180) RDF descriptors with different (un)weighted schemes."""
        result = {}
        ChargeCoordinates = read_coordinates(arc_file)
        result.update(RDF.calculate_unweight_rdf(ChargeCoordinates))
        result.update(RDF.calculate_charge_rdf(ChargeCoordinates))
        result.update(RDF.calculate_mass_rdf(mol, ChargeCoordinates)) if mol is not None else np.nan
        result.update(RDF.calculate_polarizability_rdf(ChargeCoordinates))
        result.update(RDF.calculate_sanderson_electronegativity_rdf(ChargeCoordinates))
        result.update(RDF.calculate_vdw_volume_rdf(ChargeCoordinates))
        return result
