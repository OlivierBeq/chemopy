# -*- coding: utf-8 -*-


"""3D Radial Distribution Function (RDF) descriptors."""

import math

import numpy as np
from rdkit import Chem

from .atom_property import get_relative_atomic_property
from .geo_opt import read_coordinates
from .utils import get_geometrical_distance_matrix, get_r

_beta = 100

class RDF:
    """3D Radial Distribution Function descriptors."""

    @staticmethod
    def calculate_unweight_rdf(charge_coordinates):
        """Calculate unweighted RDF descriptors."""
        R = get_r(n=30)
        temp = []
    #    ChargeCoordinates=_ReadCoordinates('temp.arc')
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
                    res = res + math.exp(-_beta * math.pow(Ri - DM[j, k], 2))
            RDFresult[f'RDFU{kkk + 1}'] = res
        return RDFresult

    @staticmethod
    def calculate_charge_rdf(charge_coordinates):
        """Calculate RDF descriptors with Charge schemes."""
        R = get_r(n=30)
        temp = []
        Charge = []
    #    ChargeCoordinates=_ReadCoordinates('temp.arc')
        for i in charge_coordinates:
            # if i[0]!='H':
            temp.append([float(i[1]), float(i[2]), float(i[3])])
            Charge.append(float(i[4]))
        DM = get_geometrical_distance_matrix(temp)
        nAT = len(temp)
        RDFresult = {}
        for kkk, Ri in enumerate(R):
            res = 0.0
            for j in range(nAT - 1):
                for k in range(j + 1, nAT):
                    res = res + Charge[j] * Charge[k] * math.exp(-_beta * math.pow(Ri - DM[j, k], 2))
            RDFresult[f'RDFC{kkk + 1}'] = res
        return RDFresult

    @staticmethod
    def calculate_mass_rdf(mol: Chem.Mol, charge_coordinates):
        """Calculate RDF descriptors with Mass schemes."""
        mass = [atom.GetMass() for atom in Chem.AddHs(mol).GetAtoms()]
        R = get_r(n=30)
        temp = []
    #    ChargeCoordinates=_ReadCoordinates('temp.arc')
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
                    res = res + mass[j] * mass[k] * math.exp(-_beta * math.pow(Ri - DM[j, k], 2))
            RDFresult[f'RDFM{kkk + 1}'] = res / 144
        return RDFresult

    @staticmethod
    def calculate_polarizability_rdf(charge_coordinates):
        """Calculate RDF descriptors with Polarizability schemes."""
        R = get_r(n=30)
        temp = []
        polarizability = []
    #    ChargeCoordinates=_ReadCoordinates('temp.arc')
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
                    res = res + polarizability[j] * polarizability[k] * math.exp(-_beta * math.pow(Ri - DM[j, k], 2))
            RDFresult[f'RDFP{kkk + 1}'] = res
        return RDFresult

    @staticmethod
    def calculate_sanderson_electronegativity_rdf(charge_coordinates):
        """Calculate RDF descriptors with Sanderson Electronegativity schemes."""
        R = get_r(n=30)
        temp = []
        EN = []
    #    ChargeCoordinates=_ReadCoordinates('temp.arc')
        for i in charge_coordinates:
            # if i[0]!='H':
            temp.append([float(i[1]), float(i[2]), float(i[3])])
            EN.append(get_relative_atomic_property(i[0], 'En'))
        DM = get_geometrical_distance_matrix(temp)
        nAT = len(temp)
        RDFresult = {}
        for kkk, Ri in enumerate(R):
            res = 0.0
            for j in range(nAT - 1):
                for k in range(j + 1, nAT):
                    res = res + EN[j] * EN[k] * math.exp(-_beta * math.pow(Ri - DM[j, k], 2))
            RDFresult[f'RDFE{kkk + 1}'] = res
        return RDFresult

    @staticmethod
    def calculate_vdw_volume_rdf(charge_coordinates):
        """Calculate RDF with atomic van der Waals volume shemes."""
        R = get_r(n=30)
        temp = []
        VDW = []
    #    ChargeCoordinates=_ReadCoordinates('temp.arc')
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
                    res = res + VDW[j] * VDW[k] * math.exp(-_beta * math.pow(Ri - DM[j, k], 2))
            RDFresult[f'RDFV{kkk + 1}'] = res
        return RDFresult

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
