# -*- coding: utf-8 -*-


"""Geometrical molecular descriptors.

Optimizes molecular structure using MOPAC.
"""

import math
from typing import List, Tuple

import numpy as np
import numpy.lib.scimath
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt

from .topology import periodicTable
from .geo_opt import read_coordinates
from .utils import get_atom_distance, get_geometrical_distance_matrix


class Geometric:
    """Geometrical molecular descriptors."""

    @staticmethod
    def _get_center_of_mass(mass_coordinates: List[Tuple[float, Tuple[float, float, float]]]) -> Tuple[float, float, float]:
        """Get the center of mass.

        :param mass_coordinates: list of atomic masses and coordinates
                                in the format [(atommass, (x,y,z)), ...]
        """
        res1 = 0.0
        res2 = 0.0
        res3 = 0.0
        temp = []
        for i in mass_coordinates:
            res1 += i[0] * i[1][0]
            res2 += i[0] * i[1][1]
            res3 += i[0] * i[1][2]
            temp.append(i[0])
        result = (res1 / sum(temp), res2 / sum(temp), res3 / sum(temp))
        return result

    @staticmethod
    def _get_geometrical_center(charge_coordinates: List[List[float]]) -> Tuple[float, float, float]:
        """Get the geometrical center.

        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        res1 = []
        res2 = []
        res3 = []
        for i in charge_coordinates:
            res1.append(float(i[1]))
            res2.append(float(i[2]))
            res3.append(float(i[3]))
        result = (np.mean(res1), np.mean(res2), np.mean(res3))
        return result

    @staticmethod
    def calculate_3D_wiener_with_h(charge_coordinates: List[List[float]]) -> float:
        """Calculate 3D Wiener index from geometrical distance matrix of a MOPAC optimized molecule (including Hs).

        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        temp = []
        for i in charge_coordinates:
            temp.append([float(i[1]), float(i[2]), float(i[3])])
        distance_matrix = get_geometrical_distance_matrix(temp)
        return np.sum(distance_matrix) / 2.0

    @staticmethod
    def calculate_3D_wiener_without_h(charge_coordinates: List[List[float]]) -> float:
        """Calculate 3D Wiener index from geometrical distance matrix of a MOPAC optimized molecule (not including Hs).

        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        temp = []
        for i in charge_coordinates:
            if i[0] != 'H':
                temp.append([float(i[1]), float(i[2]), float(i[3])])
        distance_matrix = get_geometrical_distance_matrix(temp)
        return np.sum(distance_matrix) / 2.0

    @staticmethod
    def calculate_petitjean_3D_index(charge_coordinates: List[List[float]]) -> float:
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

        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        temp = []
        for i in charge_coordinates:
            temp.append([float(i[1]), float(i[2]), float(i[3])])
        distance_matrix = get_geometrical_distance_matrix(temp)
        temp1 = np.amax(distance_matrix, axis=0)
        return max(temp1) / min(temp1) - 1.0

    @staticmethod
    def calculate_geometrical_diameter(charge_coordinates: List[List[float]]) -> float:
        """Calculate the geometrical diameter.

        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        temp = []
        for i in charge_coordinates:
            temp.append([float(i[1]), float(i[2]), float(i[3])])
        distance_matrix = get_geometrical_distance_matrix(temp)
        temp1 = np.amax(distance_matrix, axis=0)
        return max(temp1)

    @staticmethod
    def calculate_topo_electronic(charge_coordinates: List[List[float]]) -> float:
        """Calculate Topographic electronic descriptors.

        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        temp, charges = [], []
        for i in charge_coordinates:
            temp.append([float(i[1]), float(i[2]), float(i[3])])
            charges.append(float(i[4]))
        nAT = len(charge_coordinates)
        result = 0.0
        for i in range(nAT - 1):
            for j in range(i + 1, nAT):
                dis = get_atom_distance(temp[i], temp[j])
                result += np.absolute(charges[i] * charges[j]) / numpy.lib.scimath.power(dis, p=2)
        return result

    @staticmethod
    def calculate_gravitational_3D1(mol: Chem.Mol, charge_coordinates: List[List[float]]) -> float:
        """Calculate Gravitational 3D index from all atoms.

        :param mol: molecule
        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        mol = Chem.AddHs(mol)
        temp = []
        atoms = list(mol.GetAtoms())
        for i, j in enumerate(charge_coordinates):
            temp.append([periodicTable.GetAtomicWeight(atoms[i].GetAtomicNum()),
                         [float(j[1]), float(j[2]), float(j[3])]])
        nAT = len(temp)
        result = 0.0
        for i in range(nAT - 1):
            for j in range(i + 1, nAT):
                dis = get_atom_distance(temp[i][1], temp[j][1])
                result += temp[i][0] * temp[j][0] / numpy.lib.scimath.power(dis, p=2)
        return float(result) / 100

    @staticmethod
    def calculate_gravitational_3D2(mol: Chem.Mol, charge_coordinates: List[List[float]]) -> float:
        """Calculate Gravitational 3D index from bonded atoms.

        Katritzky, A.R. et al., J.Phys.Chem. 1996, 100, 10400-10407]
        :param mol: molecule
        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        raise NotImplementedError('Needs modification: calculated only on bonded atoms.')
        mol = Chem.AddHs(mol)
        temp = []
        atoms = list(mol.GetAtoms())
        for i, j in enumerate(charge_coordinates):
            temp.append([periodicTable.GetAtomicWeight(atoms[i].GetAtomicNum()),
                         [float(j[1]), float(j[2]), float(j[3])]])
        nAT = len(temp)
        result = 0.0
        for i in range(nAT - 1):
            for j in range(i + 1, nAT):
                dis = get_atom_distance(temp[i][1], temp[j][1])
                result += temp[i][0] * temp[j][0] / numpy.lib.scimath.power(dis, p=2)
        return float(result) / 100

    @staticmethod
    def calculate_radius_of_gyration(mol: Chem.Mol, charge_coordinates: List[List[float]]) -> float:
        """Calculate Radius of gyration.

        :param mol: molecule
        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        mol = Chem.AddHs(mol)
        temp = []
        atoms = list(mol.GetAtoms())
        for i, j in enumerate(charge_coordinates):
            temp.append([periodicTable.GetAtomicWeight(atoms[i].GetAtomicNum()),
                         [float(j[1]), float(j[2]), float(j[3])]])
        nAT = len(temp)
        com = Geometric._get_center_of_mass(temp)
        result = 0.0
        for i in range(nAT):
            dis = get_atom_distance(temp[i][1], com)
            result += temp[i][0] * numpy.lib.scimath.power(dis, p=2)
        return np.sqrt(float(result / ExactMolWt(mol)))

    @staticmethod
    def get_inertia_matrix(mol: Chem.Mol, charge_coordinates: List[List[float]]) -> np.array:
        """Get Inertia matrix based on atomic mass and optimized coordinates.

        :param mol: molecule
        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        mol = Chem.AddHs(mol)
        temp = []
        atoms = list(mol.GetAtoms())
        for i, j in enumerate(charge_coordinates):
            temp.append([periodicTable.GetAtomicWeight(atoms[i].GetAtomicNum()),
                         [float(j[1]), float(j[2]), float(j[3])]])
        nAT = len(temp)
        inertia_matrix = np.zeros((3, 3))
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
        inertia_matrix[0, 0] = res11
        inertia_matrix[1, 1] = res22
        inertia_matrix[2, 2] = res33
        inertia_matrix[0, 1] = res12
        inertia_matrix[0, 2] = res13
        inertia_matrix[1, 2] = res23
        inertia_matrix[1, 0] = res12
        inertia_matrix[2, 0] = res13
        inertia_matrix[2, 1] = res23
        return inertia_matrix

    @staticmethod
    def calculate_principal_moment_of_inertia(mol: Chem.Mol, charge_coordinates: List[List[float]]) -> dict:
        """Calculate X, Y and Z-principal geometric moments.

        derived from ADAPT developed by Jurs.
        :param mol: molecule
        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        inertia_matrix = Geometric.get_inertia_matrix(mol, charge_coordinates)
        ma = np.mean(inertia_matrix, axis=1)
        ms = np.std(inertia_matrix, axis=1, ddof=1)
        bb = np.ones((3, 1))
        inertia_matrix = (inertia_matrix - bb * ma.T) / (bb * ms.T)
        u, s, v = np.linalg.svd(inertia_matrix)
        res = {}
        res['IA'] = s[2]
        res['IB'] = s[1]
        res['IC'] = s[0]
        return res

    @staticmethod
    def calculate_ratio_pmi(mol: Chem.Mol, charge_coordinates: List[List[float]]) -> dict:
        """Calculate the ratio of principal moment of inertia.

        derived from ADAPT developed by Jurs.
        :param mol: molecule
        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        temp = Geometric.calculate_principal_moment_of_inertia(mol, charge_coordinates)
        res = {}
        res['IA/B'] = temp['IA'] / temp['IB']
        res['IA/C'] = temp['IA'] / temp['IC']
        res['IB/C'] = temp['IB'] / temp['IC']
        return res

    @staticmethod
    def calculate_harary_3D(charge_coordinates: List[List[float]]) -> float:
        """Calculate 3D-Harary index as the sum of all the reciprocal geometric distances.

        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        temp = []
        for i in charge_coordinates:
            temp.append([float(i[1]), float(i[2]), float(i[3])])
        distance_matrix = get_geometrical_distance_matrix(temp)
        nAT = len(temp)
        res = 0.0
        for i in range(nAT - 1):
            for j in range(i + 1, nAT):
                if distance_matrix[i, j] == 0:
                    cds = 0.0
                else:
                    cds = 1. / distance_matrix[i, j]
                res = res + cds
        return res

    @staticmethod
    def calculate_average_geometrical_distance_degree(charge_coordinates: List[List[float]]) -> float:
        """Calculate the average geometric distance degree (AGDD).

        This is the ratio between the sum of all geometric distance degrees and the atoms.
        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        temp = []
        for i in charge_coordinates:
            temp.append([float(i[1]), float(i[2]), float(i[3])])
        distance_matrix = get_geometrical_distance_matrix(temp)
        nAT = len(temp)
        res = sum(sum(distance_matrix)) / nAT
        return res

    @staticmethod
    def calculate_abs_eigenvalue_sum_on_geometric_matrix(charge_coordinates: List[List[float]]) -> float:
        """Calculate the absolute eigenvalue sum on geometry matrix (SEig).

        This is the sum of the absolute eigenvalues of the geometry matrix.
        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        temp = []
        for i in charge_coordinates:
            temp.append([float(i[1]), float(i[2]), float(i[3])])
        distance_matrix = get_geometrical_distance_matrix(temp)
        u, s, vt = np.linalg.svd(distance_matrix)
        return sum(abs(s))

    @staticmethod
    def calculate_span_r(mol: Chem.Mol, charge_coordinates: List[List[float]]) -> float:
        """Calculate the span R.

        This is defined as the radius of the smallest sphere,
        centred on the centre of mass, completely enclosing all atoms of a molecule.

        Arteca G.A. et al., Molecular Shape Descriptors in Reviews in
        Computational Chemistry - Vol. 9, K.B. Lipkowitz, D. Boyd (Eds.),
        VCH Publishers, New York (NY), pp. 191-253, 1991.
        :param mol: molecule
        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        mol = Chem.AddHs(mol)
        temp = []
        atoms = list(mol.GetAtoms())
        for i, j in enumerate(charge_coordinates):
            temp.append([periodicTable.GetAtomicWeight(atoms[i].GetAtomicNum()),
                         [float(j[1]), float(j[2]), float(j[3])]])
        masscenter = Geometric._get_center_of_mass(temp)
        res = []
        for i in temp:
            res.append(get_atom_distance(i[1], masscenter))
        return float(max(res))

    @staticmethod
    def calculate_average_span_r(mol: Chem.Mol, charge_coordinates: List[List[float]]) -> float:
        """Calculate the average span R (SPAM).

        This is the root square of the ratio of SPAN over the number of atoms.
        :param mol: molecule
        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        mol = Chem.AddHs(mol)
        temp = []
        atoms = list(mol.GetAtoms())
        for i, j in enumerate(charge_coordinates):
            temp.append([periodicTable.GetAtomicWeight(atoms[i].GetAtomicNum()),
                         [float(j[1]), float(j[2]), float(j[3])]])
        nAT = len(temp)
        masscenter = Geometric._get_center_of_mass(temp)
        res = []
        for i in temp:
            res.append(get_atom_distance(i[1], masscenter))
        return math.pow(float(max(res)) / nAT, 0.5)

    @staticmethod
    def calculate_molecular_eccentricity(mol: Chem.Mol, charge_coordinates: List[List[float]]) -> float:
        """Calculate molecular eccentricity.

        G.A. Arteca, Molecular Shape Descriptors in Reviews
        in Computational Chemistry - Vol. 9, K.B. Lipkowitz, D. Boyd (Eds.),
        VCH Publishers, New York (NY), pp. 191-253, 1991.
        :param mol: molecule
        :param charge_coordinates: Atomic coordinates and charges as read by chemopy.geo_opt.read_coordinates
        """
        inertia_matrix = Geometric.get_inertia_matrix(mol, charge_coordinates)
        u, s, v = np.linalg.svd(inertia_matrix)
        res1 = s[0]
        res3 = s[2]
        res = math.pow(res1 * res1 - res3 * res3, 1. / 2) / res1
        return res

    @staticmethod
    def get_all(mol: Chem.Mol, arc_file: str) -> dict:
        """Get all (20) geometrical descriptors.

        :param mol: the molecule
        :param arc_file: Path to MOPAC .arc file
        """
        res = {}
        charge_coordinates = read_coordinates(arc_file)
        res['W3DH'] = Geometric.calculate_3D_wiener_with_h(charge_coordinates)
        res['W3D'] = Geometric.calculate_3D_wiener_without_h(charge_coordinates)
        res['Petitj3D'] = Geometric.calculate_petitjean_3D_index(charge_coordinates)
        res['GeDi'] = Geometric.calculate_geometrical_diameter(charge_coordinates)
        res['TE1'] = Geometric.calculate_topo_electronic(charge_coordinates)
        res['grav1'] = Geometric.calculate_gravitational_3D1(mol, charge_coordinates) if mol is not None else np.nan
        # res['grav2'] = Geometric.CalculateGravitational3D2(mol, charge_coordinates)
        res['rygr'] = Geometric.calculate_radius_of_gyration(mol, charge_coordinates) if mol is not None else np.nan
        res['Harary3D'] = Geometric.calculate_harary_3D(charge_coordinates)
        res['AGDD'] = Geometric.calculate_average_geometrical_distance_degree(charge_coordinates)
        res['SEig'] = Geometric.calculate_abs_eigenvalue_sum_on_geometric_matrix(charge_coordinates)
        res['SPAN'] = Geometric.calculate_span_r(mol, charge_coordinates) if mol is not None else np.nan
        res['ASPAN'] = Geometric.calculate_average_span_r(mol, charge_coordinates) if mol is not None else np.nan
        res['MEcc'] = Geometric.calculate_molecular_eccentricity(mol, charge_coordinates) if mol is not None else np.nan
        res.update(Geometric.calculate_principal_moment_of_inertia(mol, charge_coordinates)) if mol is not None else np.nan
        res.update(Geometric.calculate_ratio_pmi(mol, charge_coordinates)) if mol is not None else np.nan
        return res
