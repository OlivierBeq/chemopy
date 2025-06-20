# -*- coding: utf-8 -*-


"""Whole Holistic Invariant Molecular (WHIM) descriptors."""

from typing import List, Tuple

from rdkit import Chem
import numpy as np

from .atom_property import get_relative_atomic_property
from .geo_opt import read_coordinates

class WHIM:
    """Whole Holistic Invariant Molecular descriptors."""

    @staticmethod
    def get_atom_coordinate_matrix(arc_file: str) -> Tuple[np.array, List[str]]:
        """Get atom coordinate matrix and elements.
    
        :param arc_file: Path to MOPAC .arc file
        :return: The 3D coordinates of the atoms of the optimized molecule
                 along with atomic symbols.
        """
        ChargeCoordinates = read_coordinates(arc_file)
        nAtom = len(ChargeCoordinates)
        CoordinateMatrix = np.zeros([nAtom, 3])
        AtomLabel = []
        for i, j in enumerate(ChargeCoordinates):
            CoordinateMatrix[i, :] = [j[1], j[2], j[3]]
            AtomLabel.append(j[0])
        return np.array(CoordinateMatrix), AtomLabel

    @staticmethod
    def x_pre_center(X: np.array) -> np.array:
        """Center the data matrix X."""
        Xdim = np.size(X, axis=0)
        Xmean = np.mean(X, axis=0)
        Xmean = np.array(Xmean)
        Xp = X - np.ones([Xdim, 1]) * Xmean
        return Xp

    @staticmethod
    def get_property_matrix(AtomLabel: List[str], proname: str = 'm') -> np.array:
        """Get the given property from an atomic symbol."""
        res = []
        for i in AtomLabel:
            res.append(get_relative_atomic_property(i, proname))
        return np.array(np.diag(res))

    @staticmethod
    def get_svd_eig(CoordinateMatrix: np.array, AtomLabel: List[str], proname: str = 'u') -> np.array:
        """Get singular values of the weighted covariance matrix."""
        nAtom, kc = CoordinateMatrix.shape
        if proname == 'u':
            weight = np.array(np.eye(nAtom))
        else:
            weight = WHIM.get_property_matrix(AtomLabel, proname)
        S = WHIM.x_pre_center(CoordinateMatrix)
        u, s, v = np.linalg.svd(S.T @ weight @ S / sum(np.diag(weight)))
        return s

    @staticmethod
    def get_whim1(CoordinateMatrix: np.array, AtomLabel: List[str], proname: str = 'u') -> float:
        """Get L1u WHIM descriptor."""
        s = WHIM.get_svd_eig(CoordinateMatrix, AtomLabel, proname)
        return s[0]

    @staticmethod
    def get_whim2(CoordinateMatrix: np.array, AtomLabel: List[str], proname: str = 'u') -> float:
        """Get L2u WHIM descriptor."""
        s = WHIM.get_svd_eig(CoordinateMatrix, AtomLabel, proname)
        return s[1]

    @staticmethod
    def get_whim3(CoordinateMatrix: np.array, AtomLabel: List[str], proname: str = 'u') -> float:
        """Get L3u WHIM descriptor."""
        s = WHIM.get_svd_eig(CoordinateMatrix, AtomLabel, proname)
        return s[2]

    @staticmethod
    def get_whim4(CoordinateMatrix: np.array, AtomLabel: List[str], proname: str = 'u') -> float:
        """Get Tu WHIM descriptor."""
        s = WHIM.get_svd_eig(CoordinateMatrix, AtomLabel, proname)
        T = sum(s)
        return T

    @staticmethod
    def get_whim5(CoordinateMatrix: np.array, AtomLabel: List[str], proname: str = 'u') -> float:
        """Get Au WHIM descriptor."""
        s = WHIM.get_svd_eig(CoordinateMatrix, AtomLabel, proname)
        A = s[0] * s[1] + s[0] * s[2] + s[1] * s[2]
        return A

    @staticmethod
    def get_whim6(CoordinateMatrix: np.array, AtomLabel: List[str], proname: str = 'u') -> float:
        """Get Vu WHIM descriptor."""
        s = WHIM.get_svd_eig(CoordinateMatrix, AtomLabel, proname)
        A = s[0] * s[1] + s[0] * s[2] + s[1] * s[2]
        T = sum(s)
        V = A + T + s[0] * s[1] * s[2]
        return V

    @staticmethod
    def get_whim7(CoordinateMatrix: np.array, AtomLabel: List[str], proname: str = 'u') -> float:
        """Get P1u WHIM descriptor."""
        s = WHIM.get_svd_eig(CoordinateMatrix, AtomLabel, proname)
        return s[0] / (s[0] + s[1] + s[2])

    @staticmethod
    def get_whim8(CoordinateMatrix: np.array, AtomLabel: List[str], proname: str = 'u'):
        """Get P2u WHIM descriptor."""
        s = WHIM.get_svd_eig(CoordinateMatrix, AtomLabel, proname)
        return s[1] / (s[0] + s[1] + s[2])

    @staticmethod
    def get_whim9(CoordinateMatrix: np.array, AtomLabel: List[str], proname: str = 'u') -> float:
        """Get Ku WHIM descriptor."""
        s = WHIM.get_svd_eig(CoordinateMatrix, AtomLabel, proname)
        res = 0.0
        for i in s:
            res = res + abs(i / sum(s) - 1 / 3.0)
        Ku = 3.0 / 4 * res
        return Ku

    @staticmethod
    def get_whim10(CoordinateMatrix: np.array, AtomLabel: List[str], proname: str = 'u') -> float:
        """Get E1u WHIM descriptor."""
        nAtom, kc = CoordinateMatrix.shape
        if proname == 'u':
            weight = np.array(np.eye(nAtom))
        else:
            weight = WHIM.get_property_matrix(AtomLabel, proname)
        S = WHIM.x_pre_center(CoordinateMatrix)
        u, s, v = np.linalg.svd(S.T @ weight @ S / sum(np.diag(weight)))
        res = np.power(s[0], 2) * nAtom / np.power(S * np.array(u[:, 0]).T, 4).sum()
        return float(res.real)

    @staticmethod
    def get_whim11(CoordinateMatrix: np.array, AtomLabel: List[str], proname: str = 'u') -> float:
        """Get E2u WHIM descriptor."""
        nAtom, kc = CoordinateMatrix.shape
        if proname == 'u':
            weight = np.array(np.eye(nAtom))
        else:
            weight = WHIM.get_property_matrix(AtomLabel, proname)
        S = WHIM.x_pre_center(CoordinateMatrix)
        u, s, v = np.linalg.svd(S.T @ weight @ S / sum(np.diag(weight)))
        res = np.power(s[1], 2) * nAtom / np.power(S * np.array(u[:, 1]).T, 4).sum()
        return float(res.real)

    @staticmethod
    def get_whim12(CoordinateMatrix: np.array, AtomLabel: List[str], proname: str = 'u') -> float:
        """Get E3u WHIM descriptor."""
        nAtom, kc = CoordinateMatrix.shape
        if proname == 'u':
            weight = np.array(np.eye(nAtom))
        else:
            weight = WHIM.get_property_matrix(AtomLabel, proname)
        S = WHIM.x_pre_center(CoordinateMatrix)
        u, s, v = np.linalg.svd(S.T @ weight @ S / sum(np.diag(weight)))
        res = np.power(s[2], 2) * nAtom / np.power(S * np.array(u[:, 2]).T, 4).sum()
        return float(res.real)

    @staticmethod
    def get_whim13(CoordinateMatrix: np.array, AtomLabel: List[str], proname: str = 'u') -> float:
        """Get Du WHIM descriptor."""
        c1 = WHIM.get_whim10(CoordinateMatrix, AtomLabel, proname)
        c2 = WHIM.get_whim11(CoordinateMatrix, AtomLabel, proname)
        c3 = WHIM.get_whim12(CoordinateMatrix, AtomLabel, proname)
        Du = c1 + c2 + c3
        return float(Du)

    @staticmethod
    def get_whim14(CoordinateMatrix: np.array, AtomLabel: List[str], proname: str = 'u') -> float:
        """Get P3u WHIM descriptor."""
        s = WHIM.get_svd_eig(CoordinateMatrix, AtomLabel, proname)
        return s[2] / (s[0] + s[1] + s[2])

    @staticmethod
    def get_whim_unweighted() -> dict:
        """Get all unweighted WHIM descriptors."""
        res = {}
        CoordinateMatrix, AtomLabel = WHIM.get_atom_coordinate_matrix()
        res['L1u'] = WHIM.get_whim1(CoordinateMatrix, AtomLabel, proname='u')
        res['L2u'] = WHIM.get_whim2(CoordinateMatrix, AtomLabel, proname='u')
        res['L3u'] = WHIM.get_whim3(CoordinateMatrix, AtomLabel, proname='u')
        res['P1u'] = WHIM.get_whim7(CoordinateMatrix, AtomLabel, proname='u')
        res['P2u'] = WHIM.get_whim8(CoordinateMatrix, AtomLabel, proname='u')
        res['P3u'] = WHIM.get_whim14(CoordinateMatrix, AtomLabel, proname='u')
        res['E1u'] = WHIM.get_whim10(CoordinateMatrix, AtomLabel, proname='u')
        res['E2u'] = WHIM.get_whim11(CoordinateMatrix, AtomLabel, proname='u')
        res['E3u'] = WHIM.get_whim12(CoordinateMatrix, AtomLabel, proname='u')
        res['Tu'] = WHIM.get_whim4(CoordinateMatrix, AtomLabel, proname='u')
        res['Au'] = WHIM.get_whim5(CoordinateMatrix, AtomLabel, proname='u')
        res['Vu'] = WHIM.get_whim6(CoordinateMatrix, AtomLabel, proname='u')
        res['Ku'] = WHIM.get_whim9(CoordinateMatrix, AtomLabel, proname='u')
        res['Du'] = WHIM.get_whim13(CoordinateMatrix, AtomLabel, proname='u')
        return res

    @staticmethod
    def get_whim_mass() -> dict:
        """Get all WHIM descriptors based on atomic mass."""
        res = {}
        CoordinateMatrix, AtomLabel = WHIM.get_atom_coordinate_matrix()
        res['L1m'] = WHIM.get_whim1(CoordinateMatrix, AtomLabel, proname='m')
        res['L2m'] = WHIM.get_whim2(CoordinateMatrix, AtomLabel, proname='m')
        res['L3m'] = WHIM.get_whim3(CoordinateMatrix, AtomLabel, proname='m')
        res['Tm'] = WHIM.get_whim4(CoordinateMatrix, AtomLabel, proname='m')
        res['Am'] = WHIM.get_whim5(CoordinateMatrix, AtomLabel, proname='m')
        res['Vm'] = WHIM.get_whim6(CoordinateMatrix, AtomLabel, proname='m')
        res['P1m'] = WHIM.get_whim7(CoordinateMatrix, AtomLabel, proname='m')
        res['P2m'] = WHIM.get_whim8(CoordinateMatrix, AtomLabel, proname='m')
        res['Km'] = WHIM.get_whim9(CoordinateMatrix, AtomLabel, proname='m')
        res['E1m'] = WHIM.get_whim10(CoordinateMatrix, AtomLabel, proname='m')
        res['E2m'] = WHIM.get_whim11(CoordinateMatrix, AtomLabel, proname='m')
        res['E3m'] = WHIM.get_whim12(CoordinateMatrix, AtomLabel, proname='m')
        res['Dm'] = WHIM.get_whim13(CoordinateMatrix, AtomLabel, proname='m')
        res['P3m'] = WHIM.get_whim14(CoordinateMatrix, AtomLabel, proname='m')
        return res

    @staticmethod
    def get_whim_sanderson_electronegativity() -> dict:
        """Get all WHIM descriptors based on Sanderson electronegativity."""
        res = {}
        CoordinateMatrix, AtomLabel = WHIM.get_atom_coordinate_matrix()
        res['L1e'] = WHIM.get_whim1(CoordinateMatrix, AtomLabel, proname='En')
        res['L2e'] = WHIM.get_whim2(CoordinateMatrix, AtomLabel, proname='En')
        res['L3e'] = WHIM.get_whim3(CoordinateMatrix, AtomLabel, proname='En')
        res['P1e'] = WHIM.get_whim7(CoordinateMatrix, AtomLabel, proname='En')
        res['P2e'] = WHIM.get_whim8(CoordinateMatrix, AtomLabel, proname='En')
        res['P3e'] = WHIM.get_whim14(CoordinateMatrix, AtomLabel, proname='En')
        res['E1e'] = WHIM.get_whim10(CoordinateMatrix, AtomLabel, proname='En')
        res['E2e'] = WHIM.get_whim11(CoordinateMatrix, AtomLabel, proname='En')
        res['E3e'] = WHIM.get_whim12(CoordinateMatrix, AtomLabel, proname='En')
        res['Te'] = WHIM.get_whim4(CoordinateMatrix, AtomLabel, proname='En')
        res['Ae'] = WHIM.get_whim5(CoordinateMatrix, AtomLabel, proname='En')
        res['Ve'] = WHIM.get_whim6(CoordinateMatrix, AtomLabel, proname='En')
        res['Ke'] = WHIM.get_whim9(CoordinateMatrix, AtomLabel, proname='En')
        res['De'] = WHIM.get_whim13(CoordinateMatrix, AtomLabel, proname='En')
        return res

    @staticmethod
    def get_whim_vdw_volume() -> dict:
        """Get all WHIM descriptors based on vdW volume."""
        res = {}
        CoordinateMatrix, AtomLabel = WHIM.get_atom_coordinate_matrix()
        res['L1v'] = WHIM.get_whim1(CoordinateMatrix, AtomLabel, proname='V')
        res['L2v'] = WHIM.get_whim2(CoordinateMatrix, AtomLabel, proname='V')
        res['L3v'] = WHIM.get_whim3(CoordinateMatrix, AtomLabel, proname='V')
        res['P1v'] = WHIM.get_whim7(CoordinateMatrix, AtomLabel, proname='V')
        res['P2v'] = WHIM.get_whim8(CoordinateMatrix, AtomLabel, proname='V')
        res['P3v'] = WHIM.get_whim14(CoordinateMatrix, AtomLabel, proname='V')
        res['E1v'] = WHIM.get_whim10(CoordinateMatrix, AtomLabel, proname='V')
        res['E2v'] = WHIM.get_whim11(CoordinateMatrix, AtomLabel, proname='V')
        res['E3v'] = WHIM.get_whim12(CoordinateMatrix, AtomLabel, proname='V')
        res['Tv'] = WHIM.get_whim4(CoordinateMatrix, AtomLabel, proname='V')
        res['Av'] = WHIM.get_whim5(CoordinateMatrix, AtomLabel, proname='V')
        res['Vv'] = WHIM.get_whim6(CoordinateMatrix, AtomLabel, proname='V')
        res['Kv'] = WHIM.get_whim9(CoordinateMatrix, AtomLabel, proname='V')
        res['Dv'] = WHIM.get_whim13(CoordinateMatrix, AtomLabel, proname='V')
        return res

    @staticmethod
    def get_whim_polarizability() -> dict:
        """Get all WHIM descriptors based on polarizability."""
        res = {}
        CoordinateMatrix, AtomLabel = WHIM.get_atom_coordinate_matrix()
        res['L1p'] = WHIM.get_whim1(CoordinateMatrix, AtomLabel, proname='alapha')
        res['L2p'] = WHIM.get_whim2(CoordinateMatrix, AtomLabel, proname='alapha')
        res['L3p'] = WHIM.get_whim3(CoordinateMatrix, AtomLabel, proname='alapha')
        res['P1p'] = WHIM.get_whim7(CoordinateMatrix, AtomLabel, proname='alapha')
        res['P2p'] = WHIM.get_whim8(CoordinateMatrix, AtomLabel, proname='alapha')
        res['P3p'] = WHIM.get_whim14(CoordinateMatrix, AtomLabel, proname='alapha')
        res['E1p'] = WHIM.get_whim10(CoordinateMatrix, AtomLabel, proname='alapha')
        res['E2p'] = WHIM.get_whim11(CoordinateMatrix, AtomLabel, proname='alapha')
        res['E3p'] = WHIM.get_whim12(CoordinateMatrix, AtomLabel, proname='alapha')
        res['Tp'] = WHIM.get_whim4(CoordinateMatrix, AtomLabel, proname='alapha')
        res['Ap'] = WHIM.get_whim5(CoordinateMatrix, AtomLabel, proname='alapha')
        res['Vp'] = WHIM.get_whim6(CoordinateMatrix, AtomLabel, proname='alapha')
        res['Kp'] = WHIM.get_whim9(CoordinateMatrix, AtomLabel, proname='alapha')
        res['Dp'] = WHIM.get_whim13(CoordinateMatrix, AtomLabel, proname='alapha')
        return res

    @staticmethod
    def get_all(mol: Chem.Mol, arc_file: str) -> dict:
        """Get all (70) WHIM descriptors.
    
        :param dir_: Path to directory containing MOPAC .arc file
        """
        res = {}
        CoordinateMatrix, AtomLabel = WHIM.get_atom_coordinate_matrix(arc_file)
        res['L1u'] = WHIM.get_whim1(CoordinateMatrix, AtomLabel, proname='u')
        res['L2u'] = WHIM.get_whim2(CoordinateMatrix, AtomLabel, proname='u')
        res['L3u'] = WHIM.get_whim3(CoordinateMatrix, AtomLabel, proname='u')
        res['P1u'] = WHIM.get_whim7(CoordinateMatrix, AtomLabel, proname='u')
        res['P2u'] = WHIM.get_whim8(CoordinateMatrix, AtomLabel, proname='u')
        res['P3u'] = WHIM.get_whim14(CoordinateMatrix, AtomLabel, proname='u')
        res['E1u'] = WHIM.get_whim10(CoordinateMatrix, AtomLabel, proname='u')
        res['E2u'] = WHIM.get_whim11(CoordinateMatrix, AtomLabel, proname='u')
        res['E3u'] = WHIM.get_whim12(CoordinateMatrix, AtomLabel, proname='u')
        res['Tu'] = WHIM.get_whim4(CoordinateMatrix, AtomLabel, proname='u')
        res['Au'] = WHIM.get_whim5(CoordinateMatrix, AtomLabel, proname='u')
        res['Vu'] = WHIM.get_whim6(CoordinateMatrix, AtomLabel, proname='u')
        res['Ku'] = WHIM.get_whim9(CoordinateMatrix, AtomLabel, proname='u')
        res['Du'] = WHIM.get_whim13(CoordinateMatrix, AtomLabel, proname='u')
        res['L1m'] = WHIM.get_whim1(CoordinateMatrix, AtomLabel, proname='m')
        res['L2m'] = WHIM.get_whim2(CoordinateMatrix, AtomLabel, proname='m')
        res['L3m'] = WHIM.get_whim3(CoordinateMatrix, AtomLabel, proname='m')
        res['P1m'] = WHIM.get_whim7(CoordinateMatrix, AtomLabel, proname='m')
        res['P2m'] = WHIM.get_whim8(CoordinateMatrix, AtomLabel, proname='m')
        res['P3m'] = WHIM.get_whim14(CoordinateMatrix, AtomLabel, proname='m')
        res['E1m'] = WHIM.get_whim10(CoordinateMatrix, AtomLabel, proname='m')
        res['E2m'] = WHIM.get_whim11(CoordinateMatrix, AtomLabel, proname='m')
        res['E3m'] = WHIM.get_whim12(CoordinateMatrix, AtomLabel, proname='m')
        res['Tm'] = WHIM.get_whim4(CoordinateMatrix, AtomLabel, proname='m')
        res['Am'] = WHIM.get_whim5(CoordinateMatrix, AtomLabel, proname='m')
        res['Vm'] = WHIM.get_whim6(CoordinateMatrix, AtomLabel, proname='m')
        res['Km'] = WHIM.get_whim9(CoordinateMatrix, AtomLabel, proname='m')
        res['Dm'] = WHIM.get_whim13(CoordinateMatrix, AtomLabel, proname='m')
        res['L1e'] = WHIM.get_whim1(CoordinateMatrix, AtomLabel, proname='En')
        res['L2e'] = WHIM.get_whim2(CoordinateMatrix, AtomLabel, proname='En')
        res['L3e'] = WHIM.get_whim3(CoordinateMatrix, AtomLabel, proname='En')
        res['P1e'] = WHIM.get_whim7(CoordinateMatrix, AtomLabel, proname='En')
        res['P2e'] = WHIM.get_whim8(CoordinateMatrix, AtomLabel, proname='En')
        res['P3e'] = WHIM.get_whim14(CoordinateMatrix, AtomLabel, proname='En')
        res['E1e'] = WHIM.get_whim10(CoordinateMatrix, AtomLabel, proname='En')
        res['E2e'] = WHIM.get_whim11(CoordinateMatrix, AtomLabel, proname='En')
        res['E3e'] = WHIM.get_whim12(CoordinateMatrix, AtomLabel, proname='En')
        res['Te'] = WHIM.get_whim4(CoordinateMatrix, AtomLabel, proname='En')
        res['Ae'] = WHIM.get_whim5(CoordinateMatrix, AtomLabel, proname='En')
        res['Ve'] = WHIM.get_whim6(CoordinateMatrix, AtomLabel, proname='En')
        res['Ke'] = WHIM.get_whim9(CoordinateMatrix, AtomLabel, proname='En')
        res['De'] = WHIM.get_whim13(CoordinateMatrix, AtomLabel, proname='En')
        res['L1v'] = WHIM.get_whim1(CoordinateMatrix, AtomLabel, proname='V')
        res['L2v'] = WHIM.get_whim2(CoordinateMatrix, AtomLabel, proname='V')
        res['L3v'] = WHIM.get_whim3(CoordinateMatrix, AtomLabel, proname='V')
        res['P1v'] = WHIM.get_whim7(CoordinateMatrix, AtomLabel, proname='V')
        res['P2v'] = WHIM.get_whim8(CoordinateMatrix, AtomLabel, proname='V')
        res['P3v'] = WHIM.get_whim14(CoordinateMatrix, AtomLabel, proname='V')
        res['E1v'] = WHIM.get_whim10(CoordinateMatrix, AtomLabel, proname='V')
        res['E2v'] = WHIM.get_whim11(CoordinateMatrix, AtomLabel, proname='V')
        res['E3v'] = WHIM.get_whim12(CoordinateMatrix, AtomLabel, proname='V')
        res['Tv'] = WHIM.get_whim4(CoordinateMatrix, AtomLabel, proname='V')
        res['Av'] = WHIM.get_whim5(CoordinateMatrix, AtomLabel, proname='V')
        res['Vv'] = WHIM.get_whim6(CoordinateMatrix, AtomLabel, proname='V')
        res['Kv'] = WHIM.get_whim9(CoordinateMatrix, AtomLabel, proname='V')
        res['Dv'] = WHIM.get_whim13(CoordinateMatrix, AtomLabel, proname='V')
        res['L1p'] = WHIM.get_whim1(CoordinateMatrix, AtomLabel, proname='alapha')
        res['L2p'] = WHIM.get_whim2(CoordinateMatrix, AtomLabel, proname='alapha')
        res['L3p'] = WHIM.get_whim3(CoordinateMatrix, AtomLabel, proname='alapha')
        res['P1p'] = WHIM.get_whim7(CoordinateMatrix, AtomLabel, proname='alapha')
        res['P2p'] = WHIM.get_whim8(CoordinateMatrix, AtomLabel, proname='alapha')
        res['P3p'] = WHIM.get_whim14(CoordinateMatrix, AtomLabel, proname='alapha')
        res['E1p'] = WHIM.get_whim10(CoordinateMatrix, AtomLabel, proname='alapha')
        res['E2p'] = WHIM.get_whim11(CoordinateMatrix, AtomLabel, proname='alapha')
        res['E3p'] = WHIM.get_whim12(CoordinateMatrix, AtomLabel, proname='alapha')
        res['Tp'] = WHIM.get_whim4(CoordinateMatrix, AtomLabel, proname='alapha')
        res['Ap'] = WHIM.get_whim5(CoordinateMatrix, AtomLabel, proname='alapha')
        res['Vp'] = WHIM.get_whim6(CoordinateMatrix, AtomLabel, proname='alapha')
        res['Kp'] = WHIM.get_whim9(CoordinateMatrix, AtomLabel, proname='alapha')
        res['Dp'] = WHIM.get_whim13(CoordinateMatrix, AtomLabel, proname='alapha')
        return res
