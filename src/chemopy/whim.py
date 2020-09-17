# -*- coding: utf-8 -*-


"""Whole Holistic Invariant Molecular (WHIM) descriptors."""

from typing import List, Tuple

import numpy as np
import scipy.linalg

from pychem.AtomProperty import GetRelativeAtomicProperty
from pychem.GeoOpt import _ReadCoordinates


def GetAtomCoordinateMatrix(arc_file: str) -> Tuple[np.matrix, List[str]]:
    """Get atom coordinate matrix and elements.

    :param arc_file: Path to MOPAC .arc file
    :return: The 3D coordinates of the atoms of the optimized molecule
             along with atomic symbols.
    """
    ChargeCoordinates = _ReadCoordinates(arc_file)
    nAtom = len(ChargeCoordinates)
    CoordinateMatrix = np.zeros([nAtom, 3])
    AtomLabel = []
    for i, j in enumerate(ChargeCoordinates):
        CoordinateMatrix[i, :] = [j[1], j[2], j[3]]
        AtomLabel.append(j[0])
    return np.matrix(CoordinateMatrix), AtomLabel


def XPreCenter(X: scipy.matrix) -> np.matrix:
    """Center the data matrix X."""
    Xdim = np.size(X, axis=0)
    Xmean = np.mean(X, axis=0)
    Xmean = np.matrix(Xmean)
    Xp = X - np.ones([Xdim, 1]) * Xmean
    return Xp


def GetPropertyMatrix(AtomLabel: List[str], proname: str = 'm') -> np.matrix:
    """Get the given property from an atomic symbol."""
    res = []
    for i in AtomLabel:
        res.append(GetRelativeAtomicProperty(i, proname))
    return np.matrix(np.diag(res))


def GetSVDEig(CoordinateMatrix: np.matrix, AtomLabel: List[str], proname: str = 'u') -> np.matrix:
    """Get singular values of the weighted covariance matrix."""
    nAtom, kc = CoordinateMatrix.shape
    if proname == 'u':
        weight = np.matrix(np.eye(nAtom))
    else:
        weight = GetPropertyMatrix(AtomLabel, proname)
    S = XPreCenter(CoordinateMatrix)
    u, s, v = np.linalg.svd(S.T * weight * S / sum(np.diag(weight)))
    return s


def GetWHIM1(CoordinateMatrix: np.matrix, AtomLabel: List[str], proname: str = 'u') -> float:
    """Get L1u WHIM descriptor."""
    s = GetSVDEig(CoordinateMatrix, AtomLabel, proname)
    return round(s[0], 3)


def GetWHIM2(CoordinateMatrix: np.matrix, AtomLabel: List[str], proname: str = 'u') -> float:
    """Get L2u WHIM descriptor."""
    s = GetSVDEig(CoordinateMatrix, AtomLabel, proname)
    return round(s[1], 3)


def GetWHIM3(CoordinateMatrix: np.matrix, AtomLabel: List[str], proname: str = 'u') -> float:
    """Get L3u WHIM descriptor."""
    s = GetSVDEig(CoordinateMatrix, AtomLabel, proname)
    return round(s[2], 3)


def GetWHIM4(CoordinateMatrix: np.matrix, AtomLabel: List[str], proname: str = 'u') -> float:
    """Get Tu WHIM descriptor."""
    s = GetSVDEig(CoordinateMatrix, AtomLabel, proname)
    T = round(sum(s), 3)
    return T


def GetWHIM5(CoordinateMatrix: np.matrix, AtomLabel: List[str], proname: str = 'u') -> float:
    """Get Au WHIM descriptor."""
    s = GetSVDEig(CoordinateMatrix, AtomLabel, proname)
    A = s[0] * s[1] + s[0] * s[2] + s[1] * s[2]
    return round(A, 3)


def GetWHIM6(CoordinateMatrix: np.matrix, AtomLabel: List[str], proname: str = 'u') -> float:
    """Get Vu WHIM descriptor."""
    s = GetSVDEig(CoordinateMatrix, AtomLabel, proname)
    A = s[0] * s[1] + s[0] * s[2] + s[1] * s[2]
    T = sum(s)
    V = A + T + s[0] * s[1] * s[2]
    return round(V, 3)


def GetWHIM7(CoordinateMatrix: np.matrix, AtomLabel: List[str], proname: str = 'u') -> float:
    """Get P1u WHIM descriptor."""
    s = GetSVDEig(CoordinateMatrix, AtomLabel, proname)
    return round(s[0] / (s[0] + s[1] + s[2]), 3)


def GetWHIM8(CoordinateMatrix: np.matrix, AtomLabel: List[str], proname: str = 'u'):
    """Get P2u WHIM descriptor."""
    s = GetSVDEig(CoordinateMatrix, AtomLabel, proname)
    return round(s[1] / (s[0] + s[1] + s[2]), 3)


def GetWHIM9(CoordinateMatrix: np.matrix, AtomLabel: List[str], proname: str = 'u') -> float:
    """Get Ku WHIM descriptor."""
    s = GetSVDEig(CoordinateMatrix, AtomLabel, proname)
    res = 0.0
    for i in s:
        res = res + abs(i / sum(s) - 1 / 3.0)
    Ku = 3.0 / 4 * res
    return round(Ku, 3)


def GetWHIM10(CoordinateMatrix: np.matrix, AtomLabel: List[str], proname: str = 'u') -> float:
    """Get E1u WHIM descriptor."""
    nAtom, kc = CoordinateMatrix.shape
    if proname == 'u':
        weight = np.matrix(np.eye(nAtom))
    else:
        weight = GetPropertyMatrix(AtomLabel, proname)
    S = XPreCenter(CoordinateMatrix)
    u, s, v = scipy.linalg.svd(S.T * weight * S / sum(np.diag(weight)))
    res = np.power(s[0], 2) * nAtom / sum(np.power(S * np.matrix(u[:, 0]).T, 4))
    return round(float(res.real), 3)


def GetWHIM11(CoordinateMatrix: np.matrix, AtomLabel: List[str], proname: str = 'u') -> float:
    """Get E2u WHIM descriptor."""
    nAtom, kc = CoordinateMatrix.shape
    if proname == 'u':
        weight = np.matrix(np.eye(nAtom))
    else:
        weight = GetPropertyMatrix(AtomLabel, proname)
    S = XPreCenter(CoordinateMatrix)
    u, s, v = scipy.linalg.svd(S.T * weight * S / sum(np.diag(weight)))
    res = np.power(s[1], 2) * nAtom / sum(np.power(S * np.matrix(u[:, 1]).T, 4))
    return round(float(res.real), 3)


def GetWHIM12(CoordinateMatrix: np.matrix, AtomLabel: List[str], proname: str = 'u') -> float:
    """Get E3u WHIM descriptor."""
    nAtom, kc = CoordinateMatrix.shape
    if proname == 'u':
        weight = np.matrix(np.eye(nAtom))
    else:
        weight = GetPropertyMatrix(AtomLabel, proname)
    S = XPreCenter(CoordinateMatrix)
    u, s, v = scipy.linalg.svd(S.T * weight * S / sum(np.diag(weight)))
    res = np.power(s[2], 2) * nAtom / sum(np.power(S * np.matrix(u[:, 2]).T, 4))
    return round(float(res.real), 3)


def GetWHIM13(CoordinateMatrix: np.matrix, AtomLabel: List[str], proname: str = 'u') -> float:
    """Get Du WHIM descriptor."""
    c1 = GetWHIM10(CoordinateMatrix, AtomLabel, proname)
    c2 = GetWHIM11(CoordinateMatrix, AtomLabel, proname)
    c3 = GetWHIM12(CoordinateMatrix, AtomLabel, proname)
    Du = c1 + c2 + c3
    return round(float(Du), 3)


def GetWHIM14(CoordinateMatrix: np.matrix, AtomLabel: List[str], proname: str = 'u') -> float:
    """Get P3u WHIM descriptor."""
    s = GetSVDEig(CoordinateMatrix, AtomLabel, proname)
    return round(s[2] / (s[0] + s[1] + s[2]), 3)


def GetWHIMUnweighted() -> dict:
    """Get all unweighted WHIM descriptors."""
    res = {}
    CoordinateMatrix, AtomLabel = GetAtomCoordinateMatrix()
    res['L1u'] = GetWHIM1(CoordinateMatrix, AtomLabel, proname='u')
    res['L2u'] = GetWHIM2(CoordinateMatrix, AtomLabel, proname='u')
    res['L3u'] = GetWHIM3(CoordinateMatrix, AtomLabel, proname='u')
    res['Tu'] = GetWHIM4(CoordinateMatrix, AtomLabel, proname='u')
    res['Au'] = GetWHIM5(CoordinateMatrix, AtomLabel, proname='u')
    res['Vu'] = GetWHIM6(CoordinateMatrix, AtomLabel, proname='u')
    res['P1u'] = GetWHIM7(CoordinateMatrix, AtomLabel, proname='u')
    res['P2u'] = GetWHIM8(CoordinateMatrix, AtomLabel, proname='u')
    res['Ku'] = GetWHIM9(CoordinateMatrix, AtomLabel, proname='u')
    res['E1u'] = GetWHIM10(CoordinateMatrix, AtomLabel, proname='u')
    res['E2u'] = GetWHIM11(CoordinateMatrix, AtomLabel, proname='u')
    res['E3u'] = GetWHIM12(CoordinateMatrix, AtomLabel, proname='u')
    res['Du'] = GetWHIM13(CoordinateMatrix, AtomLabel, proname='u')
    res['P3u'] = GetWHIM14(CoordinateMatrix, AtomLabel, proname='u')
    return res


def GetWHIMMass() -> dict:
    """Get all WHIM descriptors based on atomic mass."""
    res = {}
    CoordinateMatrix, AtomLabel = GetAtomCoordinateMatrix()
    res['L1m'] = GetWHIM1(CoordinateMatrix, AtomLabel, proname='m')
    res['L2m'] = GetWHIM2(CoordinateMatrix, AtomLabel, proname='m')
    res['L3m'] = GetWHIM3(CoordinateMatrix, AtomLabel, proname='m')
    res['Tm'] = GetWHIM4(CoordinateMatrix, AtomLabel, proname='m')
    res['Am'] = GetWHIM5(CoordinateMatrix, AtomLabel, proname='m')
    res['Vm'] = GetWHIM6(CoordinateMatrix, AtomLabel, proname='m')
    res['P1m'] = GetWHIM7(CoordinateMatrix, AtomLabel, proname='m')
    res['P2m'] = GetWHIM8(CoordinateMatrix, AtomLabel, proname='m')
    res['Km'] = GetWHIM9(CoordinateMatrix, AtomLabel, proname='m')
    res['E1m'] = GetWHIM10(CoordinateMatrix, AtomLabel, proname='m')
    res['E2m'] = GetWHIM11(CoordinateMatrix, AtomLabel, proname='m')
    res['E3m'] = GetWHIM12(CoordinateMatrix, AtomLabel, proname='m')
    res['Dm'] = GetWHIM13(CoordinateMatrix, AtomLabel, proname='m')
    res['P3m'] = GetWHIM14(CoordinateMatrix, AtomLabel, proname='m')
    return res


def GetWHIMSandersonElectronegativity() -> dict:
    """Get all WHIM descriptors based on Sanderson electronegativity."""
    res = {}
    CoordinateMatrix, AtomLabel = GetAtomCoordinateMatrix()
    res['L1e'] = GetWHIM1(CoordinateMatrix, AtomLabel, proname='En')
    res['L2e'] = GetWHIM2(CoordinateMatrix, AtomLabel, proname='En')
    res['L3e'] = GetWHIM3(CoordinateMatrix, AtomLabel, proname='En')
    res['Te'] = GetWHIM4(CoordinateMatrix, AtomLabel, proname='En')
    res['Ae'] = GetWHIM5(CoordinateMatrix, AtomLabel, proname='En')
    res['Ve'] = GetWHIM6(CoordinateMatrix, AtomLabel, proname='En')
    res['P1e'] = GetWHIM7(CoordinateMatrix, AtomLabel, proname='En')
    res['P2e'] = GetWHIM8(CoordinateMatrix, AtomLabel, proname='En')
    res['Ke'] = GetWHIM9(CoordinateMatrix, AtomLabel, proname='En')
    res['E1e'] = GetWHIM10(CoordinateMatrix, AtomLabel, proname='En')
    res['E2e'] = GetWHIM11(CoordinateMatrix, AtomLabel, proname='En')
    res['E3e'] = GetWHIM12(CoordinateMatrix, AtomLabel, proname='En')
    res['De'] = GetWHIM13(CoordinateMatrix, AtomLabel, proname='En')
    res['P3e'] = GetWHIM14(CoordinateMatrix, AtomLabel, proname='En')
    return res


def GetWHIMVDWVolume() -> dict:
    """Get all WHIM descriptors based on vdW volume."""
    res = {}
    CoordinateMatrix, AtomLabel = GetAtomCoordinateMatrix()
    res['L1v'] = GetWHIM1(CoordinateMatrix, AtomLabel, proname='V')
    res['L2v'] = GetWHIM2(CoordinateMatrix, AtomLabel, proname='V')
    res['L3v'] = GetWHIM3(CoordinateMatrix, AtomLabel, proname='V')
    res['Tv'] = GetWHIM4(CoordinateMatrix, AtomLabel, proname='V')
    res['Av'] = GetWHIM5(CoordinateMatrix, AtomLabel, proname='V')
    res['Vv'] = GetWHIM6(CoordinateMatrix, AtomLabel, proname='V')
    res['P1v'] = GetWHIM7(CoordinateMatrix, AtomLabel, proname='V')
    res['P2v'] = GetWHIM8(CoordinateMatrix, AtomLabel, proname='V')
    res['Kv'] = GetWHIM9(CoordinateMatrix, AtomLabel, proname='V')
    res['E1v'] = GetWHIM10(CoordinateMatrix, AtomLabel, proname='V')
    res['E2v'] = GetWHIM11(CoordinateMatrix, AtomLabel, proname='V')
    res['E3v'] = GetWHIM12(CoordinateMatrix, AtomLabel, proname='V')
    res['Dv'] = GetWHIM13(CoordinateMatrix, AtomLabel, proname='V')
    res['P3v'] = GetWHIM14(CoordinateMatrix, AtomLabel, proname='V')
    return res


def GetWHIMPolarizability() -> dict:
    """Get all WHIM descriptors based on polarizability."""
    res = {}
    CoordinateMatrix, AtomLabel = GetAtomCoordinateMatrix()
    res['L1p'] = GetWHIM1(CoordinateMatrix, AtomLabel, proname='alapha')
    res['L2p'] = GetWHIM2(CoordinateMatrix, AtomLabel, proname='alapha')
    res['L3p'] = GetWHIM3(CoordinateMatrix, AtomLabel, proname='alapha')
    res['Tp'] = GetWHIM4(CoordinateMatrix, AtomLabel, proname='alapha')
    res['Ap'] = GetWHIM5(CoordinateMatrix, AtomLabel, proname='alapha')
    res['Vp'] = GetWHIM6(CoordinateMatrix, AtomLabel, proname='alapha')
    res['P1p'] = GetWHIM7(CoordinateMatrix, AtomLabel, proname='alapha')
    res['P2p'] = GetWHIM8(CoordinateMatrix, AtomLabel, proname='alapha')
    res['Kp'] = GetWHIM9(CoordinateMatrix, AtomLabel, proname='alapha')
    res['E1p'] = GetWHIM10(CoordinateMatrix, AtomLabel, proname='alapha')
    res['E2p'] = GetWHIM11(CoordinateMatrix, AtomLabel, proname='alapha')
    res['E3p'] = GetWHIM12(CoordinateMatrix, AtomLabel, proname='alapha')
    res['Dp'] = GetWHIM13(CoordinateMatrix, AtomLabel, proname='alapha')
    res['P3p'] = GetWHIM14(CoordinateMatrix, AtomLabel, proname='alapha')
    return res


def GetWHIM(dir_: str) -> dict:
    """Get all (70) WHIM descriptors.

    :param dir_: Path to dir_ectory containing MOPAC .arc file
    """
    res = {}
    CoordinateMatrix, AtomLabel = GetAtomCoordinateMatrix(dir_)
    res['L1u'] = GetWHIM1(CoordinateMatrix, AtomLabel, proname='u')
    res['L2u'] = GetWHIM2(CoordinateMatrix, AtomLabel, proname='u')
    res['L3u'] = GetWHIM3(CoordinateMatrix, AtomLabel, proname='u')
    res['Tu'] = GetWHIM4(CoordinateMatrix, AtomLabel, proname='u')
    res['Au'] = GetWHIM5(CoordinateMatrix, AtomLabel, proname='u')
    res['Vu'] = GetWHIM6(CoordinateMatrix, AtomLabel, proname='u')
    res['P1u'] = GetWHIM7(CoordinateMatrix, AtomLabel, proname='u')
    res['P2u'] = GetWHIM8(CoordinateMatrix, AtomLabel, proname='u')
    res['Ku'] = GetWHIM9(CoordinateMatrix, AtomLabel, proname='u')
    res['E1u'] = GetWHIM10(CoordinateMatrix, AtomLabel, proname='u')
    res['E2u'] = GetWHIM11(CoordinateMatrix, AtomLabel, proname='u')
    res['E3u'] = GetWHIM12(CoordinateMatrix, AtomLabel, proname='u')
    res['Du'] = GetWHIM13(CoordinateMatrix, AtomLabel, proname='u')
    res['L1m'] = GetWHIM1(CoordinateMatrix, AtomLabel, proname='m')
    res['L2m'] = GetWHIM2(CoordinateMatrix, AtomLabel, proname='m')
    res['L3m'] = GetWHIM3(CoordinateMatrix, AtomLabel, proname='m')
    res['Tm'] = GetWHIM4(CoordinateMatrix, AtomLabel, proname='m')
    res['Am'] = GetWHIM5(CoordinateMatrix, AtomLabel, proname='m')
    res['Vm'] = GetWHIM6(CoordinateMatrix, AtomLabel, proname='m')
    res['P1m'] = GetWHIM7(CoordinateMatrix, AtomLabel, proname='m')
    res['P2m'] = GetWHIM8(CoordinateMatrix, AtomLabel, proname='m')
    res['Km'] = GetWHIM9(CoordinateMatrix, AtomLabel, proname='m')
    res['E1m'] = GetWHIM10(CoordinateMatrix, AtomLabel, proname='m')
    res['E2m'] = GetWHIM11(CoordinateMatrix, AtomLabel, proname='m')
    res['E3m'] = GetWHIM12(CoordinateMatrix, AtomLabel, proname='m')
    res['Dm'] = GetWHIM13(CoordinateMatrix, AtomLabel, proname='m')
    res['L1e'] = GetWHIM1(CoordinateMatrix, AtomLabel, proname='En')
    res['L2e'] = GetWHIM2(CoordinateMatrix, AtomLabel, proname='En')
    res['L3e'] = GetWHIM3(CoordinateMatrix, AtomLabel, proname='En')
    res['Te'] = GetWHIM4(CoordinateMatrix, AtomLabel, proname='En')
    res['Ae'] = GetWHIM5(CoordinateMatrix, AtomLabel, proname='En')
    res['Ve'] = GetWHIM6(CoordinateMatrix, AtomLabel, proname='En')
    res['P1e'] = GetWHIM7(CoordinateMatrix, AtomLabel, proname='En')
    res['P2e'] = GetWHIM8(CoordinateMatrix, AtomLabel, proname='En')
    res['Ke'] = GetWHIM9(CoordinateMatrix, AtomLabel, proname='En')
    res['E1e'] = GetWHIM10(CoordinateMatrix, AtomLabel, proname='En')
    res['E2e'] = GetWHIM11(CoordinateMatrix, AtomLabel, proname='En')
    res['E3e'] = GetWHIM12(CoordinateMatrix, AtomLabel, proname='En')
    res['De'] = GetWHIM13(CoordinateMatrix, AtomLabel, proname='En')
    res['L1v'] = GetWHIM1(CoordinateMatrix, AtomLabel, proname='V')
    res['L2v'] = GetWHIM2(CoordinateMatrix, AtomLabel, proname='V')
    res['L3v'] = GetWHIM3(CoordinateMatrix, AtomLabel, proname='V')
    res['Tv'] = GetWHIM4(CoordinateMatrix, AtomLabel, proname='V')
    res['Av'] = GetWHIM5(CoordinateMatrix, AtomLabel, proname='V')
    res['Vv'] = GetWHIM6(CoordinateMatrix, AtomLabel, proname='V')
    res['P1v'] = GetWHIM7(CoordinateMatrix, AtomLabel, proname='V')
    res['P2v'] = GetWHIM8(CoordinateMatrix, AtomLabel, proname='V')
    res['Kv'] = GetWHIM9(CoordinateMatrix, AtomLabel, proname='V')
    res['E1v'] = GetWHIM10(CoordinateMatrix, AtomLabel, proname='V')
    res['E2v'] = GetWHIM11(CoordinateMatrix, AtomLabel, proname='V')
    res['E3v'] = GetWHIM12(CoordinateMatrix, AtomLabel, proname='V')
    res['Dv'] = GetWHIM13(CoordinateMatrix, AtomLabel, proname='V')
    res['L1p'] = GetWHIM1(CoordinateMatrix, AtomLabel, proname='alapha')
    res['L2p'] = GetWHIM2(CoordinateMatrix, AtomLabel, proname='alapha')
    res['L3p'] = GetWHIM3(CoordinateMatrix, AtomLabel, proname='alapha')
    res['Tp'] = GetWHIM4(CoordinateMatrix, AtomLabel, proname='alapha')
    res['Ap'] = GetWHIM5(CoordinateMatrix, AtomLabel, proname='alapha')
    res['Vp'] = GetWHIM6(CoordinateMatrix, AtomLabel, proname='alapha')
    res['P1p'] = GetWHIM7(CoordinateMatrix, AtomLabel, proname='alapha')
    res['P2p'] = GetWHIM8(CoordinateMatrix, AtomLabel, proname='alapha')
    res['Kp'] = GetWHIM9(CoordinateMatrix, AtomLabel, proname='alapha')
    res['E1p'] = GetWHIM10(CoordinateMatrix, AtomLabel, proname='alapha')
    res['E2p'] = GetWHIM11(CoordinateMatrix, AtomLabel, proname='alapha')
    res['E3p'] = GetWHIM12(CoordinateMatrix, AtomLabel, proname='alapha')
    res['Dp'] = GetWHIM13(CoordinateMatrix, AtomLabel, proname='alapha')
    res['P3p'] = GetWHIM14(CoordinateMatrix, AtomLabel, proname='alapha')
    res['P3u'] = GetWHIM14(CoordinateMatrix, AtomLabel, proname='u')
    res['P3m'] = GetWHIM14(CoordinateMatrix, AtomLabel, proname='m')
    res['P3e'] = GetWHIM14(CoordinateMatrix, AtomLabel, proname='En')
    res['P3v'] = GetWHIM14(CoordinateMatrix, AtomLabel, proname='V')
    return res

#############################################################################
# if __name__ =  = "__main__":

#     from GeoOpt import GetARCFile
#     mol = 'c1ccccc1N'
#     inputmol = pybel.readstring('smi', mol)
#     dir_  =  GetARCFile(inputmol)
#     result = GetWHIMSandersonElectronegativity()
#     print(result)
#     print(len(result))
#     shutil.rmtree(dir_, ignore_errors = True)
