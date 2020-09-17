# -*- coding: utf-8 -*-


"""Molecular physical/chemical properties."""


import math
import os
import subprocess  # noqa: S404
import tempfile
from platform import architecture
from sys import platform
from typing import List, Union

from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit.Chem import MolSurf as MS

from pychem.external import ExternalToolsParser
from pychem.GeoOpt import Dispose


def CalculateMolLogP(mol: Chem.Mol) -> float:
    """Cacluate of MolLogP.

    From Wildman and G. M. Crippen JCICS _39_ 868-873 (1999).
    """
    return round(Crippen._pyMolLogP(mol), 3)


def CalculateMolLogP2(mol: Chem.Mol) -> float:
    """Cacluate MolLogP^2.

    From Wildman and G. M. Crippen JCICS _39_ 868-873 (1999).
    """
    res = Crippen._pyMolLogP(mol)
    return round(res * res, 3)


def CalculateMolMR(mol: Chem.Mol) -> float:
    """Cacluate molecular refraction.

    From Wildman and G. M. Crippen JCICS _39_ 868-873 (1999).
    """
    return round(Crippen._pyMolMR(mol), 3)


def CalculateTPSA(mol: Chem.Mol) -> float:
    """Calculate the topological polar surface area.

    From Ertl P. et al., J.Med.Chem. (2000), 43,3714-3717.
    """
    return round(MS.TPSA(mol), 3)


def _CalculateBondNumber(mol: Chem.Mol, bondtype: str = 'SINGLE') -> float:
    """Calculate number of bond of specified type.

    :param bondtype: can be SINGLE, DOUBLE, TRIPLE or AROMATIC.
    """
    i = 0
    for bond in mol.GetBonds():
        if bond.GetBondType().name == bondtype:
            i += 1
    return i


def CalculateUnsaturationIndex(mol: Chem.Mol) -> float:
    """Calculate unsaturation index."""
    nd = _CalculateBondNumber(mol, bondtype='DOUBLE')
    nt = _CalculateBondNumber(mol, bondtype='TRIPLE')
    na = _CalculateBondNumber(mol, bondtype='AROMATIC')
    res = math.log((1 + nd + nt + na), 2)
    return round(res, 3)


def CalculateHydrophilicityFactor(mol: Chem.Mol) -> float:
    """Calculate hydrophilicity factor.

    From Todeschini R. et al., SAR QSAR Environ Res (1997), 7,173-193.
    """
    nheavy = mol.GetNumHeavyAtoms()
    nc = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            nc += 1
    nhy = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 or atom.GetAtomicNum() == 8 or atom.GetAtomicNum() == 16:
            atomn = atom.GetNeighbors()
            for i in atomn:
                if i.GetAtomicNum() == 1:
                    nhy += 1
    res = ((1 + nhy) * math.log((1 + nhy), 2) + nc * (1.0 / nheavy * math.log(1.0 / nheavy, 2)) + math.sqrt(
        (nhy + 0.0) / (nheavy * nheavy))) / math.log(1.0 + nheavy)
    return round(res, 3)


def CalculateXlogP(mol: Chem.Mol) -> float:
    """Calculate XLogP.

    From Cheng T. et al., J. Chem. Inf. Model. 2007, 47, 2140-2148.
    """
    pass


def CalculateXlogP2(mol: Chem.Mol) -> float:
    """Calculate XLogP^2.

    From Cheng T. et al., J. Chem. Inf. Model. 2007, 47, 2140-2148.
    """
    pass


MolecularProperty = {'LogP': CalculateMolLogP,
                     'LogP2': CalculateMolLogP2,
                     'MR': CalculateMolMR,
                     'TPSA': CalculateTPSA,
                     'Hy': CalculateHydrophilicityFactor,
                     'UI': CalculateUnsaturationIndex,
                     }


def GetMolecularProperties(mol: Chem.Mol) -> dict:
    """Get all (6) constitutional descriptors."""
    result = {}
    for DesLabel in MolecularProperty.keys():
        result[DesLabel] = MolecularProperty[DesLabel](mol)
    return result


##########################################################
# if __name__ =='__main__':
#     smis = ['CCCC','CCCCC','CCCCCC','CC(N)C(=O)O','CC(N)C(=O)[O-].[Na+]']
#     smi5=['CCCCCC','CCC(C)CC','CC(C)CCC','CC(C)C(C)C','CCCCCN','c1ccccc1N']
#     for index, smi in enumerate(smis):
#         m = Chem.MolFromSmiles(smi)
#         print(index+1)
#         print(smi)
#         print('\t',GetMolecularProperty(m))
#     #f.close()
