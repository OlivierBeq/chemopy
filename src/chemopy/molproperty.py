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

from chemopy.external_parser import ExternalToolsParser
from chemopy.GeoOpt import Dispose


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


def CalculateXlogP(mol: Union[Chem.Mol, List[Chem.Mol]]) -> Union[float, List[float]]:
    """Calculate XLogP.

    From Cheng T. et al., J. Chem. Inf. Model. 2007, 47, 2140-2148.
    """
    # Check XlogP is properly configured
    etp = ExternalToolsParser(path=None, required_fields=['path'], skip_errors=False)
    if 'XLOGP' not in etp.tools.keys():
        raise NotImplementedError("XLOGP is not set up.")
    # Make sure to work with a list of molecules
    if isinstance(mol, Chem.Mol):
        mol = [mol]
    # Prepare temporary files
    running_dir = tempfile.mkdtemp()
    output_sd_path = os.path.join(running_dir, 'xlogp_calc.sdf')
    output_xlogp_path = os.path.join(running_dir, 'xlogp_calc.out')
    outputmol = pybel.Outputfile('sdf', output_sd_path, overwrite=True)
    # Write all molecules to SD file
    for mol_ in mol:
        pybelmol = pybel.readstring('sdf', Chem.MolToMolBlock(mol_))
        outputmol.write(pybelmol)
    outputmol.close()
    # Prepare path to executable file
    params = etp.tools['XLOGP']
    if params['path'].startswith('.'):
        path_prefix = os.path.realpath(os.path.join(os.path.dirname(__file__), params['path']))
    else:
        path_prefix = os.path.realpath(params['path'])
    exTTDB = os.path.realpath(os.path.join(path_prefix, '../parameter/DEFAULT.TTDB'))
    if platform.startswith('win32'):
        xlogp_bin = os.path.join(path_prefix, params['win_bin'])
    elif platform.startswith('linux'):
        if architecture().startswith('32'):
            xlogp_bin = os.path.join(path_prefix, params['lnx_bin'])
        else:
            xlogp_bin = os.path.join(path_prefix, params['lnx64_bin'])
    elif platform.startswith('darwin'):
        xlogp_bin = xlogp_bin = os.path.join(path_prefix, params['mac_bin'])
    else:
        Dispose(running_dir)
        raise RuntimeError(f'Platform ({platform}) not supported.')
    # Run XLOGP
    retcode = subprocess.call(f'{xlogp_bin} {output_sd_path} {output_xlogp_path} {exTTDB}', shell=False)  # noqa: S603
    if retcode:
        Dispose(running_dir)
        raise RuntimeError('XLOGP did not succeed to run properly.')
    values = []
    with open(output_xlogp_path) as xlogp:
        for line in xlogp:
            if line.startswith('Missing'):
                raise RuntimeError('Atom parameter missing.')
            raw_data = line.split()
            if not raw_data[-1].startswith('(WARN'):
                values.append(float(raw_data[-1]))
            else:
                values.append(float(raw_data[-2]))
    Dispose(running_dir)
    return values if len(values) > 1 else values[0]


def CalculateXlogP2(mol: Chem.Mol) -> float:
    """Calculate XLogP^2.

    From Cheng T. et al., J. Chem. Inf. Model. 2007, 47, 2140-2148.
    """
    values = CalculateXlogP(mol)
    if isinstance(values, float):
        return values * values
    return [value * value for value in values]


def CalculateXlogS(mol: Union[Chem.Mol, List[Chem.Mol]]) -> Union[float, List[float]]:
    """Calculate XLogS.

    From Duan B.G. et al., Acta Phys. Sin. 2012, 28(10), 2249-2257.
    """
    # Check XlogS is properly configured
    etp = ExternalToolsParser(path=None, required_fields=['path'], skip_errors=False)
    if 'XLOGS' not in etp.tools.keys():
        raise NotImplementedError("XLOGS is not set up.")
    # Make sure to work with a list of molecules
    if isinstance(mol, Chem.Mol):
        mol = [mol]
    # Prepare temporary files
    running_dir = tempfile.mkdtemp()
    output_sd_path = os.path.join(running_dir, 'xlogs_calc.sdf')
    output_xlogs_path = os.path.join(running_dir, 'xlogs_calc.out')
    outputmol = pybel.Outputfile('sdf', output_sd_path, overwrite=True)
    # Write all molecules to SD file
    for mol_ in mol:
        pybelmol = pybel.readstring('sdf', Chem.MolToMolBlock(mol_))
        outputmol.write(pybelmol)
    outputmol.close()
    # Prepare path to executable file
    params = etp.tools['XLOGS']
    if params['path'].startswith('.'):
        path_prefix = os.path.realpath(os.path.join(os.path.dirname(__file__), params['path']))
    else:
        path_prefix = os.path.realpath(params['path'])
    if platform.startswith('win32'):
        xlogs_bin = os.path.realpath(os.path.join(path_prefix, params['win_bin']))
    elif platform.startswith('linux'):
        if architecture().startswith('32'):
            xlogs_bin = os.path.realpath(os.path.join(path_prefix, params['lnx_bin']))
        else:
            xlogs_bin = os.path.realpath(os.path.join(path_prefix, params['lnx64_bin']))
    else:
        Dispose(running_dir)
        raise RuntimeError(f'Platform ({platform}) not supported.')
    # Run XLOGS
    exTTDB = os.path.realpath(os.path.join(os.path.dirname(xlogs_bin), '../parameter/DEFAULT.TTDB'))
    retcode = subprocess.call(f'{xlogs_bin} {output_sd_path} {output_xlogs_path} {exTTDB}', shell=False)  # noqa: S603
    if retcode:
        Dispose(running_dir)
        raise RuntimeError('XLOGS did not succeed to run properly.')
    values = []
    with open(output_xlogs_path) as xlogs:
        for line in xlogs:
            if line.startswith('Missing'):
                raise RuntimeError('Atom parameter missing.')
            raw_data = line.split()
            if not raw_data[-1].startswith('(WARN'):
                values.append(float(raw_data[-1]))
            else:
                values.append(float(raw_data[-2]))
    Dispose(running_dir)
    return values if len(values) > 1 else values[0]


def CalculateXlogS2(mol: Chem.Mol) -> float:
    """Calculate XLogS^2.

    From Duan B.G. et al., Acta Phys. Sin. 2012, 28(10), 2249-2257.
    """
    values = CalculateXlogS(mol)
    if isinstance(values, float):
        return values * values
    return [value * value for value in values]


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
    try:
        result['XLogP'] = CalculateXlogP(mol)
        result['XLogP2'] = CalculateXlogP2(mol)
    except RuntimeError:
        pass
    try:
        result['XLogS'] = CalculateXlogS(mol)
        result['XLogS2'] = CalculateXlogS(mol)
    except RuntimeError:
        pass
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
