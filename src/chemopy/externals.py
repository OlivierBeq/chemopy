# -*- coding: utf-8 -*-


"""External tools for molecular descriptors/fingerprints."""


import os
import re
import subprocess  # noqa: S404
import tempfile
import warnings
from platform import architecture
from sys import platform
from typing import List, Union

import descriptastorus
import mordred
import pandas as pd
from mordred import descriptors
from openbabel import pybel
from rdkit import Chem
from rdkit.DataStructs.cDataStructs import ExplicitBitVect, SparseBitVect
from scipy.io import arff

from chemopy.external_parser import ExternalToolsParser
from chemopy.GeoOpt import Dispose


def JCompoundMapper(mol: Union[Union[pybel.Molecule, Chem.Mol], List[Union[pybel.Molecule, Chem.Mol]]],
                    bits: int = 2048) -> dict:
    """Calculate molecular fingerprints with JCompoundMapper.

    Java must be installed and its path in the PATH environment variable.

    :param mol: either one or multiple molecules the fingerprints will be calculated from.
    :param bits: length of the bitstring.
    """
    # Verify configuration of JCompoundMapper
    etp = ExternalToolsParser(path=None, required_fields=['path'], skip_errors=False)
    if 'JCompoundMapper' not in etp.tools.keys():
        raise NotImplementedError("JCompoundMapper is not set up.")
    params = etp.tools['JCompoundMapper']
    jar_path = os.path.realpath(os.path.join(os.path.dirname(__file__), params['path'], params['bin']))
    # Check types
    is_list = isinstance(mol, list)
    if is_list:
        is_pybel = all(map(lambda x: isinstance(x, pybel.Molecule), mol))
        is_rdkit = all(map(lambda x: isinstance(x, Chem.Mol), mol))
        if not (is_pybel or is_rdkit):
            raise ValueError('All molecules must have the same type: either rdkit or pybel.')
    else:
        is_pybel = isinstance(mol, pybel.Molecule)
        is_rdkit = isinstance(mol, Chem.Mol)
        if not (is_pybel or is_rdkit):
            raise ValueError('molecule must either be rdkit or pybel molecule.')
    # Keep trace of molecules being processed
    out_indices, results = [], {}
    # Create temporary folder
    running_dir = tempfile.mkdtemp()
    sdf_file = os.path.join(running_dir, 'molecules.sdf')
    if is_pybel:
        outputmol = pybel.Outputfile('sdf', sdf_file, overwrite=True)
    else:  # RDKit molecule
        outputmol = Chem.rdmolfiles.SDWriter(sdf_file)
    # Write molecule(s)
    if is_list:
        for i, mol_ in enumerate(mol):
            try:
                outputmol.write(mol_)
                out_indices.append(i)
            except IOError:
                pass
    else:
        out_indices.append(0)
        outputmol.write(mol)
    outputmol.close()
    # Call JCompoundMapper on output file for each fingerprint type
    fingerprints = ['DFS', 'ASP', 'AP2D', 'AT2D', 'CATS2D', 'PHAP2POINT2D', 'PHAP3POINT2D',
                    'SHED', 'RAD2D', 'LSTAR', 'AP3D', 'AT3D', 'CATS3D', 'PHAP2POINT3D',
                    'PHAP3POINT3D', 'RAD3D']
    for fingerprint in fingerprints:
        out_file = os.path.join(running_dir, f'fp_{fingerprint}.txt')
        command = f'java -jar {jar_path} -f {sdf_file} -c {fingerprint} -ff LIBSVM_SPARSE -o {out_file} -hs {bits}'
        retcode = subprocess.call(command, shell=False)  # noqa: S603
        if retcode:  # Error occured
            warnings.warn(f'JCompoundMapper did not succeed to run properly for fingerprint {fingerprint}.')
            results[fingerprint] = [[] for _ in range(len(out_indices))]
        else:
            fp_vectors = read_libsvmsparse(out_file, bits)  # Get results
            if len(fp_vectors) != len(out_indices):  # Not all processed molecules in result
                Dispose(running_dir)
                raise RuntimeError(f'JCompoundMapper results contained {len(fp_vectors)}'
                                   f' molecules but {len(out_indices)} were expected')
            results[fingerprint] = fp_vectors
    if is_list:
        # Transform a dictionary of lists to a list of dictionaries
        results = [dict(zip(results, t)) for t in zip(*results.values())]
    Dispose(running_dir)
    return results


def read_libsvmsparse(filename: str, length: int) -> List[str]:
    """Read the content of a LIBSVM_SPARSE file.

    :param filename: LIBSVM_SPARSE file
    """
    # Read data
    with open(filename) as inputf:
        data = inputf.readlines()
    results = []
    for result in data:
        # Get all bits set
        bits_set = result.split()[1:]
        if length <= 4096:
            fp = ExplicitBitVect(length)  # Fast
        else:
            fp = SparseBitVect(length)  # Memory-efficient
        bits = [int(x.split(':')[0]) - 1 for x in bits_set]  # Get bits index
        fp.SetBitsFromList(bits)  # Set bits
        results.append(fp.ToBitString())
    return results


def BlueDesc(mol: Union[Union[pybel.Molecule, Chem.Mol], List[Union[pybel.Molecule, Chem.Mol]]]) -> dict:
    """Calculate molecular descriptors with BlueDesc.

    Java must be installed and its path in the PATH environment variable.

    :param mol: either one or multiple molecules the fingerprints will be calculated from.
    """
    # Verify configuration of BlueDesc
    etp = ExternalToolsParser(path=None, required_fields=['path'], skip_errors=False)
    if 'BlueDesc' not in etp.tools.keys():
        raise NotImplementedError("BlueDesc is not set up.")
    params = etp.tools['BlueDesc']
    jar_path = os.path.realpath(os.path.join(os.path.dirname(__file__), params['path'], params['bin']))
    # Check types
    is_list = isinstance(mol, list)
    if is_list:
        is_pybel = all(map(lambda x: isinstance(x, pybel.Molecule), mol))
        is_rdkit = all(map(lambda x: isinstance(x, Chem.Mol), mol))
        if not (is_pybel or is_rdkit):
            raise ValueError('All molecules must have the same type: either rdkit or pybel.')
    else:
        is_pybel = isinstance(mol, pybel.Molecule)
        is_rdkit = isinstance(mol, Chem.Mol)
        if not (is_pybel or is_rdkit):
            raise ValueError('molecule must either be rdkit or pybel molecule.')
    # Keep trace of molecules being processed
    out_indices, results = [], []
    # Create temporary folder
    running_dir = tempfile.mkdtemp()
    sdf_file = os.path.join(running_dir, 'molecules.sdf')
    # Write molecule(s)
    with open(sdf_file, 'w') as output:
        if is_list:
            for i, mol_ in enumerate(mol):
                try:
                    if is_pybel:
                        output.write(f'{Chem.MolToMolBlock(Chem.MolFromMolBlock(mol_), forceV3000=True)}$$$$\n')
                    else:
                        output.write(f'{Chem.MolToMolBlock(mol_, forceV3000=True)}$$$$\n')
                    out_indices.append(i)
                except IOError:
                    pass
        else:
            out_indices.append(0)
            if is_pybel:
                output.write(f'{Chem.MolToMolBlock(Chem.MolFromMolBlock(mol), forceV3000=True)}$$$$\n')
            else:
                output.write(f'{Chem.MolToMolBlock(mol, forceV3000=True)}$$$$\n')
    # Call BlueDesc descriptors on output file
    out_file = os.path.join(running_dir, 'molecules.oddescriptors.arff')
    command = f'java -jar {jar_path} -f {sdf_file} -l \'?\''
    retcode = subprocess.call(command, shell=False)  # noqa: S603
    if retcode:  # Error occured
        raise RuntimeError('BlueDesc did not succeed to run properly.')
    else:
        data = arff.loadarff(out_file)  # Get results
        Dispose(running_dir)  # Remove temp dir
        df = pd.DataFrame(data[0]).drop("'?'", axis=1)  # Remove label column
        # Rename columns
        df = df.rename(columns=lambda x: re.sub(r'joelib2.+\.((?:(?!\.).)+)$', r'\1', x))
        results = df.to_dict(orient='records')
    if len(results) != len(out_indices):
        raise RuntimeError(f'BlueDesc results contained {len(results)} but {len(out_indices)} were expected')
    return results


def Mold2(mol: Union[Union[pybel.Molecule, Chem.Mol], List[Union[pybel.Molecule, Chem.Mol]]]) -> dict:
    """Calculate molecular descriptors with Mold2.

    Java must be installed and its path in the PATH environment variable.

    :param mol: either one or multiple molecules the fingerprints will be calculated from.
    """
    # Verify configuration of Mold2
    etp = ExternalToolsParser(path=None, required_fields=['path'], skip_errors=False)
    if 'Mold2' not in etp.tools.keys():
        raise NotImplementedError("Mold2 is not set up.")
    # Check types
    is_list = isinstance(mol, list)
    if is_list:
        is_pybel = all(map(lambda x: isinstance(x, pybel.Molecule), mol))
        is_rdkit = all(map(lambda x: isinstance(x, Chem.Mol), mol))
        if not (is_pybel or is_rdkit):
            raise ValueError('All molecules must have the same type: either rdkit or pybel.')
    else:
        is_pybel = isinstance(mol, pybel.Molecule)
        is_rdkit = isinstance(mol, Chem.Mol)
        if not (is_pybel or is_rdkit):
            raise ValueError('molecule must either be rdkit or pybel molecule.')
    # Keep trace of molecules being processed
    out_indices, results = [], []
    # Create temporary folder and write v2000 SD file
    running_dir = tempfile.mkdtemp()
    sdf_file = os.path.realpath(os.path.join(running_dir, 'molecules.sdf'))
    writer = pybel.Outputfile('sdf', sdf_file, overwrite=True)
    if is_list:
        for i, mol_ in enumerate(mol):
            try:
                if is_pybel:
                    writer.write(mol_)
                else:
                    writer.write(pybel.readstring('mol', Chem.MolToMolBlock(mol_)))
                out_indices.append(i)
            except IOError:
                pass
    else:
        out_indices.append(0)
        if is_pybel:
            writer.write(mol)
            # output.write(f'{Chem.MolToMolBlock(Chem.MolFromMolBlock(mol))}$$$$\n')
        else:
            writer.write(pybel.readstring('mol', Chem.MolToMolBlock(mol)))
    writer.close()
    # Prepare path to executable file
    params = etp.tools['Mold2']
    if params['path'].startswith('.'):
        path_prefix = os.path.realpath(os.path.join(os.path.dirname(__file__), params['path']))
    else:
        path_prefix = os.path.realpath(params['path'])
    if platform.startswith('win32'):
        mold2_bin = os.path.join(path_prefix, params['win_bin'])
        log_file = 'NUL'
        echo_cmd = 'echo.'
    elif platform.startswith('linux'):
        log_file = '/dev/null'
        echo_cmd = 'echo -e \'\n\''
        if architecture()[0].startswith('32'):
            mold2_bin = os.path.join(path_prefix, params['lnx_bin'])
        else:
            mold2_bin = os.path.join(path_prefix, params['lnx64_bin'])
    else:
        Dispose(running_dir)
        raise RuntimeError(f'Platform ({platform}) not supported.')
    # Call Mold2 descriptors on output file
    out_file = os.path.join(running_dir, 'molecules.mold2descriptors.tsv')
    mold2_bin = os.path.realpath(mold2_bin)  # Ensure separators are correct
    devnull = open(os.devnull, 'wb')
    command = f'{echo_cmd} | {mold2_bin} -i {sdf_file} -o {out_file} -r {log_file}'
    _ = subprocess.check_output(command, shell=True, stderr=devnull)  # noqa: S602
    devnull.close()
    data = pd.read_table(out_file).drop('Number', axis=1)  # Get results
    Dispose(running_dir)  # Remove temp dir
    # Rename columns
    results = data.to_dict(orient='records')
    if len(results) != len(out_indices):
        raise RuntimeError(f'Mold2 results contained {len(results)} but {len(out_indices)} were expected')
    return results


def Obabel(mol: Union[Union[pybel.Molecule, Chem.Mol], List[Union[pybel.Molecule, Chem.Mol]]],
           obabel_datadir: Union[str, None] = None) -> dict:
    """Calculate molecular descriptors with OpenBabel.

    :param mol: either one or multiple molecules the fingerprints will be calculated from.
    :param obabel_datadir: path to OpenBabel data directory (contains mmff parameters for instance).
                           This usually is of the form 'Python3/share/openbabel'.
                           If None, infers the path to Python share folder from Chemopy's install
                           directory.
    """
    initialized = False
    if obabel_datadir is not None:
        obabel_datadir = os.path.realpath(obabel_datadir)
        if os.path.exists(os.path.join(obabel_datadir, 'mmff94.ff')):
            os.environ['BABEL_DATADIR'] = obabel_datadir
            initialized = True
        else:
            warnings.warn('OpenBabel data directory is erroneous.', RuntimeWarning)
    if not initialized:
        site_packages_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        python_dir = os.path.dirname(os.path.dirname(site_packages_dir))
        obabel_datadir = os.path.join(python_dir, 'share', 'openbabel')
        if os.path.exists(os.path.join(obabel_datadir, 'mmff94.ff')):
            os.environ['BABEL_DATADIR'] = obabel_datadir
        else:
            raise RuntimeError('CCould not locate OpenBabel data directory, please indicate it.')
    # Calculate descriptors and cast molecule to pybel if needed
    results = []
    descs = ['HBA1', 'HBA2', 'HBD', 'logP', 'MR', 'MW', 'nF', 'rotors', 'TPSA']
    is_list = isinstance(mol, list)
    if is_list:
        for i in range(len(mol)):
            if isinstance(mol[i], Chem.Mol):
                mol_ = pybel.readstring('mol', Chem.MolToMolBlock(mol[i]))
                results.append(mol_.calcdesc(descs))
            else:
                results.append(mol[i].calcdesc(descs))
    elif isinstance(mol, Chem.Mol):
        mol_ = pybel.readstring('mol', Chem.MolToMolBlock(mol))
        results.append(mol_.calcdesc(descs))
    else:
        results.append(mol.calcdesc(descs))
    # Rename keys to meaningful names
    for i in range(len(results)):
        results[i]['rotB'] = results[i].pop('rotors')
    return results


def Mordred_2D(mol: Union[Union[pybel.Molecule, Chem.Mol], List[Union[pybel.Molecule, Chem.Mol]]]) -> dict:
    """Calculate 2D molecular descriptors with mordred.

    :param mol: either one or multiple molecules the fingerprints will be calculated from.
    """
    mordred_calculator = mordred.Calculator(descriptors, ignore_3D=True)
    if isinstance(mol, list):
        results = []
        for i in range(len(mol)):
            result = {}
            for descriptor in mordred_calculator.descriptors:
                try:
                    if isinstance(mol[i], Chem.Mol):
                        result[str(descriptor)] = descriptor(mol[i])
                    else:
                        mol_ = Chem.MolFromMolBlock(mol[i].write('mol'))
                        result[str(descriptor)] = descriptor(mol_)
                except (ZeroDivisionError, ValueError):
                    pass
            results.append(result)
    else:
        results = {}
        for descriptor in mordred_calculator.descriptors:
            try:
                if isinstance(mol, Chem.Mol):
                    results[str(descriptor)] = descriptor(mol)
                else:
                    mol_ = Chem.MolFromMolBlock(mol.write('mol'))
                    results[str(descriptor)] = descriptor(mol_)
            except (ZeroDivisionError, ValueError):
                pass
    return results


def descriptastorus(mol: Union[Union[pybel.Molecule, Chem.Mol], List[Union[pybel.Molecule, Chem.Mol]]]) -> dict:
    """Calculate 2D molecular descriptors with descriptastorus.

    :param mol: either one or multiple molecules the fingerprints will be calculated from.
    """
    atompaircounts = descriptastorus.MakeGenerator(('atompaircounts',)).processMol(mol[0], Chem.MolToSmiles(mol[0]))[1:]
    morgan3counts = descriptastorus.MakeGenerator(('morgan3counts',)).processMol(mol[0], Chem.MolToSmiles(mol[0]))
    morganchiral3counts = descriptastorus.MakeGenerator(('morganchiral3counts',)).processMol(mol[0], Chem.MolToSmiles(mol[0]))
    morganfeature3counts = descriptastorus.MakeGenerator(('morganfeature3counts',)).processMol(mol[0], Chem.MolToSmiles(mol[0]))
    rdkit2d = descriptastorus.MakeGenerator(('rdkit2d',)).processMol(mol[0], Chem.MolToSmiles(mol[0]))
    rdkit2dnormalized = descriptastorus.MakeGenerator(('rdkit2dnormalized',)).processMol(mol[0], Chem.MolToSmiles(mol[0]))
    rdkitfpbits = descriptastorus.MakeGenerator(('rdkitfpbits',)).processMol(mol[0], Chem.MolToSmiles(mol[0]))

