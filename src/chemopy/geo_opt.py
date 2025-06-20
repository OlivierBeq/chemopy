# -*- coding: utf-8 -*-


"""Molecule geometry optimization with MOPAC."""


import os
import shutil
import subprocess  # noqa: S404
import tempfile
import warnings
from sys import platform
from typing import List, Tuple, Union, Optional

from openbabel import pybel
from rdkit import Chem
from rdkit.rdBase import BlockLogs
from rdkit.Chem import rdDistGeom

from . import vector3d
from .utils import (get_file_in_dir_from_ext, get_lastest_created_file, is_in_subdirectory_tree)



# Default OpenMOPAC configuration
MOPAC_CONFIG = {'2016': ['mopac', ['PM7', 'PM6', 'PM3', 'AM1', 'MNDO']]}


# Disable openbabel logging
pybel.ob.obErrorLog.StopLogging()


class Atom:
    """Wrapper for atomic properties."""

    def __init__(self, Coordinates: List[float]):
        """Initialize an Atom object."""
        self.pos = vector3d.Vector3D()
        self.radius = 0.0
        self.coordinates = Coordinates
        self.element = ''

    def set_coordinates(self):
        """Parse raw ARC coordinates."""
        temp = self.coordinates
        self.pos.x = float(temp[1])
        self.pos.y = float(temp[2])
        self.pos.z = float(temp[3])

    def get_coordinates(self):
        """Get coordinates of the atom."""
        self.set_coordinates()
        return self.pos

    def set_element(self):
        """Set element from raw ARC coordinates."""
        temp = self.coordinates
        self.element = temp[0]

    def get_element(self):
        """Get element."""
        self.set_element()
        return self.element

    def set_radius(self):
        """Set radius."""
        radii = {'H': 1.20, 'N': 1.55, 'Na': 2.27, 'Cu': 1.40, 'Cl': 1.75, 'C': 1.70,
                 'O': 1.52, 'I': 1.98, 'P': 1.80, 'B': 1.85, 'Br': 1.85, 'S': 1.80, 'Se': 1.90,
                 'F': 1.47, 'Fe': 1.80, 'K': 2.75, 'Mn': 1.73, 'Mg': 1.73, 'Zn': 1.39, 'Hg': 1.8,
                 'Li': 1.8, '.': 1.8}
        temp = self.get_element()
        if temp in radii.keys():
            self.radius = radii[temp]
        else:
            self.radius = radii['.']

    def get_radius(self):
        """Get the radius of the atom."""
        self.set_radius()
        return self.radius


def get_atom_class_list(coordinates: List[List[float]]) -> List[Atom]:
    """Get a list of atoms from a list of raw ARC coordinates.

    :param coordinates: raw ARC coordinates as returned by _ReadCoordinates
    """
    atoms = []
    for coords in coordinates:
        atom = Atom(coords)
        atom.set_coordinates()
        atom.set_element()
        atom.set_radius()
        atoms.append(atom)
    return atoms


def read_coordinates(arc_file: str):
    """Read coordinates and charges of atoms from a MOPAC ARC file.

    :param arc_file: Path to MOPAC .arc file
    """
    res = []
    with open(arc_file, 'r') as f:
        templine = f.readlines()
    for line in range(len(templine)):
        if templine[line][-7: -1] == "CHARGE":
            k = line
            break
    for i in templine[k + 4: len(templine) - 1]:
        temp = i.split()
        ElementCoordinate = [temp[0].strip(), temp[1].strip(),
                             temp[3].strip(), temp[5].strip(),
                             temp[-1].strip()]
        res.append(ElementCoordinate)
    return res


class MopacInputFile:
    """File holding inputs for MOPAC molecular geometry optimization."""

    periodic_table = Chem.GetPeriodicTable()

    def __init__(self, path:str, method='PM7', version='2016', opt: str=''):
        """Create an instance of a MOPAC input file."""
        self.path = path
        self.method = method
        self.version = version
        self.opt = opt

    def open(self):
        """Open the handle"""
        self.handle = open(self.path, 'w')
        self.handle.write(self.method + (' PRTCHAR ' if self.version == '2016' else ' ') + self.opt)

    def write(self, mol: Chem.Mol, conf_id: int = -1):
        """Write a molecule to the input MOPAC file.

        :param mol: RDKit molecule to write
        :param conf_id: conformer id
        """
        # If molecule is charged
        charge = Chem.GetFormalCharge(mol)
        self.handle.write(f" CHARGE={'+' if charge > 0 else '-'}{charge}" if charge != 0 else '')
        self.handle.write('\n\n\n')
        # Write atomic positions and set them all for optimization
        conf = mol.GetConformer(conf_id)
        for atom in mol.GetAtoms():
            index = atom.GetIdx()
            elem = self.periodic_table.GetElementSymbol(atom.GetAtomicNum())
            pos = conf.GetAtomPosition(index)
            self.handle.write(f'{elem: <3}{pos.x:8.5f} 1 {pos.y:8.5f} 1 {pos.z:8.5f} 1\n')

    def close(self):
        """Close the internal handle."""
        self.handle.close()

    def __enter__(self):
        self.open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


class MopacResultDir:
    """Directory holding all inputs and outputs of MOPAC."""

    def __init__(self, path: Optional[str] = None):
        """Create an instance of a MOPAC input/output directory.

        :param exists: whether the directory already exists
        """
        if path is not None:
            if not os.path.isdir(path):
                raise ValueError('Path should either be specified if the directory already exists or left blank.')
        self.already_exists = path is not None
        self.path = path

    def open(self):
        if self.path is None:
            self.path = tempfile.mkdtemp()
        return self

    def close(self):
        if not self.already_exists:
            clean_mopac_files(self.path)

    def __enter__(self):
        return self.open()

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__del__()

    def __del__(self):
        self.close()


def format_conversion(inputmol: Chem.Mol,
                      method='PM7', version='2016', opt='',
                      outfile=None, outdir=None,
                      ) -> Tuple[MopacResultDir, str]:
    """Prepare a molecule to be optimized by MOPAC.

    The preparation goes as follows:

    1. all hydrogens are made explicit,
    2. a 3D conformer is generated from RDKit's stochastic search,
       based on distance geometry and  exerimental crystallographic knowledge,
       Wang S., J. Chem. Inf. Model. (2020), 60(4),2044-2058.
    3. the conformer is converted to a first MOPAC input file with OpenBabel,
    4. semi-empirical method and MOPAC version is added to the MOPAC input file.

    :param method: MOPAC semi-empirical method to be used for molecular geometry optimization
    :param version: version of MOPAC to use
    :param opt: optimizations (e.g. tighter SCF convergence criteria)
    :param outfile: name of the output the MOPAC input file
    :param outdir: directory where to create the MOPAC input file
        If not specified, a temporary directory is created but
        THE USER MUST REMOVE IT THEMSELVES USING geo_opt.clean_mopac_files.
    :return: the directory where the MOPAC input file was created
    """
    # Ensure method and MOPAC version are compatible
    if not is_method_supported_by_mopac(method, version):
        raise ValueError(f'Method {method} is not supported by MOPAC {version}.')
    # Step 1: add Hs
    # Step 2: generate conformer
    confs = list(inputmol.GetConformers())
    if not (len(confs) > 0 and confs[-1].Is3D()):
        inputmol = Chem.AddHs(inputmol)
        rdDistGeom.EmbedMolecule(inputmol, rdDistGeom.ETKDGv3())
    else:
        inputmol = Chem.AddHs(inputmol, addCoords=True)
    running_dir = MopacResultDir().open() if outdir is None else outdir
    mpo_name = 'temp' if outfile is None else outfile
    # Step 3: create MOPAC input file
    with MopacInputFile(os.path.join(running_dir.path, f'{mpo_name}.dat'), method=method, version=version, opt=opt) as output_file:
        output_file.write(inputmol)
    return running_dir, f'{mpo_name}.dat'


def run_mopac(filename: str, version: str = '2016') -> int:
    """Run the MOPAC on a well-prepared input file.

    Parse default MOPAC config file if not read already.

    :param filename: path to the well-prepared MOPAC input file
    :param version: MOPAC version to be used
    """
    # Ensure all requirements are set
    if not is_mopac_version_available(version):
        raise ValueError(f'MOPAC version {version} is not available. Check your MOPAC config file.')
    # Run optimization
    mopac_bin = MOPAC_CONFIG[str(version)][0]
    try:
        retcode = subprocess.call(f'{mopac_bin} {filename}', shell=False,
                                  stdin=subprocess.DEVNULL,
                                  stdout=subprocess.DEVNULL,
                                  stderr=subprocess.DEVNULL)  # noqa: S603
        return retcode
    except Exception:
        return 1


def clean_mopac_files(dir_: str, force: bool = False) -> None:
    """Properly dispose of a temporary folder.

    :param force: whether to allow removing a non-temporary directory.
    """
    # If neither is dir_ in tempdir nor force set to True
    if not (is_in_subdirectory_tree(dir_, tempfile.gettempdir()) or force):
        raise PermissionError('Directory is not temporary. If you know what you '
                              'are doing force the  deletion')
    shutil.rmtree(dir_, ignore_errors=True)


def is_mopac_version_available(version: str) -> bool:
    """Return if the desired version of MOPAC can be used.

    :param version: version of MOPAC
    """
    # Is version configured?
    if version not in MOPAC_CONFIG.keys():
        return False
    # Does the executable/environment variable exist?
    if shutil.which(MOPAC_CONFIG[version][0]) is None:
        return False
    # Ensure the executable runs without trouble
    if platform.startswith('win32'):
        echo_cmd = 'echo.'
    elif platform in ('linux', 'darwin'):
        echo_cmd = r'echo -e "\n"'
    else:
        raise RuntimeError(f'Platform ({platform}) not supported.')
    retcode = subprocess.call(f'{echo_cmd} | {MOPAC_CONFIG[version][0]}', shell=True,
                              stdin=subprocess.DEVNULL,
                              stdout=subprocess.DEVNULL,
                              stderr=subprocess.DEVNULL)
    return retcode == 0


def is_method_supported_by_mopac(method: str, version: str) -> bool:
    """Return if the version of MOPAC supports a specific method.

    :param method: semi-empirical method to be applied
    :param version: version of MOPAC
    """
    if not is_mopac_version_available(version):
        raise ValueError(f'MOPAC version {version} is not available.')
    return str(method) in MOPAC_CONFIG[f'{version}'][1]


def get_arc_file(inputmol: Chem.Mol, method: str = 'PM7', version: str = '2016', opt: str= '',
                 verbose: bool = True, exit_on_fail: bool = False,
                 ) -> Union[Tuple[MopacResultDir, str], None]:
    """Optimize molecule geometry with MOPAC.

    :param inputmol: molecule to optimize
    :param method: semi-empirical method to apply
    :param version: version of MOPAC to be used
    :param opt: optimizations (e.g. tighter SCF convergence criteria)
    :param verbose: whether to print progress messages
    :param exit_on_fail: if False, if a method fails at generating
                         a structure, others are tried from most to
                         least accurate.
                         if True, return False on failure.
    :return: Tuple of (path to folder, path to arc_file) on success,
             None otherwise.
    """
    # Create proper input file
    dir_, dat_file = format_conversion(inputmol, method, version, opt)
    # Ensure dat file exists
    full_path = os.path.join(dir_.path, dat_file)
    if not os.path.isfile(full_path):
        raise FileNotFoundError('Molecule could not be prepared for MOPAC.')
    # Run MOPAC
    retcode = run_mopac(full_path, version)
    # Get generated file
    # Different versions of MOPAC handle the outputname differently
    # e.g. 7.1 appends .arc after the .dat giving a .dat.arc file
    # while 2016 replaces the .dat by .arc
    output = get_file_in_dir_from_ext(dir_.path, '.arc')
    # Success when return code is 0
    success = not retcode and len(output) > 0
    if success:
        if verbose:
            print(f'Molecule geometry was successfully optimized using {method}.')  # noqa T001
        output = output[0] if len(output) == 1 else get_lastest_created_file(filepaths=output)
        # Rename output
        curated_filename = f'{os.path.splitext(full_path)[0]}.arc'
        os.rename(output, curated_filename)
        return dir_, curated_filename
    elif exit_on_fail:
        if verbose:
            print(f'Geometry optimization failed with {method}')  # noqa T001
        return
    else:  # neither success nor exit_on_fail
        if verbose:
            print(f'Geometry optimization failed with {method}')  # noqa T001
        methods_tried = MOPAC_CONFIG[f'{version}'][1]
        # Remove the method that was just used
        methods_tried.remove(method)
        # Try all possible methods, from most to least accurate
        skip_first_time = True
        while not success:
            if not skip_first_time:
                del methods_tried[0]  # Remove method tried to allow going to the next
                if len(methods_tried) == 0:  # No method left
                    if verbose:
                        print('Molecule could not be optimized.')  # noqa T001
                    else:
                        warnings.warn('Molecule could not be optimized', RuntimeWarning)
                    return
            skip_first_time = False
            if verbose:
                print(f'Attempting optimization with {methods_tried[0]}')  # noqa T001
            # Create proper input file
            new_attempt = get_arc_file(inputmol, methods_tried[0], version, verbose, exit_on_fail=True)
            if new_attempt is not None:
                return new_attempt


def get_optimized_mol(arc_file: str = None, inputmol: Chem.Mol = None,
                      method: str = 'PM7', version: str = '2016',
                      verbose: bool = False, dispose: bool = True,
                      ) -> Union[Chem.Mol, None]:
    """Optimize molecule geometry with MOPAC and return the optimized molecule.

    In order not to optimize a molecule multiple times, an ARC file may be provided.
    The path to the MOPAC output file is determined based on the ARC file provided.

    If not already optimized, a molecule may be provided.

    :param arc_file: Path to MOPAC .arc file
                     (ignored if inputmol provided).
    :param inputmol: molecule to optimize
                     (ignored if arc_file provided).
    :param method: semi-empirical method to apply
                   (ignored if arc_file provided).
    :param version: version of MOPAC to be used
                    (ignored if arc_file provided).
    :param verbose: whether to print progress messages
                    (ignored if arc_file provided).
    :param dispose: whether to remove generated MOPAC output files
                    (ignored if arc_file provided).
    :return: optimized rdkit molecule on success, None otherwise.
    """
    if arc_file is None and inputmol is None:
        raise ValueError('Either ARC file or inputmolecule must be provided.')
    if arc_file is not None:
        mopac_out_dir = os.path.dirname(arc_file)
        mopac_out_path = get_file_in_dir_from_ext(mopac_out_dir, '.out')
        mopac_out_path = mopac_out_path[0] if len(mopac_out_path) == 1 \
            else get_lastest_created_file(filepaths=mopac_out_path)
        if not len(mopac_out_path):
            return None
        pybelmol = next(pybel.readfile('mopout', mopac_out_path))
        # Convert to RDKit
        block = BlockLogs()
        mol = Chem.MolFromMol2Block(pybelmol.write(format='mol2'))
        if mol is None:
            mol = Chem.MolFromMolBlock(pybelmol.write(format='sdf'))
        del block
        return mol
    else:
        res = get_arc_file(inputmol, method, version, verbose, False)
        if res is None:
            return None
        dir_, arc_file_ = res
        mopac_out_path = get_file_in_dir_from_ext(dir_.path, '.out')
        mopac_out_path = mopac_out_path[0] if len(mopac_out_path) == 1 \
            else get_lastest_created_file(filepaths=mopac_out_path)
        if not len(mopac_out_path):
            return None
        pybelmol = next(pybel.readfile('mopout', mopac_out_path))
        if dispose:
            dir_.close()
        # Convert to RDKit
        block = BlockLogs()
        mol = Chem.MolFromMol2Block(pybelmol.write(format='mol2'), removeHs=False)
        if mol is None:
            mol = Chem.MolFromMolBlock(pybelmol.write(format='sdf'))
        del block
        return mol
