# -*- coding: utf-8 -*-


"""Molecule geometry optimization with MOPAC."""


import os
import shutil
import subprocess  # noqa: S404
import tempfile
import warnings
from typing import List, Tuple, Union

from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import rdDistGeom

from chemopy import vector3d
from chemopy.external_parser import ExternalToolsParser
from chemopy.utils import (GetFileInDirFromExt, GetLastestCreatedFile,
                           IsInSubdirectoryTree)

MOPAC_CONFIG = {}


class Atom:
    """Wrapper for atomic properties."""

    def __init__(self, Coordinates: List[float]):
        """Initialize an Atom object."""
        self.pos = vector3d.Vector3D()
        self.radius = 0.0
        self.Coordinates = Coordinates
        self.Element = ''

    def SetCoordinates(self):
        """Parse raw ARC coordinates."""
        temp = self.Coordinates
        self.pos.x = float(temp[1])
        self.pos.y = float(temp[2])
        self.pos.z = float(temp[3])

    def GetCoordinates(self):
        """Get coordinates of the atom."""
        self.SetCoordinates()
        return self.pos

    def SetElement(self):
        """Set element from raw ARC coordinates."""
        temp = self.Coordinates
        self.Element = temp[0]

    def GetElement(self):
        """Get element."""
        self.SetElement()
        return self.Element

    def SetRadius(self):
        """Set radius."""
        radii = {'H': 1.20, 'N': 1.55, 'Na': 2.27, 'Cu': 1.40, 'Cl': 1.75, 'C': 1.70,
                 'O': 1.52, 'I': 1.98, 'P': 1.80, 'B': 1.85, 'Br': 1.85, 'S': 1.80, 'Se': 1.90,
                 'F': 1.47, 'Fe': 1.80, 'K': 2.75, 'Mn': 1.73, 'Mg': 1.73, 'Zn': 1.39, 'Hg': 1.8,
                 'Li': 1.8, '.': 1.8}
        temp = self.GetElement()
        if temp in radii.keys():
            self.radius = radii[temp]
        else:
            self.radius = radii['.']

    def GetRadius(self):
        """Get the radius of the atom."""
        self.SetRadius()
        return self.radius


def GetAtomClassList(Coordinates: List[float]) -> List[Atom]:
    """Get a list of atoms from a list of raw ARC coordinates.

    :param Coordinates: raw ARC coordinates as returned by _ReadCoordinates
    """
    Atoms = []
    for i in Coordinates:
        atom = Atom(i)
        atom.SetCoordinates()
        atom.SetElement()
        atom.SetRadius()
        Atoms.append(atom)
    return Atoms


def _ReadCoordinates(arc_file):
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


def FormatConversion(inputmol: Chem.Mol,
                     method='AM1', version='7.1',
                     outfile=None, outdir=None, sdf=False,
                     ) -> Tuple[str, str]:
    """Prepare a molecule to be optimized by MOPAC.

    The preparation goes as follows:

    1. all hydrogens are made explicit,
    2. a 3D conformer is generated from RDKit's stochastic search,
       based on distance geometry and  exerimental crystallographic knowledge,
       Wang S., J. Chem. Inf. Model. (2020), 60(4),2044-2058.
    3. the conformer is converted to a first MOPAC input file with OpenBabel,
    4. semi-empirical method and MOPAC version is added to the MOPAC input file.

    A SDF can also be written to disk from the generated conformer to ease debugging.

    :param method: MOPAC semi-empirical method to be used for molecular geometry oprimization
    :param version: version of MOPAC to use
    :param outfile: name of the output the MOPAC input file
    :param outdir: directory where to create the MOPAC input file
        If not specified, a temporary directory is created but
        THE USER MUST REMOVE IT THEMSELVES USING GeOpt.Dispose.
    :param sdf: whether a sdf file should also be created for easier debugging.
    :return: the directory where the MOPAC input file was created
    """
    # Ensure method and MOPAC version are compatible
    if not IsMethodSupportedByMOPAC(method, version):
        raise ValueError(f'Method {method} is not supported by MOPAC {version}.')
    # Step 1: add Hs
    inputmol = Chem.AddHs(inputmol)
    # Step 2: generate conformer
    rdDistGeom.EmbedMolecule(inputmol, rdDistGeom.ETKDGv3())
    pybelmol = pybel.readstring('sdf', Chem.MolToMolBlock(inputmol))
    running_dir = tempfile.mkdtemp() if outdir is None else outdir
    mpo_name = 'temp' if outfile is None else outfile
    if sdf:
        outputmol = pybel.Outputfile('sdf', os.path.join(running_dir, f'{mpo_name}.sdf'), overwrite=True)
        outputmol.write(pybelmol)
        outputmol.close()
    # Step 3: create first MOPAC input file
    outputmol = pybel.Outputfile('mop', os.path.join(running_dir, f'{mpo_name}.dat'), overwrite=True)
    outputmol.write(pybelmol)
    outputmol.close()
    # Step 4: add additional information to MOPAC input file
    with open(os.path.join(running_dir, f'{mpo_name}.dat'), 'r') as f:
        data = f.readlines()
    data[0] = method + (' PRTCHAR \n' if version == '2016' else '  \n')
    with open(os.path.join(running_dir, f'{mpo_name}.dat'), 'w') as f:
        f.write("".join(data))
    return running_dir, f'{mpo_name}.dat'


def MopacCfgParser(path: str = None) -> None:
    """Read the MOPAC configuration file and locate executables.

    :param path: path to MOPACc configuration file.
    """
    global MOPAC_CONFIG
    if path is None:
        path = os.path.join(os.path.dirname(__file__), 'MOPAC.cfg')
    etp = ExternalToolsParser(path, ['version', 'path', 'bin', 'methods'], False)
    for _, section in etp.tools.items():
        version_, path_, bin_, methods_ = (section['version'], section['path'],
                                           section['bin'], section['methods'])
        methods_ = [x.strip() for x in methods_.split()]
        # Make sure path is right
        if path_.startswith('.'):
            fullpath_ = os.path.realpath(os.path.join(os.path.dirname(__file__), path_, bin_))
        else:
            fullpath_ = os.path.realpath(os.path.join(path_, bin_))
        # Make sure executable exists
        if not os.path.isfile(fullpath_):
            warnings.warn(f'Wrong path in MOPAC config file: MOPAC(v{version_}) -> {fullpath_}', ResourceWarning)
        else:
            MOPAC_CONFIG[version_] = [fullpath_, methods_]
    if not MOPAC_CONFIG:  # Empty or corrupted config file
        raise ValueError('MOPAC config file is either empty or corrupted and needs to be edited.')


def RunMOPAC(filename: str, version: str = '7.1') -> int:
    """Run the MOPAC on a well prepared input file.

    Parse default MOPAC config file if not read already.

    :param filename: path to the well prepared MOPAC input file
    :param version: MOPAC version to be used
    """
    # Ensure all requirements are set
    if not IsMOPACVersionAvailable(version):
        raise ValueError(f'MOPAC version {version} is not available. Check your MOPAC config file.')
    # Run optimization
    mopac_bin = MOPAC_CONFIG[str(version)][0]
    try:
        retcode = subprocess.call(f'{mopac_bin} {filename}', shell=False)  # noqa: S603
        return retcode
    except Exception:
        return 1


def Dispose(dir_: str, force: bool = False) -> None:
    """Properly dispose of a temporary folder.

    :param force: whether to allow removing a non temporary directory.
    """
    # If neither is dir_ in tempdir nor force set to True
    if not (IsInSubdirectoryTree(dir_, tempfile.gettempdir()) or force):
        raise PermissionError('Directory is not temporary. If you know what you '
                              'are doing force the  deletion')
    shutil.rmtree(dir_, ignore_errors=True)


def IsMOPACVersionAvailable(version: str) -> bool:
    """Return if the dsired version of MOPAC can be used.

    :param version: version of MOPAC
    """
    if not MOPAC_CONFIG:
        MopacCfgParser()
    return str(version) in MOPAC_CONFIG.keys()


def IsMethodSupportedByMOPAC(method: str, version: str) -> bool:
    """Return if the version of MOPAC supports a specific method.

    :param method: semi-empirical method to be applied
    :param version: version of MOPAC
    """
    if not IsMOPACVersionAvailable(version):
        raise ValueError(f'MOPAC version {version} is not available. Check your MOPAC config file.')
    return str(method) in MOPAC_CONFIG[f'{version}'][1]


def GetARCFile(inputmol: Chem.Mol, method: str = 'PM3', version: str = '7.1',
               verbose: bool = True, exit_on_fail: bool = False,
               ) -> Union[Tuple[str, str], None]:
    """Optimize molecule geometry with MOPAC.

    :param inputmol: molecule to optimize
    :param method: semi-empirical method to apply
    :param version: version of MOPAC to be used
    :param verbose: whether to print progress messages
    :param exit_on_fail: if False, if a method fails at generating
                         a structure, others are tried from most to
                         least accurate.
                         if True, return False on failure.
    :return: Tuple of (path to folder, path to arc_file) on success,
             None otherwise.
    """
    # Create proper input file
    dir_, dat_file = FormatConversion(inputmol, method, version)
    # Ensure dat file exists
    full_path = os.path.join(dir_, dat_file)
    if not os.path.isfile(full_path):
        raise FileNotFoundError('Molecule could not be prepared for MOPAC.')
    # Run MOPAC
    retcode = RunMOPAC(full_path, version)
    # Get generated file
    # Different versions of MOPAC handle the outputname differently
    # e.g. 7.1 appends .arc after the .dat giving a .dat.arc file
    # while 2016 replaces the .dat by .arc
    output = GetFileInDirFromExt(dir_, '.arc')
    # Success when return code is 0
    success = not retcode and len(output) > 0
    if success:
        if verbose:
            print(f'Molecule geometry was successfully optimized using {method}.')  # noqa T001
        output = output[0] if len(output) == 1 else GetLastestCreatedFile(filepaths=output)
        # Rename output
        curated_filename = f'{os.path.splitext(full_path)[0]}.arc'
        os.rename(output, curated_filename)
        return dir_, curated_filename
    elif exit_on_fail:
        if verbose:
            print(f'Geometry optimization failed with {method}')  # noqa T001
        # if IsInSubdirectoryTree(dir_, tempfile.gettempdir()):
            # Dispose(dir_)
        return
    else:  # neither success nor exit_on_fail
        if verbose:
            print(f'Geometry optimization failed with {method}')  # noqa T001
        # if IsInSubdirectoryTree(dir_, tempfile.gettempdir()):
            # Dispose(dir_)
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
            new_attempt = GetARCFile(inputmol, methods_tried[0], version, verbose, exit_on_fail=True)
            if new_attempt is not None:
                return new_attempt


def GetOptimizedMol(arc_file: str = None, inputmol: Chem.Mol = None,
                    method: str = 'PM3', version: str = '7.1',
                    verbose: bool = False, dispose: bool = True,
                    ) -> Union[pybel.Molecule, None]:
    """Optimize molecule geometry with MOPAC and return the optimized molecule.

    In order not to optimize molecule multiple times, an ARC file may be provided.
    The path to the MOPAC output file is determined based on the one of the provided
    ARC file.

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
        mopac_out_path = GetFileInDirFromExt(mopac_out_dir, '.out')
        mopac_out_path = mopac_out_path[0] if len(mopac_out_path) == 1 \
            else GetLastestCreatedFile(filepaths=mopac_out_path)
        if not len(mopac_out_path):
            return None
        pybelmol = next(pybel.readfile('mopout', mopac_out_path))
        return Chem.MolFromMol2Block(pybelmol.write(format='mol2'))
    else:
        res = GetARCFile(inputmol, method, version, verbose, False)
        if res is None:
            return None
        dir_, arc_file_ = res
        mopac_out_path = GetFileInDirFromExt(dir_, '.out')
        mopac_out_path = mopac_out_path[0] if len(mopac_out_path) == 1 \
            else GetLastestCreatedFile(filepaths=mopac_out_path)
        if not len(mopac_out_path):
            return None
        pybelmol = next(pybel.readfile('mopout', mopac_out_path))
        if dispose:
            Dispose(dir_)
        return pybelmol
        # return Chem.MolFromMol2Block(pybelmol.write(format='mol2'), removeHs=False)


##############################################################################
# if __name__=="__main__":
#     mol='C1C=CCS1'
#     mol='SCCC(=O)N1[C@@H](CCC1)C(=O)OCC'
#     inputmol=pybel.readstring('smi',mol)
#     dir = GetARCFile(inputmol)
#     res=_ReadCoordinates(dir)
#     print(res)
#     Dispose(dir)
