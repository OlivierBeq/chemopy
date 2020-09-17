# -*- coding: utf-8 -*-


"""Utility functions for pychem."""

import math
import os
from pathlib import Path
from typing import List

import scipy
from openbabel import openbabel, pybel
from rdkit import Chem


def GetR(n: int) -> List[float]:
    """Calcuate the parameters R of the RDF equation."""
    R = []
    for i in range(2, n + 2):
        R.append(float(i * 0.5))
    return R


def GetAtomDistance(x: List[float], y: List[float]) -> float:
    """Calculate Euclidean distance between two atomic coordinates."""
    temp = [math.pow(x[0] - y[0], 2), math.pow(x[1] - y[1], 2), math.pow(x[2] - y[2], 2)]
    res = math.sqrt(sum(temp))
    return res


def GetGeometricalDistanceMatrix(CoordinateList: List[List[float]]) -> scipy.matrix:
    """Calculate distance matrix from coordinate list."""
    NAtom = len(CoordinateList)
    DistanceMatrix = scipy.zeros((NAtom, NAtom))
    for i in range(NAtom - 1):
        for j in range(i + 1, NAtom):
            DistanceMatrix[i, j] = GetAtomDistance(CoordinateList[i], CoordinateList[j])
            DistanceMatrix[j, i] = DistanceMatrix[i, j]
    return DistanceMatrix


def IsInSubdirectoryTree(path: str, parent_path: str) -> bool:
    """Test if a path is in the subdirectory tree of another.

    :param path: the path to be tested
    :param parent_path: the possible parent directory
    """
    # Get full absolute paths
    path_ = os.path.realpath(path)
    parent_path_ = os.path.realpath(parent_path)
    # Get relative path
    relpath = os.path.relpath(parent_path_, path_)
    # Remove / and \ and '..'
    relpath = relpath.replace('/', '').replace('\\', '').replace('..', '')
    return len(relpath) == 0


def GetFileInDirFromExt(dirpath: str, ext: List[str]) -> List[str]:
    """Identify the files in the current directory based on their extension.

    Does not support recursuve search of files.
    :param dirpath: path to directory
    :param ext: extensions (including .) to filter files on
    :return: list of files matching the required extension(s)
    """
    return [os.path.realpath(os.path.join(dirpath, x))
            for x in os.listdir(dirpath)
            if os.path.splitext(x)[1].lower() == ext.lower()]


def AreAllPathsAbsolute(paths: List[str]) -> bool:
    """Verify all paths are absolute."""
    return all(map(os.path.isabs, paths))


def AreAllPathsRelative(paths: List[str]) -> bool:
    """Verify all paths are relative."""
    return all(map(lambda x: not os.path.isabs(x), paths))


def GetLastestCreatedFile(dirpath: str = None, filepaths: List[str] = None, strict: bool = True) -> str:
    """Get the file that was created last.

    Search can either be executed on an entire directory or a
    list of files that can come from  different directories.

    If only filepaths is provided, all paths must be absolute.
    If both dirpath and filepaths are provided, filepaths must not be absolute paths.

    :param dirpath: path to directory
    :param filepaths: paths of files
    :param strict: whether to stop when a file from filepaths does not exist
    """
    if (dirpath, filepaths) == (None, None):
        raise ValueError('Either dirpath or filepaths must be provided.')
    if dirpath is not None and filepaths is not None:
        if not AreAllPathsRelative(filepaths):
            raise ValueError('filepaths must not be absolute when dirpath provided.')
        if not os.path.isdir(dirpath):
            raise NotADirectoryError(f'{dirpath} does not exist.')
        files = []
        for file_ in filepaths:
            full_path = os.path.realpath(os.path.join(dirpath, file_))
            if not os.path.isfile(full_path) and strict:
                raise FileNotFoundError(f'{full_path} does not exist.')
            elif os.path.isfile(full_path):
                files.append((full_path, Path(full_path).stat().st_ctime_ns))
        return sorted(files, key=lambda x: x[1], reverse=True)[0][0]
    elif dirpath is not None:
        if not os.path.isdir(dirpath):
            raise NotADirectoryError(f'{dirpath} does not exist.')
        files = []
        for file_ in os.listdir(dirpath):
            full_path = os.path.realpath(os.path.join(dirpath, file_))
            if os.path.isfile(full_path):
                files.append((full_path, Path(full_path).stat().st_ctime_ns))
        return sorted(files, key=lambda x: x[1], reverse=True)[0][0]
    else:
        if not AreAllPathsAbsolute(filepaths):
            raise ValueError('filepaths must be absolute when dirpath is not provided.')
        files = []
        for file_ in filepaths:
            if not os.path.isfile(file_) and strict:
                raise FileNotFoundError(f'{file_} does not exist.')
            elif os.path.isfile(file_):
                files.append((file_, Path(file_).stat().st_ctime_ns))
        return sorted(files, key=lambda x: x[1], reverse=True)[0][0]


def rdkit_to_openbabel_mol(mol: Chem.Mol) -> openbabel.OBMol:
    """Convert a RDKit molecule to an OpenBabel molecule.

    :param mol: RDKit molecule
    """
    obmol = openbabel.OBMol()
    # Add hydrogen atoms to complete molecule if needed
    rdkitmol = Chem.Mol(mol)
    # Perceive valence and ring information before assigning hydrogens
    rdkitmol.UpdatePropertyCache(strict=False)
    rdkitmol = Chem.AddHs(rdkitmol)
    # Kekulize molecule
    Chem.rdmolops.Kekulize(rdkitmol, clearAromaticFlags=True)
    # Add atoms
    for atom in mol.GetAtoms():
        # Create new atom and assign values
        obatom = obmol.NewAtom()
        obatom.SetAtomicNum(atom.GetAtomicNum())
        obatom.SetIsotope(atom.GetIsotope())
        obatom.SetFormalCharge(atom.GetFormalCharge())
        obatom.SetPartialCharge(atom.GetDoubleProp('_PartialCharge'))
        obatom.SetSpinMultiplicity(atom.GetNumRadicalElectrons() + 1)
    for bond in mol.GetBonds():
        obmol.AddBond(bond.GetBeginAtomIdx() + 1, bond.GetEndAtomIdx() + 1, int(bond.GetBondTypeAsDouble()))
    obmol.AssignSpinMultiplicity(True)
    return obmol


def rdkit_to_pybel_mol(mol: Chem.Mol) -> pybel.Molecule:
    """Convert a RDKit molecule to an OpenBabel.Pybel molecule.

    :param mol: RDKit molecule
    """
    return pybel.Molecule(rdkit_to_openbabel_mol(mol))


def openbabel_to_rdkit_mol(mol: openbabel.OBMol) -> Chem.Mol:
    """Convert an OpenBabel molecule to a RDKit molecule.

    :param mol: OpenBabel molecule
    """
    # Create an editable molecule
    rdkitmol = Chem.rdchem.EditableMol(Chem.rdchem.Mol())
    for obatom in openbabel.OBMolAtomIter(mol):
        # Create new atom and assign values
        atom = Chem.Atom(obatom.GetAtomicNum())
        atom.SetIsotope(obatom.GetIsotope())
        atom.SetFormalCharge(obatom.GetFormalCharge())
        atom.SetDoubleProp('_PartialCharge', obatom.GetPartialCharge())
        atom.SetNumRadicalElectrons(obatom.GetSpinMultiplicity() - 1 if obatom.GetSpinMultiplicity() != 0 else 0)
        # Add it to the current molecule
        rdkitmol.AddAtom(atom)

        orders = {1: Chem.rdchem.BondType.SINGLE,
                  2: Chem.rdchem.BondType.DOUBLE,
                  3: Chem.rdchem.BondType.TRIPLE,
                  4: Chem.rdchem.BondType.QUADRUPLE,
                  5: Chem.rdchem.BondType.QUINTUPLE,
                  1.5: Chem.rdchem.BondType.AROMATIC}
    for obbond in openbabel.OBMolBondIter(mol):
        rdkitmol.AddBond(obbond.GetBeginAtomIdx() - 1, obbond.GetEndAtomIdx() - 1, orders[obbond.GetBondOrder()])
    rdkitmol = rdkitmol.GetMol()
    Chem.SanitizeMol(rdkitmol)
    return rdkitmol


def pybel_to_rdkit_mol(mol: pybel.Molecule) -> Chem.Mol:
    """Convert an OpenBabel.Pybel molecule to a RDKit molecule.

    :param mol: OpenBabel.Pybel molecule
    """
    return openbabel_to_rdkit_mol(mol.OBMol)
