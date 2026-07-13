# -*- coding: utf-8 -*-


"""Multiple 2D molecular fingerprints."""

from functools import lru_cache
from typing import Optional

import numpy as np
from map4 import MAP4Calculator
from mhfp.encoder import MHFPEncoder
from openbabel import pybel
from rdkit import Chem, DataStructs
from rdkit.Chem import MACCSkeys
from rdkit.Chem import AllChem
from rdkit.Avalon.pyAvalonTools import GetAvalonFP

from .estate import EState

similaritymeasure = [i[0] for i in DataStructs.similarityFunctions]


@lru_cache(maxsize=None)
def _get_map4_calculator(radius: int, dimensions: int) -> MAP4Calculator:
    """Return a cached MAP4Calculator for (radius, dimensions).

    Constructing one is expensive (it builds an internal MHFPEncoder with its own
    precomputed hash permutation tables) and stateless w.r.t. the molecules it will
    later `.calculate()`, so building it once per (radius, dimensions) pair instead of
    once per molecule avoids redoing that setup work on every single call.
    """
    return MAP4Calculator(radius=radius, dimensions=dimensions, is_folded=True)


@lru_cache(maxsize=None)
def _get_key_list(prefix: str, n: int) -> list:
    """Return a cached list of `f"{prefix}{i + 1}"` names, instead of rebuilding it on every call."""
    return [f"{prefix}{i + 1}" for i in range(n)]


def _bitvect_to_dict(bv, prefix: str, n: int) -> dict:
    """Convert an RDKit bit/count vector to a dict, using a bulk C-level conversion.

    `DataStructs.ConvertToNumpyArray` fills a numpy array in one call instead of
    `bv.ToList()`, which builds a Python list one element at a time.
    """
    arr = np.zeros((n,), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(bv, arr)
    return dict(zip(_get_key_list(prefix, n), arr.tolist()))


class Fingerprint:
    """Molecular 2D fingerprints."""

    @staticmethod
    def calculate_rdk_fp(mol: Chem.Mol, nbits: int = 2048) -> dict:
        """Calculate a topological Daylight-like RDKit fingerprint.

        :param mols: the molecules
        :param nbits: Number of folded bits
        """
        bv = Chem.RDKFingerprint(mol, fpSize=nbits)
        return _bitvect_to_dict(bv, "RDK_", nbits)

    @staticmethod
    def calculate_maccs_fp(mol: Chem.Mol) -> dict:
        """Calculate the MACCS public keys (166 bits)."""
        bv = MACCSkeys.GenMACCSKeys(mol)
        return _bitvect_to_dict(bv, "MACCS_", 166)

    @staticmethod
    def calculate_fp2_fp(mol: Chem.Mol, pybel_mol=None) -> dict:
        """Calculate topological FP2 fingerprints (1024 bits)."""
        m = pybel_mol if pybel_mol is not None else pybel.readstring("smi", Chem.MolToSmiles(mol))
        on_bits = set(m.calcfp("FP2").bits)
        bits = [1 if i + 1 in on_bits else 0 for i in range(1024)]
        return dict(zip(_get_key_list("FP2_", 1024), bits))

    @staticmethod
    def calculate_fp3_fp(mol: Chem.Mol, pybel_mol=None) -> dict:
        """Calculate FP3 fingerprints (55 bits)."""
        m = pybel_mol if pybel_mol is not None else pybel.readstring("smi", Chem.MolToSmiles(mol))
        on_bits = set(m.calcfp("FP3").bits)
        bits = [1 if i + 1 in on_bits else 0 for i in range(55)]
        return dict(zip(_get_key_list("FP3_", 55), bits))

    @staticmethod
    def calculate_fp4_fp(mol: Chem.Mol, pybel_mol=None) -> dict:
        """Calculate FP4 fingerprints (307 bits)."""
        m = pybel_mol if pybel_mol is not None else pybel.readstring("smi", Chem.MolToSmiles(mol))
        on_bits = set(m.calcfp("FP4").bits)
        bits = [1 if i + 1 in on_bits else 0 for i in range(307)]
        return dict(zip(_get_key_list("FP4_", 307), bits))

    @staticmethod
    def calculate_estate_fp(mol: Chem.Mol) -> dict:
        """Calculate E-state fingerprints (79 bits)."""
        values = EState.calculate_estate_fingerprint(mol, implementation="chemopy", binary=True)
        return values

    @staticmethod
    def calculate_atompairs_fp(mol: Chem.Mol, nbits: int = 2048) -> dict:
        """Calculate atom pairs fingerprints.

        :param nbits: Number of folded bits
        """
        gen = AllChem.GetAtomPairGenerator(fpSize=nbits)
        bv = gen.GetFingerprint(mol)
        return _bitvect_to_dict(bv, "AtomPair_", nbits)

    @staticmethod
    def calculate_topological_torsion_fp(mol: Chem.Mol, nbits: int = 2048) -> dict:
        """Calculate Topological Torsion fingerprints.

        :param nbits: Number of folded bits
        """
        gen = AllChem.GetTopologicalTorsionGenerator(fpSize=nbits)
        bv = gen.GetFingerprint(mol)
        return _bitvect_to_dict(bv, "TopolTorsions_", nbits)

    @staticmethod
    def calculate_morgan_fp(mol: Chem.Mol, radius=2, nbits: int = 2048) -> dict:
        """Calculate Morgan fingerprints.

        :param radius: maximum radius of atom-centered substructures.
        :param rtype: Type of output, may either be:
                      bitstring (default), returns a binary string
                      rdkit, return the native rdkit DataStructs
                      dict, for a dict of bits turned on
        :param bits: Number of folded bits
        """
        gen = AllChem.GetMorganGenerator(fpSize=nbits)
        bv = gen.GetFingerprint(mol)
        return _bitvect_to_dict(bv, f"Morgan{radius * 2}_", nbits)

    @staticmethod
    def calculate_secfp_fp(mol: Chem.Mol, radius: int = 3, nbits: int = 2048) -> dict:
        """Calculate the folded MinHash fingerpirnt (MHFP) of molecule.

        doi: 10.1186/s13321-018-0321-8.
        :param radius: maximum radius of atom-centered substructures.
        :param nbits: Number of folded bits
        """
        bv = MHFPEncoder.secfp_from_mol(mol, length=nbits, radius=radius, rings=True, kekulize=True, min_radius=1)
        return dict(zip(_get_key_list(f"SECFP{radius * 2}_", nbits), bv))

    @staticmethod
    def calculate_minhash_atompair_fp(mol: Chem.Mol, radius: int = 2, nbits: int = 2048) -> dict:
        """Calculate the MinHash fingerprint of Atom Pairs (MAP) of molecule.

        doi: 10.1186/s13321-020-00445-4.
        :param radius: maximum radius of atom-centered substructures.
        :param nbits: Number of folded bits (ignored if rtype != 'bitstring')
        """
        mapcalc = _get_map4_calculator(radius, nbits)
        bv = mapcalc.calculate(mol)
        return dict(zip(_get_key_list(f"MAP{radius * 2}_", nbits), bv))

    @staticmethod
    def calculate_avalon_fp(mol: Chem.Mol, nbits: int = 512) -> dict:
        """Calculate an Avalon fingerprint.

        :param mols: the molecules
        :param nbits: Number of folded bits
        """
        bv = GetAvalonFP(mol, nBits=nbits)
        return _bitvect_to_dict(bv, "Avalon_", nbits)

    @staticmethod
    def get_all_fps(mol: Chem.Mol, radius: Optional[int] = None, nbits: Optional[int] = None) -> dict:
        """Calculate all fingerprints."""
        values = {}
        # FP2/FP3/FP4 each independently rebuild the same OpenBabel molecule from a SMILES
        # round-trip; build it once here and pass it through instead.
        pybel_mol = pybel.readstring("smi", Chem.MolToSmiles(mol))
        values.update(Fingerprint.calculate_fp2_fp(mol, pybel_mol=pybel_mol))
        values.update(Fingerprint.calculate_fp3_fp(mol, pybel_mol=pybel_mol))
        values.update(Fingerprint.calculate_fp4_fp(mol, pybel_mol=pybel_mol))
        for des_label, (func, supported_args) in _fp_funcs.items():
            if des_label in ("FP2", "FP3", "FP4"):
                continue
            if radius is not None and nbits is not None and "radius" in supported_args and "nbits" in supported_args:
                values.update(func(mol, radius=radius, nbits=nbits))
            elif radius is not None and "radius" in supported_args:
                values.update(func(mol, radius=radius))
            elif nbits is not None and "nbits" in supported_args:
                values.update(func(mol, nbits=nbits))
            else:
                values.update(func(mol))
        return values


_fp_funcs = {
    "FP2": (Fingerprint.calculate_fp2_fp, tuple()),
    "FP3": (Fingerprint.calculate_fp3_fp, tuple()),
    "FP4": (Fingerprint.calculate_fp4_fp, tuple()),
    "MACCS": (Fingerprint.calculate_maccs_fp, tuple()),
    "Estate": (Fingerprint.calculate_estate_fp, tuple()),
    "topological": (Fingerprint.calculate_rdk_fp, ("nbits")),
    "atompairs": (Fingerprint.calculate_atompairs_fp, ("nbits")),
    "torsions": (Fingerprint.calculate_topological_torsion_fp, ("nbits")),
    "morgan": (Fingerprint.calculate_morgan_fp, ("radius", "nbits")),
    "SECFP": (Fingerprint.calculate_secfp_fp, ("radius", "nbits")),
    "MAP": (Fingerprint.calculate_minhash_atompair_fp, ("radius", "nbits")),
    "Avalon": (Fingerprint.calculate_avalon_fp, ("nbits")),
}
