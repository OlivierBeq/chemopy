# -*- coding: utf-8 -*-


"""Multiple 2D molecular fingerprints."""


import inspect
from typing import Optional

from map4 import MAP4Calculator
from mhfp.encoder import MHFPEncoder
from openbabel import pybel
from rdkit import Chem, DataStructs
from rdkit.Chem import MACCSkeys
from rdkit.Chem.rdMolDescriptors import (GetHashedAtomPairFingerprintAsBitVect,
                                         GetMorganFingerprintAsBitVect,
                                         GetHashedTopologicalTorsionFingerprintAsBitVect)

from .estate import EState

similaritymeasure = [i[0] for i in DataStructs.similarityFunctions]

class Fingerprint:
    """Molecular 2D fingerprints."""

    @staticmethod
    def calculate_rdk_fp(mol: Chem.Mol, nbits: int = 2048) -> dict:
        """Calculate a topological Daylight-like RDKit fingerprint.

        :param mols: the molecules
        :param nbits: Number of folded bits
        """
        bv = Chem.RDKFingerprint(mol, fpSize=nbits).ToList()
        values = dict(zip([f'RDK_{i+1}' for i in range(nbits)], bv))
        return values

    @staticmethod
    def calculate_maccs_fp(mol: Chem.Mol) -> dict:
        """Calculate the MACCS public keys (166 bits)."""
        bv = MACCSkeys.GenMACCSKeys(mol).ToList()
        values = dict(zip([f'MACCS_{i + 1}' for i in range(166)], bv))
        return values

    @staticmethod
    def calculate_fp2_fp(mol: Chem.Mol) -> dict:
        """Calculate topological FP2 fingerprints (1024 bits)."""
        m = pybel.readstring('smi', Chem.MolToSmiles(mol))
        on_bits = m.calcfp('FP2').bits
        values = dict(zip([f'FP2_{i + 1}' for i in range(1024)],
                          [1 if i + 1 in on_bits else 0 for i in range(1024)]))
        return values

    @staticmethod
    def calculate_fp3_fp(mol: Chem.Mol) -> dict:
        """Calculate FP3 fingerprints (55 bits)."""
        m = pybel.readstring('smi', Chem.MolToSmiles(mol))
        on_bits = m.calcfp('FP3').bits
        values = dict(zip([f'FP3_{i + 1}' for i in range(55)],
                          [1 if i + 1 in on_bits else 0 for i in range(55)]))
        return values

    @staticmethod
    def calculate_fp4_fp(mol: Chem.Mol) -> dict:
        """Calculate FP4 fingerprints (307 bits)."""
        m = pybel.readstring('smi', Chem.MolToSmiles(mol))
        on_bits = m.calcfp('FP4').bits
        values = dict(zip([f'FP4_{i + 1}' for i in range(307)],
                          [1 if i + 1 in on_bits else 0 for i in range(307)]))
        return values

    @staticmethod
    def calculate_estate_fp(mol: Chem.Mol) -> dict:
        """Calculate E-state fingerprints (79 bits)."""
        values = EState.calculate_estate_fingerprint(mol, implementation='chemopy', binary=True)
        return values

    @staticmethod
    def calculate_atompairs_fp(mol: Chem.Mol, nbits: int = 2048) -> dict:
        """Calculate atom pairs fingerprints.

        :param nbits: Number of folded bits
        """
        bv = GetHashedAtomPairFingerprintAsBitVect(mol, nBits=nbits).ToList()
        values = dict(zip([f'AtomPair_{i + 1}' for i in range(nbits)], bv))
        return values

    @staticmethod
    def calculate_topological_torsion_fp(mol: Chem.Mol, nbits: int = 2048) -> dict:
        """Calculate Topological Torsion fingerprints.

        :param nbits: Number of folded bits
        """
        bv = GetHashedTopologicalTorsionFingerprintAsBitVect(mol, nBits=nbits)
        values = dict(zip([f'TopolTorsions_{i + 1}' for i in range(nbits)], bv))
        return values

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
        bv = GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=nbits)
        values = dict(zip([f'Morgan{radius * 2}_{i + 1}' for i in range(nbits)], bv))
        return values

    @staticmethod
    def calculate_secfp_fp(mol: Chem.Mol, radius: int = 3, nbits: int = 2048) -> dict:
        """Calculate the folded MinHash fingerpirnt (MHFP) of molecule.

        doi: 10.1186/s13321-018-0321-8.
        :param radius: maximum radius of atom-centered substructures.
        :param nbits: Number of folded bits
        """
        bv = MHFPEncoder.secfp_from_mol(mol, length=nbits, radius=radius,
                                        rings=True, kekulize=True, min_radius=1)
        values = dict(zip([f'SECFP{radius * 2}_{i + 1}' for i in range(nbits)], bv))
        return values

    @staticmethod
    def calculate_minhash_atompair_fp(mol: Chem.Mol, radius: int = 2, nbits: int = 2048) -> dict:
        """Calculate the MinHash fingerprint of Atom Pairs (MAP) of molecule.

        doi: 10.1186/s13321-020-00445-4.
        :param radius: maximum radius of atom-centered substructures.
        :param nbits: Number of folded bits (ignored if rtype != 'bitstring')
        """
        mapcalc = MAP4Calculator(radius=radius, dimensions=nbits, is_folded=True)
        bv = mapcalc.calculate(mol)
        values = dict(zip([f'MAP{radius * 2}_{i + 1}' for i in range(nbits)], bv))
        return values
    @staticmethod
    def get_all_fps(mol: Chem.Mol, radius: Optional[int] = None, nbits: Optional[int] = None) -> dict:
        """Calculate all fingerprints."""
        values = {}
        for des_label, (func, supported_args) in _fp_funcs.items():
            if radius is not None and nbits is not None and 'radius' in supported_args and 'nbits' in supported_args:
                values.update(func(mol, radius=radius, nbits=nbits))
            elif radius is not None and 'radius' in supported_args:
                values.update(func(mol, radius=radius))
            elif nbits is not None and 'nbits' in supported_args:
                values.update(func(mol, nbits=nbits))
            else:
                values.update(func(mol))
        return values


_fp_funcs = {'FP2': (Fingerprint.calculate_fp2_fp, tuple()),
             'FP3': (Fingerprint.calculate_fp3_fp, tuple()),
             'FP4': (Fingerprint.calculate_fp4_fp, tuple()),
             'MACCS': (Fingerprint.calculate_maccs_fp, tuple()),
             'Estate': (Fingerprint.calculate_estate_fp, tuple()),
             'topological': (Fingerprint.calculate_rdk_fp, ('nbits')),
             'atompairs': (Fingerprint.calculate_atompairs_fp, ('nbits')),
             'torsions': (Fingerprint.calculate_topological_torsion_fp, ('nbits')),
             'morgan': (Fingerprint.calculate_morgan_fp, ('radius', 'nbits')),
             'SECFP': (Fingerprint.calculate_secfp_fp, ('radius', 'nbits')),
             'MAP': (Fingerprint.calculate_minhash_atompair_fp, ('radius', 'nbits'))
             }
