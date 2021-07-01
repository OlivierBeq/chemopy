# -*- coding: utf-8 -*-


"""Multiple molecular fingerprints."""


from typing import Any, Tuple, Union

import numpy as np
from map4 import MAP4Calculator
from mhfp.encoder import MHFPEncoder
from openbabel import pybel
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, MACCSkeys
from rdkit.Chem.AtomPairs import Pairs, Torsions
from rdkit.Chem.Fingerprints import FingerprintMols

from chemopy.estate import CalculateEstateFingerprint as EstateFingerprint

similaritymeasure = [i[0] for i in DataStructs.similarityFunctions]


# TODO: wrap fingerprints in a dedicated class
#       result is a tuple form. The first is the number of
#       fingerprints. The second is a dict form whose keys are the
#       position which this molecule has some substructure. The third
#       is the DataStructs which is used for calculating the similarity.


def CalculateDaylightFingerprint(mol: Chem.Mol, rtype: str = 'bitstring', bits: int = 2048) -> Union[str, dict, Any]:
    """Calculate Daylight-like fingerprint or topological fingerprint.

    :param rtype: Type of output, may either be:
                  bitstring (default), returns a binary string
                  rdkit, return the native rdkit DataStructs
                  dict, for a dict of bits turned on
    :param bits: Number of folded bits (ignored if rtype != 'bitstring')
    """
    bv = FingerprintMols.FingerprintMol(mol)
    if rtype == 'rdkit':
        return bv
    elif rtype == 'dict':
        return {f'Daylight{i+1}': 1 for i in bv.GetOnBits()}
    return bv.ToBitString()


def CalculateMACCSFingerprint(mol: Chem.Mol, rtype: str = 'bitstring') -> Tuple[str, dict, Any]:
    """Calculate MACCS public keys (166 bits).

    :param rtype: Type of output, may either be:
                  bitstring (default), returns a binary string
                  rdkit, return the native rdkit DataStructs
                  dict, for a dict of bits turned on
    """
    bv = MACCSkeys.GenMACCSKeys(mol)
    if rtype == 'rdkit':
        return bv
    elif rtype == 'dict':
        return {f'MACCS{i+1}': 1 for i in bv.GetOnBits()}
    return bv.ToBitString()


def CalculateFP2Fingerprint(mol: Chem.Mol, rtype: str = 'bitstring') -> Tuple[str, dict, Any]:
    """Calculate topological FP2 fingerprints (1024 bits).

    :param rtype: Type of output, may either be:
                  bitstring (default), returns a binary string
                  rdkit, return the native rdkit DataStructs
                  dict, for a dict of bits turned on
    """
    m = pybel.readstring('smi', Chem.MolToSmiles(mol))
    temp = m.calcfp('FP2').bits
    if rtype == 'dict':
        return {f'FP2_{i}': 1 for i in temp}
    bv = DataStructs.ExplicitBitVect(1024)
    bv.SetBitsFromList([x - 1 for x in temp])
    if rtype == 'rdkit':
        return bv
    return bv.ToBitString()


def CalculateFP3Fingerprint(mol: Chem.Mol, rtype: str = 'bitstring') -> Tuple[str, dict, Any]:
    """Calculate FP3 fingerprints (55 bits).

    :param rtype: Type of output, may either be:
                  bitstring (default), returns a binary string
                  rdkit, return the native rdkit DataStructs
                  dict, for a dict of bits turned on
    """
    m = pybel.readstring('smi', Chem.MolToSmiles(mol))
    temp = m.calcfp('FP3').bits
    if rtype == 'dict':
        return {f'FP3_{i}': 1 for i in temp}
    bv = DataStructs.ExplicitBitVect(55)
    bv.SetBitsFromList([x - 1 for x in temp])
    if rtype == 'rdkit':
        return bv
    return bv.ToBitString()


def CalculateFP4Fingerprint(mol: Chem.Mol, rtype: str = 'bitstring') -> Tuple[str, dict, Any]:
    """Calculate FP4 fingerprints (307 bits).

    :param rtype: Type of output, may either be:
                  bitstring (default), returns a binary string
                  rdkit, return the native rdkit DataStructs
                  dict, for a dict of bits turned on
    """
    m = pybel.readstring('smi', Chem.MolToSmiles(mol))
    temp = m.calcfp('FP3').bits
    if rtype == 'dict':
        return {f'FP3_{i}': 1 for i in temp}
    bv = DataStructs.ExplicitBitVect(307)
    bv.SetBitsFromList([x - 1 for x in temp])
    if rtype == 'rdkit':
        return bv
    return bv.ToBitString()


def CalculateEstateFingerprint(mol: pybel.Molecule, rtype: str = 'doublestring') -> Tuple[str, dict, Any]:
    """Calculate E-state fingerprints (79 bits).

    :param rtype: Type of output, may either be:
                  doublestring (default), returns a binary string
                  dict, for a dict of bits turned on
    """
    if rtype == 'dict':
        res = {k: v for k, v in EstateFingerprint(mol).items() if v != 0}
        return res
    return ';'.join(list(EstateFingerprint(mol).values()))


def CalculateAtomPairsFingerprint(mol: Chem.Mol, rtype: str = 'countstring', bits: int = 2048) -> Tuple[str, dict, Any]:
    """Calculate atom pairs fingerprints.

    :param rtype: Type of output, may either be:
                  countstring (default), returns a binary string
                  rdkit, return the native rdkit DataStructs
                  dict, for a dict of bits turned on
    :param bits: Number of folded bits (ignored if rtype != 'countstring')
    """
    res = Pairs.GetAtomPairFingerprint(mol)
    if rtype == 'rdkit':
        return res
    counts = res.GetNonzeroElements()
    if rtype == 'dict':
        return {f'AtomPair_{k}': v for k, v in counts.items()}
    folded = np.zeros(bits)
    for k, v in counts.items():
        folded[k % bits] += v
    return ';'.join(folded.tolist())


def CalculateTopologicalTorsionFingerprint(mol: Chem.Mol, rtype: str = 'countstring',
                                           bits: int = 2048) -> Tuple[str, dict, Any]:
    """Calculate Topological Torsion fingerprints.

    :param rtype: Type of output, may either be:
                  countstring (default), returns a binary string
                  rdkit, return the native rdkit DataStructs
                  dict, for a dict of bits turned on
    :param bits: Number of folded bits (ignored if rtype != 'countstring')
    """
    res = Torsions.GetTopologicalTorsionFingerprint(mol)
    if rtype == 'rdkit':
        return res
    counts = res.GetNonzeroElements()
    if rtype == 'dict':
        return {f'TopolTorsions_{k}': v for k, v in counts.items()}
    folded = np.zeros(bits)
    for k, v in counts.items():
        folded[k % bits] += v
    return ';'.join(folded.tolist())


def CalculateMorganFingerprint(mol: Chem.Mol, radius=2, rtype: str = 'bitstring',
                               bits: int = 2048) -> Tuple[str, dict, Any]:
    """Calculate Morgan fingerprints.

    :param radius: maximum radius of atom-centered substructures.
    :param rtype: Type of output, may either be:
                  bitstring (default), returns a binary string
                  rdkit, return the native rdkit DataStructs
                  dict, for a dict of bits turned on
    :param bits: Number of folded bits (ignored if rtype != 'bitstring')
    """
    bv = AllChem.GetMorganFingerprintAsBitVect(mol, radius, bits)
    if rtype == 'rdkit':
        return bv
    elif rtype == 'dict':
        return {f'Morgan{i+1}': 1 for i in bv.GetOnBits()}
    else:
        return bv.ToBitString()


def CalculateMinHashFingerprint(mol: Chem.Mol, radius: int = 3, rtype: str = 'bitstring',
                                bits: int = 2048) -> Tuple[str, dict, Any]:
    """Calculate the MinHash Fingerprint (MHFP) of molecule.

    doi: 10.1186/s13321-018-0321-8.
    :param radius: maximum radius of atom-centered substructures.
    :param rtype: Type of output, may either be:
                  bitstring (default), returns a binary string
                  numpy, return the underlying numpy array
                  dict, for a dict of bits turned on
    :param bits: Number of folded bits (ignored if rtype != 'bitstring')
    """
    mhfp = MHFPEncoder()
    shingles = mhfp.shingling_from_mol(mol, radius, True, True, 1)
    hash_values = mhfp.hash(shingles)
    if rtype == 'numpy':
        return hash_values
    elif rtype == 'dict':
        return {x: 1 for x in hash_values.tolist()}
    else:
        folded = mhfp.fold(hash_values, bits)
        return ''.join(map(str, folded))


def CalculateMinHashAtomPairFingerprint(mol: Chem.Mol, radius: int = 2, rtype: str = 'bitstring',
                                        bits: int = 2048) -> Tuple[str, dict, Any]:
    """Calculate the MinHash Fingerprint of Atom Pairs (MAP) of molecule.

    doi: 10.1186/s13321-020-00445-4.
    :param radius: maximum radius of atom-centered substructures.
    :param rtype: Type of output, may either be:
                  bitstring (default), returns a binary string
                  numpy, return the underlying numpy array
                  dict, for a dict of bits turned on
    :param bits: Number of folded bits (ignored if rtype != 'bitstring')
    """
    if rtype == 'numpy':
        mapcalc = MAP4Calculator(radius=radius, is_folded=False)
        return mapcalc.calculate(mol)
    elif rtype == 'dict':
        mapcalc = MAP4Calculator(radius=radius, is_folded=False)
        return {x: 1 for x in mapcalc.calculate(mol).tolist()}
    else:
        folded = MAP4Calculator(radius=radius, is_folded=True, dimensions=bits).tolist()
        return ''.join(map(str, folded))


def CalculateSmilesExtendedConnectivityFingerprint(mol: Chem.Mol, radius: int = 2, rtype: str = 'bitstring',
                                                   bits: int = 2048) -> Tuple[str, dict, Any]:
    """Calculate SMILES extended connectivity fingerprint (SECFP), doi: 10.1186/s13321-018-0321-8.

    :param radius: maximum radius of atom-centered substructures.
    :param rtype: Type of output, may either be:
                  bitstring (default), returns a binary string
                  numpy, return the underlying numpy array
                  rdkit, return the native rdkit DataStructs
                  dict, for a dict of bits turned on
    :param bits: Number of folded bits (ignored if rtype != 'bitstring')
    """
    secfp = MHFPEncoder.secfp_from_mol(mol, length=bits, radius=radius, rings=True, kekulize=True, min_radius=1)
    if rtype == 'numpy':
        return secfp
    elif rtype == 'dict':
        return {x: 1 for x in secfp.tolist() if x != 0}
    bv = DataStructs.ExplicitBitVect(bits)
    bv.SetBitsFromList([x for x, y in enumerate(secfp.tolist()) if y != 0])
    if rtype == 'rdkit':
        return bv
    else:
        return bv.ToBitString()


def CalculateSimilarity(fp1, fp2, similarity="Tanimoto") -> float:
    """Calculate similarity between two molecules.

    :param fp1: DataStructs fingerprint.
    :param fp2: DataStructs fingerprint.
    """
    temp = DataStructs.similarityFunctions
    for i in temp:
        if similarity in i[0]:
            similarityfunction = i[1]
        else:
            similarityfunction = temp[0][1]
    res = similarityfunction(fp1, fp2)
    return round(res, 3)


_FingerprintFuncs = {'topological': CalculateDaylightFingerprint,
                     'Estate': CalculateEstateFingerprint,
                     'FP2': CalculateFP2Fingerprint,
                     'FP3': CalculateFP3Fingerprint,
                     'FP4': CalculateFP4Fingerprint,
                     'atompairs': CalculateAtomPairsFingerprint,
                     'torsions': CalculateTopologicalTorsionFingerprint,
                     'morgan': CalculateMorganFingerprint,
                     'MACCS': CalculateMACCSFingerprint}


################################################################
# if __name__=="__main__":
#     print("fingerprint......")
#     ms = [Chem.MolFromSmiles('CCOC=N'), Chem.MolFromSmiles('CCO')]
#     res1=CalculateTopologicalTorsionFingerprint(ms[0])
#     print(res1)
#     res2=CalculateTopologicalTorsionFingerprint(ms[1])
#     print(res2)
#     print(CalculateSimilarity(res1[2],res2[2]))
