# -*- coding: utf-8 -*-


"""Multiple molecular fingerprints."""

from openbabel import pybel
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, MACCSkeys
from rdkit.Chem.AtomPairs import Pairs, Torsions
from rdkit.Chem.Fingerprints import FingerprintMols

from pychem.estate import CalculateEstateFingerprint as EstateFingerprint

similaritymeasure = [i[0] for i in DataStructs.similarityFunctions]


# TODO: wrap fingerprints in a dedicated class
#       result is a tuple form. The first is the number of
#       fingerprints. The second is a dict form whose keys are the
#       position which this molecule has some substructure. The third
#       is the DataStructs which is used for calculating the similarity.


def CalculateDaylightFingerprint(mol: Chem.Mol):
    """Calculate Daylight-like fingerprint or topological fingerprint (2048 bits)."""
    res = {}
    NumFinger = 2048
    bv = FingerprintMols.FingerprintMol(mol)
    temp = tuple(bv.GetOnBits())
    for i in temp:
        res.update({i: 1})
    return NumFinger, res, bv


def CalculateMACCSFingerprint(mol: Chem.Mol):
    """Calculate MACCS public keys (166 bits)."""
    res = {}
    NumFinger = 166
    bv = MACCSkeys.GenMACCSKeys(mol)
    temp = tuple(bv.GetOnBits())
    for i in temp:
        res.update({i: 1})
    return NumFinger, res, bv


def CalculateFP4Fingerprint(mol: Chem.Mol):
    """Calculate FP4 fingerprints (307 bits)."""
    res = {}
    NumFinger = 307
    m = pybel.readstring('smi', Chem.MolToSmiles(mol))
    temp = m.calcfp('FP4').bits
    for i in temp:
        res.update({i: 1})
    vec = DataStructs.ExplicitBitVect(307)
    vec.SetBitsFromList([x - 1 for x in temp])
    return NumFinger, res, vec


def CalculateEstateFingerprint(mol: pybel.Molecule):
    """Calculate E-state fingerprints (79 bits)."""
    NumFinger = 79
    res = {}
    temp = EstateFingerprint(mol)
    for i in temp:
        if temp[i] > 0:
            res[i[7:]] = 1
    return NumFinger, res, temp


def CalculateAtomPairsFingerprint(mol: Chem.Mol):
    """Calculate atom pairs fingerprints."""
    res = Pairs.GetAtomPairFingerprint(mol)
    return res.GetLength(), res.GetNonzeroElements(), res


def CalculateTopologicalTorsionFingerprint(mol: Chem.Mol):
    """Calculate Topological Torsion fingerprints."""
    res = Torsions.GetTopologicalTorsionFingerprint(mol)
    return res.GetLength(), res.GetNonzeroElements(), res


def CalculateMorganFingerprint(mol: Chem.Mol, radius=2):
    """Calculate Morgan fingerprints."""
    res = AllChem.GetMorganFingerprint(mol, radius)
    return res.GetLength(), res.GetNonzeroElements(), res


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
