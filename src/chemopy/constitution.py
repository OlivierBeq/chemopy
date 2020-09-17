# -*- coding: utf-8 -*-


"""Molecular constitutional and topological indices."""

from rdkit import Chem
from rdkit.Chem import Lipinski as LPK


def CalculateHeavyMolWeight(mol: Chem.Mol) -> float:
    """Calculate molecular weight of heavy atoms."""
    MolWeight = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 1:
            MolWeight += atom.GetMass()
    return MolWeight


def CalculateAverageMolWeight(mol: Chem.Mol) -> float:
    """Calculate average molecular weight of heavy atoms."""
    MolWeight = 0
    for atom in mol.GetAtoms():
        MolWeight += atom.GetMass()
    return MolWeight / mol.GetNumAtoms()


def CalculateHydrogenNumber(mol: Chem.Mol) -> float:
    """Calculate number of Hydrogens."""
    i = 0
    Hmol = Chem.AddHs(mol)
    for atom in Hmol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            i += 1
    return i


def CalculateHalogenNumber(mol: Chem.Mol) -> float:
    """Calculate number of Halogens."""
    i = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in [9, 17, 35, 53]:
            i += 1
    return i


def CalculateHeteroNumber(mol: Chem.Mol) -> float:
    """Calculate number of Heteroatoms."""
    i = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6]:
            i += 1
    return mol.GetNumAtoms() - i


def CalculateHeavyAtomNumber(mol: Chem.Mol) -> float:
    """Calculate number of Heavy atoms."""
    return mol.GetNumHeavyAtoms()


def _CalculateElementNumber(mol: Chem.Mol, AtomicNumber=6) -> float:
    """Calculate number of atoms with specified element."""
    i = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == AtomicNumber:
            i += 1
    return i


def CalculateFluorinNumber(mol: Chem.Mol) -> float:
    """Calculate number of Fluorine atoms."""
    return _CalculateElementNumber(mol, AtomicNumber=9)


def CalculateChlorinNumber(mol: Chem.Mol) -> float:
    """Calculate number of Fluorine atoms."""
    return _CalculateElementNumber(mol, AtomicNumber=17)


def CalculateBromineNumber(mol: Chem.Mol) -> float:
    """Calculate number of Bromine atoms."""
    return _CalculateElementNumber(mol, AtomicNumber=35)


def CalculateIodineNumber(mol: Chem.Mol) -> float:
    """Calculate number of Iodine atoms."""
    return _CalculateElementNumber(mol, AtomicNumber=53)


def CalculateCarbonNumber(mol: Chem.Mol) -> float:
    """Calculate number of Carbon atoms."""
    return _CalculateElementNumber(mol, AtomicNumber=6)


def CalculatePhosphorNumber(mol: Chem.Mol) -> float:
    """Calcualtion number of Phosphor atoms."""
    return _CalculateElementNumber(mol, AtomicNumber=15)


def CalculateSulfurNumber(mol: Chem.Mol) -> float:
    """Calculate number of Sulfur atoms."""
    return _CalculateElementNumber(mol, AtomicNumber=16)


def CalculateOxygenNumber(mol: Chem.Mol) -> float:
    """Calculate number of Oxygen atoms."""
    return _CalculateElementNumber(mol, AtomicNumber=8)


def CalculateNitrogenNumber(mol: Chem.Mol) -> float:
    """Calculate number of Nitrogen atoms."""
    return _CalculateElementNumber(mol, AtomicNumber=7)


def CalculateRingNumber(mol: Chem.Mol) -> float:
    """Calculate number of rings."""
    return Chem.GetSSSR(mol)


def CalculateRotationBondNumber(mol: Chem.Mol) -> float:
    """Calculate number of rotatable bonds."""
    return LPK.NumRotatableBonds(mol)


def CalculateHdonorNumber(mol: Chem.Mol) -> float:
    """Calculate number of Hydrongen bond donors."""
    return LPK.NumHDonors(mol)


def CalculateHacceptorNumber(mol: Chem.Mol) -> float:
    """Calculate number of Hydrogen bond acceptors."""
    return LPK.NumHAcceptors(mol)


def CalculateSingleBondNumber(mol: Chem.Mol) -> float:
    """Calculate number of single bonds."""
    i = 0
    for bond in mol.GetBonds():
        if bond.GetBondType().name == 'SINGLE':
            i += 1
    return i


def CalculateDoubleBondNumber(mol: Chem.Mol) -> float:
    """Calculate number of double bonds."""
    i = 0
    for bond in mol.GetBonds():
        if bond.GetBondType().name == 'DOUBLE':
            i += 1
    return i


def CalculateTripleBondNumber(mol: Chem.Mol) -> float:
    """Calculate number of triple bonds."""
    i = 0
    for bond in mol.GetBonds():
        if bond.GetBondType().name == 'TRIPLE':
            i += 1
    return i


def CalculateAromaticBondNumber(mol: Chem.Mol) -> float:
    """Calculate number of aromatic bonds."""
    i = 0
    for bond in mol.GetBonds():
        if bond.GetBondType().name == 'AROMATIC':
            i += 1
    return i


def CalculateAllAtomNumber(mol: Chem.Mol) -> float:
    """Calculate number of all atoms."""
    return Chem.AddHs(mol).GetNumAtoms()


def _CalculatePathN(mol: Chem.Mol, PathLength=2) -> float:
    """Calculate number of path of length N."""
    return len(Chem.FindAllPathsOfLengthN(mol, PathLength, useBonds=1))


def CalculatePath1(mol: Chem.Mol) -> float:
    """Calculate number of path length of 1."""
    return _CalculatePathN(mol, 1)


def CalculatePath2(mol: Chem.Mol) -> float:
    """Calculate number of path length of 2."""
    return _CalculatePathN(mol, 2)


def CalculatePath3(mol: Chem.Mol) -> float:
    """Calculate number of path length of 3."""
    return _CalculatePathN(mol, 3)


def CalculatePath4(mol: Chem.Mol) -> float:
    """Calculate number of path length of 4."""
    return _CalculatePathN(mol, 4)


def CalculatePath5(mol: Chem.Mol) -> float:
    """Calculate number of path length of 5."""
    return _CalculatePathN(mol, 5)


def CalculatePath6(mol: Chem.Mol) -> float:
    """Calculate number of path length of 6."""
    return _CalculatePathN(mol, 6)


_constitutional = {'Weight': CalculateHeavyMolWeight,
                   'AWeight': CalculateAverageMolWeight,
                   'nH': CalculateHydrogenNumber,
                   'nHal': CalculateHalogenNumber,
                   'nHet': CalculateHeteroNumber,
                   'nHA': CalculateHeavyAtomNumber,
                   'nF': CalculateFluorinNumber,
                   'nCl': CalculateChlorinNumber,
                   'nBr': CalculateBromineNumber,
                   'nI': CalculateIodineNumber,
                   'nC': CalculateCarbonNumber,
                   'nP': CalculatePhosphorNumber,
                   'nS': CalculateOxygenNumber,
                   'nO': CalculateOxygenNumber,
                   'nN': CalculateNitrogenNumber,
                   'nRing': CalculateRingNumber,
                   'nRotB': CalculateRotationBondNumber,
                   'nHBD': CalculateHdonorNumber,
                   'nHBA': CalculateHacceptorNumber,
                   'nSBond': CalculateSingleBondNumber,
                   'nDBond': CalculateDoubleBondNumber,
                   'nAroBond': CalculateAromaticBondNumber,
                   'nTBond': CalculateTripleBondNumber,
                   'nAtom': CalculateAllAtomNumber,
                   'PathL1': CalculatePath1,
                   'PathL2': CalculatePath2,
                   'PathL3': CalculatePath3,
                   'PathL4': CalculatePath4,
                   'PathL5': CalculatePath5,
                   'PathL6': CalculatePath6}


def GetConstitutional(mol: Chem.Mol) -> dict:
    """Get all (30) constitutional descriptors."""
    result = {}
    for DesLabel in _constitutional.keys():
        result[DesLabel] = round(_constitutional[DesLabel](mol), 3)
    return result


#############################
# if __name__ =='__main__':
#     smis = ['CCCC','CCCCC','CCCCCC','CC(N)C(=O)O','CC(N)C(=O)[O-].[Na+]']
#     smi5=['CCCCCC','CCC(C)CC','CC(C)CCC','CC(C)C(C)C','CCCCCN','c1ccccc1N']
#     for index, smi in enumerate(smis):
#         m = Chem.MolFromSmiles(smi)
#         print(index+1)
#         print(smi)
#         print('\t',GetConstitutional(m))
#         print(len(GetConstitutional(m)))
