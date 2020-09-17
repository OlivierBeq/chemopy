# -*- coding: utf-8 -*-


"""Molecular topological indices."""

from typing import Iterable, List

import numpy
import scipy
from rdkit import Chem
from rdkit.Chem import GraphDescriptors as GD
from rdkit.Chem import rdchem

periodicTable = rdchem.GetPeriodicTable()


def _GetPrincipleQuantumNumber(atNum: int) -> int:
    """Get principle quantum number from atomic number."""
    if atNum <= 2:
        return 1
    elif atNum <= 10:
        return 2
    elif atNum <= 18:
        return 3
    elif atNum <= 36:
        return 4
    elif atNum <= 54:
        return 5
    elif atNum <= 86:
        return 6
    else:
        return 7


def CalculateWeiner(mol: Chem.Mol) -> float:
    """Get Weiner values of a molecule.

    Or W.
    """
    return 1.0 / 2 * sum(sum(Chem.GetDistanceMatrix(mol)))


def CalculateMeanWeiner(mol: Chem.Mol) -> float:
    """Get Mean Weiner index of a molecule.

    Or AW.
    """
    N = mol.GetNumAtoms()
    WeinerNumber = CalculateWeiner(mol)
    return 2.0 * WeinerNumber / (N * (N - 1))


def CalculateBalaban(mol: Chem.Mol) -> float:
    """Get Balaban index of a molecule.

    Or J.
    """
    adjMat = Chem.GetAdjacencyMatrix(mol)
    Distance = Chem.GetDistanceMatrix(mol)
    Nbond = mol.GetNumBonds()
    Natom = mol.GetNumAtoms()
    S = numpy.sum(Distance, axis=1)
    mu = Nbond - Natom + 1
    sumk = 0.
    for i in range(len(Distance)):
        si = S[i]
        for j in range(i, len(Distance)):
            if adjMat[i, j] == 1:
                sumk += 1. / numpy.sqrt(si * S[j])
    if mu + 1 != 0:
        J = float(Nbond) / float(mu + 1) * sumk
    else:
        J = 0
    return J


def CalculateGraphDistance(mol: Chem.Mol) -> float:
    """Get graph distance index.

    Or Tigdi.
    """
    Distance = Chem.GetDistanceMatrix(mol)
    n = int(Distance.max())
    res = 0.0
    for i in range(n):
        temp = 1. / 2 * sum(sum(Distance == i + 1))
        res = res + temp ** 2
    return numpy.log10(res)


def CalculateDiameter(mol: Chem.Mol) -> float:
    """Get largest value of the distance matrix.

    Or diametert.
    From Petitjean, M. J. Chem. Inf. Comput. Sci. 1992, 32, 4, 331-337.
    """
    Distance = Chem.GetDistanceMatrix(mol)
    return Distance.max()


def CalculateRadius(mol: Chem.Mol) -> float:
    """Get radius based on topology.

    Or radiust.
    From Petitjean, M. J. Chem. Inf. Comput. Sci. 1992, 32, 4, 331-337.
    """
    Distance = Chem.GetDistanceMatrix(mol)
    temp = []
    for i in Distance:
        temp.append(max(i))
    return min(temp)


def CalculatePetitjean(mol: Chem.Mol) -> float:
    """Get Petitjean based on topology.

    Or petitjeant.
    """
    diameter = CalculateDiameter(mol)
    radius = CalculateRadius(mol)
    return 1 - radius / float(diameter)


def CalculateXuIndex(mol: Chem.Mol) -> float:
    """Get Xu index.

    Or Xu.
    """
    nAT = mol.GetNumAtoms()
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    Distance = Chem.GetDistanceMatrix(mol)
    sigma = scipy.sum(Distance, axis=1)
    temp1 = 0.0
    temp2 = 0.0
    for i in range(nAT):
        temp1 = temp1 + deltas[i] * ((sigma[i]) ** 2)
        temp2 = temp2 + deltas[i] * (sigma[i])
    Xu = numpy.sqrt(nAT) * numpy.log(temp1 / temp2)
    return Xu


def CalculateGutmanTopo(mol: Chem.Mol) -> float:
    """Get Gutman molecular topological simple vertex index.

    Or GMTI.
    """
    nAT = mol.GetNumAtoms()
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    Distance = Chem.GetDistanceMatrix(mol)
    res = 0.0
    for i in range(nAT):
        for j in range(i + 1, nAT):
            res = res + deltas[i] * deltas[j] * Distance[i, j]
    return numpy.log10(res)


def CalculatePolarityNumber(mol: Chem.Mol) -> float:
    """Get Polarity number.

    Or Pol.
    """
    Distance = Chem.GetDistanceMatrix(mol)
    res = 1. / 2 * sum(sum(Distance == 3))
    return res


def CalculatePoglianiIndex(mol: Chem.Mol) -> float:
    """Get Poglicani index.

    Or DZ.
    From Pogliani L. J.Phys.Chem. (1996), 100,18065-18077.
    """
    res = 0.0
    for atom in mol.GetAtoms():
        n = atom.GetAtomicNum()
        nV = periodicTable.GetNOuterElecs(n)
        mP = _GetPrincipleQuantumNumber(n)
        res = res + (nV + 0.0) / mP
    return res


def CalculateIpc(mol: Chem.Mol) -> float:
    """Get Bonchev-Trinajstic complexity index.

    Or Ipc.
    From Bonchev D. & Trinajstic N., J. Chem. Phys. (1977) 67,4517-4533.
    """
    return numpy.log10(GD.Ipc(mol))


def CalculateBertzCT(mol: Chem.Mol) -> float:
    """Get Bertz complexity index.

    Or BertzCT.
    From Bertz S. H., J. Am. Chem. Soc. (1981) 103,3599-3601.
    """
    return numpy.log10(GD.BertzCT(mol))


def CalculateHarary(mol: Chem.Mol) -> float:
    """Get Harary number.

    Or Thara.
    """
    Distance = numpy.array(Chem.GetDistanceMatrix(mol), 'd')
    return 1.0 / 2 * (sum(1.0 / Distance[Distance != 0]))


def CalculateSchiultz(mol: Chem.Mol) -> float:
    """Get Schiultz number.

    Or Tsch.
    """
    Distance = numpy.array(Chem.GetDistanceMatrix(mol), 'd')
    Adjacent = numpy.array(Chem.GetAdjacencyMatrix(mol), 'd')
    VertexDegree = sum(Adjacent)
    return sum(scipy.dot((Distance + Adjacent), VertexDegree))


def CalculateZagreb1(mol: Chem.Mol) -> float:
    """Get Zagreb index with order 1.

    Or ZM1.
    """
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    return sum(numpy.array(deltas) ** 2)


def CalculateZagreb2(mol: Chem.Mol) -> float:
    """Get Zagreb index with order 2.

    Or ZM2.
    """
    ke = [x.GetBeginAtom().GetDegree() * x.GetEndAtom().GetDegree() for x in mol.GetBonds()]
    return sum(ke)


def CalculateMZagreb1(mol: Chem.Mol) -> float:
    """Get Modified Zagreb index with order 1.

    Or MZM1.
    """
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    while 0 in deltas:
        deltas.remove(0)
    deltas = numpy.array(deltas, 'd')
    res = sum((1. / deltas) ** 2)
    return res


def CalculateMZagreb2(mol: Chem.Mol) -> float:
    """Get Modified Zagreb index with order 2.

    Or MZM2.
    """
    cc = [x.GetBeginAtom().GetDegree() * x.GetEndAtom().GetDegree() for x in mol.GetBonds()]
    if len(cc) == 0:
        return 0.0
    while 0 in cc:
        cc.remove(0)
    cc = numpy.array(cc, 'd')
    res = sum((1. / cc) ** 2)
    return res


def CalculateQuadratic(mol: Chem.Mol) -> float:
    """Get Quadratic index.

    Or Qindex.
    """
    M = CalculateZagreb1(mol)
    N = mol.GetNumAtoms()
    return 3 - 2 * N + M / 2.0


def CalculatePlatt(mol: Chem.Mol) -> float:
    """Get Platt number.

    Or Platt.
    """
    cc = [x.GetBeginAtom().GetDegree() + x.GetEndAtom().GetDegree() - 2 for x in mol.GetBonds()]
    return sum(cc)


def CalculateSimpleTopoIndex(mol: Chem.Mol) -> float:
    """Get the logarithm of the simple topological index.

    Or Sito.
    From Narumi H., MATCH (Comm. Math. Comp. Chem.), (1987), 22,195-207.
    """
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    while 0 in deltas:
        deltas.remove(0)
    deltas = numpy.array(deltas, 'd')
    res = numpy.prod(deltas)
    return numpy.log(res)


def CalculateHarmonicTopoIndex(mol: Chem.Mol) -> float:
    """Get harmonic topological index.

    Or Hato.
    From Narumi H., MATCH (Comm. Math. Comp. Chem.), (1987), 22,195-207.
    """
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    while 0 in deltas:
        deltas.remove(0)
    if len(deltas) == 0:
        return 0.0
    deltas = numpy.array(deltas, 'd')
    nAtoms = mol.GetNumAtoms()
    res = nAtoms / sum(1. / deltas)
    return res


def CalculateGeometricTopoIndex(mol: Chem.Mol) -> float:
    """Get Geometric topological index.

    Or Geto.
    From Narumi H., MATCH (Comm. Math. Comp. Chem.), (1987), 22,195-207.
    """
    nAtoms = mol.GetNumAtoms()
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    while 0 in deltas:
        deltas.remove(0)
    if len(deltas) == 0:
        return 0.0
    deltas = numpy.array(deltas, 'd')
    temp = numpy.prod(deltas)
    res = numpy.power(temp, 1. / nAtoms)
    return res


def CalculateArithmeticTopoIndex(mol: Chem.Mol) -> float:
    """Get Arithmetic topological index.

    Or Arto.
    From Narumi H., MATCH (Comm. Math. Comp. Chem.), (1987), 22,195-207.
    """
    nAtoms = mol.GetNumAtoms()
    nBonds = mol.GetNumBonds()
    res = 2. * nBonds / nAtoms
    return res


def CalculateMolSizeTotalInf(mol: Chem.Mol) -> float:
    """Get total information index on molecular size.

    Or ISIZ.
    """
    Hmol = Chem.AddHs(mol)
    nAT = Hmol.GetNumAtoms()
    ISIZ = nAT * numpy.log2(nAT)
    return ISIZ


def CalculateAtomCompTotalInf(mol: Chem.Mol) -> float:
    """Ge total information index on atomic composition.

    Or TIAC.
    """
    Hmol = Chem.AddHs(mol)
    nAtoms = Hmol.GetNumAtoms()
    IC = []
    for i in range(nAtoms):
        at = Hmol.GetAtomWithIdx(i)
        IC.append(at.GetAtomicNum())
    Unique = numpy.unique(IC)
    NAtomType = len(Unique)
    res = 0.0
    for i in range(NAtomType):
        cc = IC.count(Unique[i])
        res += cc * numpy.log2(cc)
    if nAtoms != 0:
        return nAtoms * numpy.log2(nAtoms) - res
    else:
        return 0.0


def CalculateDistanceEqualityTotalInf(mol: Chem.Mol) -> float:
    """Get total information index on distance equality.

    Or DET.
    """
    Distance = Chem.GetDistanceMatrix(mol)
    nAT = mol.GetNumAtoms()
    n = 1. / 2 * nAT ** 2 - nAT
    DisType = int(Distance.max())
    res = 0.0
    for i in range(DisType):
        cc = 1. / 2 * sum(sum(Distance == i + 1))
        res += cc * numpy.log2(cc)
    return n * numpy.log2(n) - res


def _CalculateEntropy(Probability: Iterable[float]) -> float:
    """Calculate entropy (Information content) of given probability."""
    res = 0.0
    for i in Probability:
        if i != 0:
            res = res - i * numpy.log2(i)
    return res


def CalculateDistanceEqualityMeanInf(mol: Chem.Mol) -> float:
    """Get the mean information index on distance equality.

    Or IDE.
    """
    Distance = Chem.GetDistanceMatrix(mol)
    nAT = mol.GetNumAtoms()
    n = 1. / 2 * nAT ** 2 - nAT
    DisType = int(Distance.max())
    res = 0.0
    cc = numpy.zeros(DisType, numpy.float)
    for i in range(DisType):
        cc[i] = 1. / 2 * sum(sum(Distance == i + 1))
    res = _CalculateEntropy(cc / n)
    return res


def CalculateVertexEqualityTotalInf(mol: Chem.Mol) -> float:
    """Get the total information index on vertex equality.

    Or IVDE.
    """
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    res = 0.0
    while 0 in deltas:
        deltas.remove()
    for i in range(max(deltas)):
        cc = deltas.count(i + 1)
        if cc == 0:
            res = res
        else:
            res += cc * numpy.log2(cc)
    n = len(deltas)
    return n * numpy.log2(n) - res


def _HKDeltas(mol: Chem.Mol, skipHs: bool = True) -> List[float]:
    """Calculate Kier & Hall valence delta-values for molecular connectivity.

    From Kier L. and Hall L., J. Pharm. Sci. (1983), 72(10),1170-1173.
    """
    global periodicTable
    res = []
    for atom in mol.GetAtoms():
        n = atom.GetAtomicNum()
        if n > 1:
            nV = periodicTable.GetNOuterElecs(n)
            nHs = atom.GetTotalNumHs()
            if n < 10:
                res.append(float(nV - nHs))
            else:
                res.append(float(nV - nHs) / float(n - nV - 1))
        elif not skipHs:
            res.append(0.0)
    return res


def CalculateSimpleTopovIndex(mol: Chem.Mol) -> float:
    """Get the logarithm of the simple topological index.

    Or Sitov.
    From Narumi H., MATCH (Comm. Math. Comp. Chem.), (1987), 22,195-207.

    Kier and Hall's valence delta-values are used in place of atom degrees.
    From Kier L. and Hall L., J. Pharm. Sci. (1983), 72(10),1170-1173.
    """
    deltas = _HKDeltas(mol, skipHs=0)
    while 0 in deltas:
        deltas.remove(0)
    deltas = numpy.array(deltas, 'd')
    res = numpy.prod(deltas)
    return numpy.log(res)


def CalculateHarmonicTopovIndex(mol: Chem.Mol) -> float:
    """Get harmonic topological index.

    Or Hatov.
    From Narumi H., MATCH (Comm. Math. Comp. Chem.), (1987), 22,195-207.

    Kier and Hall's valence delta-values are used in place of atom degrees.
    From Kier L. and Hall L., J. Pharm. Sci. (1983), 72(10),1170-1173.
    """
    deltas = _HKDeltas(mol, skipHs=0)
    while 0 in deltas:
        deltas.remove(0)
    deltas = numpy.array(deltas, 'd')
    nAtoms = mol.GetNumAtoms()
    res = nAtoms / sum(1. / deltas)
    return res


def CalculateGeometricTopovIndex(mol: Chem.Mol) -> float:
    """Get Geometric topological index.

    Or Getov.
    From Narumi H., MATCH (Comm. Math. Comp. Chem.), (1987), 22,195-207.

    Kier and Hall's valence delta-values are used in place of atom degrees.
    From Kier L. and Hall L., J. Pharm. Sci. (1983), 72(10),1170-1173.
    """
    nAtoms = mol.GetNumAtoms()
    deltas = _HKDeltas(mol, skipHs=0)
    while 0 in deltas:
        deltas.remove(0)
    deltas = numpy.array(deltas, 'd')
    temp = numpy.prod(deltas)
    res = numpy.power(temp, 1. / nAtoms)
    return res


def CalculateGravitationalTopoIndex(mol: Chem.Mol) -> float:
    """Get Gravitational topological index based on topological distance.

    Or Gravto
    From Katritzky, A. J. Phys. Chem., (1996), 100,10400-10407.
    """
    nAT = mol.GetNumAtoms()
    Distance = Chem.GetDistanceMatrix(mol)
    res = 0.0
    Atom = mol.GetAtoms()
    for i in range(nAT - 1):
        for j in range(i + 1, nAT):
            temp = Atom[i].GetMass() * Atom[j].GetMass()
            res = res + temp / numpy.power(Distance[i][j], 2)
    return res / 100


def CalculateGutmanVTopo(mol: Chem.Mol) -> float:
    """Get molecular topological index based on valence vertex degree.

    Or GMTIV.
    From Gutman,I. J. Chem. Inf. Comput. Sci., (1994), 34,1037-1039.
    """
    nAT = mol.GetNumAtoms()
    deltas = _HKDeltas(mol)
    Distance = Chem.GetDistanceMatrix(mol)
    res = 0.0
    for i in range(nAT):
        for j in range(i + 1, nAT):
            res = res + deltas[i] * deltas[j] * Distance[i, j]
    return numpy.log10(res)


_Topology = {'W': CalculateWeiner,
             'AW': CalculateMeanWeiner,
             'J': CalculateBalaban,
             'Tigdi': CalculateGraphDistance,
             'Xu': CalculateXuIndex,
             'GMTI': CalculateGutmanTopo,
             'Pol': CalculatePolarityNumber,
             'DZ': CalculatePoglianiIndex,
             'Ipc': CalculateIpc,
             'BertzCT': CalculateBertzCT,
             'Thara': CalculateHarary,
             'Tsch': CalculateSchiultz,
             'ZM1': CalculateZagreb1,
             'ZM2': CalculateZagreb2,
             'MZM1': CalculateMZagreb1,
             'MZM2': CalculateMZagreb2,
             'Qindex': CalculateQuadratic,
             'Platt': CalculatePlatt,
             'diametert': CalculateDiameter,
             'radiust': CalculateRadius,
             'petitjeant': CalculatePetitjean,
             'Sito': CalculateSimpleTopoIndex,
             'Hato': CalculateHarmonicTopoIndex,
             'Geto': CalculateGeometricTopoIndex,
             'Arto': CalculateArithmeticTopoIndex,
             'ISIZ': CalculateMolSizeTotalInf,
             'TIAC': CalculateAtomCompTotalInf,
             'IDET': CalculateDistanceEqualityTotalInf,
             'IDE': CalculateDistanceEqualityMeanInf,
             'IVDE': CalculateVertexEqualityTotalInf,
             'Sitov': CalculateSimpleTopovIndex,
             'Hatov': CalculateHarmonicTopovIndex,
             'Gravto': CalculateGravitationalTopoIndex,
             'Getov': CalculateGeometricTopovIndex,
             'GMTIV': CalculateGutmanVTopo,
             }


def GetTopology(mol: Chem.Mol) -> float:
    """Get all (35) constitutional descriptors."""
    result = {}
    for DesLabel in _Topology.keys():
        result[DesLabel] = round(_Topology[DesLabel](mol), 3)
    return result

#####################################################################
# if __name__ =='__main__':
#     smis = ['CCCC','CCCCC','CCCCCC','CC(N)C(=O)O','CC(N)C(=O)[O-]']
#     for index, smi in enumerate(smis):
#         m = Chem.MolFromSmiles(smi)
#         print(index+1)
#         print(smi)
#         print('\t',GetTopology(m))
#         print('\t',len(GetTopology(m)))
