# -*- coding: utf-8 -*-


"""Molecular topological indices."""

from typing import Iterable, List

import numpy as np
from rdkit import Chem
from rdkit.Chem import GraphDescriptors as GD
from rdkit.Chem import rdchem

periodicTable = rdchem.GetPeriodicTable()


class Topology:
    """Topological indices."""

    @staticmethod
    def _get_principle_quantum_number(atNum: int) -> int:
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

    @staticmethod
    def calculate_weiner(mol: Chem.Mol) -> float:
        """Get Weiner values of a molecule.
    
        Or W.
        """
        return 1.0 / 2 * sum(sum(Chem.GetDistanceMatrix(mol)))

    @staticmethod
    def calculate_mean_weiner(mol: Chem.Mol) -> float:
        """Get Mean Weiner index of a molecule.
    
        Or AW.
        """
        N = mol.GetNumAtoms()
        WeinerNumber = Topology.calculate_weiner(mol)
        return 2.0 * WeinerNumber / (N * (N - 1))

    @staticmethod
    def calculate_balaban(mol: Chem.Mol) -> float:
        """Get Balaban index of a molecule.
    
        Or J.
        """
        adjMat = Chem.GetAdjacencyMatrix(mol)
        Distance = Chem.GetDistanceMatrix(mol)
        Nbond = mol.GetNumBonds()
        Natom = mol.GetNumAtoms()
        S = np.sum(Distance, axis=1)
        mu = Nbond - Natom + 1
        sumk = 0.
        for i in range(len(Distance)):
            si = S[i]
            for j in range(i, len(Distance)):
                if adjMat[i, j] == 1:
                    sumk += 1. / np.sqrt(si * S[j])
        if mu + 1 != 0:
            J = float(Nbond) / float(mu + 1) * sumk
        else:
            J = 0
        return J

    @staticmethod
    def calculate_graph_distance(mol: Chem.Mol) -> float:
        """Get graph distance index.
    
        Or Tigdi.
        """
        Distance = Chem.GetDistanceMatrix(mol)
        n = int(Distance.max())
        res = 0.0
        for i in range(n):
            temp = 1. / 2 * sum(sum(Distance == i + 1))
            res = res + temp ** 2
        return np.log10(res)

    @staticmethod
    def calculate_diameter(mol: Chem.Mol) -> float:
        """Get largest value of the distance matrix.
    
        Or diametert.
        From Petitjean, M. J. Chem. Inf. Comput. Sci. 1992, 32, 4, 331-337.
        """
        Distance = Chem.GetDistanceMatrix(mol)
        return Distance.max()

    @staticmethod
    def calculate_radius(mol: Chem.Mol) -> float:
        """Get radius based on topology.
    
        Or radiust.
        From Petitjean, M. J. Chem. Inf. Comput. Sci. 1992, 32, 4, 331-337.
        """
        Distance = Chem.GetDistanceMatrix(mol)
        temp = []
        for i in Distance:
            temp.append(max(i))
        return min(temp)

    @staticmethod
    def calculate_petitjean(mol: Chem.Mol) -> float:
        """Get Petitjean based on topology.
    
        Or petitjeant.
        """
        diameter = Topology.calculate_diameter(mol)
        radius = Topology.calculate_radius(mol)
        return (float(diameter) - radius) / radius

    @staticmethod
    def calculate_xu_index(mol: Chem.Mol) -> float:
        """Get Xu index.
    
        Or Xu.
        """
        nAT = mol.GetNumAtoms()
        deltas = [x.GetDegree() for x in mol.GetAtoms()]
        Distance = Chem.GetDistanceMatrix(mol)
        sigma = np.sum(Distance, axis=1)
        temp1 = 0.0
        temp2 = 0.0
        for i in range(nAT):
            temp1 = temp1 + deltas[i] * ((sigma[i]) ** 2)
            temp2 = temp2 + deltas[i] * (sigma[i])
        Xu = np.sqrt(nAT) * np.log(temp1 / temp2)
        return Xu

    @staticmethod
    def calculate_gutman_topo(mol: Chem.Mol) -> float:
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
        return np.log10(res)

    @staticmethod
    def calculate_polarity_number(mol: Chem.Mol) -> float:
        """Get Polarity number.
    
        Or Pol.
        """
        Distance = Chem.GetDistanceMatrix(mol)
        res = 1. / 2 * sum(sum(Distance == 3))
        return res

    @staticmethod
    def calculate_pogliani_index(mol: Chem.Mol) -> float:
        """Get Poglicani index.
    
        Or DZ.
        From Pogliani L. J.Phys.Chem. (1996), 100,18065-18077.
        """
        res = 0.0
        for atom in mol.GetAtoms():
            n = atom.GetAtomicNum()
            nV = periodicTable.GetNOuterElecs(n)
            mP = Topology._get_principle_quantum_number(n)
            res = res + (nV + 0.0) / mP
        return res

    @staticmethod
    def calculate_ipc(mol: Chem.Mol) -> float:
        """Get Bonchev-Trinajstic complexity index.
    
        Or Ipc.
        From Bonchev D. & Trinajstic N., J. Chem. Phys. (1977) 67,4517-4533.
        """
        return np.log10(GD.Ipc(mol))

    @staticmethod
    def calculate_bertzct(mol: Chem.Mol) -> float:
        """Get Bertz complexity index.
    
        Or BertzCT.
        From Bertz S. H., J. Am. Chem. Soc. (1981) 103,3599-3601.
        """
        return np.log10(GD.BertzCT(mol))

    @staticmethod
    def calculate_harary(mol: Chem.Mol) -> float:
        """Get Harary number.
    
        Or Thara.
        """
        Distance = np.array(Chem.GetDistanceMatrix(mol), 'd')
        return 1.0 / 2 * (sum(1.0 / Distance[Distance != 0]))

    @staticmethod
    def calculate_schiultz(mol: Chem.Mol) -> float:
        """Get Schiultz number.
    
        Or Tsch.
        """
        Distance = np.array(Chem.GetDistanceMatrix(mol), 'd')
        Adjacent = np.array(Chem.GetAdjacencyMatrix(mol), 'd')
        VertexDegree = sum(Adjacent)
        return sum(np.dot((Distance + Adjacent), VertexDegree))

    @staticmethod
    def calculate_zagreb1(mol: Chem.Mol) -> float:
        """Get Zagreb index with order 1.
    
        Or ZM1.
        """
        deltas = [x.GetDegree() for x in mol.GetAtoms()]
        return sum(np.array(deltas) ** 2)

    @staticmethod
    def calculate_zagreb2(mol: Chem.Mol) -> float:
        """Get Zagreb index with order 2.
    
        Or ZM2.
        """
        ke = [x.GetBeginAtom().GetDegree() * x.GetEndAtom().GetDegree() for x in mol.GetBonds()]
        return sum(ke)

    @staticmethod
    def calculate_mzagreb1(mol: Chem.Mol) -> float:
        """Get Modified Zagreb index with order 1.
    
        Or MZM1.
        """
        deltas = [x.GetDegree() for x in mol.GetAtoms()]
        while 0 in deltas:
            deltas.remove(0)
        deltas = np.array(deltas, 'd')
        res = sum((1. / deltas) ** 2)
        return res

    @staticmethod
    def calculate_mzagreb2(mol: Chem.Mol) -> float:
        """Get Modified Zagreb index with order 2.
    
        Or MZM2.
        """
        cc = [x.GetBeginAtom().GetDegree() * x.GetEndAtom().GetDegree() for x in mol.GetBonds()]
        if len(cc) == 0:
            return 0.0
        while 0 in cc:
            cc.remove(0)
        cc = np.array(cc, 'd')
        res = sum((1. / cc) ** 2)
        return res

    @staticmethod
    def calculate_quadratic(mol: Chem.Mol) -> float:
        """Get Quadratic index.
    
        Or Qindex.
        """
        M = Topology.calculate_zagreb1(mol)
        N = mol.GetNumAtoms()
        return 3 - 2 * N + M / 2.0

    @staticmethod
    def calculate_platt(mol: Chem.Mol) -> float:
        """Get Platt number.
    
        Or Platt.
        """
        cc = [x.GetBeginAtom().GetDegree() + x.GetEndAtom().GetDegree() - 2 for x in mol.GetBonds()]
        return sum(cc)

    @staticmethod
    def calculate_simple_topo_index(mol: Chem.Mol) -> float:
        """Get the logarithm of the simple topological index.
    
        Or Sito.
        From Narumi H., MATCH (Comm. Math. Comp. Chem.), (1987), 22,195-207.
        """
        deltas = [x.GetDegree() for x in mol.GetAtoms()]
        while 0 in deltas:
            deltas.remove(0)
        deltas = np.array(deltas, 'd')
        res = np.prod(deltas)
        return np.log(res)

    @staticmethod
    def calculate_harmonic_topo_index(mol: Chem.Mol) -> float:
        """Get harmonic topological index.
    
        Or Hato.
        From Narumi H., MATCH (Comm. Math. Comp. Chem.), (1987), 22,195-207.
        """
        deltas = [x.GetDegree() for x in mol.GetAtoms()]
        while 0 in deltas:
            deltas.remove(0)
        if len(deltas) == 0:
            return 0.0
        deltas = np.array(deltas, 'd')
        nAtoms = mol.GetNumAtoms()
        res = nAtoms / sum(1. / deltas)
        return res

    @staticmethod
    def calculate_geometric_topo_index(mol: Chem.Mol) -> float:
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
        deltas = np.array(deltas, 'd')
        temp = np.prod(deltas)
        res = np.power(temp, 1. / nAtoms)
        return res

    @staticmethod
    def calculate_arithmetic_topo_index(mol: Chem.Mol) -> float:
        """Get Arithmetic topological index.
    
        Or Arto.
        From Narumi H., MATCH (Comm. Math. Comp. Chem.), (1987), 22,195-207.
        """
        nAtoms = mol.GetNumAtoms()
        nBonds = mol.GetNumBonds()
        res = 2. * nBonds / nAtoms
        return res

    @staticmethod
    def calculate_mol_size_total_inf(mol: Chem.Mol) -> float:
        """Get total information index on molecular size.
    
        Or ISIZ.
        """
        Hmol = Chem.AddHs(mol)
        nAT = Hmol.GetNumAtoms()
        ISIZ = nAT * np.log2(nAT)
        return ISIZ

    @staticmethod
    def calculate_atom_comp_total_inf(mol: Chem.Mol) -> float:
        """Ge total information index on atomic composition.
    
        Or TIAC.
        """
        Hmol = Chem.AddHs(mol)
        nAtoms = Hmol.GetNumAtoms()
        IC = []
        for i in range(nAtoms):
            at = Hmol.GetAtomWithIdx(i)
            IC.append(at.GetAtomicNum())
        Unique = np.unique(IC)
        NAtomType = len(Unique)
        res = 0.0
        for i in range(NAtomType):
            cc = IC.count(Unique[i])
            res += cc * np.log2(cc)
        if nAtoms != 0:
            return nAtoms * np.log2(nAtoms) - res
        else:
            return 0.0

    @staticmethod
    def calculate_distance_equality_total_inf(mol: Chem.Mol) -> float:
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
            res += cc * np.log2(cc)
        return n * np.log2(n) - res if n > 0 else np.NaN

    @staticmethod
    def _calculate_entropy(Probability: Iterable[float]) -> float:
        """calculate_ entropy (Information content) of given probability."""
        res = 0.0
        for i in Probability:
            if i != 0:
                res = res - i * np.log2(i)
        return res

    @staticmethod
    def calculate_distance_equality_mean_inf(mol: Chem.Mol) -> float:
        """Get the mean information index on distance equality.
    
        Or IDE.
        """
        Distance = Chem.GetDistanceMatrix(mol)
        nAT = mol.GetNumAtoms()
        n = 1. / 2 * nAT ** 2 - nAT
        DisType = int(Distance.max())
        res = 0.0
        cc = np.zeros(DisType, dtype=float)
        for i in range(DisType):
            cc[i] = 1. / 2 * sum(sum(Distance == i + 1))
        res = Topology._calculate_entropy(cc / n)
        return res

    @staticmethod
    def calculate_vertex_equality_total_inf(mol: Chem.Mol) -> float:
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
                res += cc * np.log2(cc)
        n = len(deltas)
        return n * np.log2(n) - res

    @staticmethod
    def _hall_kier_deltas(mol: Chem.Mol, skipHs: bool = True) -> List[float]:
        """calculate_ Kier & Hall valence delta-values for molecular connectivity.
    
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

    @staticmethod
    def calculate_simple_topo_vindex(mol: Chem.Mol) -> float:
        """Get the logarithm of the simple topological index.
    
        Or Sitov.
        From Narumi H., MATCH (Comm. Math. Comp. Chem.), (1987), 22,195-207.
    
        Kier and Hall's valence delta-values are used in place of atom degrees.
        From Kier L. and Hall L., J. Pharm. Sci. (1983), 72(10),1170-1173.
        """
        deltas = Topology._hall_kier_deltas(mol, skipHs=0)
        while 0 in deltas:
            deltas.remove(0)
        deltas = np.array(deltas, 'd')
        res = np.prod(deltas)
        return np.log(res)

    @staticmethod
    def calculate_harmonic_topo_vindex(mol: Chem.Mol) -> float:
        """Get harmonic topological index.
    
        Or Hatov.
        From Narumi H., MATCH (Comm. Math. Comp. Chem.), (1987), 22,195-207.
    
        Kier and Hall's valence delta-values are used in place of atom degrees.
        From Kier L. and Hall L., J. Pharm. Sci. (1983), 72(10),1170-1173.
        """
        deltas = Topology._hall_kier_deltas(mol, skipHs=0)
        while 0 in deltas:
            deltas.remove(0)
        deltas = np.array(deltas, 'd')
        nAtoms = mol.GetNumAtoms()
        res = nAtoms / sum(1. / deltas)
        return res

    @staticmethod
    def calculate_geometric_topo_vindex(mol: Chem.Mol) -> float:
        """Get Geometric topological index.
    
        Or Getov.
        From Narumi H., MATCH (Comm. Math. Comp. Chem.), (1987), 22,195-207.
    
        Kier and Hall's valence delta-values are used in place of atom degrees.
        From Kier L. and Hall L., J. Pharm. Sci. (1983), 72(10),1170-1173.
        """
        nAtoms = mol.GetNumAtoms()
        deltas = Topology._hall_kier_deltas(mol, skipHs=0)
        while 0 in deltas:
            deltas.remove(0)
        deltas = np.array(deltas, 'd')
        temp = np.prod(deltas)
        res = np.power(temp, 1. / nAtoms)
        return res

    @staticmethod
    def calculate_gravitational_topo_index(mol: Chem.Mol) -> float:
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
                res = res + temp / np.power(Distance[i][j], 2)
        return res / 100

    @staticmethod
    def calculate_gutman_vtopo(mol: Chem.Mol) -> float:
        """Get molecular topological index based on valence vertex degree.
    
        Or GMTIV.
        From Gutman,I. J. Chem. Inf. Comput. Sci., (1994), 34,1037-1039.
        """
        mol_ = Chem.RemoveHs(mol)
        nAT = mol_.GetNumAtoms()
        deltas = Topology._hall_kier_deltas(mol_)
        Distance = Chem.GetDistanceMatrix(mol_)
        res = 0.0
        for i in range(nAT):
            for j in range(i + 1, nAT):
                res += deltas[i] * deltas[j] * Distance[i, j]
        return np.log10(res)
    
    def get_all(mol: Chem.Mol) -> dict:
        """Get all (35) constitutional descriptors."""
        result = {}
        for DesLabel in _topology.keys():
            result[DesLabel] = _topology[DesLabel](mol)
        return result

_topology = {'W': Topology.calculate_weiner,
             'AW': Topology.calculate_mean_weiner,
             'J': Topology.calculate_balaban,
             'Thara': Topology.calculate_harary,
             'Tsch': Topology.calculate_schiultz,
             'Tigdi': Topology.calculate_graph_distance,
             'Platt': Topology.calculate_platt,
             'Xu': Topology.calculate_xu_index,
             'Pol': Topology.calculate_polarity_number,
             'DZ': Topology.calculate_pogliani_index,
             'Ipc': Topology.calculate_ipc,
             'BertzCT': Topology.calculate_bertzct,
             'GMTI': Topology.calculate_gutman_topo,
             'ZM1': Topology.calculate_zagreb1,
             'ZM2': Topology.calculate_zagreb2,
             'MZM1': Topology.calculate_mzagreb1,
             'MZM2': Topology.calculate_mzagreb2,
             'Qindex': Topology.calculate_quadratic,
             'diametert': Topology.calculate_diameter,
             'radiust': Topology.calculate_radius,
             'petitjeant': Topology.calculate_petitjean,
             'Sito': Topology.calculate_simple_topo_index,
             'Hato': Topology.calculate_harmonic_topo_index,
             'Geto': Topology.calculate_geometric_topo_index,
             'Arto': Topology.calculate_arithmetic_topo_index,
             'ISIZ': Topology.calculate_mol_size_total_inf,
             'TIAC': Topology.calculate_atom_comp_total_inf,
             'IDET': Topology.calculate_distance_equality_total_inf,
             'IDE': Topology.calculate_distance_equality_mean_inf,
             'IVDE': Topology.calculate_vertex_equality_total_inf,
             'Sitov': Topology.calculate_simple_topo_vindex,
             'Hatov': Topology.calculate_harmonic_topo_vindex,
             'Getov': Topology.calculate_geometric_topo_vindex,
             'Gravto': Topology.calculate_gravitational_topo_index,
             'GMTIV': Topology.calculate_gutman_vtopo,
             }
