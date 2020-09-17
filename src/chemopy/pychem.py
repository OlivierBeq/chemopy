# -*- coding: utf-8 -*-


"""Main classes for the calculation of molecular descriptors."""

from typing import Dict, Tuple

from openbabel import pybel
from rdkit import Chem

from pychem import (basak, bcut, charge, connectivity, constitution, cpsa,
                    estate, fingerprint, geary, geometric, getmol, kappa, moe,
                    molproperty, moran, moreaubroto, morse, quanchem, rdf,
                    topology, whim)
from pychem.GeoOpt import Dispose, GetARCFile

FingerprintName = list(fingerprint._FingerprintFuncs.keys())
# ['topological','Estate','FP4','atompairs','torsions',
#                      'morgan','MACCS']

# TODO: 1) merge both 2D and 3D in a unique class keeping PyChem3D's obmol and rdmol duality
#       2) create a MOPAC opti method that generates an ARC file
#       3) force the use of with to dispose ARC file


class PyChem2D:
    """A PyDrug class used for computing drug descriptors."""

    def __init__(self):
        """Initialise a PyChem2D object.

        :noindex:
        """
        pass

    def ReadMolFromMolFile(self, filename: str = "") -> None:
        """Read a molecule from SDF or MOL input file."""
        self.mol = Chem.MolFromMolFile(filename)

    def ReadMolFromSmiles(self, smi: str = "") -> None:
        """Read a molecule from SMILES string."""
        self.mol = Chem.MolFromSmiles(smi.strip())

    def ReadMolFromInchi(self, inchi: str = "") -> None:
        """Read a molecule from InChI string."""
        self.mol = Chem.MolFromInchi(inchi.strip())

    def ReadMolFromMol(self, filename: str = "") -> None:
        """Read a molecule from a MOL input file."""
        self.mol = Chem.MolFromMolFile(filename)

    def GetMolFromNCBI(self, ID: str = "") -> None:
        """Get a molecule by NCBI id.

        :param ID: CID NCBI compound identifier (e.g., 2244).
        """
        self.mol = getmol.GetMolFromNCBI(cid=ID)

    def GetMolFromEBI(self, ID: str = "") -> None:
        """Get a molecule by EBI id.

        :param ID: ChEBI or ChEMBL compound identifier.
        """
        self.mol = getmol.GetMolFromEBI(ID)

    def GetMolFromCAS(self, ID: str = "") -> None:
        """Get a molecule by CAS id.

        :param ID: CAS compound identifier (e.g., 50-29-3).
        """
        self.mol = getmol.GetMolFromCAS(casid=ID)

    def GetMolFromKegg(self, ID: str = "") -> None:
        """Get a molecule by kegg id.

        :param ID: KEGG compound identifier (e.g., D02176).
        """
        self.mol = getmol.GetMolFromKegg(kid=ID)

    def GetMolFromDrugbank(self, ID: str = "") -> None:
        """Get a molecule by drugbank id.

        :param ID: Drugbank compound identifier (e.g. DB00133)
        """
        self.mol = getmol.GetMolFromDrugbank(dbid=ID)

    def GetKappa(self) -> dict:
        """Calculate all 7 kappa descriptors."""
        res = kappa.GetKappa(self.mol)
        return res

    def GetCharge(self) -> dict:
        """Calculate all 25 charge descriptors."""
        res = charge.GetCharge(self.mol)
        return res

    def GetConnectivity(self) -> dict:
        """Calculate all 44 connectivity descriptors."""
        res = connectivity.GetConnectivity(self.mol)
        return res

    def GetConstitution(self) -> dict:
        """Calculate all 30 constitutional descriptors."""
        res = constitution.GetConstitutional(self.mol)
        return res

    def GetEstateDescriptors(self) -> dict:
        """Calculate 316 estate descriptors."""
        res = estate.GetEstateDescriptors(self.mol)
        return res

    def GetGeary(self) -> dict:
        """Calculate all 32 Geary autocorrelation descriptors."""
        res = geary.GetGearyAuto(self.mol)
        return res

    def GetMOE(self) -> dict:
        """Calculate all 60 MOE-type descriptors."""
        res = moe.GetMOE(self.mol)
        return res

    def GetMolProperties(self) -> dict:
        """Calculate all 6 molecular properties."""
        res = molproperty.GetMolecularProperties(self.mol)
        return res

    def GetMoran(self) -> dict:
        """Calculate all 32 Moran autocorrealtion descriptors."""
        res = moran.GetMoranAuto(self.mol)
        return res

    def GetMoreauBroto(self) -> dict:
        """Calculate all 32 Moreau-Broto autocorrelation descriptors."""
        res = moreaubroto.GetMoreauBrotoAuto(self.mol)
        return res

    def GetTopology(self) -> dict:
        """Calculate all 35 topological descriptors."""
        res = topology.GetTopology(self.mol)
        return res

    def GetBcut(self) -> dict:
        """Calculate all 64 bcut/burden descriptors."""
        res = bcut.GetBurden(self.mol)
        return res

    def GetBasak(self) -> dict:
        """Calculate all 21 Basak information descriptors."""
        res = basak.Getbasak(self.mol)
        return res

    def GetAllFingerprints(self, radius: int = 2) -> Dict[str, Tuple[int, dict]]:
        """Calculate all molecular fingerprints.

        :param radius: maximum radius of Morgan fingeprints.
        """
        fps = {}
        for FPName in FingerprintName:
            temp = fingerprint._FingerprintFuncs[FPName]
            if FPName == 'morgan':
                size, bits, _ = temp(self.mol, radius)
            else:
                size, bits, _ = temp(self.mol)
            fps.update({FPName: (size, bits)})
        fps.update(estate.GetEstateFingerprints(self.mol))
        return fps

    def GetFingerprint(self, FPName: str = 'topological') -> Tuple[int, dict]:
        """Calculate molecular fingerprint.

        :param FPName: fingerprint name as in pychem.FingerprintName
        """
        if FPName in FingerprintName:
            size, bits, _ = fingerprint._FingerprintFuncs[FPName](self.mol)
            return size, bits
        raise NotImplementedError(f'Fingerprint type {FPName} is not available.')

    def GetAllDescriptors(self) -> dict:
        """Calculate all 633 2D descriptors."""
        res = {}
        res.update(self.GetKappa())
        res.update(self.GetCharge())
        res.update(self.GetConnectivity())
        res.update(self.GetConstitution())
        res.update(self.GetEstateDescriptors())
        res.update(self.GetGeary())
        res.update(self.GetMOE())
        res.update(self.GetMoran())
        res.update(self.GetMoreauBroto())
        res.update(self.GetTopology())
        res.update(self.GetMolProperties())
        res.update(self.GetBasak())
        res.update(self.GetBcut())
        return res


class PyChem3D:
    """PyDrug class used for computing drug descriptors."""

    def __init__(self):
        """Initialise a PyChem3D object.

        :noindex:
        """
        pass

    def ReadMol(self, molstr: str = "", molformat: str = 'smi') -> None:
        """Read a molecular input string.

        :param molstr: input molecular string
        :param molformat: 3-letters code for openbabel supported format
        """
        self.mol = pybel.readstring(molformat, molstr)
        self.rdmol = Chem.MolFromMolBlock(self.mol.write(format='sdf'))

    def GetMolFromNCBI(self, ID: str = "") -> None:
        """Get a molecule by NCBI id.

        :param ID: CID NCBI compound identifier (e.g., 2244).
        """
        self.rdmol = getmol.GetMolFromNCBI(cid=ID)
        self.mol = pybel.readstring('sdf', Chem.MolToMolBlock(self.rdmol))

    def GetMolFromEBI(self, ID: str = "") -> None:
        """Get a molecule by EBI id.

        :param ID: ChEBI or ChEMBL compound identifier.
        """
        self.rdmol = getmol.GetMolFromEBI(ID)
        self.mol = pybel.readstring('sdf', Chem.MolToMolBlock(self.rdmol))

    def GetMolFromCAS(self, ID="") -> None:
        """Get a molecule by CAS id.

        :param ID: CAS compound identifier (e.g., 50-29-3).
        """
        self.rdmol = getmol.GetMolFromCAS(casid=ID)
        self.mol = pybel.readstring('sdf', Chem.MolToMolBlock(self.rdmol))

    def GetMolFromKegg(self, ID: str = "") -> None:
        """Get a molecule by kegg id.

        :param ID: KEGG compound identifier (e.g., D02176).
        """
        self.rdmol = getmol.GetMolFromKegg(kid=ID)
        self.mol = pybel.readstring('sdf', Chem.MolToMolBlock(self.rdmol))

    def GetMolFromDrugbank(self, ID: str = "") -> None:
        """Get a molecule by drugbank id.

        :param ID: Drugbank compound identifier (e.g. DB00133)
        """
        self.rdmol = getmol.GetMolFromDrugbank(dbid=ID)
        self.mol = pybel.readstring('sdf', Chem.MolToMolBlock(self.rdmol))

    def GetGeometric(self) -> dict:
        """Calculate all 12 geometric descriptors."""
        GetARCFile(self.mol)
        res = geometric.GetGeometric(self.mol)
        return res

    def GetMoRSE(self) -> dict:
        """Calculate all 210 3D MoRSE descriptors."""
        GetARCFile(self.mol)
        res = morse.GetMoRSE(self.mol)
        return res
    
    def GetRDF(self) -> dict:
        """Calculate all 180 3D RDF descriptors."""
        GetARCFile(self.mol)
        res = rdf.GetRDF(self.mol)
        return res

    def GetWHIM(self) -> dict:
        """Calculate all 70 WHIM descriptors."""
        GetARCFile(self.mol)
        res = whim.GetWHIM()
        return res

    def GetCPSA(self) -> dict:
        """Calculate all 30 CPSA descriptors."""
        GetARCFile(self.mol)
        res = cpsa.GetCPSA()
        return res

    def GetAllDescriptor(self) -> dict:
        """Calculate all 502 3D descriptors."""
        res = {}
        GetARCFile(self.mol)
        res.update(cpsa.GetCPSA())
        res.update(rdf.GetRDF(self.mol))
        res.update(whim.GetWHIM())
        res.update(morse.GetMoRSE(self.mol))
        res.update(geometric.GetGeometric(self.mol))
        return res


# if __name__ == "__main__":
#     drugclass = PyChem2D()
#     drugclass.ReadMolFromSmile("CCC1(c2ccccc2)C(=O)N(C)C(=N1)O")
#     print(drugclass.GetKappa())
#     print(len(drugclass.GetTopology()))
#     print(len(drugclass.GetBasak()))
#     print(len(drugclass.GetKappa()))
#     print(len(drugclass.GetConnectivity()))
#     print(len(drugclass.GetConstitution()))
#     print(len(drugclass.GetMoran()))
#     print(len(drugclass.GetMOE()))
#     print(len(drugclass.GetGeary()))
#     print(len(drugclass.GetMolProperty()))
#     print(len(drugclass.GetBcut()))
#     print(len(drugclass.GetEstate()))
#     print(len(drugclass.GetMoreauBroto()))
#     print(len(drugclass.GetCharge()))
#     print(len(drugclass.GetAllDescriptor()))
#     print(drugclass.GetAllDescriptor())
# #    print(drugclass.GetMolFromDrugbank(ID="DB00133"))
#     res = drugclass.GetFingerprint(FPName='Estate')
#     print(res)

#     ###############################
#     drug = PyChem3D()
#     molsmi = drug.GetMolFromDrugbank("DB00133")
#     print(molsmi)
#     drug.ReadMol(molsmi, 'smi')
#     print(len(drug.GetGeometric()))
#     print(len(drug.GetCPSA()))
#     print(len(drug.GetMoRSE()))
#     print(len(drug.GetRDF()))
#     print(len(drug.GetWHIM()))
#     print(drug.GetAllDescriptor())
#     print(len(drug.GetAllDescriptor()))
