# -*- coding: utf-8 -*-


"""Molecule format parsers."""

import re
import urllib.request

import defusedxml.ElementTree as ET
from openbabel import pybel
from rdkit import Chem


def ReadMolFromFile(self, filename: str = "") -> Chem.Mol:
    """Read molecular SDF or MOL file format."""
    mol = Chem.MolFromMolFile(filename)
    return mol


def ReadMolFromSmile(smi: str = "") -> Chem.Mol:
    """Read molecular SMILES string."""
    mol = Chem.MolFromSmiles(smi.strip())
    return mol


def ReadMolFromInchi(inchi: str = "") -> Chem.Mol:
    """Read molecular InChI string."""
    temp = pybel.readstring("inchi", inchi)
    smi = temp.write("smi")
    mol = Chem.MolFromSmiles(smi.strip())
    return mol


def ReadMolFromMol(filename: str = "") -> Chem.Mol:
    """Read MOL format file."""
    mol = Chem.MolFromMolFile(filename)
    return mol


def ReadMol(molstructure: str, molformat: str = 'smi') -> Chem.Mol:
    """Read molecular text of the specified format.

    :param molstructure: molecular text
    :param molformat: 3-letters code for openbabel supported format
    """
    mol = pybel.readstring(molformat, molstructure)
    return mol


def GetMolFromCAS(casid: str = "") -> Chem.Mol:
    """Download molecule from ChemNet by CAS ID."""
    casid = casid.strip()
    request = urllib.request.Request(f'http://www.chemnet.com/cas/supplier.cgi?terms={casid}&l=&exact=dict')
    with urllib.request.urlopen(request) as response:
        temp = [line.decode("utf-8") for line in response.readlines()]
    for i in temp:
        if re.findall('InChI=', i) == ['InChI=']:
            k = i.split('    <td align="left">')
            kk = k[1].split('</td>\r\n')
            if kk[0][0:5] == "InChI":
                res = kk[0]
            else:
                res = "None"
    m = Chem.MolFromInchi(res.strip())
    return m


def GetMolFromEBI(chid: str = "") -> Chem.Mol:
    """Donwload molecule from ChEBI or ChEMBL using ChEBI or ChEMBL id."""
    chid = chid.strip().upper()
    if chid.startswith('CHEBI'):
        request = urllib.request.Request(f'https://www.ebi.ac.uk/chebi/saveStructure.do?sdf=true&chebiId={chid}')
        with urllib.request.urlopen(request) as response:
            sdf = response.read()
        if not len(sdf):
            raise Exception(f'Not a valid ChEBI ID: {chid}')
    elif chid.startswith('CHEMBL'):
        request = urllib.request.Request(f'https://www.ebi.ac.uk/chembl/api/data/molecule?chembl_id={chid}')
        with urllib.request.urlopen(request) as response:
            xml = response.read()
        xml_tree = ET.fromstring(xml)
        structure = xml_tree.findall('./molecules/molecule/molecule_structures/molfile')
        if not len(structure):
            raise Exception(f'Not a valid ChEMBL ID: {chid}')
        sdf = structure[0].text
    else:
        raise Exception('Valid ID starts with CHEBI: or CHEMBL')
    m = Chem.MolFromMolBlock(sdf)
    return m


def GetMolFromNCBI(cid: str = "") -> Chem.Mol:
    """Download molecule from PubChem using PubChem CID."""
    cid = cid.strip()
    cid = cid.upper().replace('CID', '') if 'CID' in cid.upper() else cid
    request = urllib.request.Request(f'http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid={cid}&disopt=SaveSDF')
    with urllib.request.urlopen(request) as response:
        m = Chem.MolFromMolBlock(response.read())
    return m


def GetMolFromDrugbank(dbid: str = "") -> Chem.Mol:
    """Download molecule from DrugBank using DrugBank id."""
    dbid = dbid.strip()
    request = urllib.request.Request(f'http://www.drugbank.ca/drugs/{dbid}.sdf')
    with urllib.request.urlopen(request) as response:
        m = Chem.MolFromMolBlock(response.read())
    return m


def GetMolFromKegg(kid: str = ""):
    """Download molecule from KEGG using KEGG id."""
    ID = str(kid)
    request = urllib.request.Request(f'http://www.genome.jp/dbget-bin/www_bget?-f+m+drug{ID}')
    with urllib.request.urlopen(request) as response:  # nosec: S310
        m = Chem.MolFromMolBlock(response.read())
    return m


############
# if __name__=="__main__":
#     print("Downloading......")
#     temp=GetMolFromCAS(casid="50-12-4")
#     print(temp)
#     temp=GetMolFromNCBI(cid="2244")
#     print(temp)
#     temp=GetMolFromDrugbank(dbid="DB00133")
#     print(temp)
#     temp=GetMolFromKegg(kid="D02176")
#     print(temp)
