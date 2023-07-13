# -*- coding: utf-8 -*-


"""Molecule format parsers."""

import re
import urllib.request

from rdkit import Chem


class MolFrom:
    """Obtain molecular structures from common resources."""

    @staticmethod
    def CAS(cas_id: str = "") -> Chem.Mol:
        """Download molecule from ChemNet by CAS ID."""
        casid = cas_id.strip()
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

    @staticmethod
    def EBI(chebi_id: str = "") -> Chem.Mol:
        """Donwload molecule from ChEBI or ChEMBL using ChEBI or ChEMBL id."""
        chid = chebi_id.strip().upper()
        if chid.startswith('CHEBI'):
            request = urllib.request.Request(f'https://www.ebi.ac.uk/chebi/saveStructure.do?sdf=true&chebiId={chid}')
            with urllib.request.urlopen(request) as response:
                sdf = response.read()
            if not len(sdf):
                raise Exception(f'Not a valid ChEBI ID: {chid}')
        elif chid.startswith('CHEMBL'):
            request = urllib.request.Request(f'https://www.ebi.ac.uk/chembl/api/data/molecule?chembl_id={chid}&format=sdf')
            with urllib.request.urlopen(request) as response:
                sdf = response.read().decode()
            # xml_tree = ET.fromstring(xml)
            # structure = xml_tree.findall('./molecules/molecule/molecule_structures/molfile')
            if not len(sdf):
                raise Exception(f'Not a valid ChEMBL ID: {chid}')
            # sdf = structure[0].text
        else:
            raise Exception('Valid ID starts with CHEBI: or CHEMBL')
        m = Chem.MolFromMolBlock(sdf)
        return m

    @staticmethod
    def NCBI(cid: str = "") -> Chem.Mol:
        """Download molecule from PubChem using PubChem CID."""
        cid = cid.strip()
        cid = cid.upper().replace('CID', '') if 'CID' in cid.upper() else cid
        request = urllib.request.Request(f'http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid={cid}&disopt=SaveSDF')
        with urllib.request.urlopen(request) as response:
            m = Chem.MolFromMolBlock(response.read())
        return m

    @staticmethod
    def DrugBank(drugbank_id: str = "") -> Chem.Mol:
        """Download molecule from DrugBank using DrugBank id."""
        dbid = drugbank_id.strip()
        request = urllib.request.Request(f'https://go.drugbank.com/structures/small_molecule_drugs/{dbid}.sdf',
                                         headers={"User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_10_1) "
                                                                "AppleWebKit/537.36 (KHTML, like Gecko) "
                                                                "Chrome/39.0.2171.95 "
                                                                "Safari/537.36"})
        with urllib.request.urlopen(request) as response:
            m = Chem.MolFromMolBlock(response.read())
        return m

    @staticmethod
    def KEGG(kegg_id: str = ""):
        """Download molecule from KEGG using KEGG id."""
        ID = str(kegg_id)
        request = urllib.request.Request(f'http://www.genome.jp/dbget-bin/www_bget?-f+m+drug+{ID}')
        with urllib.request.urlopen(request) as response:  # nosec: S310
            m = Chem.MolFromMolBlock(response.read())
        return m
