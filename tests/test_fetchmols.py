# -*- coding: utf-8 -*-

"""Tests for molecule retrieval from online resources (network calls are mocked)."""

import unittest
from contextlib import contextmanager
from unittest import mock

from rdkit import Chem

from chemopy import MolFrom

_ETHANOL_MOLBLOCK = Chem.MolToMolBlock(Chem.MolFromSmiles("CCO"))

_CAS_RESPONSE_LINES = [
    b"<html>\r\n",
    b'    <td align="left">InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3</td>\r\n',
    b"</html>\r\n",
]


def _fake_urlopen(read_value: bytes = None, readlines_value=None):
    """Build a context-manager-compatible fake for urllib.request.urlopen."""

    @contextmanager
    def _cm(*args, **kwargs):
        response = mock.Mock()
        if read_value is not None:
            response.read.return_value = read_value
        if readlines_value is not None:
            response.readlines.return_value = readlines_value
        yield response

    return _cm


class TestMolFromCAS(unittest.TestCase):
    def test_cas_valid_response(self):
        with mock.patch("urllib.request.urlopen", side_effect=_fake_urlopen(readlines_value=_CAS_RESPONSE_LINES)):
            mol = MolFrom.CAS("64-17-5")
        self.assertIsNotNone(mol)
        self.assertNotEqual(Chem.MolToSmiles(mol), "")

    def test_cas_no_inchi_in_response_raises_clear_error(self):
        with mock.patch("urllib.request.urlopen", side_effect=_fake_urlopen(readlines_value=[b"<html></html>\r\n"])):
            with self.assertRaises(ValueError):
                MolFrom.CAS("64-17-5")


class TestMolFromEBI(unittest.TestCase):
    def test_ebi_chebi_id(self):
        with mock.patch("urllib.request.urlopen", side_effect=_fake_urlopen(read_value=_ETHANOL_MOLBLOCK.encode())):
            mol = MolFrom.EBI("CHEBI:16236")
        self.assertIsNotNone(mol)
        self.assertNotEqual(Chem.MolToSmiles(mol), "")

    def test_ebi_chembl_id(self):
        with mock.patch("urllib.request.urlopen", side_effect=_fake_urlopen(read_value=_ETHANOL_MOLBLOCK.encode())):
            mol = MolFrom.EBI("CHEMBL545")
        self.assertIsNotNone(mol)

    def test_ebi_invalid_id_raises(self):
        with self.assertRaises(Exception):
            MolFrom.EBI("NOTVALID")

    def test_ebi_empty_response_raises(self):
        with mock.patch("urllib.request.urlopen", side_effect=_fake_urlopen(read_value=b"")):
            with self.assertRaises(Exception):
                MolFrom.EBI("CHEBI:16236")


class TestMolFromNCBI(unittest.TestCase):
    def test_ncbi(self):
        with mock.patch("urllib.request.urlopen", side_effect=_fake_urlopen(read_value=_ETHANOL_MOLBLOCK.encode())):
            mol = MolFrom.NCBI("CID702")
        self.assertIsNotNone(mol)
        self.assertNotEqual(Chem.MolToSmiles(mol), "")


class TestMolFromDrugBank(unittest.TestCase):
    def test_drugbank(self):
        with mock.patch("urllib.request.urlopen", side_effect=_fake_urlopen(read_value=_ETHANOL_MOLBLOCK.encode())):
            mol = MolFrom.DrugBank("DB00898")
        self.assertIsNotNone(mol)
        self.assertNotEqual(Chem.MolToSmiles(mol), "")


class TestMolFromKEGG(unittest.TestCase):
    def test_kegg(self):
        with mock.patch("urllib.request.urlopen", side_effect=_fake_urlopen(read_value=_ETHANOL_MOLBLOCK.encode())):
            mol = MolFrom.KEGG("D00013")
        self.assertIsNotNone(mol)
        self.assertNotEqual(Chem.MolToSmiles(mol), "")


if __name__ == "__main__":
    unittest.main()
