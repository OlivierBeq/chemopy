# -*- coding: utf-8 -*-

"""Unit tests for chemopy."""

from unittest import TestCase

# from openbabel import pybel
from rdkit import Chem

from chemopy import PyChem3D
from tests.constants import MOL_FILE


class TestPyChem3D(TestCase):
    """Tests for PyChem3D."""

    def setUp(self):
        """Create PyChem3D molecule."""
        self.mol1 = PyChem3D()
        self.mol2 = PyChem3D()

    def test_ReadMol(self):
        """Test reading the sample MOL file."""
        self.mol1.ReadMol(MOL_FILE)
        self.assertEqual(self.mol1.mol.GetNumHeavyAtoms(), self.num_heavy_atoms)
        self.assertEqual(Chem.rdMolDescriptors.CalcMolFormula(self.mol1.mol), self.formula)
        self.molecules.append(self.mol1)

    def test_GetMolFromNCBI(self):
        """Test GetConnectivity."""
        pass

    def test_GetMolFromEBI(self):
        """Test GetConnectivity."""
        pass

    def test_GetMolFromCAS(self):
        """Test GetConnectivity."""
        pass

    def test_GetMolFromKegg(self):
        """Test GetConnectivity."""
        pass

    def test_GetMolFromDrugbank(self):
        """Test GetConnectivity."""
        pass

    def test_SetMOPACParameters(self):
        """Test GetConnectivity."""
        pass

    def test_GetGeometric(self):
        """Test GetConnectivity."""
        pass

    def test_GetMoRSE(self):
        """Test GetConnectivity."""
        pass

    def test_GetRDF(self):
        """Test GetConnectivity."""
        pass

    def test_GetWHIM(self):
        """Test GetConnectivity."""
        pass

    def test_GetCPSA(self):
        """Test GetConnectivity."""
        pass

    def test_GetQuantumChemistry(self):
        """Test GetConnectivity."""
        pass

    def test_GetAllDescriptors(self):
        """Test GetConnectivity."""
        pass
