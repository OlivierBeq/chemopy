# -*- coding: utf-8 -*-

"""Tests for EState SMARTS atom-type definitions and matching."""

import unittest
from unittest import mock

from rdkit import Chem

from chemopy import atom_types


class TestBuildPatterns(unittest.TestCase):
    """Tests build_patterns() SMARTS compilation."""

    def setUp(self):
        self.addCleanup(atom_types.build_patterns)  # restore the default patterns afterwards

    def test_build_patterns_default(self):
        atom_types.build_patterns()
        self.assertIsNotNone(atom_types.esPatterns)
        self.assertEqual(len(atom_types.esPatterns), len(atom_types._raw_d))
        self.assertTrue(all(entry is not None for entry in atom_types.esPatterns))

    def test_invalid_pattern_warns_with_user_warning(self):
        with mock.patch.object(Chem, "MolFromSmarts", side_effect=ValueError("bad SMARTS")):
            with self.assertWarns(UserWarning):
                atom_types.build_patterns([("bogus", "not a real smarts")])


class TestTypeAtoms(unittest.TestCase):
    """Tests EState atom type assignment."""

    def test_type_atoms_ethanol(self):
        mol = Chem.MolFromSmiles("CCO")
        res = atom_types.type_atoms(mol)
        self.assertEqual(len(res), mol.GetNumAtoms())
        self.assertTrue(any(len(types) > 0 for types in res))

    def test_get_atom_label(self):
        mol = Chem.MolFromSmiles("CCO")
        res = atom_types.get_atom_label(mol)
        self.assertEqual(len(res), len(atom_types._raw_d))


if __name__ == "__main__":
    unittest.main()
