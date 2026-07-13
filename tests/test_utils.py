# -*- coding: utf-8 -*-

"""Tests for chemopy.utils helper functions."""

import os
import tempfile
import unittest

from rdkit import Chem
from rdkit.Chem import AllChem

from chemopy import utils


class TestGetR(unittest.TestCase):
    def test_get_r(self):
        self.assertEqual(utils.get_r(3), [1.0, 1.5, 2.0])


class TestDistances(unittest.TestCase):
    def test_get_atom_distance(self):
        self.assertAlmostEqual(utils.get_atom_distance([0, 0, 0], [3, 4, 0]), 5.0)

    def test_get_geometrical_distance_matrix(self):
        coords = [[0, 0, 0], [3, 4, 0], [0, 0, 0]]
        mat = utils.get_geometrical_distance_matrix(coords)
        self.assertEqual(mat.shape, (3, 3))
        self.assertAlmostEqual(mat[0, 1], 5.0)
        self.assertAlmostEqual(mat[1, 0], 5.0)
        self.assertEqual(mat[0, 0], 0.0)


class TestPathHelpers(unittest.TestCase):
    def test_are_all_paths_absolute(self):
        self.assertTrue(utils.are_all_paths_absolute(["/a/b", "/c/d"]))
        self.assertFalse(utils.are_all_paths_absolute(["/a/b", "c/d"]))

    def test_are_all_paths_relative(self):
        self.assertTrue(utils.are_all_paths_relative(["a/b", "c/d"]))
        self.assertFalse(utils.are_all_paths_relative(["/a/b", "c/d"]))

    def test_is_in_subdirectory_tree_true(self):
        with tempfile.TemporaryDirectory() as parent:
            child = os.path.join(parent, "sub")
            os.makedirs(child)
            self.assertTrue(utils.is_in_subdirectory_tree(child, parent))

    def test_is_in_subdirectory_tree_false(self):
        with tempfile.TemporaryDirectory() as d1, tempfile.TemporaryDirectory() as d2:
            self.assertFalse(utils.is_in_subdirectory_tree(d1, d2))


class TestFileHelpers(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)
        for name in ("a.arc", "b.arc", "c.txt"):
            with open(os.path.join(self.tmpdir.name, name), "w") as fh:
                fh.write("x")

    def test_get_file_in_dir_from_ext(self):
        found = utils.get_file_in_dir_from_ext(self.tmpdir.name, ".arc")
        self.assertEqual(len(found), 2)
        self.assertTrue(all(f.endswith(".arc") for f in found))

    def test_get_lastest_created_file_from_dirpath(self):
        latest = utils.get_lastest_created_file(dirpath=self.tmpdir.name)
        self.assertTrue(os.path.isfile(latest))

    def test_get_lastest_created_file_from_absolute_filepaths(self):
        paths = [os.path.join(self.tmpdir.name, "a.arc"), os.path.join(self.tmpdir.name, "b.arc")]
        latest = utils.get_lastest_created_file(filepaths=paths)
        self.assertIn(latest, paths)

    def test_get_lastest_created_file_from_dirpath_and_relative_filepaths(self):
        latest = utils.get_lastest_created_file(dirpath=self.tmpdir.name, filepaths=["a.arc", "b.arc"])
        self.assertTrue(os.path.isfile(latest))

    def test_get_lastest_created_file_neither_arg_raises(self):
        with self.assertRaises(ValueError):
            utils.get_lastest_created_file()

    def test_get_lastest_created_file_missing_dir_raises(self):
        with self.assertRaises(NotADirectoryError):
            utils.get_lastest_created_file(dirpath="/no/such/dir")

    def test_get_lastest_created_file_strict_missing_file_raises(self):
        with self.assertRaises(FileNotFoundError):
            utils.get_lastest_created_file(dirpath=self.tmpdir.name, filepaths=["missing.arc"])

    def test_get_lastest_created_file_non_strict_missing_file_ignored(self):
        latest = utils.get_lastest_created_file(
            dirpath=self.tmpdir.name, filepaths=["missing.arc", "a.arc"], strict=False
        )
        self.assertTrue(os.path.isfile(latest))


class TestMolConversions(unittest.TestCase):
    def test_rdkit_to_pybel_and_back_without_charges(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        pybel_mol = utils.rdkit_to_pybel_mol(mol)
        self.assertIsNotNone(pybel_mol)
        back = utils.pybel_to_rdkit_mol(pybel_mol)
        self.assertIsNotNone(back)
        self.assertEqual(back.GetNumAtoms(), mol.GetNumAtoms())

    def test_rdkit_to_pybel_and_back_with_gasteiger_charges(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        AllChem.ComputeGasteigerCharges(mol)
        pybel_mol = utils.rdkit_to_pybel_mol(mol)
        back = utils.pybel_to_rdkit_mol(pybel_mol)
        self.assertEqual(back.GetNumAtoms(), mol.GetNumAtoms())
        self.assertTrue(back.GetAtomWithIdx(0).HasProp("_GasteigerCharge"))

    def test_needs_hydrogens_true(self):
        mol = Chem.MolFromSmiles("CCO")
        self.assertTrue(utils.needs_hydrogens(mol))

    def test_needs_hydrogens_false(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        self.assertFalse(utils.needs_hydrogens(mol))


if __name__ == "__main__":
    unittest.main()
