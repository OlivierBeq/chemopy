# -*- coding: utf-8 -*-

"""Tests for details about molecular descriptors/fingerprints."""

import unittest

from chemopy import ChemoPy

from .constants import SMALL_MOLECULES_3D


class TestDetailsOfDescriptorsAndFingerprints(unittest.TestCase):
    """Tests molecular descriptors."""

    @classmethod
    def setUpClass(cls) -> None:
        """Create the molecular descriptor calculator (once for the whole class: MOPAC is slow)."""
        cls.cmp = ChemoPy(ignore_3D=False)
        cls.desc_names = cls.cmp.calculate([SMALL_MOLECULES_3D[0]], show_banner=False).columns.tolist()

    def test_each_descriptor_details(self):
        """Test the details for each descriptor available."""
        for desc_name in self.desc_names:
            details = self.cmp.get_details(desc_name)
            self.assertIsNotNone(details)
            self.assertFalse(details.empty)
            self.assertEqual(details.shape[0], 1)
            self.assertEqual(details.Name.item(), desc_name)

    def test_all_descriptor_details(self):
        """Test the details of all descriptors/fingerprints available."""
        details = self.cmp.get_details()
        self.assertIsNotNone(details)
        self.assertFalse(details.empty)
        # Includes the 12 additional fingerprints
        self.assertEqual(details.shape[0], 1196)
