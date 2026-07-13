# -*- coding: utf-8 -*-

"""Tests for chemopy.radius (van der Waals radii lookup table)."""

import unittest

from chemopy import radius


class TestRadii(unittest.TestCase):
    def test_known_elements(self):
        self.assertEqual(radius.radii["H"], 1.20)
        self.assertEqual(radius.radii["C"], 1.70)
        self.assertEqual(radius.radii["O"], 1.52)

    def test_fallback_entry_present(self):
        self.assertIn(".", radius.radii)


if __name__ == "__main__":
    unittest.main()
