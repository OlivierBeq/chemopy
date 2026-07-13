# -*- coding: utf-8 -*-

"""Tests that ChemoPy.calculate engages exactly the number of cores requested."""

import multiprocessing
import unittest
import warnings
from unittest import mock

from chemopy import chemopy
from chemopy.chemopy import ChemoPy, _pin_worker

from .constants import SMALL_MOLECULES_3D


class TestNjobsResolution(unittest.TestCase):
    """Tests the njobs<0 (joblib-style "relative to all cores") resolution and clipping."""

    def test_negative_njobs_resolves_relative_to_cpu_count(self):
        with mock.patch("os.cpu_count", return_value=8):
            cmp = ChemoPy(ignore_3D=True)
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                values = cmp.calculate([SMALL_MOLECULES_3D[0]], show_banner=False, njobs=-1)
            self.assertEqual(values.shape[0], 1)
            self.assertFalse(any(issubclass(w.category, UserWarning) for w in caught))

    def test_njobs_all_but_one_core(self):
        with mock.patch("os.cpu_count", return_value=8):
            cmp = ChemoPy(ignore_3D=True)
            # njobs=-2 should resolve to 8 - 2 + 1 = 7 workers: just verify it runs without
            # oversubscription warning (7 <= 8).
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                cmp.calculate([SMALL_MOLECULES_3D[0]], show_banner=False, njobs=-2)
            self.assertFalse(any(issubclass(w.category, UserWarning) for w in caught))

    def test_njobs_exceeding_cpu_count_is_clipped_with_warning(self):
        with mock.patch("os.cpu_count", return_value=4):
            cmp = ChemoPy(ignore_3D=True)
            with self.assertWarns(UserWarning):
                values = cmp.calculate([SMALL_MOLECULES_3D[0]], show_banner=False, njobs=100)
            self.assertEqual(values.shape[0], 1)


class TestWorkerCorePinning(unittest.TestCase):
    """Tests that each pool worker is assigned a unique, permanent core id."""

    def setUp(self):
        # _pin_worker mutates a module-level global; restore it so other tests
        # in the same process (e.g. single-process `calculate()` calls) aren't
        # affected by a leftover fake core assignment.
        self._original_core_id = chemopy._WORKER_CORE_ID
        self.addCleanup(setattr, chemopy, "_WORKER_CORE_ID", self._original_core_id)

    def test_pin_worker_assigns_sequential_distinct_ids(self):
        counter = multiprocessing.Value("i", 0)
        seen = []
        for _ in range(5):
            _pin_worker(counter)
            seen.append(chemopy._WORKER_CORE_ID)
        self.assertEqual(seen, [0, 1, 2, 3, 4])
        self.assertEqual(len(set(seen)), 5)

    def test_default_worker_core_id_is_none_outside_a_pool(self):
        # In the parent process (not a pool worker), no core has been assigned.
        self.assertIsNone(self._original_core_id)


if __name__ == "__main__":
    unittest.main()
