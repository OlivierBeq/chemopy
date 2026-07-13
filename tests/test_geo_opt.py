# -*- coding: utf-8 -*-

"""Tests for MOPAC geometry optimization and core-usage correctness."""

import unittest
from unittest import mock

from chemopy import geo_opt


class TestMopacAvailability(unittest.TestCase):
    """Tests MOPAC configuration lookups."""

    def test_is_mopac_version_available_true(self):
        self.assertTrue(geo_opt.is_mopac_version_available("2016"))

    def test_is_mopac_version_available_unknown_version(self):
        self.assertFalse(geo_opt.is_mopac_version_available("1999"))

    def test_is_method_supported_by_mopac(self):
        self.assertTrue(geo_opt.is_method_supported_by_mopac("PM7", "2016"))
        self.assertFalse(geo_opt.is_method_supported_by_mopac("BOGUS", "2016"))

    def test_is_method_supported_by_mopac_unknown_version(self):
        with self.assertRaises(ValueError):
            geo_opt.is_method_supported_by_mopac("PM7", "1999")


class TestRunMopacCoreUsage(unittest.TestCase):
    """Tests that run_mopac engages exactly the requested number of distinct cores."""

    def setUp(self):
        self.calls = []
        patcher_call = mock.patch("subprocess.call", side_effect=self._fake_call)
        patcher_avail = mock.patch.object(geo_opt, "is_mopac_version_available", return_value=True)
        patcher_cpu = mock.patch("multiprocessing.cpu_count", return_value=8)
        self.addCleanup(patcher_call.stop)
        self.addCleanup(patcher_avail.stop)
        self.addCleanup(patcher_cpu.stop)
        patcher_call.start()
        patcher_avail.start()
        patcher_cpu.start()

    def _fake_call(self, cmd, **kwargs):
        self.calls.append(cmd)
        return 0

    def _last_core_ids(self):
        # precmd is 'taskset --cpu-list <ids> mopac <file>'
        return self.calls[-1].split()[2].split(",")

    def test_default_single_core(self):
        geo_opt.run_mopac("dummy.dat")
        self.assertEqual(len(self._last_core_ids()), 1)

    def test_distinct_cores_without_replacement(self):
        geo_opt.run_mopac("dummy.dat", n_jobs=4)
        ids = self._last_core_ids()
        self.assertEqual(len(ids), 4)
        self.assertEqual(len(set(ids)), 4)

    def test_explicit_affinity_is_honored_exactly(self):
        geo_opt.run_mopac("dummy.dat", affinity=[3])
        self.assertEqual(self._last_core_ids(), ["3"])

    def test_oversubscription_is_clipped_with_warning(self):
        with self.assertWarns(UserWarning):
            geo_opt.run_mopac("dummy.dat", n_jobs=20)
        ids = self._last_core_ids()
        self.assertEqual(len(ids), 8)
        self.assertEqual(len(set(ids)), 8)


class TestGetArcFileAffinity(unittest.TestCase):
    """Tests that get_arc_file forwards affinity down to run_mopac."""

    def test_affinity_is_forwarded_to_run_mopac(self):
        with (
            mock.patch.object(geo_opt, "format_conversion") as mock_format,
            mock.patch.object(geo_opt, "run_mopac", return_value=0) as mock_run,
            mock.patch("os.path.isfile", return_value=True),
            mock.patch.object(geo_opt, "get_file_in_dir_from_ext", return_value=["out.arc"]),
            mock.patch("os.rename"),
        ):
            mock_dir = mock.Mock()
            mock_dir.path = "/tmp/fake"
            mock_format.return_value = (mock_dir, "temp.dat")
            geo_opt.get_arc_file(mock.Mock(), n_jobs=3, affinity=[5, 6, 7], verbose=False)
            mock_run.assert_called_once()
            self.assertEqual(mock_run.call_args.kwargs.get("affinity"), [5, 6, 7])
            self.assertEqual(mock_run.call_args.kwargs.get("n_jobs"), 3)


if __name__ == "__main__":
    unittest.main()
