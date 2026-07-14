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
    """Tests that run_mopac engages exactly the requested number of distinct cores.

    `geo_opt.platform` (captured once at import time via `from sys import platform`) is
    mocked directly rather than relying on `sys.platform`, so each of the three platform
    branches (linux/win32/darwin) is exercised deterministically regardless of which OS
    actually runs this test file -- e.g. the Windows-specific bitmask/`/WAIT` logic and
    macOS's "no pinning available" fallback are both verified from a single Linux CI run.
    """

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

    def _last_core_ids(self, plat):
        """Extract the core ids the last run_mopac() call was pinned to, as a list of ints.

        precmd is either 'taskset --cpu-list <ids> mopac <file>' (Linux) or
        'START /WAIT /AFFINITY <hexmask> mopac <file>' (Windows).
        """
        cmd = self.calls[-1]
        if plat == "win32":
            mask = int(cmd.split()[3], 16)
            return [i for i in range(64) if mask & (1 << i)]
        return [int(x) for x in cmd.split()[2].split(",")]

    def test_default_single_core_linux(self):
        with mock.patch.object(geo_opt, "platform", "linux"):
            geo_opt.run_mopac("dummy.dat")
        self.assertEqual(len(self._last_core_ids("linux")), 1)

    def test_distinct_cores_without_replacement_linux(self):
        with mock.patch.object(geo_opt, "platform", "linux"):
            geo_opt.run_mopac("dummy.dat", n_jobs=4)
        ids = self._last_core_ids("linux")
        self.assertEqual(len(ids), 4)
        self.assertEqual(len(set(ids)), 4)

    def test_explicit_affinity_is_honored_exactly_linux(self):
        with mock.patch.object(geo_opt, "platform", "linux"):
            geo_opt.run_mopac("dummy.dat", affinity=[3])
        self.assertEqual(self._last_core_ids("linux"), [3])

    def test_oversubscription_is_clipped_with_warning_linux(self):
        with mock.patch.object(geo_opt, "platform", "linux"):
            with self.assertWarns(UserWarning):
                geo_opt.run_mopac("dummy.dat", n_jobs=20)
        ids = self._last_core_ids("linux")
        self.assertEqual(len(ids), 8)
        self.assertEqual(len(set(ids)), 8)

    def test_default_single_core_windows(self):
        with mock.patch.object(geo_opt, "platform", "win32"):
            geo_opt.run_mopac("dummy.dat")
        self.assertEqual(len(self._last_core_ids("win32")), 1)
        self.assertIn("/WAIT", self.calls[-1])

    def test_distinct_cores_without_replacement_windows(self):
        with mock.patch.object(geo_opt, "platform", "win32"):
            geo_opt.run_mopac("dummy.dat", n_jobs=4)
        ids = self._last_core_ids("win32")
        self.assertEqual(len(ids), 4)
        self.assertEqual(len(set(ids)), 4)

    def test_explicit_affinity_is_honored_exactly_windows(self):
        with mock.patch.object(geo_opt, "platform", "win32"):
            geo_opt.run_mopac("dummy.dat", affinity=[3])
        self.assertEqual(self._last_core_ids("win32"), [3])

    def test_oversubscription_is_clipped_with_warning_windows(self):
        with mock.patch.object(geo_opt, "platform", "win32"):
            with self.assertWarns(UserWarning):
                geo_opt.run_mopac("dummy.dat", n_jobs=20)
        ids = self._last_core_ids("win32")
        self.assertEqual(len(ids), 8)
        self.assertEqual(len(set(ids)), 8)

    def test_runs_unpinned_with_warning_on_macos(self):
        with mock.patch.object(geo_opt, "platform", "darwin"):
            with self.assertWarns(UserWarning):
                geo_opt.run_mopac("dummy.dat", affinity=[3])
        self.assertNotIn("taskset", self.calls[-1])
        self.assertNotIn("AFFINITY", self.calls[-1])
        self.assertEqual(self.calls[-1].strip(), "mopac dummy.dat")

    def test_unsupported_platform_raises(self):
        with mock.patch.object(geo_opt, "platform", "aix"):
            with self.assertRaises(RuntimeError):
                geo_opt.run_mopac("dummy.dat")


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
