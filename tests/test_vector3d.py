# -*- coding: utf-8 -*-

"""Tests for the Vector3D class and related helper functions."""

import math
import unittest

from chemopy import vector3d
from chemopy.vector3d import Vector3D


class TestIsNearZero(unittest.TestCase):
    def test_default_epsilon(self):
        self.assertTrue(vector3d.is_near_zero(1e-9))
        self.assertFalse(vector3d.is_near_zero(1.0))

    def test_custom_epsilon(self):
        self.assertTrue(vector3d.is_near_zero(0.5, epsilon=1.0))
        self.assertFalse(vector3d.is_near_zero(2.0, epsilon=1.0))


class TestVector3D(unittest.TestCase):
    def test_construction_defaults(self):
        v = Vector3D()
        self.assertEqual((v.x, v.y, v.z), (0.0, 0.0, 0.0))

    def test_add_sub(self):
        a = Vector3D(1, 2, 3)
        b = Vector3D(4, 5, 6)
        self.assertEqual(a + b, Vector3D(5, 7, 9))
        self.assertEqual(b - a, Vector3D(3, 3, 3))

    def test_hadamard_mul_and_pow(self):
        a = Vector3D(1, 2, 3)
        b = Vector3D(2, 2, 2)
        self.assertEqual(a * b, Vector3D(2, 4, 6))
        self.assertEqual(a**2, Vector3D(1, 4, 9))

    def test_neg_pos(self):
        a = Vector3D(1, -2, 3)
        self.assertEqual(-a, Vector3D(-1, 2, -3))
        self.assertEqual(+a, a)

    def test_eq_ne(self):
        a = Vector3D(1, 2, 3)
        b = Vector3D(1, 2, 3 + 1e-9)
        c = Vector3D(1, 2, 4)
        self.assertEqual(a, b)
        self.assertNotEqual(a, c)
        self.assertTrue(a != c)

    def test_str_repr(self):
        a = Vector3D(1, 2, 3)
        self.assertEqual(str(a), "(1.00, 2.00, 3.00)")
        self.assertIn("Vector3D", repr(a))

    def test_copy_is_independent(self):
        a = Vector3D(1, 2, 3)
        b = a.copy()
        b.x = 100
        self.assertEqual(a.x, 1)

    def test_dot_and_cross(self):
        a = Vector3D(1, 0, 0)
        b = Vector3D(0, 1, 0)
        self.assertEqual(a.dot(b), 0)
        self.assertEqual(a.cross(b), Vector3D(0, 0, 1))

    def test_length_and_length_sq(self):
        a = Vector3D(3, 4, 0)
        self.assertEqual(a.length_sq(), 25)
        self.assertEqual(a.length(), 5)

    def test_scale_and_normalize(self):
        a = Vector3D(1, 0, 0)
        a.scale(3)
        self.assertEqual(a, Vector3D(3, 0, 0))
        a.normalize()
        self.assertAlmostEqual(a.length(), 1.0)

    def test_scaled_vec_and_normal_vec(self):
        a = Vector3D(2, 0, 0)
        self.assertEqual(a.scaled_vec(2), Vector3D(4, 0, 0))
        self.assertEqual(a.normal_vec(), Vector3D(1, 0, 0))
        # original vector is untouched by scaled_vec/normal_vec
        self.assertEqual(a, Vector3D(2, 0, 0))

    def test_parallel_and_perpendicular_vec(self):
        a = Vector3D(1, 1, 0)
        axis = Vector3D(1, 0, 0)
        parallel = a.parallel_vec(axis)
        self.assertEqual(parallel, Vector3D(1, 0, 0))
        perpendicular = a.perpendicular_vec(axis)
        self.assertEqual(perpendicular, Vector3D(0, 1, 0))

    def test_parallel_vec_zero_length_axis(self):
        a = Vector3D(1, 1, 0)
        axis = Vector3D(0, 0, 0)
        self.assertEqual(a.parallel_vec(axis), a)

    def test_diff_length_sq_and_diff_length(self):
        a = Vector3D(0, 0, 0)
        b = Vector3D(3, 4, 0)
        self.assertEqual(a.diff_length_sq(b), 25)
        self.assertEqual(a.diff_length(b), 5)

    def test_to_numpy_row_and_col(self):
        a = Vector3D(1, 2, 3)
        row = a.to_numpy_row()
        self.assertEqual(list(row), [1.0, 2.0, 3.0])

    def test_from_numpy(self):
        a = Vector3D.from_numpy([1, 2, 3])
        self.assertEqual(a, Vector3D(1, 2, 3))

    def test_from_numpy_too_many_values_raises(self):
        with self.assertRaises(ValueError):
            Vector3D.from_numpy([1, 2, 3, 4])


class TestPosDistance(unittest.TestCase):
    def test_pos_distance(self):
        a = Vector3D(0, 0, 0)
        b = Vector3D(3, 4, 0)
        self.assertEqual(vector3d.pos_distance(a, b), 5)

    def test_pos_distance_sq(self):
        a = Vector3D(0, 0, 0)
        b = Vector3D(3, 4, 0)
        self.assertEqual(vector3d.pos_distance_sq(a, b), 25)


if __name__ == "__main__":
    unittest.main()
