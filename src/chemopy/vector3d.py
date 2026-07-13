# -*- coding: utf-8 -*-


"""Class to manipulate 3D vectors."""

from __future__ import annotations

import math

import numpy as np

RAD2DEG = 180.0 / math.pi
DEG2RAD = math.pi / 180.0
SMALL = 1e-6


def is_near_zero(x: float, epsilon: float = None) -> bool:
    """Return if x is lower than epsilon.

    :param x: small number to test
    :param epsilon: extremelly small number x needs to be smaller than
                    (default: 1e-6)
    """
    return abs(x) < SMALL if epsilon is None else abs(x) < epsilon


class Vector3D:
    """Class holding 3D vector coordinates."""

    def __init__(self, x: float = 0.0, y: float = 0.0, z: float = 0.0) -> None:
        """Initialize a Vector3D."""
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    def __add__(self, rhs: Vector3D) -> Vector3D:
        """Return the resulting vector from adding another vector to the current one.

        :param rhs: other Vector3D.
        """
        return Vector3D(rhs.x + self.x, rhs.y + self.y, rhs.z + self.z)

    def __sub__(self, rhs: Vector3D) -> Vector3D:
        """Return the resulting vector from substracting another vector from the current one.

        :param rhs: other Vector3D.
        """
        return Vector3D(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)

    def __mul__(self, rhs: Vector3D) -> Vector3D:
        """Return the Hadamard product of the two vectors."""
        return Vector3D(self.x * rhs.x, self.y * rhs.y, self.z * rhs.z)

    def __pow__(self, n: int) -> Vector3D:
        """Execute n times the Hadamard product."""
        temp = self.copy()
        acc = Vector3D(1.0, 1.0, 1.0)
        for _ in range(n):
            acc *= temp
        return acc

    def __neg__(self) -> Vector3D:
        """Return the opposite vector of the current one."""
        return Vector3D(-self.x, -self.y, -self.z)

    def __pos__(self) -> Vector3D:
        """Return a copy of the current vector."""
        return Vector3D(self.x, self.y, self.z)

    def __eq__(self, rhs: Vector3D) -> bool:
        """Test for close equality between two vectors.

        :param rhs: other Vector3D.
        """
        return is_near_zero(self.x - rhs.x) and is_near_zero(self.y - rhs.y) and is_near_zero(self.z - rhs.z)

    def __ne__(self, rhs: Vector3D) -> bool:
        """Test for inequality between two vectors.

        :param rhs: other Vector3D.
        """
        return not (self == rhs)

    def __str__(self) -> str:
        """Pretty print representation."""
        return f"({self.x:.2f}, {self.y:.2f}, {self.z:.2f})"

    def __repr__(self) -> str:
        """Unambiguous representation."""
        return f"Vector3D({self.x:.2f}, {self.y:.2f}, {self.z:.2f})"

    def copy(self) -> Vector3D:
        """Copy the current vector."""
        return Vector3D(self.x, self.y, self.z)

    def dot(self, rhs: Vector3D) -> Vector3D:
        """Return the dot product of the vector with another one."""
        return self.x * rhs.x + self.y * rhs.y + self.z * rhs.z

    def cross(self, rhs: Vector3D) -> Vector3D:
        """Return the cross product of the vector with another one."""
        return Vector3D(
            self.y * rhs.z - self.z * rhs.y, self.z * rhs.x - self.x * rhs.z, self.x * rhs.y - self.y * rhs.x
        )

    def length_sq(self) -> float:
        """Return the square value of the euclidean norm."""
        return self.x * self.x + self.y * self.y + self.z * self.z

    def length(self) -> float:
        """Return the euclidean norm."""
        return math.sqrt(self.length_sq())

    def scale(self, scale: float) -> Vector3D:
        """Scale the vector by a scalar.

        :param scale: value to scale each coordinate by.
        """
        self.x *= scale
        self.y *= scale
        self.z *= scale

    def normalize(self) -> Vector3D:
        """Scale the vector to unit length."""
        self.scale(1.0 / self.length())

    def scaled_vec(self, scale: float) -> Vector3D:
        """Return a copy scaled to specified length.

        :param scale: value to scale each coordinate by.
        """
        v = self.copy()
        v.scale(scale)
        return v

    def normal_vec(self) -> Vector3D:
        """Return a copy scaled to unit length."""
        return self.scaled_vec(1.0 / self.length())

    def parallel_vec(self, axis: Vector3D) -> Vector3D:
        """Return a vector parallel to the current one using the specified axis.

        :param axis: axis along which to search for the parallel vector
        """
        axis_len = axis.length()
        if is_near_zero(axis_len):
            result = self
        else:
            result = axis.scaled_vec(self.dot(axis) / axis.length() / axis.length())
        return result

    def perpendicular_vec(self, axis) -> Vector3D:
        """Return a vector perpendicular to the current one using the specified axis.

        :param axis: axis along which to search for the perpendicular vector
        """
        return self - self.parallel_vec(axis)

    def transform(self, matrix) -> Vector3D:
        """Apply transformation from a matrix to the vector.

        :param matrix: transforamtion matrix
        """
        x = matrix.elem00 * self.x + matrix.elem10 * self.y + matrix.elem20 * self.z + matrix.elem30
        y = matrix.elem01 * self.x + matrix.elem11 * self.y + matrix.elem21 * self.z + matrix.elem31
        z = matrix.elem02 * self.x + matrix.elem12 * self.y + matrix.elem22 * self.z + matrix.elem32
        self.x, self.y, self.z = x, y, z

    def diff_length_sq(self, rhs: Vector3D) -> float:
        """Return the squared euclidean norm of the difference of the two vectors.

        :param rhs: other Vector3D.
        """
        diff = self - rhs
        return diff.dot(diff)

    def diff_length(self, rhs: Vector3D) -> float:
        """Return the squared euclidean norm of the difference of the two vectors.

        :param rhs: other Vector3D.
        """
        return math.sqrt(self.diff_length_sq(rhs))

    def to_numpy_row(self):
        """Transform to a numpy matrix row."""
        return np.array([self.x, self.y, self.z], dtype=float)

    def to_numpy_col(self):
        """Transform to a numpy matrix column."""
        return self.to_numpy_row().T

    @staticmethod
    def from_numpy(array: np.array) -> Vector3D:
        """Instanciate a Vector3D from numpy array."""
        if len(array) > 3:
            raise ValueError(f"numpy.array({array}) has more values than Vector3D can hold.")
        return Vector3D(array[0], array[1], array[2])


def pos_distance(p1, p2):
    """Return the euclidean norm of the difference of the two vectors."""
    return math.sqrt(pos_distance_sq(p2, p1))


def pos_distance_sq(p1, p2):
    """Return the squared euclidean norm of the difference of the two vectors."""
    x = p1.x - p2.x
    y = p1.y - p2.y
    z = p1.z - p2.z
    return x * x + y * y + z * z
