# -*- coding: utf-8 -*-


"""Class to manipulate 3D vectors."""


from __future__ import annotations

import math
import random
from typing import Iterable

import numpy as np

RAD2DEG = 180.0 / math.pi
DEG2RAD = math.pi / 180.0
SMALL = 1E-6


def is_near_zero(x: float, epsilon: float = None) -> bool:
    """Return if x is lower than epsilon.

    :param x: small number to test
    :param epsilon: extremelly small number x needs to be smaller than
                    (default: 1e-6)
    """
    return x < SMALL if epsilon is None else x < epsilon


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
        temp, acc = self.copy(), self.copy()
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
        return (is_near_zero(self.x - rhs.x) and is_near_zero(self.y - rhs.y) and is_near_zero(self.z - rhs.z))

    def __ne__(self, rhs: Vector3D) -> bool:
        """Test for inequality between two vectors.

        :param rhs: other Vector3D.
        """
        return not (self == rhs)

    def __str__(self) -> str:
        """Pretty print representation."""
        return f'({self.x:.2f}, {self.y:.2f}, {self.z:.2f})'

    def __repr__(self) -> str:
        """Unambiguous representation."""
        return f'Vector3D({self.x:.2f}, {self.y:.2f}, {self.z:.2f})'

    def copy(self) -> Vector3D:
        """Copy the current vector."""
        return Vector3D(self.x, self.y, self.z)

    def dot(self, rhs: Vector3D) -> Vector3D:
        """Return the dot product of the vector with another one."""
        return self.x * rhs.x + self.y * rhs.y + self.z * rhs.z

    def cross(self, rhs: Vector3D) -> Vector3D:
        """Return the cross product of the vector with another one."""
        return Vector3D(self.y * rhs.z - self.z * rhs.y,
                        self.z * rhs.x - self.x * rhs.z,
                        self.x * rhs.y - self.y * rhs.x)

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
        x = matrix.elem00 * self.x + \
            matrix.elem10 * self.y + \
            matrix.elem20 * self.z + \
            matrix.elem30
        y = matrix.elem01 * self.x + \
            matrix.elem11 * self.y + \
            matrix.elem21 * self.z + \
            matrix.elem31
        z = matrix.elem02 * self.x + \
            matrix.elem12 * self.y + \
            matrix.elem22 * self.z + \
            matrix.elem32
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
        return np.matrix([self.x, self.y, self.z], dtype=float)

    def to_numpy_col(self):
        """Transform to a numpy matrix column."""
        return self.to_numpy_row().T

    @staticmethod
    def from_numpy(array: np.array) -> Vector3D:
        """Instanciate a Vector3D from numpy array."""
        if len(array) > 3:
            raise ValueError(f'numpy.array({array}) has more values than Vector3D can hold.')
        return Vector3D(array[0], array[1], array[2])


def normalize_angle(angle: float):
    """Normalize angle to be between -pi and +pi.

    :param angle: angle in radians
    """
    while abs(angle) > math.pi:
        if angle > math.pi:
            angle -= math.pi * 2
        if angle < -math.pi:
            angle += 2 * math.pi
    if is_near_zero(abs(angle + math.pi)):
        angle = math.pi
    return angle


def angle_diff(angle1: float, angle2: float, radians: bool = True) -> float:
    """Calculate the difference between two angles.

    :param angle1: angle (in radians unless specified)
    :param angle2: angle (in radians unless specified)
    """
    if not radians:
        angle1 = angle1 * math.pi / 180
        angle2 = angle2 * math.pi / 180
    norm_angle1 = normalize_angle(angle1)
    norm_angle2 = normalize_angle(angle2)
    diff = normalize_angle(norm_angle1 - norm_angle2)
    return diff if radians else radians / math.pi * 180


def pos_distance(p1, p2):
    """Return the euclidean norm of the difference of the two vectors."""
    return math.sqrt(pos_distance_sq(p2, p1))


def pos_distance_sq(p1, p2):
    """Return the squared euclidean norm of the difference of the two vectors."""
    x = p1.x - p2.x
    y = p1.y - p2.y
    z = p1.z - p2.z
    return x * x + y * y + z * z


def vec_angle(a, b):
    """Get the angle between the two vectors."""
    a_len = a.length()
    b_len = b.length()
    if a_len * b_len < 1E-6:
        return 0.0
    c = a.dot(b) / a_len / b_len
    if c >= 1.0:
        return 0.0
    elif c <= -1.0:
        return math.pi
    else:
        return math.acos(c)


def pos_angle(p1, p2, p3):
    """Get the angle from three points."""
    return vec_angle(p1 - p2, p3 - p2)


def vec_dihedral(a, axis, c):
    """Get the dihedral angle between a and c along the axis."""
    ap = a.perpendicular_vec(axis)
    cp = c.perpendicular_vec(axis)
    angle = vec_angle(ap, cp)
    if ap.cross(cp).dot(axis) > 0:
        angle = -angle
    return angle


def pos_dihedral(p1, p2, p3, p4):
    """Get the dihedral angle from four points."""
    return vec_dihedral(p1 - p2, p2 - p3, p4 - p3)


def RotatedPos(theta, anchor, center, pos):
    """Rotate the position by theta around the center."""
    return Rotation(center - anchor, theta, center).transform_vec(pos)


# def ProjectedPos(length, angle, dihedral, p1, p2, p3):
#     """Project."""
#     norm = plane_normal(p1, p2, p3)
#     axis = p3 - p2
#     vec_diff = axis.scaled_vec(-length)
#     vec_diff = RotationAtOrigin(norm, -angle).transform_vec(vec_diff)
#     vec_diff = RotationAtOrigin(axis, dihedral).transform_vec(vec_diff)
    # return p3 + vec_diff


class Matrix3D(np.matrix):
    """Object to manipulate three dimensional matrices.

    Wrapper around numpy matrix object with shape (4, 4).
    """

    def __new__(cls, values: Iterable[float] = None, order='C') -> None:
        """Create a new Matrix3D object.

        :param values: values of the matrix
        :param order: index order of provided values
                      'C': column-first indexing (first items of values will make first row)
                      'F': row-first indexing (first items of values will make first column)
        """
        if values is None:
            self = np.zeros((4, 4), dtype=float)
        else:
            len_ = len(list(values))
            if not len_ == 16:
                raise ValueError(f'Matrix3D cannot be instanciated with {len_} values')
            self = np.matrix(values, dtype=float).reshape((4, 4), order=order)
        return self.view(cls)

    def transform_vector(self, vector: Vector3D) -> Vector3D:
        """Apply the transformation to the vector.

        :param vector: vector to be transformed
        """
        return Vector3D.from_numpy(self[:3, :3].T * vector.to_numpy_col() + self[3, :3].T)


def RotationAtOrigin(axis, theta):
    """Create a matrix for a rotation around the origin.

    :param axis: axis around which the rotation should be performed
    :param theta: angle of the rotation
    """
    v = axis.normal_vec()
    c = math.cos(float(theta))
    s = math.sin(float(theta))
    t = 1.0 - c
    m = Matrix3D([[t * v.x * v.x + c,
                   t * v.x * v.y + v.z * s,
                   t * v.x * v.z - v.y * s,
                   0],
                  [t * v.y * v.x - v.z * s,
                   t * v.y * v.y + c,
                   t * v.y * v.z + v.x * s,
                   0],
                  [t * v.z * v.x + v.y * s,
                   t * v.z * v.y - v.x * s,
                   t * v.z * v.z + c,
                   0],
                  [0, 0, 0, 0]])
    return m


def Translation(vector):
    """Create a translation matrix.

    :param vector: vector of translation
    """
    m = Matrix3D()
    m[3, :3] = [vector.x, vector.y, vector.z]
    return m


def Rotation(axis, theta, center):
    """Create the matrix to rotate around an axis at center."""
    rot = RotationAtOrigin(axis, theta)
    trans = Translation(center - rot.transform_vec(center))
    return trans * rot


def Superposition3(ref1, ref2, ref3, mov1, mov2, mov3):
    """Superimpose."""
    mov_diff = mov2 - mov1
    ref_diff = ref2 - ref1

    m1 = Matrix3D()
    if math.fabs(vec_angle(mov_diff, ref_diff)) < SMALL:
        m1 = Translation(ref1 - mov1)
    else:
        axis = mov_diff.cross(ref_diff)
        torsion = vec_dihedral(ref_diff, axis, mov_diff)
        rot = RotationAtOrigin(axis, torsion)
        trans = Translation(ref2 - rot.transform_vec(mov2))
        m1 = trans * rot

    mov_diff = ref2 - m1.transform_vec(mov3)
    ref_diff = ref2 - ref3

    m = Matrix3D()
    if math.fabs(vec_angle(mov_diff, ref_diff)) < SMALL:
        m = m1
    else:
        axis = ref2 - ref1
        torsion = vec_dihedral(ref_diff, axis, mov_diff)
        m2 = RotationAtOrigin(axis, torsion)
        m3 = Translation(ref2 - m2.transform_vec(ref2))
        m = m3 * m2 * m1

    return m


def RandomVec():
    """Get a random vector."""
    return Vector3D(random.uniform(-100, 100),  # noqa: S311
                    random.uniform(-100, 100),  # noqa: S311
                    random.uniform(-100, 100))  # noqa: S311


def RandomOriginRotation():
    """Get a random rotation around the origin."""
    axis = RandomVec()
    angle = random.uniform(-math.pi / 2, math.pi / 2)  # noqa: S311
    return RotationAtOrigin(axis, angle)


def RandomTransform():
    """Get a random transformation matrix."""
    axis = RandomVec()
    angle = random.uniform(-math.pi / 2, math.pi / 2)  # noqa: S311
    center = RandomVec()
    return Rotation(axis, angle, center)
