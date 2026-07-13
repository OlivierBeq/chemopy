# -*- coding: utf-8 -*-


"""Routines to calculate the Accessible Surface Area of a set of atoms.

The algorithm is adapted from the Rose lab's chasa.py, which uses
the dot density technique found in:

Shrake, A., and J. A. Rupley. "Environment and Exposure to Solvent
of Protein Atoms. Lysozyme and Insulin." JMB (1973) 79:351-371.
"""

import math
from typing import List

import numpy as np
from numba import njit

from .geo_opt import Atom
from .vector3d import Vector3D, pos_distance


def generate_sphere_points(n: int) -> List[Vector3D]:
    """Distribute points on a sphere using the Golden Section Spiral algorithm.

    :param n: number of points to generate along the sphere surface.
    """
    points = []
    inc = math.pi * (3 - math.sqrt(5))
    offset = 2 / float(n)
    for k in range(int(n)):
        y = k * offset - 1 + (offset / 2)
        r = math.sqrt(1 - y * y)
        phi = k * inc
        points.append(Vector3D(math.cos(phi) * r, y, math.sin(phi) * r))
        # points.append([math.cos(phi) * r, y, math.sin(phi) * r])
    return points


def find_neighbor_indices(atoms: List[Atom], probe: float, k: int) -> List[int]:
    """Return indices of atoms within probe distance to atom k.

    If another atom u is found to be distant from atom k by less than
    the sum of its radius, the radius of atom k, and the diameter of
    the probe, then the probe cannot fit in between those and atom u.
    Hnce atom u is considered a neighbour of atom k.

    :param atoms: list of atoms
    :param probe: radius of probe
    :param k: atom from which the search is performed.
    """
    neighbor_indices = []
    atom_k = atoms[k]
    radius = atom_k.radius + 2 * probe
    indices = list(range(k))
    indices.extend(range(k + 1, len(atoms)))
    for i in indices:
        atom_i = atoms[i]
        dist = pos_distance(atom_k.pos, atom_i.pos)
        if dist < radius + atom_i.radius:
            neighbor_indices.append(i)
    return neighbor_indices


@njit(cache=True)
def _calculate_asa_jit(positions: np.ndarray, radii: np.ndarray, sphere_points: np.ndarray, probe: float) -> np.ndarray:
    """JIT-compiled core of the Shrake-Rupley accessible surface area algorithm.

    Numerically identical to the reference nested-loop implementation (same neighbor
    search radius, same "start from last-known-overlapping-neighbor" cycling order and
    early break on first burial), just operating on plain arrays so numba can compile it,
    instead of looping over Vector3D/Atom Python objects.
    """
    n_atoms = positions.shape[0]
    n_points = sphere_points.shape[0]
    const = 4.0 * np.pi / n_points
    areas = np.empty(n_atoms, dtype=np.float64)
    for i in range(n_atoms):
        # Neighbor search: same widened radius (2 * probe) as the reference find_neighbor_indices.
        search_radius = radii[i] + 2 * probe
        neighbor_indices = np.empty(n_atoms - 1, dtype=np.int64)
        n_neighbor = 0
        for j in range(n_atoms):
            if j == i:
                continue
            dx = positions[i, 0] - positions[j, 0]
            dy = positions[i, 1] - positions[j, 1]
            dz = positions[i, 2] - positions[j, 2]
            dist = math.sqrt(dx * dx + dy * dy + dz * dz)
            if dist < search_radius + radii[j]:
                neighbor_indices[n_neighbor] = j
                n_neighbor += 1

        radius = probe + radii[i]
        j_closest_neighbor = 0
        n_accessible_point = 0
        for p in range(n_points):
            test_x = sphere_points[p, 0] * radius + positions[i, 0]
            test_y = sphere_points[p, 1] * radius + positions[i, 1]
            test_z = sphere_points[p, 2] * radius + positions[i, 2]

            is_accessible = True
            for c in range(n_neighbor):
                j = (j_closest_neighbor + c) % n_neighbor
                jj = neighbor_indices[j]
                r = radii[jj] + probe
                dx = positions[jj, 0] - test_x
                dy = positions[jj, 1] - test_y
                dz = positions[jj, 2] - test_z
                diff_sq = dx * dx + dy * dy + dz * dz
                if diff_sq < r * r:
                    j_closest_neighbor = j
                    is_accessible = False
                    break
            if is_accessible:
                n_accessible_point += 1

        areas[i] = const * n_accessible_point * radius * radius
    return areas


def calculate_asa(atoms: List[Atom], probe: float, n_sphere_point: int = 960) -> List[float]:
    """Return the accessible surface area of the atoms.

    :param atoms: list of atoms to get the surface area of
    :param probe: radius of the probe
    :param n_sphere_point: number of evenly distributed
                           points along the atoms' surface
    """
    sphere_points = generate_sphere_points(n_sphere_point)
    sphere_points_arr = np.array([[p.x, p.y, p.z] for p in sphere_points], dtype=np.float64)
    positions = np.array([[atom.pos.x, atom.pos.y, atom.pos.z] for atom in atoms], dtype=np.float64)
    radii = np.array([atom.radius for atom in atoms], dtype=np.float64)
    areas = _calculate_asa_jit(positions, radii, sphere_points_arr, probe)
    return list(areas)
