# -*- coding: utf-8 -*-


"""Routines to calculate the Accessible Surface Area of a set of atoms.

The algorithm is adapted from the Rose lab's chasa.py, which uses
the dot density technique found in:

Shrake, A., and J. A. Rupley. "Environment and Exposure to Solvent
of Protein Atoms. Lysozyme and Insulin." JMB (1973) 79:351-371.
"""


import math
from typing import List

from chemopy.GeoOpt import Atom
from chemopy.vector3d import Vector3D, pos_distance, pos_distance_sq

# TODO: change list[x, y, z] to Vector3D


def generate_sphere_points(n: int) -> List[float]:
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
        points.append([math.cos(phi) * r, y, math.sin(phi) * r])
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
    radius = atom_k.radius + probe + probe
    indices = list(range(k))
    indices.extend(range(k + 1, len(atoms)))
    for i in indices:
        atom_i = atoms[i]
        dist = pos_distance(atom_k.pos, atom_i.pos)
        if dist < radius + atom_i.radius:
            neighbor_indices.append(i)
    return neighbor_indices


def calculate_asa(atoms: List[Atom], probe: float, n_sphere_point: int = 960) -> List[float]:
    """Return the accessible surface area of the atoms.

    :param atoms: list of atoms to get the surface area of
    :param probe: radius of the probe
    :param n_sphere_point: number of evenly distributed
                           points along the atoms' surface
    """
    sphere_points = generate_sphere_points(n_sphere_point)

    const = 4.0 * math.pi / len(sphere_points)
    test_point = Vector3D()
    areas = []
    for i, atom_i in enumerate(atoms):

        neighbor_indices = find_neighbor_indices(atoms, probe, i)
        n_neighbor = len(neighbor_indices)
        j_closest_neighbor = 0
        radius = probe + atom_i.radius

        n_accessible_point = 0
        for point in sphere_points:
            is_accessible = True

            test_point.x = point[0] * radius + atom_i.pos.x
            test_point.y = point[1] * radius + atom_i.pos.y
            test_point.z = point[2] * radius + atom_i.pos.z

            cycled_indices = list(range(j_closest_neighbor, n_neighbor))
            cycled_indices.extend(range(j_closest_neighbor))

            for j in cycled_indices:
                atom_j = atoms[neighbor_indices[j]]
                r = atom_j.radius + probe
                diff_sq = pos_distance_sq(atom_j.pos, test_point)
                if diff_sq < r * r:
                    j_closest_neighbor = j
                    is_accessible = False
                    break
            if is_accessible:
                n_accessible_point += 1

        area = const * n_accessible_point * radius * radius
        areas.append(area)
    return areas


# if __name__ == "__main__":
#    main()
#    print(generate_sphere_points(10))
