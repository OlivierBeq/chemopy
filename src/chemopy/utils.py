# -*- coding: utf-8 -*-


"""Utility functions for pychem."""

import math
import os
from pathlib import Path
from typing import List

import scipy
from openbabel import openbabel, pybel
from rdkit import Chem


def GetR(n: int) -> List[float]:
    """Calcuate the parameters R of the RDF equation."""
    R = []
    for i in range(2, n + 2):
        R.append(float(i * 0.5))
    return R


def GetAtomDistance(x: List[float], y: List[float]) -> float:
    """Calculate Euclidean distance between two atomic coordinates."""
    temp = [math.pow(x[0] - y[0], 2), math.pow(x[1] - y[1], 2), math.pow(x[2] - y[2], 2)]
    res = math.sqrt(sum(temp))
    return res


def GetGeometricalDistanceMatrix(CoordinateList: List[List[float]]) -> scipy.matrix:
    """Calculate distance matrix from coordinate list."""
    NAtom = len(CoordinateList)
    DistanceMatrix = scipy.zeros((NAtom, NAtom))
    for i in range(NAtom - 1):
        for j in range(i + 1, NAtom):
            DistanceMatrix[i, j] = GetAtomDistance(CoordinateList[i], CoordinateList[j])
            DistanceMatrix[j, i] = DistanceMatrix[i, j]
    return DistanceMatrix

