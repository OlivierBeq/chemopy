# -*- coding: utf-8 -*-

#####
# $Id: atom_types.py 997 2009-02-25 06:12:43Z glandrum $
#
# Copyright (C) 2002-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#

"""SMARTS definitions and calculators for EState atom types.

Defined in: Hall and Kier JCICS _35_ 1039-1045 (1995)  Table 1
"""

import warnings
from typing import List, Tuple

import numpy as np
from rdkit import Chem

_raw_d = [
    ('sLi', '[LiD1]-*'),

    ('ssBe', '[BeD2](-*)-*'),
    ('ssssBe', '[BeD4](-*)(-*)(-*)-*'),

    ('ssBH', '[BD2H](-*)-*'),
    ('sssB', '[BD3](-*)(-*)-*'),
    ('ssssB', '[BD4](-*)(-*)(-*)-*'),

    ('sCH3', '[CD1H3]-*'),
    ('dCH2', '[CD1H2]=*'),
    ('ssCH2', '[CD2H2](-*)-*'),
    ('tCH', '[CD1H]#*'),
    ('dsCH', '[CD2H](=*)-*'),
    ('aaCH', '[C,c;D2H](:*):*'),
    ('sssCH', '[CD3H](-*)(-*)-*'),
    ('ddC', '[CD2H0](=*)=*'),
    ('tsC', '[CD2H0](#*)-*'),
    ('dssC', '[CD3H0](=*)(-*)-*'),
    ('aasC', '[C,c;D3H0](:*)(:*)-*'),
    ('aaaC', '[C,c;D3H0](:*)(:*):*'),
    ('ssssC', '[CD4H0](-*)(-*)(-*)-*'),

    ('sNH3', '[ND1H3]-*'),
    ('sNH2', '[ND1H2]-*'),
    ('ssNH2', '[ND2H2](-*)-*'),
    ('dNH', '[ND1H]=*'),
    ('ssNH', '[ND2H](-*)-*'),
    ('aaNH', '[N,nD2H](:*):*'),
    ('tN', '[ND1H0]#*'),
    ('sssNH', '[ND3H](-*)(-*)-*'),
    ('dsN', '[ND2H0](=*)-*'),
    ('aaN', '[N,nD2H0](:*):*'),
    ('sssN', '[ND3H0](-*)(-*)-*'),
    ('ddsN', '[ND3H0](~[OD1H0])(~[OD1H0])-,:*'),  # mod
    ('aasN', '[N,nD3H0](:*)(:*)-,:*'),              # mod
    ('ssssN', '[ND4H0](-*)(-*)(-*)-*'),

    ('sOH', '[OD1H]-*'),
    ('dO', '[OD1H0]=*'),
    ('ssO', '[OD2H0](-*)-*'),
    ('aaO', '[O,oD2H0](:*):*'),

    ('sF', '[FD1]-*'),

    ('sSiH3', '[SiD1H3]-*'),
    ('ssSiH2', '[SiD2H2](-*)-*'),
    ('sssSiH', '[SiD3H1](-*)(-*)-*'),
    ('ssssSi', '[SiD4H0](-*)(-*)(-*)-*'),

    ('sPH2', '[PD1H2]-*'),
    ('ssPH', '[PD2H1](-*)-*'),
    ('sssP', '[PD3H0](-*)(-*)-*'),
    ('dsssP', '[PD4H0](=*)(-*)(-*)-*'),
    ('sssssP', '[PD5H0](-*)(-*)(-*)(-*)-*'),

    ('sSH', '[SD1H1]-*'),
    ('dS', '[SD1H0]=*'),
    ('ssS', '[SD2H0](-*)-*'),
    ('aaS', '[S,sD2H0](:*):*'),
    ('dssS', '[SD3H0](=*)(-*)-*'),
    ('ddssS', '[SD4H0](~[OD1H0])(~[OD1H0])(-*)-*'),  # mod

    ('sCl', '[ClD1]-*'),

    ('sGeH3', '[GeD1H3](-*)'),
    ('ssGeH2', '[GeD2H2](-*)-*'),
    ('sssGeH', '[GeD3H1](-*)(-*)-*'),
    ('ssssGe', '[GeD4H0](-*)(-*)(-*)-*'),

    ('sAsH2', '[AsD1H2]-*'),
    ('ssAsH', '[AsD2H1](-*)-*'),
    ('sssAs', '[AsD3H0](-*)(-*)-*'),
    ('sssdAs', '[AsD4H0](=*)(-*)(-*)-*'),
    ('sssssAs', '[AsD5H0](-*)(-*)(-*)(-*)-*'),

    ('sSeH', '[SeD1H1]-*'),
    ('dSe', '[SeD1H0]=*'),
    ('ssSe', '[SeD2H0](-*)-*'),
    ('aaSe', '[SeD2H0](:*):*'),
    ('dssSe', '[SeD3H0](=*)(-*)-*'),
    ('ddssSe', '[SeD4H0](=*)(=*)(-*)-*'),

    ('sBr', '[BrD1]-*'),

    ('sSnH3', '[SnD1H3]-*'),
    ('ssSnH2', '[SnD2H2](-*)-*'),
    ('sssSnH', '[SnD3H1](-*)(-*)-*'),
    ('ssssSn', '[SnD4H0](-*)(-*)(-*)-*'),

    ('sI', '[ID1]-*'),

    ('sPbH3', '[PbD1H3]-*'),
    ('ssPbH2', '[PbD2H2](-*)-*'),
    ('sssPbH', '[PbD3H1](-*)(-*)-*'),
    ('ssssPb', '[PbD4H0](-*)(-*)(-*)-*'),
]


esPatterns = None


def build_patterns(rawV: List[Tuple[str]] = None) -> None:
    """Create SMARTS molecule from patters."""
    global esPatterns, _raw_d
    if rawV is None:
        rawV = _raw_d

    esPatterns = [None] * len(rawV)
    for i, (name, sma) in enumerate(rawV):
        try:
            patt = Chem.MolFromSmarts(sma)
        except Exception:
            warnings.warn(f'Problems with pattern {sma} (name: {name}), skipped.', ResourceWarning)
        else:
            esPatterns[i] = name, patt


def type_atoms(mol) -> List[Tuple[int]]:
    """Assign each atom to an EState type.

    :return AtomType match.
    #  list of tuples (atoms can possibly match multiple patterns) with atom types
    """
    if esPatterns is None:
        build_patterns()
    nAtoms = mol.GetNumAtoms()
    res = [None] * nAtoms
    for name, patt in esPatterns:
        matches = mol.GetSubstructMatches(patt, uniquify=0)
        for match in matches:
            idx = match[0]
            if res[idx] is None:
                res[idx] = [name]
            elif name not in res[idx]:
                res[idx].append(name)
    for i, v in enumerate(res):
        if v is not None:
            res[i] = tuple(v)
        else:
            res[i] = ()
    return res


def get_atom_label(mol: Chem.Mol) -> List[float]:
    """Calculate the atom index in a molecule for the above given atom types."""
    if esPatterns is None:
        build_patterns()
    res = []
    for _, patt in esPatterns:
        matches = mol.GetSubstructMatches(patt, uniquify=0)
        cc = []
        for match in matches:
            cc.append(match[0])
        bb = list(np.unique(np.array(cc)))
        res.append(bb)
    return res
