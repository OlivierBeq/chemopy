# -*- coding: utf-8 -*-


"""Main classes for the calculation of molecular descriptors."""

import os
from copy import deepcopy
from typing import List, Optional

import more_itertools
import numpy as np
import pandas as pd
from bounded_pool_executor import BoundedProcessPoolExecutor
from rdkit import Chem

from . import getmol, geo_opt
from .basak import Basak
from .bcut import BCUT
from .charge import Charge
from .connectivity import Connectivity
from .constitution import Constitution
from .cpsa import CPSA
from .estate import EState
from .fingerprint import Fingerprint
from .fingerprint3d import Fingerprint3D
from .geary import Geary
from .geometric import Geometric
from .kappa import Kappa
from .moe import MOE
from .molproperty import MolecularProperties
from .moran import Moran
from .moreaubroto import MoreauBroto
from .morse import MoRSE
from .quanchem import QuantumChemistry
from .rdf import RDF
from .topology import Topology
from .whim import WHIM


class ChemoPy:
    """Molecular descriptor calculator."""

    def __init__(self, ignore_3D: bool = True, include_fps: bool = False):
        """Instantiate a molecular descriptor calculator.

        :param ignore_3D: Avoid calculating 3D molecular descriptors (default: True).
        If `False`, then molecules must have 3D conformers.
        :param include_fps: Should 2D fingerprints be also calculated (default: False).
        """
        self.ignore_3D = ignore_3D
        self.include_fps = include_fps

    def calculate(self, mols: List[Chem.Mol], show_banner: bool = True, njobs: int = 1,
                  chunksize: Optional[int] = 1000) -> pd.DataFrame:
        """Calculate molecular descriptors and/or fingerprints.

        :param mols: RDKit molecules for which descriptors/fingerprints should be calculated
        :param show_banner: If True, show notice on ChemoPy's usage
        :param njobs: number of concurrent processes
        :param chunksize: number of molecules to be processed by a process; ignored if njobs is 1
        :return: a pandas DataFrame containing all ChemoPy descriptor/fingerprint values
        """
        if show_banner:
            self._show_banner()
        if njobs < 0:
            njobs = os.cpu_count() - njobs + 1
        # Parallelize should need be
        if njobs > 1:
            with BoundedProcessPoolExecutor(max_workers=njobs) as worker:
                futures = [worker.submit(self._calculate, list(chunk))
                           for chunk in more_itertools.batched(mols, chunksize)
                           ]
            return pd.concat([future.result()
                              for future in futures]
                             ).reset_index(drop=True).fillna(0)
        # Single process
        return self._calculate(list(mols))

    def _show_banner(self):
        """Print info message for citing."""
        print("""ChemoPy computes 1178 descriptor values (632 2D and 546 3D)
of 16 different chemical feature groups along with 11 2D and 1 3D fingerprint.
It relies the semi-empirical quantum chemistry program MOPAC to compute all
of its 3D molecular descriptors and fingerprint.

###################################

Should you publish results based on the ChemoPy descriptors, please cite:

ChemoPy: freely available python package for computational biology and chemoinformatics
Dong-Sheng Cao, Qing-Song Xu, Qian-Nan Hu and Yi-Zeng Liang
Bioinformatics 2013, 29(8), 1092–1094
DOI: 10.1093/bioinformatics/btt105

###################################

""")

    def _calculate(self, mols: List[Chem.Mol]) -> pd.DataFrame:
        """Calculate ChemoPy descriptors/fingerprints on one process.

        :param mols: RDKit molecules for which chemoPy descriptors should be calculated
        :return: a pandas DataFrame containing all chemoPy descriptor values
        """
        descs_2D = [(Constitution, 'get_all'), (Topology, 'get_all'), (Connectivity, 'get_all'), (Kappa, 'get_all'), (Basak, 'get_all'),
                    (EState, 'get_all_fps'), (EState, 'get_all_descriptors'), (BCUT, 'get_all'), (MoreauBroto, 'get_all'),
                    (Moran, 'get_all'), (Geary, 'get_all'), (Charge, 'get_all'), (MolecularProperties, 'get_all'), (MOE, 'get_all')]
        descs_3D = []
        # Sanity check for 3D descriptors
        if not self.ignore_3D:
            for mol in mols:
                confs = list(mol.GetConformers())
                if not (len(confs) > 0 and confs[-1].Is3D()):
                    raise ValueError('Cannot calculate 3D descriptors for a conformer-less molecule')
            descs_3D.extend([(Geometric, 'get_all'), (CPSA, 'get_all') , (WHIM, 'get_all'),
                             (MoRSE, 'get_all'), (RDF, 'get_all'), (QuantumChemistry, 'get_all')])
        # Include fingerprints?
        if self.include_fps:
            descs_2D.append((Fingerprint, 'get_all_fps'))
            # Include 3D fingerprint(s)
            if not self.ignore_3D:
                descs_3D.append((Fingerprint3D, 'calculate_e3fp'))
        # Calculate descriptors
        all_values = []
        skipped_mols = []
        for i, mol in enumerate(mols):
            mol_values = {}
            # 2D descriptors
            for desc_2D, fn in descs_2D:
                try:
                    tmp = getattr(desc_2D, fn)(mol)
                    mol_values.update(tmp)
                except Exception as e:
                    raise e
                    pass
            # 3D descriptors
            if len(descs_3D):
                try:
                    dir_, arc_file = geo_opt.get_arc_file(mol, verbose=False)
                    mol3d = geo_opt.get_optimized_mol(arc_file=arc_file)
                    for desc_3D, fn in descs_3D:
                        try:
                            mol_values.update(getattr(desc_3D, fn)(mol=mol3d, arc_file=arc_file))
                        except Exception as e:
                            raise e
                            pass
                    # Dispose of MOPAC files
                    try:
                        dir_.close()
                    except Exception as e:
                        raise e
                        pass
                except Exception as e:
                    raise e
                    # Add molecule index to those skipped
                    skipped_mols.append(i)
            all_values.append(mol_values)
        # Transform to pandas dataframe
        results = pd.DataFrame(all_values)
        if results.empty:
            return results
        # Insert lines for skipped molecules
        if len(skipped_mols):
            results = pd.DataFrame(np.insert(results.values, skipped_mols,
                                             values=[np.NaN] * len(results.columns),
                                             axis=0),
                                   columns=results.columns)
        return results

    def _multiproc_calculate(self, mols: List[Chem.Mol]) -> pd.DataFrame:
        """Calculate ChemoPy descriptors/fingerprints in thread-safe manner.

        :param mols: RDKit molecules for which PaDEL descriptors should be calculated.
         Only the last conformer of molecules is considered.
        :return: a pandas DataFrame containing all PaDEL desciptor values and the path to the temp dir to be removed
        """
        # Copy self instance to make thread safe
        cmp = deepcopy(self)
        # Run copy
        result = cmp.calculate(mols, show_banner=False, njobs=1)
        return result

    @staticmethod
    def get_details(desc_name: Optional[str] = None) -> pd.DataFrame:
        """Obtain details about one/all descriptor(s)/fingerprint(s).

        :param desc_name: Name of descriptor/fingerprint (default: None).
        If `None`, obtain details for all.
        :return: The details about the descriptor(s)/fingerprint(s)
        """
        details = pd.read_json(os.path.abspath(os.path.join(__file__, os.pardir, 'descs.json')), orient='index')
        if desc_name is not None:
            if desc_name not in details.Name.tolist():
                raise ValueError(f'descriptor name {desc_name} is not available')
            details = details[details.Name == desc_name]
        return details
