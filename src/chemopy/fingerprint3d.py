# -*- coding: utf-8 -*-


"""Multiple 3D molecular fingerprints."""

from typing import Union, List, Optional

from e3fp.fingerprint.fprinter import Fingerprinter
from rdkit import Chem

from .geo_opt import get_optimized_mol


class Fingerprint3D:
    """Molecular 3D fingerprints."""

    @staticmethod
    def calculate_e3fp(mol: Union[Chem.Mol, List[Chem.Mol]], arc_file: str = None, nbits: int = 2048, level: int = 5,
                       radius_multiplier: float = 1.718, stereo: bool = True, counts: bool = False,
                       include_disconnected: bool = True, rdkit_invariants: bool = False,
                       exclude_floating: bool = True, remove_duplicate_substructs: bool = True) -> Union[dict, List[dict]]:
        """Calculate extended 3D fingerprint(s).

        :param mol: MOPAC optimized molecule (result of `geo_opt.get_optimized_mol`)
        :param arc_file: Ignored
        :param nbits: number of bits in the folded fingerprint
        :param level: maximum number of iterations for fingerprint generation
        :param radius_multiplier: ratio by which to increase shell size
        :param stereo: differentiate based on stereography
        :param counts: should a count fingerprint be generated
        :param include_disconnected: include disconnected atoms from hashes and substructure
        :param rdkit_invariants: use the atom invariants used by RDKit for its Morgan fingerprint
        :param exclude_floating: exclude atoms with no bonds
        :param remove_duplicate_substructs: drop duplicate hashes in the fingerprint
        :return: a dictionary if only one molecule is provided, otherwise a list of dicts
        """
        fper = Fingerprinter(bits=nbits, level=level,radius_multiplier=radius_multiplier,
                             stereo=stereo, counts=counts, include_disconnected=include_disconnected,
                             rdkit_invariants=rdkit_invariants, exclude_floating=exclude_floating,
                             remove_duplicate_substructs=remove_duplicate_substructs)
        # Convenience lambda function
        to_dense_fp = lambda bits, size: [1 if i in bits else 0 for i in range(size)]
        # If unique molecule transform to list
        if not isinstance(mol, list):
            mol = [mol]
        # Calculate fingerprint
        fps = []
        for mol_ in mol:
            # Ensure conformers are available, optimize otherwise
            confs = list(mol_.GetConformers())
            if not (len(confs) > 0 and confs[-1].Is3D()):
                mol_ = get_optimized_mol(inputmol=mol_)
            fper.run(mol=mol_)
            fp = to_dense_fp(fper.get_fingerprint_at_level().indices, nbits)
            values = dict(zip([f'E3FP_{i + 1}' for i in range(nbits)],
                              fp))
            fps.append(values)
        # Return only one value if input was not a list
        if len(mol) == 1:
            return fps[0]
        return fps

    _fp_funcs = {'E3FP': (calculate_e3fp, ('nbits')),
                 }

    @staticmethod
    def get_all_fps(mol: Chem.Mol, nbits: Optional[int] = None) -> dict:
        """Calculate all fingerprints."""
        values = {}
        for des_label, (func, supported_args) in Fingerprint3D._fp_funcs.items():
            if nbits is not None and 'nbits' in supported_args:
                values.update(func(mol, nbits=nbits))
            else:
                values.update(func(mol))
        return values
