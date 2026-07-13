# -*- coding: utf-8 -*-

"""Tests for individual descriptor methods that get_all()/get_all_descriptors() no longer
call internally (they were optimized to avoid redundant recomputation), but which remain
public API and must still work correctly when called directly."""

import unittest

from rdkit import Chem

from chemopy.basak import Basak
from chemopy.charge import Charge
from chemopy.estate import EState

_MOL = Chem.AddHs(Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O"))


class TestBasakIndividualMethods(unittest.TestCase):
    def test_ic_sic_cic_0_through_6(self):
        for order in range(7):
            ic = getattr(Basak, f"calculate_basak_ic{order}")(_MOL)
            sic = getattr(Basak, f"calculate_basak_sic{order}")(_MOL)
            cic = getattr(Basak, f"calculate_basak_cic{order}")(_MOL)
            self.assertIsInstance(ic, float)
            self.assertIsInstance(sic, float)
            self.assertIsInstance(cic, float)

    def test_individual_matches_get_all(self):
        all_values = Basak.get_all(_MOL)
        for order in range(7):
            ic = getattr(Basak, f"calculate_basak_ic{order}")(_MOL)
            self.assertAlmostEqual(ic, all_values[f"IC{order}"], places=8)


class TestChargeIndividualMethods(unittest.TestCase):
    def test_all_individual_charge_methods(self):
        methods = [
            "calculate_h_max_pcharge",
            "calculate_c_max_pcharge",
            "calculate_n_max_pcharge",
            "calculate_o_max_pcharge",
            "calculate_h_max_ncharge",
            "calculate_c_max_ncharge",
            "calculate_n_max_ncharge",
            "calculate_o_max_ncharge",
            "calculate_all_max_pcharge",
            "calculate_all_max_ncharge",
            "calculate_h_sum_square_charge",
            "calculate_c_sum_square_charge",
            "calculate_n_sum_square_charge",
            "calculate_o_sum_square_charge",
            "calculate_all_sum_square_charge",
            "calculate_mean_pcharge",
            "calculate_total_pcharge",
            "calculate_mean_ncharge",
            "calculate_total_ncharge",
            "calculate_mean_absolute_charge",
            "calculate_total_absolute_charge",
            "calculate_relative_ncharge",
            "calculate_relative_pcharge",
            "calculate_submol_polarity_param",
            "calculate_local_dipole_index",
        ]
        for method in methods:
            value = getattr(Charge, method)(_MOL)
            self.assertIsInstance(value, (int, float))

    def test_individual_matches_get_all(self):
        all_values = Charge.get_all(_MOL)
        self.assertAlmostEqual(Charge.calculate_c_max_pcharge(_MOL), all_values["QCmax"], places=6)
        self.assertAlmostEqual(Charge.calculate_all_max_pcharge(_MOL), all_values["Qmax"], places=6)
        self.assertAlmostEqual(Charge.calculate_local_dipole_index(_MOL), all_values["LDI"], places=6)


class TestEStateIndividualMethods(unittest.TestCase):
    def test_all_individual_estate_methods(self):
        methods = [
            "calculate_heavy_atom_estate",
            "calculate_c_atom_estate",
            "calculate_halogen_estate",
            "calculate_hetero_estate",
            "calculate_average_estate",
            "calculate_max_estate",
            "calculate_min_estate",
            "calculate_diff_max_min_estate",
        ]
        for method in methods:
            value = getattr(EState, method)(_MOL)
            self.assertIsInstance(value, float)

    def test_individual_matches_get_all_descriptors(self):
        all_values = EState.get_all_descriptors(_MOL)
        self.assertAlmostEqual(EState.calculate_heavy_atom_estate(_MOL), all_values["Shev"], places=6)
        self.assertAlmostEqual(EState.calculate_c_atom_estate(_MOL), all_values["Scar"], places=6)

    def test_max_min_atom_type_estate_and_fingerprint(self):
        max_es = EState.calculate_max_atom_type_estate(_MOL)
        min_es = EState.calculate_min_atom_type_estate(_MOL)
        fp = EState.calculate_estate_fingerprint(_MOL)
        all_fps = EState.get_all_fps(_MOL)
        for k, v in max_es.items():
            self.assertAlmostEqual(all_fps[k], v, places=6)
        for k, v in min_es.items():
            self.assertAlmostEqual(all_fps[k], v, places=6)
        for k, v in fp.items():
            self.assertAlmostEqual(all_fps[k], v, places=6)


if __name__ == "__main__":
    unittest.main()
