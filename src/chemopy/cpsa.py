# -*- coding: utf-8 -*-


"""Charged partial surface area (CPSA) descriptors."""

from rdkit import Chem

from . import asa
from .geo_opt import get_atom_class_list, read_coordinates


class CPSA:
    """Charge partial surface area descriptors."""

    @staticmethod
    def get_charge_sa(arc_file, radius_probe=1.5, n_sphere_point=960):
        """Get atom symbol, charge and partial solvent-accessible surface areas for all atoms.

        :param arc_file: Path to MOPAC .arc file
        :param RadiusProbe: radius of the probe used to calculate SASA
        :param n_sphere_point: number of points per atom to calculate SASA
        """
        charge_coordinates = read_coordinates(arc_file)
        atoms = get_atom_class_list(charge_coordinates)
        fasa = asa.calculate_asa(atoms, radius_probe, n_sphere_point)
        res = []
        for i in range(len(fasa)):
            res.append([charge_coordinates[i][0], charge_coordinates[i][4], fasa[i]])
        return res
    
    @staticmethod
    def calculate_asa(charge_sa):
        """calculate solvent-accessible surface area."""
        res = 0.0
        for i in charge_sa:
            res = res + i[2]
        return res

    @staticmethod
    def calculate_msa(arc_file):
        """calculate molecular surface areas.

        :param arc_file: Path to MOPAC .arc file
        """
        charge_sa = CPSA.get_charge_sa(arc_file, radius_probe =0, n_sphere_point=960)
        res = 0.0
        for i in charge_sa:
            res = res + i[2]
        return res

    @staticmethod
    def calculate_pnsa1(charge_sa):
        """Calculate partial negative area."""
        res = 0.0
        for i in charge_sa:
            if float(i[1]) < 0:
                res = res + i[2]
        return res

    @staticmethod
    def calculate_ppsa1(charge_sa):
        """Calculate partial positive area."""
        res = 0.0
        for i in charge_sa:
            if float(i[1]) > 0:
                res = res + i[2]
        return res

    @staticmethod
    def calculate_pnsa2(charge_sa):
        """Calculate total charge weighted negative surface area."""
        temp1, temp2 = 0.0, 0.0
        for i in charge_sa:
            if float(i[1]) < 0:
                temp1 += float(i[1])
                temp2 += i[2]
        res = temp1 * temp2
        return res

    @staticmethod
    def calculate_ppsa2(charge_sa):
        """Calculate total charge weighted positive surface area."""
        temp1, temp2 = 0.0, 0.0
        for i in charge_sa:
            if float(i[1]) > 0:
                temp1 += float(i[1])
                temp2 += i[2]
        res = temp1 * temp2
        return res

    @staticmethod
    def calculate_pnsa3(charge_sa):
        """Calculate atom charge weighted negative surface area."""
        res = 0.0
        for i in charge_sa:
            if float(i[1]) < 0:
                res += float(i[1]) * i[2]
        return res

    @staticmethod
    def calculate_ppsa3(charge_sa):
        """Calculate atom charge weighted positive surface area."""
        res = 0.0
        for i in charge_sa:
            if float(i[1]) > 0:
                res += float(i[1]) * i[2]
        return res

    @staticmethod
    def calculate_dpsa1(charge_sa):
        """Calculate difference in charged partial surface area."""
        return CPSA.calculate_ppsa1(charge_sa) - CPSA.calculate_pnsa1(charge_sa)

    @staticmethod
    def calculate_dpsa2(charge_sa):
        """Calculate difference in total charge weighted partial surface area."""
        return CPSA.calculate_ppsa2(charge_sa) - CPSA.calculate_pnsa2(charge_sa)

    @staticmethod
    def calculate_dpsa3(charge_sa):
        """Calculate difference in atomic charge weighted surface area."""
        return CPSA.calculate_ppsa3(charge_sa) - CPSA.calculate_pnsa3(charge_sa)

    @staticmethod
    def calculate_fnsa1(charge_sa):
        """Calculate fractional charged partial negative surface area."""
        temp = 0.0
        for i in charge_sa:
            temp += i[2]
        return CPSA.calculate_pnsa1(charge_sa) / temp

    @staticmethod
    def calculate_fnsa2(charge_sa):
        """Calculate fractional charged total negative surface area."""
        temp = 0.0
        for i in charge_sa:
            temp += i[2]
        return CPSA.calculate_pnsa2(charge_sa) / temp

    @staticmethod
    def calculate_fnsa3(charge_sa):
        """Calculate fractional charged atom negative surface area."""
        temp = 0.0
        for i in charge_sa:
            temp += i[2]
        return CPSA.calculate_pnsa3(charge_sa) / temp

    @staticmethod
    def calculate_fpsa1(charge_sa):
        """Calculate fractional charged partial positive surface area."""
        temp = 0.0
        for i in charge_sa:
            temp += i[2]
        return CPSA.calculate_ppsa1(charge_sa) / temp

    @staticmethod
    def calculate_fpsa2(charge_sa):
        """Calculate fractional charged total positive surface area."""
        temp = 0.0
        for i in charge_sa:
            temp += i[2]
        return CPSA.calculate_ppsa2(charge_sa) / temp

    @staticmethod
    def calculate_fpsa3(charge_sa):
        """Calculate fractional charged atom positive surface area."""
        temp = 0.0
        for i in charge_sa:
            temp += i[2]
        return CPSA.calculate_ppsa3(charge_sa) / temp

    @staticmethod
    def calculate_wnsa1(charge_sa):
        """Calculate surface weighted charged partial negative surface area."""
        temp = 0.0
        for i in charge_sa:
            temp += i[2]
        return CPSA.calculate_pnsa1(charge_sa) * temp / 1000

    @staticmethod
    def calculate_wnsa2(charge_sa) -> float:
        """Calculate surface weighted charged total negative surface area."""
        temp = 0.0
        for i in charge_sa:
            temp += i[2]
        return CPSA.calculate_pnsa2(charge_sa) * temp / 1000

    @staticmethod
    def calculate_wnsa3(charge_sa) -> float:
        """Calculate surface weighted charged atom negative surface area."""
        temp = 0.0
        for i in charge_sa:
            temp += i[2]
        return CPSA.calculate_pnsa3(charge_sa) * temp / 1000

    @staticmethod
    def calculate_wpsa1(charge_sa) -> float:
        """Calculate surface weighted charged partial positive surface area."""
        temp = 0.0
        for i in charge_sa:
            temp += i[2]
        return CPSA.calculate_ppsa1(charge_sa) * temp / 1000

    @staticmethod
    def calculate_wpsa2(charge_sa) -> float:
        """Calculate surface weighted charged total positive surface area."""
        temp = 0.0
        for i in charge_sa:
            temp += i[2]
        return CPSA.calculate_ppsa2(charge_sa) * temp / 1000

    @staticmethod
    def calculate_wpsa3(charge_sa) -> float:
        """Calculate surface weighted charged atom positive surface area."""
        temp = 0.0
        for i in charge_sa:
            temp += i[2]
        return CPSA.calculate_ppsa3(charge_sa) * temp / 1000

    @staticmethod
    def calculate_tasa(charge_sa) -> float:
        """Calculate total apolar (hydrophobic) surface area."""
        res = 0.0
        for i in charge_sa:
            if abs(float(i[1])) < 0.2:
                res += i[2]
        return res

    @staticmethod
    def calculate_tpsa(charge_sa) -> float:
        """Calculate total polar surface area."""
        res = 0.0
        for i in charge_sa:
            if abs(float(i[1])) >= 0.2:
                res += i[2]
        return res

    @staticmethod
    def calculate_ratiotatp(charge_sa) -> float:
        """Calculate ratio between TASA and TPSA (FrTATP)."""
        res = 0.0
        if CPSA.calculate_tpsa(charge_sa) == 0:
            return res
        else:
            return CPSA.calculate_tasa(charge_sa) / CPSA.calculate_tpsa(charge_sa)

    @staticmethod
    def calculate_rasa(charge_sa) -> float:
        """Calculate relative hydrophobic surface area."""
        temp = 0.0
        for i in charge_sa:
            temp += i[2]
        return CPSA.calculate_tasa(charge_sa) / temp

    @staticmethod
    def calculate_rpsa(charge_sa) -> float:
        """Calculate relative polar surface area."""
        temp = 0.0
        for i in charge_sa:
            temp += i[2]
        return CPSA.calculate_tpsa(charge_sa) / temp

    @staticmethod
    def calculate_rncs(charge_sa) -> float:
        """Calculate relative negative charge surface area."""
        charge = []
        for i in charge_sa:
            charge.append(float(i[1]))
        temp = []
        for i in charge_sa:
            temp.append(i[2])
        RNCG = min(charge) / sum(i for i in charge if i < 0)
        return temp[charge.index(min(charge))] / RNCG

    @staticmethod
    def calculate_rpcs(charge_sa) -> float:
        """Calculate relative positive charge surface area."""
        charge = []
        for i in charge_sa:
            charge.append(float(i[1]))
        temp = []
        for i in charge_sa:
            temp.append(i[2])
        RPCG = max(charge) / sum(i for i in charge if i > 0)
        return temp[charge.index(min(charge))] / RPCG

    @staticmethod
    def get_all(mol: Chem.Mol, arc_file: str) -> dict:
        """Get all CPSA descriptors.

        :param arc_file: Path to MOPAC .arc file
        """
        res = {}
        charge_sa = CPSA.get_charge_sa(arc_file, radius_probe=1.5, n_sphere_point=5000)
        res['SASA'] = CPSA.calculate_asa(charge_sa)
        res['MSA'] = CPSA.calculate_msa(arc_file)
        res['PNSA1'] = CPSA.calculate_pnsa1(charge_sa)
        res['PPSA1'] = CPSA.calculate_ppsa1(charge_sa)
        res['PNSA2'] = CPSA.calculate_pnsa2(charge_sa)
        res['PPSA2'] = CPSA.calculate_ppsa2(charge_sa)
        res['PNSA3'] = CPSA.calculate_pnsa3(charge_sa)
        res['PPSA3'] = CPSA.calculate_ppsa3(charge_sa)
        res['DPSA1'] = CPSA.calculate_dpsa1(charge_sa)
        res['DPSA2'] = CPSA.calculate_dpsa2(charge_sa)
        res['DPSA3'] = CPSA.calculate_dpsa3(charge_sa)
        res['FNSA1'] = CPSA.calculate_fnsa1(charge_sa)
        res['FNSA2'] = CPSA.calculate_fnsa2(charge_sa)
        res['FNSA3'] = CPSA.calculate_fnsa3(charge_sa)
        res['FPSA1'] = CPSA.calculate_fpsa1(charge_sa)
        res['FPSA2'] = CPSA.calculate_fpsa2(charge_sa)
        res['FPSA3'] = CPSA.calculate_fpsa3(charge_sa)
        res['WNSA1'] = CPSA.calculate_wnsa1(charge_sa)
        res['WNSA2'] = CPSA.calculate_wnsa2(charge_sa)
        res['WNSA3'] = CPSA.calculate_wnsa3(charge_sa)
        res['WPSA1'] = CPSA.calculate_wpsa1(charge_sa)
        res['WPSA2'] = CPSA.calculate_wpsa2(charge_sa)
        res['WPSA3'] = CPSA.calculate_wpsa3(charge_sa)
        res['TASA'] = CPSA.calculate_tasa(charge_sa)
        res['TotPSA'] = CPSA.calculate_tpsa(charge_sa)
        res['FrTATP'] = CPSA.calculate_ratiotatp(charge_sa)
        res['RASA'] = CPSA.calculate_rasa(charge_sa)
        res['RPSA'] = CPSA.calculate_rpsa(charge_sa)
        res['RNCS'] = CPSA.calculate_rncs(charge_sa)
        res['RPCS'] = CPSA.calculate_rpcs(charge_sa)
        return res
