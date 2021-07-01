# -*- coding: utf-8 -*-

"""Unit tests for chemopy.chemopy2D."""

from unittest import TestCase

from rdkit import Chem

from chemopy import PyChem2D
from tests.constants import MOL_FILE


class TestPyChem2D(TestCase):
    """Tests for PyChem2D."""

    def setUp(self):
        """Create PyChem2D molecules for Docetaxel."""
        self.molecules = []
        self.num_heavy_atoms = 58
        self.formula = 'C43H53NO14'

    def test_ReadMolFromMolFile(self):
        """Test reading the sample MOL file."""
        molecule = PyChem2D()
        molecule.ReadMolFromMolFile(MOL_FILE)
        self.assertEqual(molecule.mol.GetNumHeavyAtoms(), self.num_heavy_atoms)
        self.assertEqual(Chem.rdMolDescriptors.CalcMolFormula(molecule.mol), self.formula)
        self.molecules.append(molecule)

    def test_ReadMolFromSmiles(self):
        """Test reading from sample SMILES."""
        molecule = PyChem2D()
        molecule.ReadMolFromSmile('CC1=C2[C@H](C(=O)[C@@]3([C@H](C[C@@H]4'
                                  '[C@]([C@H]3[C@@H]([C@@](C2(C)C)(C[C@@H]1OC(=O)[C@@H]([C@H]'
                                  '(C5=CC=CC=C5)NC(=O)OC(C)(C)C)O)O)OC(=O)C6=CC=CC=C6)(CO4)OC(=O)C)O)C)O')
        self.assertEqual(molecule.mol.GetNumHeavyAtoms(), self.num_heavy_atoms)
        self.assertEqual(Chem.rdMolDescriptors.CalcMolFormula(molecule.mol), self.formula)
        self.molecules.append(molecule)

    def test_ReadMolFromInchi(self):
        """Test reading from sample InChI."""
        molecule = PyChem2D()
        molecule.ReadMolFromInchi('InChI=1S/C43H53NO14/c1-22-26'
                                  '(55-37(51)32(48)30(24-15-11-9-12-16-24)44-38(52)58-39(3,4)5)'
                                  '20-43(53)35(56-36(50)25-17-13-10-14-18-25)33-41(8,34(49)31(47)'
                                  '29(22)40(43,6)7)27(46)19-28-42(33,21-54-28)57-23(2)45/h9-18,'
                                  '26-28,30-33,35,46-48,53H,19-21H2,1-8H3,(H,44,52)/t26-,27-,28+,'
                                  '30-,31+,32+,33-,35-,41+,42-,43+/m0/s1')
        self.assertEqual(molecule.mol.GetNumHeavyAtoms(), self.num_heavy_atoms)
        self.assertEqual(Chem.rdMolDescriptors.CalcMolFormula(molecule.mol), self.formula)
        self.molecules.append(molecule)

    def test_GetMolFromNCBI(self):
        """Test reading from sample NCBI ID."""
        molecule = PyChem2D()
        # Ensure TimeoutError is consistent
        passed, tryouts = False, 0
        while not passed and tryouts < 5:
            try:
                molecule.GetMolFromNCBI('CID148124')
                passed = True
            except TimeoutError:
                tryouts += 1
        self.assertLess(tryouts, 5, 'Too many Timeout exceptions. Verify your internet connection.')
        self.assertEqual(molecule.mol.GetNumHeavyAtoms(), self.num_heavy_atoms)
        self.assertEqual(Chem.rdMolDescriptors.CalcMolFormula(molecule.mol), self.formula)
        self.molecules.append(molecule)

    def test_ChEBI_GetMolFromEBI(self):
        """Test reading from sample CHEBI ID."""
        molecule = PyChem2D()
        # Ensure TimeoutError is consistent
        passed, tryouts = False, 0
        while not passed and tryouts < 5:
            try:
                molecule.GetMolFromNCBI('CHEBI:4672')
                passed = True
            except TimeoutError:
                tryouts += 1
        self.assertLess(tryouts, 5, 'Too many Timeout exceptions. Verify your internet connection.')
        self.assertEqual(molecule.mol.GetNumHeavyAtoms(), self.num_heavy_atoms)
        self.assertEqual(Chem.rdMolDescriptors.CalcMolFormula(molecule.mol), self.formula)
        self.molecules.append(molecule)

    def test_ChEMBL_GetMolFromEBI(self):
        """Test reading from sample ChEMBL ID."""
        molecule = PyChem2D()
        # Ensure TimeoutError is consistent
        passed, tryouts = False, 0
        while not passed and tryouts < 5:
            try:
                molecule.GetMolFromNCBI('CHEMBL92')
                passed = True
            except TimeoutError:
                tryouts += 1
        self.assertLess(tryouts, 5, 'Too many Timeout exceptions. Verify your internet connection.')
        self.assertEqual(molecule.mol.GetNumHeavyAtoms(), self.num_heavy_atoms)
        self.assertEqual(Chem.rdMolDescriptors.CalcMolFormula(molecule.mol), self.formula)
        self.molecules.append(molecule)

    def test_GetMolFromCAS(self):
        """Test reading from sample CAS ID."""
        molecule = PyChem2D()
        # Ensure TimeoutError is consistent
        passed, tryouts = False, 0
        while not passed and tryouts < 5:
            try:
                molecule.GetMolFromCAS('114977-28-5')
                passed = True
            except TimeoutError:
                tryouts += 1
        self.assertLess(tryouts, 5, 'Too many Timeout exceptions. Verify your internet connection.')
        self.assertEqual(molecule.mol.GetNumHeavyAtoms(), self.num_heavy_atoms)
        self.assertEqual(Chem.rdMolDescriptors.CalcMolFormula(molecule.mol), self.formula)
        self.molecules.append(molecule)

    def test_GetMolFromKegg(self):
        """Test reading from sample KEGG ID."""
        molecule = PyChem2D()
        # Ensure TimeoutError is consistent
        passed, tryouts = False, 0
        while not passed and tryouts < 5:
            try:
                molecule.GetMolFromCAS('C11231')
                passed = True
            except TimeoutError:
                tryouts += 1
        self.assertLess(tryouts, 5, 'Too many Timeout exceptions. Verify your internet connection.')
        self.assertEqual(molecule.mol.GetNumHeavyAtoms(), self.num_heavy_atoms)
        self.assertEqual(Chem.rdMolDescriptors.CalcMolFormula(molecule.mol), self.formula)
        self.molecules.append(molecule)

    def test_GetMolFromDrugbank(self):
        """Test reading from sample DrugBank ID."""
        molecule = PyChem2D()
        # Ensure TimeoutError is consistent
        passed, tryouts = False, 0
        while not passed and tryouts < 5:
            try:
                molecule.GetMolFromCAS('DB01248')
                passed = True
            except TimeoutError:
                tryouts += 1
        self.assertLess(tryouts, 5, 'Too many Timeout exceptions. Verify your internet connection.')
        self.assertEqual(molecule.mol.GetNumHeavyAtoms(), self.num_heavy_atoms)
        self.assertEqual(Chem.rdMolDescriptors.CalcMolFormula(molecule.mol), self.formula)
        self.molecules.append(molecule)

    def test_molecules_equality(self):
        """Test if all sample molecules are the same."""
        equalities = []
        for i in range(len(self.molecules)):
            for j in range(i + 1, len(self.molecules)):
                mol1, mol2 = self.molecules[i], self.molecules[j]
                equalities.append(mol1.HasSubstrctureMatch(mol2) and mol2.HasSubstrctureMatch(mol1))
        self.assertTrue(all(equalities))

    def test_GetKappa(self):
        """Test GetKappa."""
        for i in range(len(self.molecules)):
            kappas_k, kappas_v = list(zip(*sorted(self.molecules[i].GetKappa().items())))
            self.assertEqual(len(kappas_k), 7)
            self.assertListEqual(kappas_k, ['kappa1', 'kappa2', 'kappa3', 'kappam1', 'kappam2', 'kappam3', 'phi'])
            for i, value in enumerate([47.478, 17.875, 8.8, 42.923, 15.172, 7.241, 11.228]):
                self.assertAlmostEqual(kappas_v[i], value, 3)

    def test_GetCharge(self):
        """Test GetCharge."""
        for i in range(len(self.molecules)):
            charge_k, charge_v = list(zip(*sorted(self.molecules[i].GetCharge().items())))
            self.assertEqual(len(charge_k), 25)
            self.assertListEqual(charge_k, ['LDI', 'Mac', 'Mnc', 'Mpc', 'QCmax', 'QCmin', 'QCss', 'QHmax', 'QHmin',
                                            'QHss', 'QNmax', 'QNmin', 'QNss', 'QOmax', 'QOmin', 'QOss', 'Qass',
                                            'Qmax', 'Qmin', 'Rnc', 'Rpc', 'SPP', 'Tac', 'Tnc', 'Tpc'])
            for i, value in enumerate([0.249, 0.112, -0.178, 0.082, 0.408, -0.062, 0.747, 0.212, 0.024,
                                       0.312, -0.311, -0.311, 0.097, -0.225, -0.456, 1.869, 3.025, 0.408,
                                       -0.456, 0.073, 0.066, 0.864, 12.442, -6.221, 6.221]):
                self.assertAlmostEqual(charge_v[i], value, 3)

    def test_GetConnectivity(self):
        """Test GetConnectivity."""
        for i in range(len(self.molecules)):
            connectivity_k, connectivity_v = list(zip(*sorted(self.molecules[i].GetConnectivity().items())))
            self.assertEqual(len(connectivity_k), 44)
            self.assertListEqual(connectivity_k, ['Chi0', 'Chi1', 'Chi10', 'Chi2', 'Chi3', 'Chi3c', 'Chi3ch', 'Chi4',
                                                  'Chi4c', 'Chi4ch', 'Chi4pc', 'Chi5', 'Chi5ch', 'Chi6', 'Chi6ch',
                                                  'Chi7', 'Chi8', 'Chi9', 'Chiv0', 'Chiv1', 'Chiv10', 'Chiv2', 'Chiv3',
                                                  'Chiv3c', 'Chiv3ch', 'Chiv4', 'Chiv4c', 'Chiv4ch', 'Chiv4pc', 'Chiv5',
                                                  'Chiv5ch', 'Chiv6', 'Chiv6ch', 'Chiv7', 'Chiv8', 'Chiv9', 'dchi0',
                                                  'dchi1', 'dchi2', 'dchi3', 'dchi4', 'knotp', 'knotpv', 'mChi1'])
            for i, value in enumerate([42.75, 26.99, 3.302, 27.641, 22.626, 7.33, 0.0, 19.181, 0.78, 0.144, 15.544,
                                       15.825, 0.0, 11.818, 0.272, 9.232, 6.455, 4.577, 33.885, 19.335, 0.975,
                                       17.574, 12.567, 4.392, 0.0, 9.524, 0.506, 0.083, 7.921, 7.026, 0.0, 5.007,
                                       0.124, 3.436, 2.224, 1.456, 8.865, 7.656, 10.067, 10.059, 9.656, 8.213,
                                       3.529, 0.428]):
                self.assertAlmostEqual(connectivity_v[i], value, 3)

    def test_GetConstitution(self):
        """Test GetConstitution."""
        for i in range(len(self.molecules)):
            constitution_k, constitution_v = list(zip(*sorted(self.molecules[i].GetConstitution().items())))
            self.assertEqual(len(constitution_k), 44)
            self.assertListEqual(constitution_k, ['AWeight', 'PathL1', 'PathL2', 'PathL3', 'PathL4', 'PathL5',
                                                  'PathL6', 'Weight', 'nAroBond', 'nAtom', 'nBr', 'nC', 'nCl',
                                                  'nDBond', 'nF', 'nH', 'nHA', 'nHBA', 'nHBD', 'nHal', 'nHet',
                                                  'nI', 'nN', 'nO', 'nP', 'nRing', 'nRotB', 'nS', 'nSBond',
                                                  'nTBond'])
            for i, value in enumerate([13.008, 63, 100, 140, 196, 272, 348, 754.466, 12, 111, 0, 43, 0, 6, 0,
                                       53, 58, 14, 5, 0, 43, 0, 1, 14, 0, 6, 8, 14, 45, 0]):
                self.assertAlmostEqual(constitution_v[i], value, 3)

    def test_GetEstateDescriptors(self):
        """Test GetEstateDescriptors."""
        for i in range(len(self.molecules)):
            estate_k, estate_v = list(zip(*sorted(self.molecules[i].GetEstateDescriptors().items())))
            self.assertEqual(len(estate_k), 8)
            self.assertListEqual(estate_k, ['DS', 'Save', 'Scar', 'Shal', 'Shet', 'Shev', 'Smax', 'Smin'])
            for i, value in enumerate([17.299, 2.611, 1.436586, 0.0, 149.98, 151.417, 14.945, -2.353]):
                self.assertAlmostEqual(estate_v[i], value, 3)

    def test_GetGeary(self):
        """Test GetGeary."""
        for i in range(len(self.molecules)):
            geary_k, geary_v = list(zip(*sorted(self.molecules[i].GetGeary().items())))
            self.assertEqual(len(geary_k), 32)
            self.assertListEqual(geary_k, ['GATSe1', 'GATSe2', 'GATSe3', 'GATSe4', 'GATSe5', 'GATSe6',
                                           'GATSe7', 'GATSe8', 'GATSm1', 'GATSm2', 'GATSm3', 'GATSm4',
                                           'GATSm5', 'GATSm6', 'GATSm7', 'GATSm8', 'GATSp1', 'GATSp2',
                                           'GATSp3', 'GATSp4', 'GATSp5', 'GATSp6', 'GATSp7', 'GATSp8',
                                           'GATSv1', 'GATSv2', 'GATSv3', 'GATSv4', 'GATSv5', 'GATSv6',
                                           'GATSv7', 'GATSv8'])
            for i, value in enumerate([0.83, 0.822, 0.936, 1.251, 1.174, 1.04, 1.008, 0.954, 0.83, 0.821,
                                       0.937, 1.251, 1.174, 1.04, 1.008, 0.954, 0.839, 0.815, 0.942, 1.251,
                                       1.167, 1.035, 1.008, 0.952, 0.836, 0.817, 0.94, 1.251, 1.17, 1.038,
                                       1.008, 0.953]):
                self.assertAlmostEqual(geary_v[i], value, 3)

    def test_GetMOE(self):
        """Test GetMOE."""
        for i in range(len(self.molecules)):
            moe_k, moe_v = list(zip(*sorted(self.molecules[i].GetMOE().items())))
            self.assertEqual(len(moe_k), 59)
            self.assertListEqual(moe_k, ['EstateVSA0', 'EstateVSA1', 'EstateVSA10', 'EstateVSA2', 'EstateVSA3',
                                         'EstateVSA4', 'EstateVSA5', 'EstateVSA6', 'EstateVSA7', 'EstateVSA8',
                                         'EstateVSA9', 'LabuteASA', 'MRVSA0', 'MRVSA1', 'MRVSA2', 'MRVSA3',
                                         'MRVSA4', 'MRVSA5', 'MRVSA6', 'MRVSA7', 'MRVSA8', 'MRVSA9', 'PEOEVSA0',
                                         'PEOEVSA1', 'PEOEVSA10', 'PEOEVSA11', 'PEOEVSA12', 'PEOEVSA13',
                                         'PEOEVSA2', 'PEOEVSA3', 'PEOEVSA4', 'PEOEVSA5', 'PEOEVSA6', 'PEOEVSA7',
                                         'PEOEVSA8', 'PEOEVSA9', 'TPSA1', 'VSAEstate0', 'VSAEstate1',
                                         'VSAEstate2', 'VSAEstate3', 'VSAEstate4', 'VSAEstate5', 'VSAEstate6',
                                         'VSAEstate7', 'VSAEstate8', 'VSAEstate9', 'slogPVSA0', 'slogPVSA1',
                                         'slogPVSA10', 'slogPVSA11', 'slogPVSA2', 'slogPVSA3', 'slogPVSA4',
                                         'slogPVSA5', 'slogPVSA6', 'slogPVSA7', 'slogPVSA8', 'slogPVSA9'])
            for i, value in enumerate([112.422, 29.737, 0.0, 5.563, 6.924, 25.98, 83.15, 0.0, 5.317, 23.684,
                                       44.399, 336.191, 68.083, 0.0, 5.317, 16.748, 127.7, 6.607, 82.937,
                                       0.0, 0.0, 29.784, 49.427, 9.589, 17.488, 0.0, 0.0, 24.001, 14.384,
                                       0.0, 0.0, 62.378, 63.461, 25.18, 35.649, 35.618, 224.45, 29.508, 68.91,
                                       51.562, -8.136, -5.363, 14.627, -11.482, 11.79, 0.0, 0.0, 5.317,
                                       110.244, 0.0, 0.0, 38.068, 16.748, 90.195, 71.81, 0.0, 0.0, 0.0, 4.795]):
                self.assertAlmostEqual(moe_v[i], value, 3)

    def test_GetMolProperties(self):
        """Test GetMolPropertie."""
        for i in range(len(self.molecules)):
            molprop_k, molprop_v = list(zip(*sorted(self.molecules[i].GetMolProperties().items())))
            self.assertEqual(len(molprop_k), 10)
            self.assertListEqual(molprop_k, ['Hy', 'LogP', 'LogP2', 'MR', 'TPSA', 'UI', 'XLogP', 'XLogP2'])
            for i, value in enumerate([-1.065, 3.26, 10.625, 203.604, 224.45, 4.248, 2.81, 7.89610, -5.07, -5.07]):
                self.assertAlmostEqual(molprop_v[i], value, 3)

    def test_GetMoran(self):
        """Test GetMoran."""
        for i in range(len(self.molecules)):
            moran_k, moran_v = list(zip(*sorted(self.molecules[i].GetMoran().items())))
            self.assertEqual(len(moran_k), 32)
            self.assertListEqual(moran_k, ['MATSe0', 'MATSe1', 'MATSe2', 'MATSe3', 'MATSe4', 'MATSe5', 'MATSe6',
                                           'MATSe7', 'MATSm0', 'MATSm1', 'MATSm2', 'MATSm3', 'MATSm4', 'MATSm5',
                                           'MATSm6', 'MATSm7', 'MATSp0', 'MATSp1', 'MATSp2', 'MATSp3', 'MATSp4',
                                           'MATSp5', 'MATSp6', 'MATSp7', 'MATSv0', 'MATSv1', 'MATSv2', 'MATSv3',
                                           'MATSv4', 'MATSv5', 'MATSv6', 'MATSv7'])
            for i, value in enumerate([-0.092, 0.048, -0.025, -0.216, -0.009, 0.053, 0.039, 0.018, -0.092, 0.049,
                                       -0.025, -0.215, -0.009, 0.052, 0.038, 0.018, -0.096, 0.06, -0.026, -0.212,
                                       -0.014, 0.047, 0.032, 0.018, -0.095, 0.056, -0.026, -0.213, -0.013, 0.049,
                                       0.034, 0.018]):
                self.assertAlmostEqual(moran_v[i], value, 3)

    def test_GetMoreauBroto(self):
        """Test GetMoreauBroto."""
        for i in range(len(self.molecules)):
            moreaubroto_k, moreaubroto_v = list(zip(*sorted(self.molecules[i].GetMoreauBroto().items())))
            self.assertEqual(len(moreaubroto_k), 32)
            self.assertListEqual(moreaubroto_k, ['ATSe0', 'ATSe1', 'ATSe2', 'ATSe3', 'ATSe4', 'ATSe5', 'ATSe6',
                                                 'ATSe7', 'ATSm0', 'ATSm1', 'ATSm2', 'ATSm3', 'ATSm4', 'ATSm5',
                                                 'ATSm6', 'ATSm7', 'ATSp0', 'ATSp1', 'ATSp2', 'ATSp3', 'ATSp4',
                                                 'ATSp5', 'ATSp6', 'ATSp7', 'ATSv0', 'ATSv1', 'ATSv2', 'ATSv3',
                                                 'ATSv4', 'ATSv5', 'ATSv6', 'ATSv7'])
            for i, value in enumerate([4.256, 4.727, 4.945, 5.123, 5.224, 5.234, 5.193, 5.066, 4.258, 4.729, 4.947,
                                       5.126, 5.226, 5.237, 5.195, 5.068, 3.968, 4.354, 4.537, 4.611, 4.655, 4.722,
                                       4.708, 4.627, 3.991, 4.382, 4.569, 4.654, 4.699, 4.76, 4.744, 4.66]):
                self.assertAlmostEqual(moreaubroto_v[i], value, 3)

    def test_GetTopology(self):
        """Test GetTopology."""
        for i in range(len(self.molecules)):
            topology_k, topology_v = list(zip(*sorted(self.molecules[i].GetTopology().items())))
            self.assertEqual(len(topology_k), 35)
            self.assertListEqual(topology_k, ['AW', 'Arto', 'BertzCT', 'DZ', 'GMTI', 'GMTIV', 'Geto', 'Getov',
                                              'Gravto', 'Hato', 'Hatov', 'IDE', 'IDET', 'ISIZ', 'IVDE', 'Ipc',
                                              'J', 'MZM1', 'MZM2', 'Platt', 'Pol', 'Qindex', 'Sito', 'Sitov',
                                              'TIAC', 'Thara', 'Tigdi', 'Tsch', 'W', 'Xu', 'ZM1', 'ZM2',
                                              'diametert', 'petitjeant', 'radiust'])
            for i, value in enumerate([7.615, 2.172, 3.296, 130.5, 4.739, 5.192, 1.951, 3.175, 213.804, 1.736,
                                       2.671, 3.989, 6169.42, 754.18, 108.469, 12.166, 1.415, 23.951, 2.782,
                                       200, 121.0, 50.0, 38.778, 67.004, 163.968, 337.788, 5.279, 52970.0,
                                       12587.0, 46.303, 326, 403, 18.0, 1.0, 9.0]):
                self.assertAlmostEqual(topology_v[i], value, 3)

    def test_GetBcut(self):
        """Test GetBcut."""
        for i in range(len(self.molecules)):
            bcut_k, bcut_v = list(zip(*sorted(self.molecules[i].GetBcut().items())))
            self.assertEqual(len(bcut_k), 64)
            self.assertListEqual(bcut_k, ['bcute1', 'bcute10', 'bcute11', 'bcute12', 'bcute13', 'bcute14',
                                          'bcute15', 'bcute16', 'bcute2', 'bcute3', 'bcute4', 'bcute5',
                                          'bcute6', 'bcute7', 'bcute8', 'bcute9', 'bcutm1', 'bcutm10',
                                          'bcutm11', 'bcutm12', 'bcutm13', 'bcutm14', 'bcutm15', 'bcutm16',
                                          'bcutm2', 'bcutm3', 'bcutm4', 'bcutm5', 'bcutm6', 'bcutm7',
                                          'bcutm8', 'bcutm9', 'bcutp1', 'bcutp10', 'bcutp11', 'bcutp12',
                                          'bcutp13', 'bcutp14', 'bcutp15', 'bcutp16', 'bcutp2', 'bcutp3',
                                          'bcutp4', 'bcutp5', 'bcutp6', 'bcutp7', 'bcutp8', 'bcutp9',
                                          'bcutv1', 'bcutv10', 'bcutv11', 'bcutv12', 'bcutv13', 'bcutv14',
                                          'bcutv15', 'bcutv16', 'bcutv2', 'bcutv3', 'bcutv4', 'bcutv5',
                                          'bcutv6', 'bcutv7', 'bcutv8', 'bcutv9'])
            for i, value in enumerate([4.075, 1.87, 1.839, 1.815, 1.591, 1.463, 1.446, 1.351, 3.892, 3.85, 3.847,
                                       3.649, 3.554, 3.523, 3.421, 1.984, 4.025, 1.975, 1.943, 1.938, 1.814, 1.622,
                                       1.608, 1.56, 3.836, 3.796, 3.793, 3.567, 3.496, 3.455, 3.355, 2.078, 3.976,
                                       1.97, 1.944, 1.918, 1.793, 1.681, 1.656, 1.558, 3.822, 3.789, 3.744, 3.472,
                                       3.37, 3.343, 3.245, 2.102, 3.975, 1.974, 1.95, 1.926, 1.804, 1.682, 1.661,
                                       1.564, 3.818, 3.785, 3.742, 3.468, 3.371, 3.344, 3.244, 2.106]):
                self.assertAlmostEqual(bcut_v[i], value, 3)

    def test_GetBasak(self):
        """Test GetBasak."""
        for i in range(len(self.molecules)):
            basak_k, basak_v = list(zip(*sorted(self.molecules[i].GetBasak().items())))
            self.assertEqual(len(basak_k), 21)
            self.assertListEqual(basak_k, ['CIC0', 'CIC1', 'CIC2', 'CIC3', 'CIC4', 'CIC5', 'CIC6', 'IC0', 'IC1',
                                           'IC2', 'IC3', 'IC4', 'IC5', 'IC6', 'SIC0', 'SIC1', 'SIC2', 'SIC3',
                                           'SIC4', 'SIC5', 'SIC6'])
            for i, value in enumerate([5.317, 3.593, 1.944, 1.059, 0.802, 0.73, 0.73, 1.477, 3.202, 4.85, 5.735,
                                       5.992, 6.064, 6.064, 0.217, 0.471, 0.714, 0.844, 0.882, 0.893, 0.893]):
                self.assertAlmostEqual(basak_v[i], value, 3)

    def test_GetAllFingerprints(self):
        """Test GetAllFingerprints."""
        pass

    def test_GetFingerprint(self):
        """Test GetFingerprint."""
        pass

    def test_GetAllDescriptors(self):
        """Test GetAllDescriptors."""
        pass
