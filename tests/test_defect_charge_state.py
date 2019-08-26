import unittest

from py_sc_fermi.defect_charge_state import DefectChargeState

class TestDefectChargeStateInit(unittest.TestCase):

    def test_defect_charge_state_is_initialised(self):
        charge = 1.0
        energy = 123.4
        deg = 2
        defect_charge_state = DefectChargeState(charge=charge, energy=energy, deg=deg)
        self.assertEqual( defect_charge_state.charge, charge )
        self.assertEqual( defect_charge_state.energy, energy )
        self.assertEqual( defect_charge_state.deg, deg )
        self.assertEqual( defect_charge_state.fixed_concentration, False )

class test_defect_charge_state(unittest.TestCase):

    def setUp(self):
        charge = 1.0
        energy = 0.1234
        deg = 2
        self.defect_charge_state = DefectChargeState(charge=charge, energy=energy, deg=deg)
   
    def test_get_formation_energy(self):
        e_fermi = 1.2
        formation_energy = self.defect_charge_state.get_formation_energy( e_fermi )
        self.assertEqual( formation_energy, 0.1234 + ( 1.0*1.2 ) )

    def test_get_concentration(self):
        e_fermi = 1.2
        temperature = 298.0 
        conc = self.defect_charge_state.get_concentration( e_fermi=e_fermi,
            temperature=temperature )
        self.assertEqual( conc, 8.311985602942568e-23 )
 
if __name__ == '__main__':
    unittest.main()
