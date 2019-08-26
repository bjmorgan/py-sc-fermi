import unittest

from py_sc_fermi.defect_charge_state import DefectChargeState

class TestDefectChargeState(unittest.TestCase):

    def test_defect_charge_state_is_initialised(self):
        charge = 1.0
        energy = 123.4
        deg = 2
        defect_charge_state = DefectChargeState(charge=charge, energy=energy, deg=deg)
        self.assertEqual( defect_charge_state.charge, charge )
        self.assertEqual( defect_charge_state.energy, energy )
        self.assertEqual( defect_charge_state.deg, deg )
        self.assertEqual( defect_charge_state.fixed_concentration, False )

if __name__ == '__main__':
    unittest.main()
