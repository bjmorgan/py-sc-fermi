import unittest
from py_sc_fermi.defect_charge_state import DefectChargeState


class TestDefectChargeStateInit(unittest.TestCase):
    def test_defect_charge_state_is_initialised(self):
        charge = 1.0
        energy = 123.4
        degeneracy = 2
        defect_charge_state = DefectChargeState(
            charge=charge, energy=energy, degeneracy=degeneracy
        )
        self.assertEqual(defect_charge_state._charge, charge)
        self.assertEqual(defect_charge_state._energy, energy)
        self.assertEqual(defect_charge_state._degeneracy, degeneracy)
        self.assertEqual(defect_charge_state.fixed_concentration, None)

    def test_bad_energy_and_concentration(self):
        with self.assertRaises(ValueError):
            DefectChargeState(1, None, None)


class TestDefectChargeState(unittest.TestCase):
    def setUp(self):
        charge = 1.0
        energy = 0.1234
        degeneracy = 2
        self.defect_charge_state = DefectChargeState(
            charge=charge, energy=energy, degeneracy=degeneracy
        )

    def test_charge_property(self):
        self.assertEqual(
            self.defect_charge_state.charge, self.defect_charge_state._charge
        )

    def test_energy_property(self):
        self.assertEqual(
            self.defect_charge_state.energy, self.defect_charge_state._energy
        )

    def test_degeneracy_property(self):
        self.assertEqual(
            self.defect_charge_state.degeneracy, self.defect_charge_state._degeneracy
        )

    def test_fix_concentration(self):
        self.assertEqual(self.defect_charge_state.fixed_concentration, None)
        self.defect_charge_state.fix_concentration(1)
        self.assertEqual(self.defect_charge_state.fixed_concentration, 1)

    def test_get_formation_energy(self):
        e_fermi = 1.2
        formation_energy = self.defect_charge_state.get_formation_energy(e_fermi)
        self.assertEqual(formation_energy, 0.1234 + (1.0 * 1.2))

    def test_get_formation_energy_raises(self):
        with self.assertRaises(ValueError):
            self.defect_charge_state._energy = None
            self.defect_charge_state.get_formation_energy(0.1234)

    def test_get_concentration(self):
        e_fermi = 1.2
        temperature = 298.0
        conc = self.defect_charge_state.get_concentration(
            e_fermi=e_fermi, temperature=temperature
        )
        self.assertEqual(conc, 8.311501552630706e-23)

    def test_get_concentration_with_fixed_concentration(self):
        e_fermi = 1.2
        temperature = 298.0
        self.defect_charge_state.fix_concentration(1.0)
        conc = self.defect_charge_state.get_concentration(
            e_fermi=e_fermi, temperature=temperature
        )
        self.assertEqual(conc, 1.0)

    def test_defect_charge_state_from_dict(self):
        dictionary = {"degeneracy": 2, "energy": 0.1234, "charge": 1}
        defect_charge_state = DefectChargeState.from_dict(dictionary)
        self.assertEqual(defect_charge_state.degeneracy, 2)
        self.assertEqual(defect_charge_state.energy, 0.1234)
        self.assertEqual(defect_charge_state.charge, 1)
        self.assertEqual(defect_charge_state.fixed_concentration, None)

    def test_defect_charge_state_from_dict_with_fixed_concentration(self):
        dictionary = {
            "degeneracy": 2,
            "energy": 0.1234,
            "charge": 1,
            "fixed_concentration": 0.1234,
        }
        defect_charge_state = DefectChargeState.from_dict(dictionary)
        self.assertEqual(defect_charge_state.degeneracy, 2)
        self.assertEqual(defect_charge_state.charge, 1)
        self.assertEqual(defect_charge_state.fixed_concentration, 0.1234)

    def test_defect_system_as_dict(self):
        dictionary = self.defect_charge_state.as_dict()
        self.assertEqual(dictionary["degeneracy"],2)
        self.assertEqual(dictionary["energy"],0.1234)
        self.assertEqual(dictionary["charge"],1)

    def test_defect_system_as_dict_fixed_concentration(self):
        self.defect_charge_state.fix_concentration(1)
        dictionary = self.defect_charge_state.as_dict()
        self.assertEqual(dictionary["degeneracy"],2)
        self.assertEqual(dictionary["energy"],0.1234)
        self.assertEqual(dictionary["charge"],1)
        self.assertEqual(dictionary["fixed_concentration"],1)

    def test_defect_charge_state_from_dict_warns(self):
        dictionary = {
            "degeneracy": 2,
            "energy": 0.1234,
            "charge": 1,
            "fixed_concentration": 0.1234,
            "foo": "bar"
        }
        with self.assertWarns(UserWarning):
            DefectChargeState.from_dict(dictionary)

    def test_defect_charge_state_from_string(self):
        string = "1 0.1234 2"
        defect_charge_state = DefectChargeState.from_string(string)
        self.assertEqual(defect_charge_state.degeneracy, 2)
        self.assertEqual(defect_charge_state.energy, 0.1234)
        self.assertEqual(defect_charge_state.charge, 1)
        self.assertEqual(defect_charge_state.fixed_concentration, None)

    def test_defect_charge_state_from_string_with_fixed_concentration(self):
        string = "V_O 1 0.1234"
        defect_charge_state = DefectChargeState.from_string(
            string, frozen=True, volume=100
        )
        self.assertEqual(defect_charge_state.fixed_concentration, 1.234e-23)
        self.assertEqual(defect_charge_state.charge, 1)

    def test_defect_charge_state_from_string_raises(self):
        string = "V_O 1 0.1234"
        with self.assertRaises(ValueError):
            DefectChargeState.from_string(string, frozen=True, volume=None)

    def test__repr__(self):
        self.assertEqual(
            str(self.defect_charge_state), "q=+1.0, e=0.1234, deg=2",
        )


if __name__ == "__main__":
    unittest.main()
