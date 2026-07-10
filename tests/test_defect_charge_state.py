import unittest

import numpy as np
import yaml

from py_sc_fermi.defect_charge_state import DefectChargeState
from py_sc_fermi.warnings import UnrecognisedKeyWarning


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

    def test_init_raises_error_on_none_energy_and_concentration(self):
        with self.assertRaises(ValueError):
            DefectChargeState(1, None, None)
            
    def test_init_raises_error_on_invalid_degeneracy(self):
        with self.assertRaises(ValueError):
            DefectChargeState(charge=1, energy=0.5, degeneracy=0)
        with self.assertRaises(ValueError):
            DefectChargeState(charge=1, energy=0.5, degeneracy=-1)


class TestDefectChargeStateChargeProperty(unittest.TestCase):
    def setUp(self):
        charge = 1
        energy = 0.1234
        degeneracy = 2
        self.defect_charge_state = DefectChargeState(
            charge=charge, energy=energy, degeneracy=degeneracy
        )

    def test_charge_property(self):
        self.assertEqual(
            self.defect_charge_state.charge, self.defect_charge_state._charge
        )


class TestDefectChargeStateEnergyProperty(unittest.TestCase):
    def setUp(self):
        charge = 1
        energy = 0.1234
        degeneracy = 2
        self.defect_charge_state = DefectChargeState(
            charge=charge, energy=energy, degeneracy=degeneracy
        )

    def test_energy_property(self):
        self.assertEqual(
            self.defect_charge_state.energy, self.defect_charge_state._energy
        )


class TestDefectChargeStateDegeneracyProperty(unittest.TestCase):
    def setUp(self):
        charge = 1
        energy = 0.1234
        degeneracy = 2
        self.defect_charge_state = DefectChargeState(
            charge=charge, energy=energy, degeneracy=degeneracy
        )

    def test_degeneracy_property(self):
        self.assertEqual(
            self.defect_charge_state.degeneracy, self.defect_charge_state._degeneracy
        )


class TestDefectChargeStateFixConcentration(unittest.TestCase):
    def setUp(self):
        charge = 1
        energy = 0.1234
        degeneracy = 2
        self.defect_charge_state = DefectChargeState(
            charge=charge, energy=energy, degeneracy=degeneracy
        )

    def test_fix_concentration(self):
        self.assertEqual(self.defect_charge_state.fixed_concentration, None)
        self.defect_charge_state.fix_concentration(1)
        self.assertEqual(self.defect_charge_state.fixed_concentration, 1)


class TestDefectChargeStateGetFormationEnergy(unittest.TestCase):
    def setUp(self):
        charge = 1
        energy = 0.1234
        degeneracy = 2
        self.defect_charge_state = DefectChargeState(
            charge=charge, energy=energy, degeneracy=degeneracy
        )

    def test_get_formation_energy(self):
        e_fermi = 1.2
        formation_energy = self.defect_charge_state.get_formation_energy(e_fermi)
        self.assertAlmostEqual(formation_energy, 0.1234 + (1.0 * 1.2), places=4)

    def test_get_formation_energy_raises(self):
        with self.assertRaises(ValueError):
            self.defect_charge_state._energy = None
            self.defect_charge_state.get_formation_energy(0.1234)


class TestDefectChargeStateGetConcentration(unittest.TestCase):
    def setUp(self):
        charge = 1
        energy = 0.1234
        degeneracy = 2
        self.defect_charge_state = DefectChargeState(
            charge=charge, energy=energy, degeneracy=degeneracy
        )

    def test_get_concentration(self):
        e_fermi = 1.2
        temperature = 298.0
        conc = self.defect_charge_state.get_concentration(
            e_fermi=e_fermi, temperature=temperature
        )
        self.assertAlmostEqual(conc, 8.311501552630706e-23, places=25)

    def test_get_concentration_with_fixed_concentration(self):
        e_fermi = 1.2
        temperature = 298.0
        self.defect_charge_state.fix_concentration(1.0)
        conc = self.defect_charge_state.get_concentration(
            e_fermi=e_fermi, temperature=temperature
        )
        self.assertEqual(conc, 1.0)


class TestDefectChargeStateDictionaryOperations(unittest.TestCase):
    def setUp(self):
        charge = 1
        energy = 0.1234
        degeneracy = 2
        self.defect_charge_state = DefectChargeState(
            charge=charge, energy=energy, degeneracy=degeneracy
        )

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
        self.assertEqual(dictionary["degeneracy"], 2)
        self.assertEqual(dictionary["energy"], 0.1234)
        self.assertEqual(dictionary["charge"], 1)

    def test_defect_system_as_dict_fixed_concentration(self):
        self.defect_charge_state.fix_concentration(1)
        dictionary = self.defect_charge_state.as_dict()
        self.assertEqual(dictionary["degeneracy"], 2)
        self.assertEqual(dictionary["energy"], 0.1234)
        self.assertEqual(dictionary["charge"], 1)
        self.assertEqual(dictionary["fixed_concentration"], 1)

    def test_as_dict_emits_native_floats_for_numpy_values(self):
        cs = DefectChargeState(charge=1, energy=np.float64(0.5), degeneracy=2)
        cs.fix_concentration(np.float64(1e20))
        dictionary = cs.as_dict()
        self.assertIs(type(dictionary["energy"]), float)
        self.assertIs(type(dictionary["fixed_concentration"]), float)
        yaml.safe_dump(dictionary)

    def test_as_dict_preserves_none_energy(self):
        cs = DefectChargeState(charge=1, fixed_concentration=1e20)
        dictionary = cs.as_dict()
        self.assertIsNone(dictionary["energy"])
        yaml.safe_dump(dictionary)

    def test_defect_charge_state_from_dict_warns_unrecognised_key_category(self):
        """from_dict emits UnrecognisedKeyWarning (not a bare UserWarning) for
        keys it does not recognise."""
        dictionary = {
            "degeneracy": 2,
            "energy": 0.1234,
            "charge": 1,
            "foo": "bar",
        }
        with self.assertWarns(UnrecognisedKeyWarning):
            DefectChargeState.from_dict(dictionary)


    def test_name_defaults_to_charge_string(self):
        self.assertEqual(DefectChargeState(charge=2, energy=0.5).name, "q+2")
        self.assertEqual(DefectChargeState(charge=-1, energy=0.5).name, "q-1")
        self.assertEqual(DefectChargeState(charge=0, energy=0.5).name, "q+0")

    def test_explicit_name_overrides_default(self):
        cs = DefectChargeState(charge=2, energy=1.2, name="V_O_2+")
        self.assertEqual(cs.name, "V_O_2+")

    def test_name_round_trips_through_dict(self):
        cs = DefectChargeState(charge=2, energy=1.2, degeneracy=1, name="V_O_2+")
        d = cs.as_dict()
        self.assertEqual(d["name"], "V_O_2+")
        cs2 = DefectChargeState.from_dict(d)
        self.assertEqual(cs2.name, "V_O_2+")

    def test_as_dict_omits_defaulted_name(self):
        cs = DefectChargeState(charge=1, energy=0.5)
        self.assertNotIn("name", cs.as_dict())

    def test_from_dict_without_name_uses_charge_default(self):
        cs = DefectChargeState.from_dict({"charge": 1, "energy": 0.5, "degeneracy": 1})
        self.assertEqual(cs.name, "q+1")

    def test_round_trip_preserves_energy_of_fixed_concentration_state(self):
        cs = DefectChargeState(charge=1, energy=0.7, degeneracy=1, name="V_O_frozen")
        cs.fix_concentration(0.01)
        reloaded = DefectChargeState.from_dict(cs.as_dict())
        self.assertEqual(reloaded.energy, 0.7)
        self.assertEqual(reloaded.fixed_concentration, 0.01)
        self.assertEqual(reloaded.name, "V_O_frozen")


class TestDefectChargeStateRepr(unittest.TestCase):

    def setUp(self):
        charge = 1
        energy = 0.1234
        degeneracy = 2
        self.defect_charge_state = DefectChargeState(
            charge=charge, energy=energy, degeneracy=degeneracy
        )

    def test_repr(self):
        self.assertEqual(
            str(self.defect_charge_state), "q=+1, e=0.1234, deg=2",
        )

    def test_repr_with_name(self):
        cs = DefectChargeState(charge=1, energy=0.1234, degeneracy=2, name="V_O_1+")
        self.assertEqual(str(cs), "q=+1, e=0.1234, deg=2, name=V_O_1+")


if __name__ == "__main__":
    unittest.main()
