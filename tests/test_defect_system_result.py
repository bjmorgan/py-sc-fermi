import unittest
from dataclasses import FrozenInstanceError

from py_sc_fermi.defect_system_result import DefectSystemResult


class DefectSystemResultViewsTestCase(unittest.TestCase):
    def _result(self, volume=100.0, label=None):
        return DefectSystemResult(
            temperature=300.0,
            fermi_energy=0.5,
            volume=volume,
            label=label,
            p0_per_cell=2.0,
            n0_per_cell=4.0,
            charge_state_concentrations_per_cell={
                "V_O": {"q+2": 1.0, "q0": 3.0},
                "Ti_i": {"q+4": 5.0},
            },
        )

    def test_per_cell_fields_stored_verbatim(self):
        r = self._result()
        self.assertEqual(r.p0_per_cell, 2.0)
        self.assertEqual(r.n0_per_cell, 4.0)
        self.assertEqual(
            r.charge_state_concentrations_per_cell["V_O"], {"q+2": 1.0, "q0": 3.0}
        )

    def test_cm3_views_scale_per_cell_by_1e24_over_volume(self):
        r = self._result(volume=100.0)
        scale = 1e24 / 100.0
        self.assertEqual(r.p0, 2.0 * scale)
        self.assertEqual(r.n0, 4.0 * scale)
        self.assertEqual(
            r.charge_state_concentrations,
            {
                "V_O": {"q+2": 1.0 * scale, "q0": 3.0 * scale},
                "Ti_i": {"q+4": 5.0 * scale},
            },
        )

    def test_concentrations_per_cell_are_species_totals(self):
        r = self._result()
        self.assertEqual(r.concentrations_per_cell, {"V_O": 4.0, "Ti_i": 5.0})

    def test_concentrations_equal_summed_breakdown_in_cm3(self):
        r = self._result(volume=100.0)
        scale = 1e24 / 100.0
        self.assertEqual(
            r.concentrations, {"V_O": 4.0 * scale, "Ti_i": 5.0 * scale}
        )

    def test_frozen(self):
        r = self._result()
        with self.assertRaises(FrozenInstanceError):
            r.fermi_energy = 0.9

    def test_equal_by_value(self):
        self.assertEqual(self._result(), self._result())

    def test_unhashable(self):
        with self.assertRaises(TypeError):
            hash(self._result())


if __name__ == "__main__":
    unittest.main()
