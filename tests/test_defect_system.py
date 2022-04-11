import unittest
from unittest.mock import Mock, patch

from py_sc_fermi.defect_system import DefectSystem
from py_sc_fermi.defect_species import DefectSpecies
from py_sc_fermi.defect_charge_state import DefectChargeState
from py_sc_fermi.dos import DOS

input_string = "1\n12\n0.1\n298\n1\nv_O 1 1\n 1 1 1\n1\nO_i 1e+22\n1\nO_i 1 1e+22\n"


class TestDefectSystemInit(unittest.TestCase):
    def test_defect_system_is_initialised(self):
        volume = 100
        mock_defect_species = [Mock(spec=DefectSpecies), Mock(spec=DefectSpecies)]
        dos = Mock(spec=DOS)
        temperature = 298
        defect_system = DefectSystem(
            defect_species=mock_defect_species,
            volume=volume,
            dos=dos,
            temperature=temperature,
        )
        self.assertEqual(defect_system.volume, volume)
        self.assertEqual(defect_system.dos, dos)
        self.assertEqual(defect_system.temperature, temperature)
        self.assertEqual(defect_system.defect_species[0], mock_defect_species[0])
        self.assertEqual(defect_system.defect_species[1], mock_defect_species[1])


class TestDefectSystem(unittest.TestCase):
    def setUp(self):
        volume = 100
        mock_defect_species = [Mock(spec=DefectSpecies), Mock(spec=DefectSpecies)]
        mock_defect_species[0].name = "v_O"
        mock_defect_species[1].name = "O_i"
        dos = Mock(spec=DOS)
        dos.spin_polarised = True
        dos._nelect = 12
        dos._egap = 0.1
        temperature = 298
        self.defect_system = DefectSystem(
            defect_species=mock_defect_species,
            volume=volume,
            dos=dos,
            temperature=temperature,
        )

    def test_defect_species_by_name(self):
        self.assertEqual(
            self.defect_system.defect_species_by_name("v_O"),
            self.defect_system.defect_species[0],
        )

    def test_total_defect_charge_contributions(self):
        self.defect_system.defect_species[0].defect_charge_contributions = Mock(
            return_value=[1, 2]
        )
        self.defect_system.defect_species[1].defect_charge_contributions = Mock(
            return_value=[2, 1]
        )
        self.assertEqual(
            self.defect_system.total_defect_charge_contributions(1), (3, 3)
        )

    def test_q_tot(self):
        self.defect_system.dos.carrier_concentrations = Mock(return_value=(1, 1))
        self.defect_system.total_defect_charge_contributions = Mock(return_value=(1, 1))
        self.assertEqual(self.defect_system.q_tot(2), 0)

    def test_as_dict(self):
        self.defect_system.get_sc_fermi = Mock(return_value=[1, {}])
        self.defect_system.dos.carrier_concentrations = Mock(return_value=(1, 1))
        self.defect_system.defect_species[0].get_concentration = Mock(return_value=1)
        self.defect_system.defect_species[1].get_concentration = Mock(return_value=1)
        self.defect_system.defect_species[0].name = "v_O"
        self.defect_system.defect_species[1].name = "O_i"
        volume = self.defect_system.volume
        self.assertEqual(
            self.defect_system.as_dict(),
            {
                "Fermi Energy": 1,
                "p0": 1 / volume * 1e24,
                "n0": 1 / volume * 1e24,
                "O_i": 1 / volume * 1e24,
                "v_O": 1 / volume * 1e24,
            },
        )

    def test__collect_defect_species_with_fixed_chg_states(self):
        toy_defect_species = [DefectChargeState(1, 1, 1, 1)]
        self.defect_system.defect_species[
            0
        ].fixed_conc_charge_states = toy_defect_species
        self.defect_system.defect_species[1].fixed_conc_charge_states = {}
        self.assertEqual(
            self.defect_system._collect_defect_species_with_fixed_chg_states(),
            {"v_O": toy_defect_species},
        )

    def test__get_input_string(self):

        self.defect_system.defect_species[0].name = "v_O"
        self.defect_system.defect_species[0].nsites = 1
        self.defect_system.defect_species[0].variable_conc_charge_states = {
            1: Mock(spec=DefectChargeState)
        }
        self.defect_system.defect_species[0].variable_conc_charge_states[1].energy = 1
        self.defect_system.defect_species[0].variable_conc_charge_states[
            1
        ].degeneracy = 1
        self.defect_system.defect_species[0]._fixed_concentration = None

        self.defect_system.defect_species[1].name = "O_i"
        self.defect_system.defect_species[1].nsites = 1
        self.defect_system.defect_species[1].variable_conc_charge_states = {}
        self.defect_system.defect_species[1]._fixed_concentration = 1

        self.defect_system._collect_defect_species_with_fixed_chg_states = {
            "O_i": [DefectChargeState(1, 1, 1, 1)]
        }
        self.assertEqual(self.defect_system._get_input_string(), input_string)


if __name__ == "__main__":
    unittest.main()
