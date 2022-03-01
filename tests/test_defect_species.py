import unittest
from unittest.mock import Mock, patch

from py_sc_fermi.defect_species import DefectSpecies
from py_sc_fermi.defect_charge_state import DefectChargeState, FrozenDefectChargeState


class TestDefectSpeciesInit(unittest.TestCase):
    def test_defect_species_is_initialised(self):
        name = "foo"
        nsites = 2
        mock_charge_states = [
            Mock(spec=DefectChargeState),
            Mock(spec=DefectChargeState),
        ]
        mock_charge_states[0].charge = 0
        mock_charge_states[1].charge = 1
        defect_species = DefectSpecies(
            name=name, nsites=nsites, charge_states=mock_charge_states
        )
        self.assertEqual(defect_species._name, name)
        self.assertEqual(defect_species._nsites, nsites)
        self.assertEqual(defect_species._charge_states[0], mock_charge_states[0])
        self.assertEqual(defect_species._charge_states[1], mock_charge_states[1])
        self.assertEqual(defect_species._fixed_concentration, None)


class TestDefectSpecies(unittest.TestCase):
    def setUp(self):
        name = "V_O"
        nsites = 2
        mock_charge_states = [
            Mock(spec=DefectChargeState),
            Mock(spec=DefectChargeState),
        ]
        mock_charge_states[0].charge = 0
        mock_charge_states[1].charge = 1
        mock_charge_states[0].concentration_is_fixed = False
        mock_charge_states[1].concentration_is_fixed = False
        self.defect_species = DefectSpecies(
            name=name, nsites=nsites, charge_states=mock_charge_states
        )

    def test_name_property(self):
        self.assertEqual(self.defect_species.name, self.defect_species._name)

    def test_nsites_property(self):
        self.assertEqual(self.defect_species.nsites, self.defect_species._nsites)

    def test_charge_states_property(self):
        self.assertEqual(
            self.defect_species.charge_states, self.defect_species._charge_states
        )

    def test_fixed_concentration_property(self):
        self.defect_species._fixed_concentration = 0.1234
        self.assertEqual(
            self.defect_species.fixed_concentration,
            self.defect_species._fixed_concentration,
        )

    def test_charge_states_by_formation_energy(self):
        self.defect_species.charge_states[0].get_formation_energy = Mock(
            return_value=0.3
        )
        self.defect_species.charge_states[1].get_formation_energy = Mock(
            return_value=0.1
        )
        self.defect_species.variable_conc_charge_states = Mock(
            return_value={
                0: self.defect_species.charge_states[0],
                1: self.defect_species.charge_states[1],
            }
        )
        sorted_charge_states = self.defect_species.charge_states_by_formation_energy(
            e_fermi=0.0
        )
        self.assertEqual(sorted_charge_states[0], self.defect_species.charge_states[1])
        self.assertEqual(sorted_charge_states[1], self.defect_species.charge_states[0])

    def test_charge_states_by_formation_energy_with_frozen_charge_state(self):
        name = "V_O"
        nsites = 2
        mock_charge_states = [
            Mock(spec=DefectChargeState),
            Mock(spec=FrozenDefectChargeState),
        ]
        mock_charge_states[0].charge = 0
        mock_charge_states[1].charge = 1
        defect_species = DefectSpecies(
            name=name, nsites=nsites, charge_states=mock_charge_states
        )
        defect_species.variable_conc_charge_states = Mock(
            return_value={0: defect_species.charge_states[0]}
        )
        defect_species.charge_states[0].get_formation_energy = Mock(return_value=0.3)
        sorted_charge_states = defect_species.charge_states_by_formation_energy(
            e_fermi=0.0
        )
        self.assertEqual(sorted_charge_states, [defect_species.charge_states[0]])

    def test_get_formation_energies(self):
        self.defect_species.charge_states[0].get_formation_energy = Mock(
            return_value=0.3
        )
        self.defect_species.charge_states[1].get_formation_energy = Mock(
            return_value=0.1
        )
        self.defect_species.variable_conc_charge_states = Mock(
            return_value={
                0: self.defect_species.charge_states[0],
                1: self.defect_species.charge_states[1],
            }
        )
        formation_energies_dict = self.defect_species.get_formation_energies(0.0)
        self.assertEqual(formation_energies_dict, {0: 0.3, 1: 0.1})


if __name__ == "__main__":
    unittest.main()
