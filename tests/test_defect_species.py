import unittest
from unittest.mock import Mock, patch

from numpy.testing import assert_almost_equal, assert_equal

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

    def test_fix_concentration(self):
        assert self.defect_species.fixed_concentration == None
        self.defect_species.fix_concentration(0.1234)
        self.assertEqual(self.defect_species.fixed_concentration,0.1234)

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

    def test_min_energy_charge_state(self):

        charge_state_1 = DefectChargeState(-1, 0.1, 1)
        charge_state_2 = DefectChargeState(-2, 0.2, 1)
        defect = DefectSpecies('v_O', 1, [charge_state_1, charge_state_2])
        self.assertEqual(defect.min_energy_charge_state(0), defect.charge_states[-1])
        self.assertEqual(defect.min_energy_charge_state(2), defect.charge_states[-2])

    def test_get_concentrations(self):

        charge_state_1 = DefectChargeState(-1, 1.1, 1)
        charge_state_2 = DefectChargeState(-2, 1.2, 1)
        defect = DefectSpecies('v_O', 1, [charge_state_1, charge_state_2])
        
        assert_almost_equal(defect.get_concentration(0.2, 300), 3.7117030892903665e-14)
        defect.fix_concentration(0.1234)
        assert_almost_equal(defect.get_concentration(0.2, 300), 0.1234)

        #self.assertEqual(defect.min_energy_charge_state(0), defect.charge_states[-1])
        #self.assertEqual(defect.min_energy_charge_state(2), defect.charge_states[-2])


if __name__ == "__main__":
    unittest.main()
