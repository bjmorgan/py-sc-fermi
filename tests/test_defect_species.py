import unittest
from unittest.mock import Mock, PropertyMock, patch

from copy import deepcopy

from numpy.testing import assert_almost_equal, assert_equal

from py_sc_fermi.defect_species import DefectSpecies
from py_sc_fermi.defect_charge_state import DefectChargeState


class TestDefectSpeciesInit(unittest.TestCase):
    def test_defect_species_is_initialised(self):
        name = "foo"
        nsites = 2
        mock_charge_states = {
            0 : Mock(spec=DefectChargeState),
            1 : Mock(spec=DefectChargeState),
        }
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
        mock_charge_states = {
            0 : Mock(spec=DefectChargeState),
            1 : Mock(spec=DefectChargeState),
            2 : Mock(spec=DefectChargeState),
        }
        mock_charge_states[0].charge = 0
        mock_charge_states[1].charge = 1
        mock_charge_states[2].charge = 2
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
        self.assertEqual(self.defect_species.fixed_concentration, None)
        self.defect_species.fix_concentration(0.1234)
        self.assertEqual(self.defect_species.fixed_concentration, 0.1234)

    def test_charge_states_by_formation_energy(self):
        self.defect_species.charge_states[0].get_formation_energy = Mock(
            return_value=0.3
        )
        self.defect_species.charge_states[1].get_formation_energy = Mock(
            return_value=0.1
        )
        self.defect_species.charge_states[2].get_formation_energy = Mock(
            return_value=0.5
        )
        self.defect_species.variable_conc_charge_states = Mock(
            return_value={
                0: self.defect_species.charge_states[0],
                1: self.defect_species.charge_states[1],
                2: self.defect_species.charge_states[2],
            }
        )
        sorted_charge_states = self.defect_species.charge_states_by_formation_energy(
            e_fermi=0.0
        )
        self.assertEqual(sorted_charge_states[0], self.defect_species.charge_states[1])
        self.assertEqual(sorted_charge_states[1], self.defect_species.charge_states[0])
        self.assertEqual(sorted_charge_states[2], self.defect_species.charge_states[2])

    def test_charge_states_by_formation_energy_with_frozen_charge_state(self):
        self.defect_species.charge_states[0].get_formation_energy = Mock(
            return_value=0.3
        )
        self.defect_species.charge_states[1].fixed_concentration = Mock(
            return_value=0.1234
        )
        self.defect_species.charge_states[2].get_formation_energy = Mock(
            return_value=0.5
        )
        self.defect_species.variable_conc_charge_states = Mock(
            return_value={
                0: self.defect_species.charge_states[0],
                2: self.defect_species.charge_states[2],
            }
        )
        sorted_charge_states = self.defect_species.charge_states_by_formation_energy(
            e_fermi=0.0
        )
        self.assertEqual(sorted_charge_states[0], self.defect_species.charge_states[0])
        self.assertEqual(sorted_charge_states[1], self.defect_species.charge_states[2])

    def test_get_formation_energies(self):
        self.defect_species.charge_states[0].get_formation_energy = Mock(
            return_value=0.3
        )
        self.defect_species.charge_states[1].get_formation_energy = Mock(
            return_value=0.1
        )
        self.defect_species.charge_states[2].get_formation_energy = Mock(
            return_value=0.5
        )
        self.defect_species.variable_conc_charge_states = Mock(
            return_value={
                0: self.defect_species.charge_states[0],
                1: self.defect_species.charge_states[1],
                2: self.defect_species.charge_states[2],
            }
        )
        formation_energies_dict = self.defect_species.get_formation_energies(0.0)
        self.assertEqual(formation_energies_dict, {0: 0.3, 1: 0.1, 2: 0.5})

    def test_min_energy_charge_state(self):
        self.defect_species.charge_states[0].get_formation_energy = Mock(
            return_value=0.3
        )
        self.defect_species.charge_states[1].get_formation_energy = Mock(
            return_value=0.1
        )
        self.defect_species.charge_states[2].get_formation_energy = Mock(
            return_value=0.5
        )
        self.defect_species.variable_conc_charge_states = Mock(
            return_value={
                0: self.defect_species.charge_states[0],
                1: self.defect_species.charge_states[1],
                2: self.defect_species.charge_states[2],
            }
        )
        self.assertEqual(
            self.defect_species.min_energy_charge_state(0),
            self.defect_species.charge_states[1],
        )

    def test_get_concentrations(self):
        with patch(
            "py_sc_fermi.defect_species.DefectSpecies.fixed_concentration",
            new_callable=PropertyMock,
        ) as mock_fixed_concentration:
            mock_fixed_concentration.return_value = 0.1234
            self.assertEqual(self.defect_species.get_concentration(1.5, 298), 0.1234)

        self.defect_species.charge_state_concentrations = Mock(
            return_value={0: 0.1234, 1: 0.1234, 2: 0.1234}
        )
        self.assertEqual(self.defect_species.get_concentration(1.5, 298), 0.1234 * 3)

    def test_get_transition_level_and_energy(self):
        self.defect_species.get_formation_energies = Mock(return_value={0: 1, 1: 0})
        self.assertEqual(
            self.defect_species.get_transition_level_and_energy(0, 1), (1, 1)
        )

    def test_fixed_concentration_charge_states(self):
        self.defect_species.charge_states[0].fixed_concentration = Mock(
            return_value=0.1234
        )
        self.defect_species.charge_states[1].fixed_concentration = None
        self.defect_species.charge_states[2].fixed_concentration = None
        self.assertEqual(
            self.defect_species.fixed_conc_charge_states(),
            {0: self.defect_species.charge_states[0]},
        )

    def test_variable_concentration_charge_states(self):
        self.defect_species.charge_states[0].fixed_concentration = Mock(
            return_value=0.1234
        )
        self.defect_species.charge_states[1].fixed_concentration = Mock(
            return_value=0.1234
        )
        self.defect_species.charge_states[2].fixed_concentration = None
        self.assertEqual(
            self.defect_species.variable_conc_charge_states(),
            {2: self.defect_species.charge_states[2]},
        )

    def test_charge_state_concentrations(self):
        self.defect_species.charge_states[0].fixed_concentration = None
        self.defect_species.charge_states[1].fixed_concentration = None
        self.defect_species.charge_states[2].fixed_concentration = None
        self.defect_species.charge_states[0].get_concentration = Mock(
            return_value=0.1234
        )
        self.defect_species.charge_states[1].get_concentration = Mock(
            return_value=0.1234
        )
        self.defect_species.charge_states[2].get_concentration = Mock(
            return_value=0.1234
        )
        self.assertEqual(
            self.defect_species.charge_state_concentrations(1.5, 298),
            {
                0: 0.1234 * self.defect_species.nsites,
                1: 0.1234 * self.defect_species.nsites,
                2: 0.1234 * self.defect_species.nsites,
            },
        )

    def test_charge_state_concentrations_with_fixed_concentration(self):
        self.defect_species.charge_states[0].fixed_concentration = None
        self.defect_species.charge_states[1].fixed_concentration = None
        self.defect_species.charge_states[2].fixed_concentration = None
        self.defect_species.charge_states[0].get_concentration = Mock(
            return_value=0.1234/3
        )
        self.defect_species.charge_states[1].get_concentration = Mock(
            return_value=0.1234/3
        )
        self.defect_species.charge_states[2].get_concentration = Mock(
            return_value=0.1234/3
        )
        self.defect_species._fixed_concentration = 0.1234
        self.assertEqual(
            self.defect_species.charge_state_concentrations(1.5, 298),
            {
                0: 0.04113333333333333,
                1: 0.04113333333333333,
                2: 0.04113333333333333,
            },
        )

    def test_defect_charge_contributions(self):
        self.defect_species.charge_state_concentrations = Mock(return_value={1: 0.1234})
        self.assertEqual(
            self.defect_species.defect_charge_contributions(1.5, 298), (0.1234, 0)
        )
        self.defect_species.charge_state_concentrations = Mock(
            return_value={-1: 0.1234}
        )
        self.assertEqual(
            self.defect_species.defect_charge_contributions(1.5, 298), (0, 0.1234)
        )

    def test_tl_profile(self):
        # TODO: ideally, this test should more directly check the
        # funtionality of this method
        charge_state_1 = DefectChargeState(0, energy=2, degeneracy=1)
        charge_state_2 = DefectChargeState(2, energy=-1, degeneracy=1)
        defect = DefectSpecies("foo", 1, {0 : charge_state_1, 2 :charge_state_2})
        assert_equal(defect.tl_profile(0, 5), [[0, -1], [1.5, 2], [5, 2]])

    def test__repr__(self):
        self.defect_species._charge_states = {2 : DefectChargeState(2, energy=-1, degeneracy=1)}
        self.defect_species._fixed_concentration = 100
        self.assertEqual(
            str(self.defect_species),
            "\nV_O, nsites=2\nfixed [c] = 100\n  q=+2, e=-1, deg=1\n"
        )

    def test_from_dict(self):
        d = {"V_O" : {
            "nsites": 2,
            "charge_states": {1 : {"formation_energy": 0, "degeneracy": 1}}}}
        self.assertEqual(DefectSpecies.from_dict(d).name, "V_O")
        self.assertEqual(DefectSpecies.from_dict(d).nsites, 2)
        self.assertEqual(DefectSpecies.from_dict(d).charge_states[1].charge, 1)
        self.assertEqual(DefectSpecies.from_dict(d).charge_states[1].energy, 0)

    def test_from_dict_with_fixed_concentration(self):
        d = {"V_O" : {
            "nsites": 2,
            "charge_states": {1 : {"fixed_concentration": 100, "degeneracy": 1}},
            "fixed_concentration": 100}}
        self.assertEqual(DefectSpecies.from_dict(d, volume=1e24).name, "V_O")
        self.assertEqual(DefectSpecies.from_dict(d, volume=1e24).nsites, 2)
        self.assertEqual(DefectSpecies.from_dict(d, volume = 1e24).fixed_concentration, 100)
        self.assertEqual(DefectSpecies.from_dict(d, volume = 1e24).charge_states[1].fixed_concentration, 100)

    def test__from_string(self):
        string = "V_O 1 2\n2 2 2"
        string = string.splitlines()
        self.assertEqual(DefectSpecies._from_list_of_strings(deepcopy(string)).name, "V_O")
        self.assertEqual(DefectSpecies._from_list_of_strings(deepcopy(string)).nsites, 2)



if __name__ == "__main__":
    unittest.main()
