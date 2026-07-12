import unittest
from unittest.mock import Mock, PropertyMock, patch

import numpy as np
from scipy.special import logsumexp

from py_sc_fermi.defect_charge_state import DefectChargeState
from py_sc_fermi.defect_species import DefectSpecies, kboltz


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
        
    def test_defect_species_is_initialised_with_fixed_concentration(self):
        name = "foo"
        nsites = 2
        fixed_concentration = 0.1234
        mock_charge_states = [
            Mock(spec=DefectChargeState),
            Mock(spec=DefectChargeState),
        ]
        mock_charge_states[0].charge = 0
        mock_charge_states[1].charge = 1
        defect_species = DefectSpecies(
            name=name,
            nsites=nsites,
            charge_states=mock_charge_states,
            fixed_concentration=fixed_concentration
        )
        self.assertEqual(defect_species._name, name)
        self.assertEqual(defect_species._nsites, nsites)
        self.assertEqual(defect_species._charge_states[0], mock_charge_states[0])
        self.assertEqual(defect_species._charge_states[1], mock_charge_states[1])
        self.assertEqual(defect_species._fixed_concentration, fixed_concentration)

    def test_defect_species_raises_with_empty_charge_states(self):
        with self.assertRaises(ValueError) as cm:
            DefectSpecies(name="V_O", nsites=2, charge_states=[])
        self.assertIn("V_O", str(cm.exception))

    def test_defect_species_raises_with_non_positive_nsites(self):
        cs = Mock(spec=DefectChargeState)
        cs.charge = 0
        for bad in (0, -1):
            with self.assertRaises(ValueError) as cm:
                DefectSpecies(name="V_O", nsites=bad, charge_states=[cs])
            self.assertIn("V_O", str(cm.exception))

    def test_construction_consumes_charge_states_once(self):
        # A single-pass iterable must not be exhausted by validation before
        # the states are stored.
        states = (DefectChargeState(charge=q, energy=1.0) for q in (0, 1))
        species = DefectSpecies(name="V_O", nsites=1, charge_states=states)
        self.assertEqual([cs.charge for cs in species.charge_states], [0, 1])

    def test_duplicate_explicit_names_raise(self):
        cs_a = DefectChargeState(charge=1, energy=1.0, name="dup")
        cs_b = DefectChargeState(charge=2, energy=2.0, name="dup")
        with self.assertRaises(ValueError):
            DefectSpecies(name="V_O", nsites=1, charge_states=[cs_a, cs_b])

    def test_unnamed_metastable_states_raise(self):
        # Two unnamed states sharing a charge collide on the default name.
        cs_a = DefectChargeState(charge=1, energy=1.0)
        cs_b = DefectChargeState(charge=1, energy=1.4)
        with self.assertRaisesRegex(ValueError, "q\\+1"):
            DefectSpecies(name="V_O", nsites=1, charge_states=[cs_a, cs_b])

    def test_named_metastable_states_are_accepted(self):
        cs_a = DefectChargeState(charge=1, energy=1.0, name="V_O_tet")
        cs_b = DefectChargeState(charge=1, energy=1.4, name="V_O_oct")
        species = DefectSpecies(name="V_O", nsites=1, charge_states=[cs_a, cs_b])
        self.assertEqual([cs.name for cs in species.charge_states], ["V_O_tet", "V_O_oct"])

    def test_explicit_name_colliding_with_default_raises(self):
        cs_a = DefectChargeState(charge=1, energy=1.0)              # name -> "q+1"
        cs_b = DefectChargeState(charge=2, energy=2.0, name="q+1")  # collides
        with self.assertRaises(ValueError):
            DefectSpecies(name="V_O", nsites=1, charge_states=[cs_a, cs_b])

    def test_charge_state_by_name(self):
        cs_a = DefectChargeState(charge=0, energy=1.0)
        cs_b = DefectChargeState(charge=2, energy=0.5, name="V_O_2+")
        species = DefectSpecies(name="V_O", nsites=1, charge_states=[cs_a, cs_b])
        self.assertIs(species.charge_state_by_name("q+0"), cs_a)
        self.assertIs(species.charge_state_by_name("V_O_2+"), cs_b)

    def test_charge_state_by_name_raises_on_unknown(self):
        cs_a = DefectChargeState(charge=0, energy=1.0)
        species = DefectSpecies(name="V_O", nsites=1, charge_states=[cs_a])
        with self.assertRaisesRegex(ValueError, "available: q\\+0"):
            species.charge_state_by_name("nope")

    def test_charge_states_collection_is_owned(self):
        # The species owns its charge-state collection: neither mutating the
        # constructor list afterwards nor the returned collection can bypass
        # the name-uniqueness check.
        states = [DefectChargeState(charge=0, energy=1.0)]
        species = DefectSpecies(name="V_O", nsites=1, charge_states=states)
        states.append(DefectChargeState(charge=0, energy=1.4))
        self.assertEqual(len(species.charge_states), 1)
        with self.assertRaises(AttributeError):
            species.charge_states.append(DefectChargeState(charge=0, energy=1.4))

    def test_from_dict_with_duplicate_charge_state_names_raises(self):
        d = {
            "name": "V_O",
            "nsites": 1,
            "charge_states": [
                {"charge": 1, "energy": 1.0, "degeneracy": 1},
                {"charge": 1, "energy": 1.4, "degeneracy": 1},
            ],
        }
        with self.assertRaisesRegex(ValueError, "q\\+1"):
            DefectSpecies.from_dict(d)

    def test_init_rejects_invalid_fixed_concentration(self):
        for bad in (-0.1, float("nan"), float("inf")):
            with self.assertRaisesRegex(ValueError, r"'V_O'.*finite and non-negative"):
                DefectSpecies(
                    "V_O", nsites=1,
                    charge_states=[DefectChargeState(charge=0, energy=1.0)],
                    fixed_concentration=bad,
                )

    def test_init_accepts_zero_fixed_concentration(self):
        species = DefectSpecies(
            "V_O", nsites=1,
            charge_states=[DefectChargeState(charge=0, energy=1.0)],
            fixed_concentration=0.0,
        )
        self.assertEqual(species.fixed_concentration, 0.0)


class TestDefectSpecies(unittest.TestCase):
    def setUp(self):
        name = "V_O"
        nsites = 2
        mock_charge_states = [
            Mock(spec=DefectChargeState),
            Mock(spec=DefectChargeState),
            Mock(spec=DefectChargeState),
        ]
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
            self.defect_species.charge_states, tuple(self.defect_species._charge_states)
        )

    def test_fixed_concentration_property(self):
        self.defect_species._fixed_concentration = 0.1234
        self.assertEqual(
            self.defect_species.fixed_concentration,
            self.defect_species._fixed_concentration,
        )

    def test__fix_concentration(self):
        self.assertEqual(self.defect_species.fixed_concentration, None)
        self.defect_species._fix_concentration(0.1234)
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
            return_value=[
                self.defect_species.charge_states[0],
                self.defect_species.charge_states[1],
                self.defect_species.charge_states[2],
            ]
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
            return_value=[
                self.defect_species.charge_states[0],
                self.defect_species.charge_states[2],
            ]
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
            return_value=[
                self.defect_species.charge_states[0],
                self.defect_species.charge_states[1],
                self.defect_species.charge_states[2],
            ]
        )
        formation_energies_dict = self.defect_species.get_formation_energies(0.0)
        self.assertEqual(formation_energies_dict, {0: 0.3, 1: 0.1, 2: 0.5})

    def test_get_formation_energies_with_metastable_charge_state(self):
        # two q=+1 states (a ground state at 0.5 eV and a metastable form at
        # 0.9 eV) must collapse to the lower-energy (0.5 eV) value, not
        # whichever is listed last.
        cs_low = DefectChargeState(charge=1, energy=0.5, degeneracy=1, name="V_O_1+_ground")
        cs_high = DefectChargeState(charge=1, energy=0.9, degeneracy=1, name="V_O_1+_metastable")
        cs_other = DefectChargeState(charge=0, energy=0.0, degeneracy=1)
        defect = DefectSpecies(
            name="V_O", nsites=1, charge_states=[cs_low, cs_high, cs_other]
        )
        formation_energies = defect.get_formation_energies(0.0)
        self.assertEqual(formation_energies[1], 0.5)
        self.assertEqual(formation_energies[0], 0.0)

    def test_get_formation_energies_at_finite_temperature(self):
        cs_a = DefectChargeState(charge=1, energy=0.5, degeneracy=1, name="V_O_1+_a")
        cs_b = DefectChargeState(charge=1, energy=0.6, degeneracy=2, name="V_O_1+_b")
        defect = DefectSpecies(name="V_O", nsites=1, charge_states=[cs_a, cs_b])

        temperature = 300.0
        kT = kboltz * temperature
        expected = -kT * logsumexp([np.log(1) - 0.5 / kT, np.log(2) - 0.6 / kT])

        formation_energies = defect.get_formation_energies(0.0, temperature=temperature)
        self.assertAlmostEqual(formation_energies[1], expected, places=10)
        # the Boltzmann-weighted energy of two states is below the energy of
        # the lowest individual state
        self.assertLess(formation_energies[1], 0.5)

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
            return_value=[
                self.defect_species.charge_states[0],
                self.defect_species.charge_states[1],
                self.defect_species.charge_states[2],
            ]
        )
        self.assertEqual(
            self.defect_species.min_energy_charge_state(0),
            self.defect_species.charge_states[1],
        )

    def test_min_energy_charge_state_raises_with_no_variable_charge_states(self):
        cs_0 = DefectChargeState(charge=0, fixed_concentration=0.5, degeneracy=1)
        cs_1 = DefectChargeState(charge=1, fixed_concentration=0.3, degeneracy=1)
        defect = DefectSpecies(name="V_O", nsites=1, charge_states=[cs_0, cs_1])
        with self.assertRaises(ValueError) as cm:
            defect.min_energy_charge_state(e_fermi=0.0)
        self.assertIn("V_O", str(cm.exception))

    def test_effective_formation_energy_at_zero_temperature_matches_min_energy_state(
        self,
    ):
        defect = DefectSpecies(
            "V_O",
            nsites=1,
            charge_states=[
                DefectChargeState(charge=1, energy=0.5, degeneracy=1, name="V_O_1+_ground"),
                DefectChargeState(charge=1, energy=0.9, degeneracy=1, name="V_O_1+_meta"),
                DefectChargeState(charge=0, energy=2.0, degeneracy=1),
            ],
        )
        e_fermi = 0.3
        expected = defect.min_energy_charge_state(e_fermi).get_formation_energy(
            e_fermi
        )
        self.assertAlmostEqual(
            defect.effective_formation_energy(e_fermi), expected, places=10
        )

    def test_effective_formation_energy_at_finite_temperature(self):
        cs_a = DefectChargeState(charge=1, energy=0.5, degeneracy=1)
        cs_b = DefectChargeState(charge=0, energy=0.6, degeneracy=2)
        defect = DefectSpecies("V_O", nsites=1, charge_states=[cs_a, cs_b])

        e_fermi = 0.3
        temperature = 300.0
        kT = kboltz * temperature
        expected = -kT * logsumexp(
            [
                np.log(cs_a.degeneracy) - cs_a.get_formation_energy(e_fermi) / kT,
                np.log(cs_b.degeneracy) - cs_b.get_formation_energy(e_fermi) / kT,
            ]
        )

        result = defect.effective_formation_energy(e_fermi, temperature=temperature)
        self.assertAlmostEqual(result, expected, places=10)
        # the Boltzmann sum over multiple states is below the lowest
        # individual formation energy
        self.assertLess(
            result,
            min(
                cs_a.get_formation_energy(e_fermi), cs_b.get_formation_energy(e_fermi)
            ),
        )

    def test_effective_formation_energy_raises_with_no_variable_charge_states(self):
        cs_0 = DefectChargeState(charge=0, fixed_concentration=0.5, degeneracy=1)
        cs_1 = DefectChargeState(charge=1, fixed_concentration=0.3, degeneracy=1)
        defect = DefectSpecies(name="V_O", nsites=1, charge_states=[cs_0, cs_1])
        with self.assertRaises(ValueError) as cm:
            defect.effective_formation_energy(e_fermi=0.0)
        self.assertIn("V_O", str(cm.exception))

    def test_get_concentrations(self):
        with patch(
            "py_sc_fermi.defect_species.DefectSpecies.fixed_concentration",
            new_callable=PropertyMock,
        ) as mock_fixed_concentration:
            mock_fixed_concentration.return_value = 0.1234
            self.assertEqual(self.defect_species.get_concentration(1.5, 298), 0.1234)

        self.defect_species.charge_state_concentrations = Mock(
            return_value=[(Mock(spec=DefectChargeState), 0.1234) for _ in range(3)]
        )
        self.assertEqual(self.defect_species.get_concentration(1.5, 298), 0.1234 * 3)

    def test_get_transition_level_and_energy(self):
        self.defect_species.get_formation_energies = Mock(return_value={0: 1, 1: 0})
        self.assertEqual(
            self.defect_species.get_transition_level_and_energy(0, 1), (1, 1)
        )

    def test_fixed_concentration_charge_states(self):
        self.defect_species.charge_states[0].fixed_concentration = 0.1234
        self.defect_species.charge_states[1].fixed_concentration = None
        self.defect_species.charge_states[2].fixed_concentration = None
        self.assertEqual(
            self.defect_species.fixed_conc_charge_states(),
            [self.defect_species.charge_states[0]],
        )

    def test_charge_states_and_defect_species_have_fixed_concentration(self):
        self.defect_species._fix_concentration(concentration=0.1234*3)
        for cs in self.defect_species.charge_states:
            cs.fixed_concentration = 0.1234
            cs.get_concentration = Mock(return_value=0.1234)
        self.assertEqual(
            self.defect_species.fixed_conc_charge_states(),
            [
                self.defect_species.charge_states[0],
                self.defect_species.charge_states[1],
                self.defect_species.charge_states[2]
            ],
        )
        print(self.defect_species.charge_state_concentrations(e_fermi=1.5, temperature=298.0))

    def test_variable_concentration_charge_states(self):
        self.defect_species.charge_states[0].fixed_concentration = 0.1234
        self.defect_species.charge_states[1].fixed_concentration = 0.1234
        self.defect_species.charge_states[2].fixed_concentration = None
        self.assertEqual(
            self.defect_species.variable_conc_charge_states(),
            [self.defect_species.charge_states[2]],
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
            [
                (self.defect_species.charge_states[0], 0.1234 * self.defect_species.nsites),
                (self.defect_species.charge_states[1], 0.1234 * self.defect_species.nsites),
                (self.defect_species.charge_states[2], 0.1234 * self.defect_species.nsites),
            ],
        )

    def test_charge_state_concentrations_with_fixed_concentration(self):
        self.defect_species.charge_states[0].fixed_concentration = None
        self.defect_species.charge_states[1].fixed_concentration = None
        self.defect_species.charge_states[2].fixed_concentration = None
        self.defect_species.charge_states[0].get_formation_energy = Mock(return_value=0.0)
        self.defect_species.charge_states[1].get_formation_energy = Mock(return_value=0.0)
        self.defect_species.charge_states[2].get_formation_energy = Mock(return_value=0.0)
        self.defect_species.charge_states[0].get_concentration = Mock(return_value=0.1)
        self.defect_species.charge_states[1].get_concentration = Mock(return_value=0.1)
        self.defect_species.charge_states[2].get_concentration = Mock(return_value=0.1)
        self.defect_species.charge_states[0].degeneracy = 1
        self.defect_species.charge_states[1].degeneracy = 1
        self.defect_species.charge_states[2].degeneracy = 1
        self.defect_species._fixed_concentration = 0.1234

        result = self.defect_species.charge_state_concentrations(1.5, 298)

        # With equal formation energies and degeneracies, should be equal distribution
        expected = 0.1234 / 3
        for _, conc in result:
            self.assertAlmostEqual(conc, expected, places=10)

    def test_defect_charge_contributions(self):
        cs_positive = Mock(spec=DefectChargeState)
        cs_positive.charge = 1
        self.defect_species.charge_state_concentrations = Mock(
            return_value=[(cs_positive, 0.1234)]
        )
        self.assertEqual(
            self.defect_species.defect_charge_contributions(1.5, 298), (0.1234, 0)
        )
        cs_negative = Mock(spec=DefectChargeState)
        cs_negative.charge = -1
        self.defect_species.charge_state_concentrations = Mock(
            return_value=[(cs_negative, 0.1234)]
        )
        self.assertEqual(
            self.defect_species.defect_charge_contributions(1.5, 298), (0, 0.1234)
        )

    def test_tl_profile(self):
        # Updated test to check the functionality of this method more directly
        charge_state_1 = DefectChargeState(0, energy=2, degeneracy=1)
        charge_state_2 = DefectChargeState(2, energy=-1, degeneracy=1)
        defect = DefectSpecies("foo", 1, [charge_state_1, charge_state_2])
        tl_profile = defect.tl_profile(0, 5)
        self.assertEqual(tl_profile[0][0], 0)
        self.assertEqual(tl_profile[0][1], -1)
        self.assertEqual(tl_profile[1][0], 1.5)
        self.assertEqual(tl_profile[1][1], 2)
        self.assertEqual(tl_profile[2][0], 5)
        self.assertEqual(tl_profile[2][1], 2)

    def test_tl_profile_with_metastable_charge_state(self):
        # q=+1 has a ground state (0.5 eV) and a higher-energy metastable
        # form (0.9 eV). The resulting profile must reflect the ground state
        # and must not depend on the order the charge states are listed in.
        defect_a = DefectSpecies(
            "foo",
            1,
            charge_states=[
                DefectChargeState(charge=1, energy=0.5, degeneracy=1, name="foo_1+_ground"),
                DefectChargeState(charge=1, energy=0.9, degeneracy=1, name="foo_1+_meta"),
                DefectChargeState(charge=0, energy=2.0, degeneracy=1),
                DefectChargeState(charge=-1, energy=2.0, degeneracy=1),
            ],
        )
        defect_b = DefectSpecies(
            "foo",
            1,
            charge_states=[
                DefectChargeState(charge=1, energy=0.9, degeneracy=1, name="foo_1+_meta"),
                DefectChargeState(charge=1, energy=0.5, degeneracy=1, name="foo_1+_ground"),
                DefectChargeState(charge=0, energy=2.0, degeneracy=1),
                DefectChargeState(charge=-1, energy=2.0, degeneracy=1),
            ],
        )
        profile_a = defect_a.tl_profile(0, 5)
        profile_b = defect_b.tl_profile(0, 5)
        np.testing.assert_array_almost_equal(profile_a, profile_b)
        # the q=+1 endpoint reflects the lower-energy (0.5 eV) state
        self.assertAlmostEqual(profile_a[0][1], 0.5)

    def test__repr__(self):
        self.defect_species._charge_states = [
            DefectChargeState(2, energy=-1, degeneracy=1)
        ]
        self.defect_species._fixed_concentration = 100
        self.assertEqual(
            str(self.defect_species),
            "\nV_O, nsites=2\nfixed [c] = 100\n  q=+2, e=-1, deg=1\n",
        )

    def test_from_dict(self):
        d = {
            "name": "V_O",
            "nsites": 2,
            "charge_states": [{"charge": 1, "energy": 0, "degeneracy": 1}],
        }
        self.assertEqual(DefectSpecies.from_dict(d).name, "V_O")
        self.assertEqual(DefectSpecies.from_dict(d).nsites, 2)
        self.assertEqual(DefectSpecies.from_dict(d).charge_states[0].charge, 1)
        self.assertEqual(DefectSpecies.from_dict(d).charge_states[0].energy, 0)

    def test_from_dict_with_fixed_concentration(self):
        d = {
            "name": "V_O",
            "nsites": 2,
            "charge_states": [
                {"charge": 1, "fixed_concentration": 100, "degeneracy": 1}
            ],
            "fixed_concentration": 100,
        }
        self.assertEqual(DefectSpecies.from_dict(d).name, "V_O")
        self.assertEqual(DefectSpecies.from_dict(d).nsites, 2)
        self.assertEqual(DefectSpecies.from_dict(d).fixed_concentration, 100)
        self.assertEqual(
            DefectSpecies.from_dict(d).charge_states[0].fixed_concentration,
            100,
        )

    def test_as_dict(self):
        # Setup for the test:
        mock_charge_states = [
            Mock(spec=DefectChargeState),
            Mock(spec=DefectChargeState),
            Mock(spec=DefectChargeState),
        ]
        mock_charge_states[0].as_dict.return_value = {"charge": 0}
        mock_charge_states[1].as_dict.return_value = {"charge": 1}
        mock_charge_states[2].as_dict.return_value = {"charge": 2}

        self.defect_species._charge_states = mock_charge_states
        self.defect_species._name = "v_O"
        self.defect_species._nsites = 2
        self.defect_species._fixed_concentration = 0.1234

        # Call the method and get the result
        result = self.defect_species.as_dict()

        # Expected result
        expected_result = {
            "name": "v_O",
            "nsites": 2,
            "charge_states": [{"charge": 0}, {"charge": 1}, {"charge": 2}],
            "fixed_concentration": 0.1234,
        }

        # Verify the result
        self.assertEqual(result, expected_result)

    def test_charge_state_concentrations_all_fixed(self):
        """When all charge states have fixed concentrations, return those values."""
        cs_0 = DefectChargeState(charge=0, fixed_concentration=0.5, degeneracy=1)
        cs_1 = DefectChargeState(charge=1, fixed_concentration=0.3, degeneracy=1)

        defect = DefectSpecies(
            name="test",
            nsites=1,
            charge_states=[cs_0, cs_1]
        )

        result = dict(
            (cs.charge, conc)
            for cs, conc in defect.charge_state_concentrations(e_fermi=0.0, temperature=298)
        )

        self.assertEqual(result[0], 0.5)
        self.assertEqual(result[1], 0.3)

    def test_charge_state_concentrations_with_fixed_concentration_zero_variable_concs(self):
        """Fixed concentration scaling should not produce NaN when variable
        concentrations underflow."""
        cs_0 = DefectChargeState(charge=0, energy=1.0, degeneracy=1)
        cs_minus1 = DefectChargeState(charge=-1, energy=2.0, degeneracy=1)

        # Mock get_formation_energy to return very large values (simulating underflow)
        cs_0.get_formation_energy = Mock(return_value=1000.0)
        cs_minus1.get_formation_energy = Mock(return_value=1000.0)

        defect = DefectSpecies(
            name="v_Na",
            nsites=1,
            charge_states=[cs_0, cs_minus1],
            fixed_concentration=1e-5,
        )

        result = defect.charge_state_concentrations(e_fermi=0.0, temperature=300)

        for cs, conc in result:
            self.assertFalse(np.isnan(conc), f"Concentration for charge {cs.charge} is NaN")

        # Total should equal fixed concentration
        total = sum(conc for _, conc in result)
        self.assertAlmostEqual(total, 1e-5, places=15)
            
    def test_charge_state_concentrations_raises_if_fixed_exceeds_total(self):
        """Should raise ValueError if fixed charge state concentrations exceed species total."""
        cs_0 = DefectChargeState(charge=0, energy=0.5, degeneracy=1)
        cs_1 = DefectChargeState(charge=1, energy=0.6, degeneracy=1, fixed_concentration=0.5)
        
        defect = DefectSpecies(
            name="test",
            nsites=1,
            charge_states=[cs_0, cs_1],
            fixed_concentration=0.1,  # Less than fixed charge state concentration
        )

        with self.assertRaises(ValueError):
            defect.charge_state_concentrations(e_fermi=0.0, temperature=298)

    def test_charge_state_concentrations_raises_if_below_total_with_no_variable(self):
        """Should raise ValueError if all charge states are fixed and sum to
        less than the species total, with no variable state to make up the
        difference."""
        cs_0 = DefectChargeState(charge=0, fixed_concentration=3.0)
        defect = DefectSpecies(
            name="test", nsites=10, charge_states=[cs_0], fixed_concentration=5.0
        )

        with self.assertRaises(ValueError):
            defect.charge_state_concentrations(e_fermi=0.0, temperature=298)


class TestDefectSpeciesHashability(unittest.TestCase):
    def test_defect_species_is_hashable(self):
        # Load-bearing: ElementPool accepts DefectSpecies as member keys, and the
        # solver keys dict[DefectSpecies, float] / set[DefectSpecies] on them.
        # Adding a value __eq__ without a matching __hash__ would break those;
        # this fails loudly here rather than mysteriously inside a pool dict.
        cs = Mock(spec=DefectChargeState)
        cs.charge = 0
        sp = DefectSpecies(name="V_O", nsites=2, charge_states=[cs])
        self.assertEqual(hash(sp), hash(sp))
        self.assertEqual({sp: 1}[sp], 1)


if __name__ == "__main__":
    unittest.main()
