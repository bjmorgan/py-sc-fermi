import json
import os
import unittest
import warnings
from io import StringIO
from unittest.mock import Mock, patch

import numpy as np
import yaml
from scipy.constants import physical_constants

from py_sc_fermi.defect_charge_state import DefectChargeState
from py_sc_fermi.defect_species import DefectSpecies
from py_sc_fermi.defect_system import (
    DefectSystem,
    DefectSystemFactory,
    DiluteLimitWarning,
)
from py_sc_fermi.dos import DOS
from py_sc_fermi.element_pools import ElementPoolError

input_string = "1\n12\n0.1\n298\n1\nv_O 1 1\n 1 1 1\n1\nO_i 1e+22\n1\nO_i 1 1e+22\n"
input_string_spin = (
    "0\n12\n0.1\n298\n1\nv_O 1 1\n 1 1 1\n1\nO_i 1e+22\n1\nO_i 1 1e+22\n"
)
test_data_dir = "dummy_inputs/"
test_yaml_filename = os.path.join(
    os.path.dirname(__file__), test_data_dir, "defect_system.yaml"
)
test_exception_yaml_filename = os.path.join(
    os.path.dirname(__file__), test_data_dir, "bad_yaml.yaml"
)
test_vasprun_filename = os.path.join(
    os.path.dirname(__file__), test_data_dir, "vasprun_nsp.xml"
)

kboltz = physical_constants["Boltzmann constant in eV/K"][0]



class TestDefectSystemInit(unittest.TestCase):
    def test_defect_system_is_initialised(self):
        volume = 100
        mock_defect_species = [Mock(spec=DefectSpecies), Mock(spec=DefectSpecies)]
        mock_defect_species[0].name = "v_O"
        mock_defect_species[1].name = "O_i"
        mock_defect_species[0].charge_states = []
        mock_defect_species[1].charge_states = []
        mock_defect_species[0].fixed_concentration = None
        mock_defect_species[1].fixed_concentration = None
        mock_defect_species[0].nsites = 1
        mock_defect_species[1].nsites = 1
        dos = Mock(spec=DOS)
        temperature = 298
        defect_system = DefectSystem(
            defect_species=mock_defect_species,
            volume=volume,
            dos=dos,
            temperature=temperature,
            convergence_tolerance=1e-6,
        )
        self.assertEqual(defect_system.volume, volume)
        self.assertEqual(defect_system.dos, dos)
        self.assertEqual(defect_system.temperature, temperature)
        # `defect_species` is deep-copied, so these are independent copies,
        # not the same objects as `mock_defect_species`.
        self.assertEqual(defect_system.defect_species[0].name, "v_O")
        self.assertEqual(defect_system.defect_species[1].name, "O_i")


class TestDefectSystem(unittest.TestCase):
    def setUp(self):
        volume = 100
        mock_defect_species = [Mock(spec=DefectSpecies), Mock(spec=DefectSpecies)]
        mock_defect_species[0].name = "v_O"
        mock_defect_species[1].name = "O_i"
        mock_defect_species[0].charge_states = []
        mock_defect_species[1].charge_states = []
        mock_defect_species[0].fixed_concentration = None
        mock_defect_species[1].fixed_concentration = None
        mock_defect_species[0].nsites = 1
        mock_defect_species[1].nsites = 1
        dos = Mock(spec=DOS)
        dos.spin_polarised = True
        dos._nelect = 12
        dos.bandgap = 0.1
        dos._bandgap = 0.1
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

    def test_defect_species_by_name_raises_for_unknown_name(self):
        with self.assertRaisesRegex(
            ValueError, "no defect species named 'Xx'; available: v_O, O_i"
        ):
            self.defect_system.defect_species_by_name("Xx")

    def test_defect_species_names(self):
        self.assertEqual(self.defect_system.defect_species_names, ["v_O", "O_i"])

    def test_total_defect_charge_contributions(self):
        cs_pos = DefectChargeState(charge=1, fixed_concentration=2)
        cs_neg = DefectChargeState(charge=-1, fixed_concentration=3)
        cs_neutral = DefectChargeState(charge=0, fixed_concentration=5)
        self.defect_system._global_defect_concs = Mock(
            return_value={cs_pos: 2.0, cs_neg: 3.0, cs_neutral: 5.0}
        )
        lhs, rhs = self.defect_system.total_defect_charge_contributions(1)
        self.assertEqual(lhs, 2.0)
        self.assertEqual(rhs, 3.0)

    def test_q_tot(self):
        self.defect_system.dos.carrier_concentrations = Mock(return_value=(1, 1))
        self.defect_system.total_defect_charge_contributions = Mock(return_value=(1, 1))
        self.assertEqual(self.defect_system.q_tot(2), 0)

    def test_as_dict(self):
        self.defect_system.dos = DOS.from_vasprun(test_vasprun_filename, nelect=12)
        defect_dict = self.defect_system.as_dict()
        self.assertEqual(defect_dict["volume"], 100)
        self.assertEqual(defect_dict["temperature"], 298)

    def test_from_dict(self):
        dictionary = {
            "volume": 100,
            "temperature": 100,
            "convergence_tolerance": 1,
            "defect_species": [{
                "name": "V_O",
                "nsites": 2,
                "charge_states": [{"charge": 1, "energy": 0, "degeneracy": 1}],
            }],
            "dos": {
                "dos": np.ones(101),
                "edos": np.linspace(-10.0, 10.0, 101),
                "bandgap": 3.0,
                "nelect": 10,
                "spin_pol": False
            }
        }
        defect_system = self.defect_system.from_dict(dictionary)
        self.assertEqual(defect_system.volume, 100)
        self.assertEqual(defect_system.temperature, 100)
        self.assertEqual(defect_system.convergence_tolerance, 1)

    def test_get_sc_fermi(self):
        self.defect_system.dos.emin = Mock(return_value=0)
        self.defect_system.dos.emax = Mock(return_value=1)
        self.defect_system.q_tot = lambda e_fermi: 0.5 - e_fermi
        
        e_fermi, residual = self.defect_system.get_sc_fermi()
        
        self.assertAlmostEqual(e_fermi, 0.5, places=10)
        self.assertAlmostEqual(residual, 0.0, places=10)

    def test_get_sc_fermi_raises_when_no_solution_positive(self):
        """RuntimeError when q_tot is always positive (no root)."""
        self.defect_system.dos.emin = Mock(return_value=0)
        self.defect_system.dos.emax = Mock(return_value=1)
        self.defect_system.q_tot = lambda e_fermi: 0.1
        with self.assertRaises(RuntimeError):
            self.defect_system.get_sc_fermi()
    
    def test_get_sc_fermi_raises_when_no_solution_negative(self):
        """RuntimeError when q_tot is always negative (no root)."""
        self.defect_system.dos.emin = Mock(return_value=0)
        self.defect_system.dos.emax = Mock(return_value=1)
        self.defect_system.q_tot = lambda e_fermi: -0.1
        with self.assertRaises(RuntimeError):
            self.defect_system.get_sc_fermi()

    def test_get_sc_fermi_no_solution_error_carries_brentq_message(self):
        """The no-solution RuntimeError surfaces brentq's own ValueError
        message rather than discarding it."""
        self.defect_system.dos.emin = Mock(return_value=0)
        self.defect_system.dos.emax = Mock(return_value=1)
        self.defect_system.q_tot = lambda e_fermi: 0.1
        with self.assertRaises(RuntimeError) as ctx:
            self.defect_system.get_sc_fermi()
        self.assertIsInstance(ctx.exception.__cause__, ValueError)
        self.assertIn(str(ctx.exception.__cause__), str(ctx.exception))

    def test_get_sc_fermi_finds_charge_neutral_fermi_energy(self):
        """Solver should find Fermi energy where total charge is zero."""
        dos = DOS(
            dos=np.ones(101),
            edos=np.linspace(-10.0, 10.0, 101),
            bandgap=3.0,
            nelect=10,
        )
        charge_state = DefectChargeState(charge=1, energy=0.5, degeneracy=1)
        defect_species = DefectSpecies(
            name="test_defect",
            nsites=1,
            charge_states=[charge_state],
        )
        defect_system = DefectSystem(
            dos=dos,
            volume=100,
            temperature=300,
            defect_species=[defect_species],
        )
        
        e_fermi, residual = defect_system.get_sc_fermi()
        
        # Verify charge neutrality at converged Fermi energy
        q_tot = defect_system.q_tot(e_fermi)
        self.assertAlmostEqual(q_tot, 0.0, places=10)
        self.assertLess(residual,  1e-10)

    def test_get_sc_fermi_with_extreme_formation_energy_remains_finite(self):
        """Site-exclusion statistics bound every charge state's concentration
        by its group's `nsites`, so q_tot can no longer overflow even for a
        very low formation energy at the edge of the DOS window -- brentq
        brackets directly to the finite interior root."""
        dos = DOS(
            dos=np.ones(101),
            edos=np.linspace(-10.0, 10.0, 101),
            bandgap=3.0,
            nelect=10,
        )
        charge_state = DefectChargeState(charge=-2, energy=0.0, degeneracy=1)
        defect_species = DefectSpecies(
            name="deep_acceptor",
            nsites=1,
            charge_states=[charge_state],
        )
        defect_system = DefectSystem(
            dos=dos,
            volume=100,
            temperature=300,
            defect_species=[defect_species],
        )

        self.assertTrue(np.isfinite(defect_system.q_tot(dos.emax())))
        e_fermi, residual = defect_system.get_sc_fermi()
        self.assertTrue(np.isfinite(e_fermi))
        self.assertLess(residual, 1e-10)

    def test_get_transition_levels(self):
        self.defect_system.defect_species_by_name("v_O").tl_profile = Mock(
            return_value=[[1, 2], [1, 2]]
        )
        self.defect_system.defect_species_by_name("O_i").tl_profile = Mock(
            return_value=[[1, 2], [1, 2]]
        )
        self.assertEqual(
            self.defect_system.get_transition_levels(),
            {"v_O": [[1, 1], [2, 2]], "O_i": [[1, 1], [2, 2]]},
        )

    def test_concentration_dict(self):
        cs_v_O = DefectChargeState(charge=1, fixed_concentration=1)
        cs_O_i = DefectChargeState(charge=-1, fixed_concentration=1)
        self.defect_system.get_sc_fermi = Mock(return_value=[1, {}])
        self.defect_system.dos.carrier_concentrations = Mock(return_value=(1, 1))
        self.defect_system.defect_species[0].get_concentration = Mock(return_value=1)
        self.defect_system.defect_species[1].get_concentration = Mock(return_value=1)
        self.defect_system.defect_species[0].charge_state_concentrations = Mock(
            return_value=[(cs_v_O, 1)]
        )
        self.defect_system.defect_species[1].charge_state_concentrations = Mock(
            return_value=[(cs_O_i, 1)]
        )
        self.defect_system.defect_species[0].charge_states = [cs_v_O]
        self.defect_system.defect_species[1].charge_states = [cs_O_i]
        self.defect_system.defect_species[0].fixed_concentration = None
        self.defect_system.defect_species[1].fixed_concentration = None
        self.defect_system.defect_species[0].nsites = 1
        self.defect_system.defect_species[1].nsites = 1
        self.defect_system.defect_species[0].name = "v_O"
        self.defect_system.defect_species[1].name = "O_i"
        self.defect_system.volume = 100

        expected_dict = {
            "Fermi Energy": 1.0,
            "p0": 1.0e22,
            "n0": 1.0e22,
            "v_O": 1.0e22,
            "O_i": 1.0e22
        }
        result_dict = self.defect_system.concentration_dict()
        self.assertEqual(result_dict, expected_dict)

        expected_decomposed_dict = {
            "Fermi Energy": 1.0,
            "p0": 1.0e22,
            "n0": 1.0e22,
            "v_O": {1: 1.0e22},
            "O_i": {-1: 1.0e22}
        }
        result_decomposed_dict = self.defect_system.concentration_dict(decomposed=True)
        self.assertEqual(result_decomposed_dict, expected_decomposed_dict)



    def test__repr__(self):
        self.defect_system.defect_species = []
        self.defect_system.dos.nelect = 100
        self.defect_system.dos.bandgap = 0.1
        self.assertEqual(
            str(self.defect_system).strip(),
            "DefectSystem\n"
            "  bandgap:     0.1 eV    nelect: 100\n"
            "  volume:      100 Å³    temperature: 298 K\n"
            "\n"
            "  0 defect species:",
        )


class TestDefectSystemSitePools(unittest.TestCase):
    def setUp(self):
        self.dos = DOS(
            dos=np.ones(101),
            edos=np.linspace(-5.0, 5.0, 101),
            bandgap=2.0,
            nelect=10,
        )
        self.species_a = DefectSpecies(
            "A",
            nsites=5,
            charge_states=[
                DefectChargeState(charge=0, energy=1.0, degeneracy=1),
                DefectChargeState(charge=1, energy=1.2, degeneracy=2),
            ],
        )
        self.species_b = DefectSpecies(
            "B",
            nsites=1,
            charge_states=[
                DefectChargeState(charge=0, energy=0.8, degeneracy=1),
                DefectChargeState(charge=-1, energy=1.5, degeneracy=2),
            ],
        )

    def test_no_pools_matches_dilute_limit_for_small_weights(self):
        system = DefectSystem(
            defect_species=[self.species_a, self.species_b],
            dos=self.dos,
            volume=100,
            temperature=300,
        )
        e_fermi = 1.0
        exclusion_total = sum(system._global_defect_concs(e_fermi).values())
        dilute_total = sum(
            conc
            for sp in (self.species_a, self.species_b)
            for _, conc in sp.charge_state_concentrations(e_fermi, 300)
        )
        self.assertAlmostEqual(exclusion_total, dilute_total, places=8)

    def test_unpooled_species_saturates_at_own_nsites(self):
        species = DefectSpecies(
            "C",
            nsites=3,
            charge_states=[DefectChargeState(charge=0, energy=-2.0, degeneracy=1)],
        )
        system = DefectSystem(
            defect_species=[species],
            dos=self.dos,
            volume=100,
            temperature=300,
        )
        concs = system._global_defect_concs(1.0)
        total = sum(concs[cs] for cs in system.defect_species[0].charge_states)
        self.assertLessEqual(total, species.nsites)
        self.assertGreater(total, 0.99 * species.nsites)

    def test_pooled_species_share_sites_within_pool_capacity(self):
        n_pool = 10.0
        system = DefectSystem(
            defect_species=[self.species_a, self.species_b],
            dos=self.dos,
            volume=100,
            temperature=300,
            site_pools={"shared": (n_pool, [self.species_a, self.species_b])},
        )
        concs = system._global_defect_concs(1.0)
        total_occupied = sum(
            concs[cs]
            for sp in system.defect_species
            for cs in sp.charge_states
        )
        self.assertGreater(total_occupied, 0.0)
        self.assertLessEqual(total_occupied, n_pool)

    def test_pool_references_are_normalised_to_names(self):
        system = DefectSystem(
            defect_species=[self.species_a, self.species_b],
            dos=self.dos,
            volume=100,
            temperature=300,
            site_pools={"shared": (10.0, [self.species_a, "B"])},
        )
        self.assertEqual(system.site_pools, {"shared": (10.0, ["A", "B"])})

    def test_pool_raises_when_fixed_concentrations_exceed_site_count(self):
        # Over-budget fixed concentrations are a static constraint violation,
        # rejected at construction rather than surfacing from the solve.
        self.species_a.charge_states[0].fix_concentration(20.0)
        with self.assertRaises(ValueError):
            DefectSystem(
                defect_species=[self.species_a],
                dos=self.dos,
                volume=100,
                temperature=300,
                site_pools={"shared": (10.0, [self.species_a])},
            )

    def test_pooled_species_honours_species_level_fixed_concentration(self):
        self.species_a.fix_concentration(3.0)
        no_pool_system = DefectSystem(
            defect_species=[self.species_a],
            dos=self.dos,
            volume=100,
            temperature=300,
        )
        pooled_system = DefectSystem(
            defect_species=[self.species_a, self.species_b],
            dos=self.dos,
            volume=100,
            temperature=300,
            site_pools={"shared": (10.0, [self.species_a, self.species_b])},
        )
        no_pool_total = sum(
            no_pool_system._global_defect_concs(1.0)[cs]
            for cs in no_pool_system.defect_species[0].charge_states
        )
        pooled_total = sum(
            pooled_system._global_defect_concs(1.0)[cs]
            for cs in pooled_system.defect_species[0].charge_states
        )
        self.assertAlmostEqual(no_pool_total, 3.0, places=8)
        self.assertAlmostEqual(pooled_total, 3.0, places=8)

    def test_pool_raises_when_species_level_fixed_concentration_exceeds_site_count(
        self,
    ):
        self.species_a.fix_concentration(20.0)
        with self.assertRaises(ValueError):
            DefectSystem(
                defect_species=[self.species_a],
                dos=self.dos,
                volume=100,
                temperature=300,
                site_pools={"shared": (10.0, [self.species_a])},
            )

    def test_unpooled_species_raises_when_fixed_concentration_exceeds_nsites(self):
        # An unpooled species' implicit group budget is its own nsites; the
        # error names the species and the budget-versus-occupancy.
        self.species_a.fix_concentration(20.0)  # nsites=5
        with self.assertRaisesRegex(ValueError, r"'A' has 5 .*occupy 20"):
            DefectSystem(
                defect_species=[self.species_a],
                dos=self.dos,
                volume=100,
                temperature=300,
            )

    def test_species_fixed_below_its_fixed_charge_states_raises_at_construction(self):
        # A species fixed below the sum of its own fixed charge states is a
        # static inconsistency, independent of the Fermi level; reject it at
        # construction rather than wrapped as a Fermi-window error mid-solve.
        species = DefectSpecies(
            "F",
            nsites=10,
            charge_states=[DefectChargeState(charge=0, fixed_concentration=5.0)],
            fixed_concentration=1.0,
        )
        with self.assertRaisesRegex(ValueError, r"'F' is fixed at 1.*require 5"):
            DefectSystem(
                defect_species=[species],
                dos=self.dos,
                volume=100,
                temperature=300,
            )

    def test_fixed_concentration_equal_to_sites_is_accepted(self):
        # Exactly saturating the sites (n_free == 0) is valid: the strict ">"
        # budget check must admit it and the solve must succeed.
        species = DefectSpecies(
            "F",
            nsites=2,
            charge_states=[DefectChargeState(charge=0, energy=-0.5, degeneracy=1)],
            fixed_concentration=2.0,
        )
        system = DefectSystem(
            defect_species=[species],
            dos=self.dos,
            volume=100,
            temperature=300,
        )
        total = sum(
            system._global_defect_concs(system.get_sc_fermi()[0])[cs]
            for cs in system.defect_species[0].charge_states
        )
        self.assertAlmostEqual(total, 2.0, places=8)

    def test_pool_raises_when_members_jointly_exceed_site_count(self):
        # Each member fits the pool alone -- A across its two fixed charge
        # states, B in one -- but together they exceed the shared budget.
        self.species_a.charge_states[0].fix_concentration(4.0)
        self.species_a.charge_states[1].fix_concentration(4.0)
        self.species_b.charge_states[0].fix_concentration(4.0)
        with self.assertRaises(ValueError):
            DefectSystem(
                defect_species=[self.species_a, self.species_b],
                dos=self.dos,
                volume=100,
                temperature=300,
                site_pools={"shared": (10.0, [self.species_a, self.species_b])},
            )

    def test_species_fixed_above_all_fixed_charge_states_raises_at_construction(self):
        # A species pinned at a total its all-fixed charge states cannot reach,
        # with no variable charge state to make up the difference, is a static
        # inconsistency; reject it at construction rather than silently
        # under-reporting the total.
        species = DefectSpecies(
            "F",
            nsites=10,
            charge_states=[DefectChargeState(charge=0, fixed_concentration=3.0)],
            fixed_concentration=5.0,
        )
        with self.assertRaisesRegex(
            ValueError, r"'F' is fixed at 5.*sum to 3.*no variable"
        ):
            DefectSystem(
                defect_species=[species],
                dos=self.dos,
                volume=100,
                temperature=300,
            )

    def test_species_fixed_equal_to_all_fixed_charge_states_is_accepted(self):
        # All charge states fixed and summing to the species total -- within
        # floating-point noise, since 0.1 + 0.2 != 0.3 in binary -- is
        # consistent: construct AND solve must both succeed, exercising the
        # isclose tolerance through the whole path, not just construction.
        species = DefectSpecies(
            "F",
            nsites=10,
            charge_states=[
                DefectChargeState(charge=0, fixed_concentration=0.1),
                DefectChargeState(charge=0, fixed_concentration=0.2),
            ],
            fixed_concentration=0.3,
        )
        system = DefectSystem(
            defect_species=[species],
            dos=self.dos,
            volume=100,
            temperature=300,
        )
        total = sum(
            system._global_defect_concs(system.get_sc_fermi()[0])[cs]
            for cs in system.defect_species[0].charge_states
        )
        self.assertAlmostEqual(total, 0.3, places=8)

    def test_species_fixed_above_fixed_charge_states_with_variable_succeeds(self):
        # With a variable charge state to absorb the remainder, a species
        # total above its fixed charge states is fine.
        species = DefectSpecies(
            "F",
            nsites=10,
            charge_states=[
                DefectChargeState(charge=0, fixed_concentration=3.0),
                DefectChargeState(charge=0, energy=-0.5, degeneracy=1),
            ],
            fixed_concentration=5.0,
        )
        system = DefectSystem(
            defect_species=[species],
            dos=self.dos,
            volume=100,
            temperature=300,
        )
        total = sum(
            system._global_defect_concs(system.get_sc_fermi()[0])[cs]
            for cs in system.defect_species[0].charge_states
        )
        self.assertAlmostEqual(total, 5.0, places=8)

    def test_repr_lists_site_pools_by_name(self):
        system = DefectSystem(
            defect_species=[self.species_a, self.species_b],
            dos=self.dos,
            volume=100,
            temperature=300,
            site_pools={"shared": (10.0, [self.species_a, "B"])},
            element_pools={"dE": (1.0, [("A", 2.0)])},
        )
        self.assertIn("shared: 10 sites", repr(system))
        self.assertIn("[A, B]", repr(system))
        self.assertIn("dE: 1 per cell", repr(system))
        self.assertIn("A ×2", repr(system))


class TestDefectSystemSitePercentages(unittest.TestCase):
    def setUp(self):
        self.dos = DOS(
            dos=np.ones(101),
            edos=np.linspace(-5.0, 5.0, 101),
            bandgap=2.0,
            nelect=10,
        )

    def test_near_saturation_species_reports_at_most_100_percent(self):
        # Reproduction: a low-energy single species saturates its sites.
        # The unbounded dilute formula reported ~2.5e10%; the solved
        # concentration caps it at nsites, i.e. <= 100%.
        species = DefectSpecies(
            "R",
            nsites=2,
            charge_states=[DefectChargeState(charge=0, energy=-0.5, degeneracy=1)],
        )
        system = DefectSystem(
            defect_species=[species],
            dos=self.dos,
            volume=100,
            temperature=300,
        )
        # Two-sided: a saturating species fills almost all its sites, so the
        # occupancy is both bounded by 100% (the bug) and near it (not zero).
        pct = system.site_percentages()["R"]
        self.assertLessEqual(pct, 100.0)
        self.assertGreater(pct, 99.0)

    def test_site_percentages_agree_with_concentration_dict(self):
        # The documented contract: site_percentages reports the same solved
        # concentrations as concentration_dict, just as an occupancy fraction.
        # Comparing the two public methods pins that contract without
        # duplicating the internal formula.
        species = DefectSpecies(
            "R",
            nsites=2,
            charge_states=[
                DefectChargeState(charge=0, energy=-0.5, degeneracy=1),
                DefectChargeState(charge=1, energy=0.3, degeneracy=2),
            ],
        )
        system = DefectSystem(
            defect_species=[species],
            dos=self.dos,
            volume=100,
            temperature=300,
        )
        per_cell = system.concentration_dict(decomposed=False, per_volume=False)
        expected = per_cell["R"] / species.nsites * 100
        self.assertAlmostEqual(system.site_percentages()["R"], expected, places=8)

    def test_pooled_species_percentages_bounded_by_pool(self):
        # Two species sharing a site pool: each <= 100% and, since the pool
        # size is the shared denominator, their percentages sum to <= 100%.
        species_a = DefectSpecies(
            "A",
            nsites=5,
            charge_states=[DefectChargeState(charge=0, energy=-0.5, degeneracy=1)],
        )
        species_b = DefectSpecies(
            "B",
            nsites=5,
            charge_states=[DefectChargeState(charge=0, energy=-0.3, degeneracy=1)],
        )
        system = DefectSystem(
            defect_species=[species_a, species_b],
            dos=self.dos,
            volume=100,
            temperature=300,
            site_pools={"shared": (4.0, [species_a, species_b])},
        )
        pct = system.site_percentages()
        self.assertLessEqual(pct["A"], 100.0)
        self.assertLessEqual(pct["B"], 100.0)
        self.assertLessEqual(pct["A"] + pct["B"], 100.0)

    def test_dilute_regime_matches_old_expression(self):
        # A high formation energy is firmly in the dilute limit, where site
        # exclusion (c = n_free * w / (1 + w)) and the old dilute expression
        # (c = nsites * w) agree. Compared as a ratio because the absolute
        # percentage is astronomically small at this energy.
        species = DefectSpecies(
            "D",
            nsites=2,
            charge_states=[DefectChargeState(charge=0, energy=1.5, degeneracy=1)],
        )
        system = DefectSystem(
            defect_species=[species],
            dos=self.dos,
            volume=100,
            temperature=300,
        )
        e_fermi = system.get_sc_fermi()[0]
        old_expression = species.get_concentration(e_fermi, 300) / species.nsites * 100
        new_pct = system.site_percentages()["D"]
        self.assertAlmostEqual(new_pct / old_expression, 1.0, places=6)

    def test_fixed_concentration_within_sites_reports_faithful_ratio(self):
        # A fixed concentration takes a distinct path (_build_group writes it
        # straight into the concentrations); in budget it is reported
        # faithfully as fixed / nsites, e.g. 3 of 4 sites -> 75%.
        species = DefectSpecies(
            "F",
            nsites=4,
            charge_states=[DefectChargeState(charge=0, energy=-0.5, degeneracy=1)],
            fixed_concentration=3.0,
        )
        system = DefectSystem(
            defect_species=[species],
            dos=self.dos,
            volume=100,
            temperature=300,
        )
        self.assertAlmostEqual(system.site_percentages()["F"], 75.0, places=8)

    def test_mixed_pooled_and_unpooled_species_use_their_own_denominators(self):
        # A pooled species is divided by the pool size; an unpooled species in
        # the same system is divided by its own nsites.
        pooled = DefectSpecies(
            "P",
            nsites=5,
            charge_states=[DefectChargeState(charge=0, energy=-0.5, degeneracy=1)],
        )
        unpooled = DefectSpecies(
            "U",
            nsites=3,
            charge_states=[DefectChargeState(charge=0, energy=-0.5, degeneracy=1)],
        )
        system = DefectSystem(
            defect_species=[pooled, unpooled],
            dos=self.dos,
            volume=100,
            temperature=300,
            site_pools={"shared": (4.0, [pooled])},
        )
        per_cell = system.concentration_dict(decomposed=False, per_volume=False)
        pct = system.site_percentages()
        # P's own nsites is 5, but as a pool member its denominator is the
        # pool size 4.0; U falls back to its own nsites 3.
        self.assertAlmostEqual(pct["P"], per_cell["P"] / 4.0 * 100, places=8)
        self.assertAlmostEqual(pct["U"], per_cell["U"] / 3.0 * 100, places=8)
        self.assertLessEqual(pct["P"], 100.0)
        self.assertLessEqual(pct["U"], 100.0)


class TestDefectSystemPoolValidation(unittest.TestCase):
    def setUp(self):
        self.dos = DOS(
            dos=np.ones(101),
            edos=np.linspace(-5.0, 5.0, 101),
            bandgap=2.0,
            nelect=10,
        )
        self.species_a = DefectSpecies(
            "A",
            nsites=5,
            charge_states=[DefectChargeState(charge=0, energy=1.0, degeneracy=1)],
        )
        self.species_b = DefectSpecies(
            "B",
            nsites=1,
            charge_states=[DefectChargeState(charge=0, energy=0.8, degeneracy=1)],
        )

    def _system(self, **kwargs):
        return DefectSystem(
            defect_species=[self.species_a, self.species_b],
            dos=self.dos,
            volume=100,
            temperature=300,
            **kwargs,
        )

    def test_duplicate_roster_names_raise(self):
        twin = DefectSpecies(
            "A",
            nsites=2,
            charge_states=[DefectChargeState(charge=0, energy=1.5, degeneracy=1)],
        )
        with self.assertRaisesRegex(ValueError, "duplicate names: A"):
            DefectSystem(
                defect_species=[self.species_a, twin],
                dos=self.dos,
                volume=100,
                temperature=300,
            )

    def test_site_pool_referencing_unknown_species_raises(self):
        with self.assertRaisesRegex(
            ValueError, "site pool 'shared' references species not in defect_species: C"
        ):
            self._system(site_pools={"shared": (10.0, ["A", "C"])})

    def test_element_pool_referencing_unknown_species_raises(self):
        with self.assertRaisesRegex(
            ValueError, "element pool 'X' references species not in defect_species: C"
        ):
            self._system(element_pools={"X": (1.0, [("C", 1.0)])})

    def test_site_pool_listing_a_species_twice_raises(self):
        with self.assertRaisesRegex(
            ValueError, "site pool 'shared' lists species more than once: A"
        ):
            self._system(site_pools={"shared": (10.0, [self.species_a, "A"])})

    def test_element_pool_listing_a_species_twice_raises(self):
        with self.assertRaisesRegex(
            ValueError, "element pool 'X' lists species more than once: A"
        ):
            self._system(element_pools={"X": (1.0, [("A", 1.0), ("A", -1.0)])})

    def test_species_in_two_site_pools_raises(self):
        with self.assertRaisesRegex(
            ValueError, "species 'A' appears in site pools 'p1' and 'p2'"
        ):
            self._system(site_pools={"p1": (5.0, ["A"]), "p2": (5.0, ["A", "B"])})

    def test_species_may_appear_in_multiple_element_pools(self):
        system = self._system(
            element_pools={
                "X": (0.1, [("A", 1.0)]),
                "Y": (0.1, [("A", 1.0), ("B", 1.0)]),
            }
        )
        self.assertEqual(set(system.element_pools), {"X", "Y"})


class TestDefectSystemElementPools(unittest.TestCase):
    def setUp(self):
        self.dos = DOS(
            dos=np.ones(101),
            edos=np.linspace(-5.0, 5.0, 101),
            bandgap=2.0,
            nelect=10,
        )
        self.species = DefectSpecies(
            "Mg_Zn",
            nsites=10,
            charge_states=[
                DefectChargeState(charge=0, energy=1.0, degeneracy=1),
                DefectChargeState(charge=1, energy=1.3, degeneracy=2),
            ],
        )

    def test_element_pool_drives_total_content_to_target(self):
        target = 5.0
        system = DefectSystem(
            defect_species=[self.species],
            dos=self.dos,
            volume=100,
            temperature=300,
            element_pools={"Mg": (target, [(self.species, 1.0)])},
        )
        concs = system._global_defect_concs(1.0)
        total_mg = sum(concs[cs] for cs in system.defect_species[0].charge_states)
        self.assertAlmostEqual(total_mg, target, places=6)

    def test_element_pool_can_reference_species_by_name(self):
        target = 5.0
        system = DefectSystem(
            defect_species=[self.species],
            dos=self.dos,
            volume=100,
            temperature=300,
            element_pools={"Mg": (target, [("Mg_Zn", 1.0)])},
        )
        concs = system._global_defect_concs(1.0)
        total_mg = sum(concs[cs] for cs in system.defect_species[0].charge_states)
        self.assertAlmostEqual(total_mg, target, places=6)

    def test_element_pool_references_are_normalised_to_names(self):
        o_i = DefectSpecies(
            "O_i",
            nsites=4,
            charge_states=[DefectChargeState(charge=0, energy=1.0, degeneracy=1)],
        )
        v_o = DefectSpecies(
            "V_O",
            nsites=4,
            charge_states=[DefectChargeState(charge=0, energy=1.0, degeneracy=1)],
        )
        # Mixed object/name spelling, two members, asymmetric stoichiometries:
        # pins that normalisation preserves both order and the name-stoichiometry
        # pairing (a swap would invert the balance constraint).
        system = DefectSystem(
            defect_species=[o_i, v_o],
            dos=self.dos,
            volume=100,
            temperature=300,
            element_pools={"dO": (0.0, [(o_i, 2.0), ("V_O", -3.0)])},
        )
        self.assertEqual(
            system.element_pools, {"dO": (0.0, [("O_i", 2.0), ("V_O", -3.0)])}
        )

    def test_element_pool_leaves_fixed_charge_states_unscaled(self):
        fixed_value = 2.0
        self.species.charge_states[0].fix_concentration(fixed_value)
        target = 5.0
        system = DefectSystem(
            defect_species=[self.species],
            dos=self.dos,
            volume=100,
            temperature=300,
            element_pools={"Mg": (target, [(self.species, 1.0)])},
        )
        concs = system._global_defect_concs(1.0)
        self.assertEqual(concs[system.defect_species[0].charge_states[0]], fixed_value)
        total_mg = sum(concs[cs] for cs in system.defect_species[0].charge_states)
        self.assertAlmostEqual(total_mg, target, places=6)

    def test_element_pool_raises_when_fixed_states_exceed_target(self):
        self.species.charge_states[0].fix_concentration(10.0)
        system = DefectSystem(
            defect_species=[self.species],
            dos=self.dos,
            volume=100,
            temperature=300,
            element_pools={"Mg": (5.0, [(self.species, 1.0)])},
        )
        with self.assertRaises(ValueError):
            system._global_defect_concs(1.0)

    def test_element_pool_honours_species_level_fixed_concentration(self):
        fixed_value = 3.0
        species_a = DefectSpecies(
            "A", nsites=10, fixed_concentration=fixed_value,
            charge_states=[DefectChargeState(charge=0, energy=1.0, degeneracy=1)],
        )
        species_b = DefectSpecies(
            "B", nsites=10,
            charge_states=[DefectChargeState(charge=0, energy=1.0, degeneracy=1)],
        )
        target = 10.0
        system = DefectSystem(
            defect_species=[species_a, species_b],
            dos=self.dos,
            volume=100,
            temperature=300,
            element_pools={"X": (target, [(species_a, 1.0), (species_b, 2.0)])},
        )
        concs = system._global_defect_concs(1.0)
        c_a = concs[system.defect_species[0].charge_states[0]]
        c_b = concs[system.defect_species[1].charge_states[0]]
        self.assertEqual(c_a, fixed_value)
        self.assertAlmostEqual(1.0 * c_a + 2.0 * c_b, target, places=6)

    def test_element_pool_raises_when_species_level_fixed_concentration_exceeds_target(
        self,
    ):
        species_a = DefectSpecies(
            "A", nsites=10, fixed_concentration=3.0,
            charge_states=[DefectChargeState(charge=0, energy=1.0, degeneracy=1)],
        )
        species_b = DefectSpecies(
            "B", nsites=10,
            charge_states=[DefectChargeState(charge=0, energy=1.0, degeneracy=1)],
        )
        system = DefectSystem(
            defect_species=[species_a, species_b],
            dos=self.dos,
            volume=100,
            temperature=300,
            element_pools={"X": (2.0, [(species_a, 1.0), (species_b, 2.0)])},
        )
        with self.assertRaises(ValueError):
            system._global_defect_concs(1.0)

    def test_element_pool_preserves_activity_ratio_across_targets(self):
        species_a = DefectSpecies(
            "A", nsites=20,
            charge_states=[DefectChargeState(charge=0, energy=0.0, degeneracy=1)],
        )
        species_b = DefectSpecies(
            "B", nsites=20,
            charge_states=[DefectChargeState(charge=0, energy=0.0, degeneracy=1)],
        )
        ratios = []
        for target in (6.0, 20.0):
            system = DefectSystem(
                defect_species=[species_a, species_b],
                dos=self.dos,
                volume=100,
                temperature=300,
                element_pools={"X": (target, [(species_a, 1.0), (species_b, 2.0)])},
            )
            concs = system._global_defect_concs(1.0)
            c_a = concs[system.defect_species[0].charge_states[0]]
            c_b = concs[system.defect_species[1].charge_states[0]]
            self.assertAlmostEqual(1.0 * c_a + 2.0 * c_b, target, places=6)
            # c_i / (N_i_free - c_i) = w_i * lambda**s_i exactly, so this
            # combination of activities is independent of the target.
            activity_a = c_a / (species_a.nsites - c_a)
            activity_b = c_b / (species_b.nsites - c_b)
            ratios.append(activity_b / activity_a**2)
        self.assertAlmostEqual(ratios[0] / ratios[1], 1.0, delta=1e-12)

    def test_two_coupled_element_pools_hit_both_targets(self):
        cs_mgo = DefectChargeState(charge=0, energy=0.0, degeneracy=1)
        cs_mgi = DefectChargeState(charge=0, energy=0.0, degeneracy=1)
        cs_oi = DefectChargeState(charge=0, energy=0.0, degeneracy=1)
        sp_mgo = DefectSpecies("Mg_O", nsites=20, charge_states=[cs_mgo])
        sp_mgi = DefectSpecies("Mg_i", nsites=20, charge_states=[cs_mgi])
        sp_oi = DefectSpecies("O_i", nsites=20, charge_states=[cs_oi])
        system = DefectSystem(
            defect_species=[sp_mgo, sp_mgi, sp_oi],
            dos=self.dos,
            volume=100,
            temperature=300,
            element_pools={
                "Mg": (8.0, [(sp_mgo, 1.0), (sp_mgi, 1.0)]),
                "O": (5.0, [(sp_mgo, 1.0), (sp_oi, 1.0)]),
            },
        )
        concs = system._global_defect_concs(1.0)
        copied_cs_mgo = system.defect_species[0].charge_states[0]
        copied_cs_mgi = system.defect_species[1].charge_states[0]
        copied_cs_oi = system.defect_species[2].charge_states[0]
        total_mg = concs[copied_cs_mgo] + concs[copied_cs_mgi]
        total_o = concs[copied_cs_mgo] + concs[copied_cs_oi]
        self.assertAlmostEqual(total_mg / 8.0, 1.0, delta=1e-8)
        self.assertAlmostEqual(total_o / 5.0, 1.0, delta=1e-8)

    def test_element_pool_raises_when_target_exceeds_site_capacity(self):
        species = DefectSpecies(
            "Mg_Zn", nsites=2,
            charge_states=[DefectChargeState(charge=0, energy=1.0, degeneracy=1)],
        )
        system = DefectSystem(
            defect_species=[species],
            dos=self.dos,
            volume=100,
            temperature=300,
            element_pools={"Mg": (5.0, [(species, 1.0)])},
        )
        with self.assertRaises(ValueError) as ctx:
            system._global_defect_concs(1.0)
        self.assertIn("Mg", str(ctx.exception))

    def test_element_pool_raises_when_coupled_targets_jointly_infeasible(self):
        """Two pools supplied only by the same species cannot demand
        different contents: each target is achievable on its own, but the
        joint solve has no solution and must raise, naming both pools."""
        species = DefectSpecies(
            "S", nsites=20,
            charge_states=[DefectChargeState(charge=0, energy=0.0, degeneracy=1)],
        )
        system = DefectSystem(
            defect_species=[species],
            dos=self.dos,
            volume=100,
            temperature=300,
            element_pools={
                "X": (5.0, [("S", 1.0)]),
                "Y": (1.0, [("S", 1.0)]),
            },
        )
        with self.assertRaises(ElementPoolError) as ctx:
            system._global_defect_concs(1.0)
        self.assertIn("X", str(ctx.exception))
        self.assertIn("Y", str(ctx.exception))


class TestDefectSystemElementPoolConvergence(unittest.TestCase):
    """Convergence of the element-pool chemical-potential solve across
    concentration regimes: from saturated sites and O(1) occupancies down
    to ~1e-18 defects per unit cell, including underflowed Boltzmann
    weights and zero net-content targets."""

    def setUp(self):
        self.dos = DOS(
            dos=np.ones(101),
            edos=np.linspace(-5.0, 5.0, 101),
            bandgap=2.0,
            nelect=10,
        )

    def coupled_system(self, energy):
        """Two coupled element pools over species P and Q, each with a
        single neutral charge state at `energy`. With targets X = 4u and
        Y = 2u for u = exp(-energy / kT), the exact solution is
        c_P = c_Q = 2u."""
        u = np.exp(-energy / (kboltz * 300.0))
        sp_p = DefectSpecies(
            "P", nsites=1,
            charge_states=[DefectChargeState(charge=0, energy=energy, degeneracy=1)],
        )
        sp_q = DefectSpecies(
            "Q", nsites=1,
            charge_states=[DefectChargeState(charge=0, energy=energy, degeneracy=1)],
        )
        system = DefectSystem(
            defect_species=[sp_p, sp_q],
            dos=self.dos,
            volume=100,
            temperature=300,
            element_pools={
                "X": (4 * u, [("P", 1.0), ("Q", 1.0)]),
                "Y": (2 * u, [("Q", 1.0)]),
            },
        )
        return system, 4 * u, 2 * u

    @staticmethod
    def species_content(system, concs):
        """Total concentration per species, keyed by species name."""
        return {
            sp.name: sum(concs[cs] for cs in sp.charge_states)
            for sp in system.defect_species
        }

    def test_single_element_pool_hits_dilute_target(self):
        """One pool, one element, dilute target (~1e-17 per cell): the
        modal calling pattern."""
        u = np.exp(-1.0 / (kboltz * 300.0))
        target = 2 * u
        sp = DefectSpecies(
            "P", nsites=1,
            charge_states=[DefectChargeState(charge=0, energy=1.0, degeneracy=1)],
        )
        system = DefectSystem(
            defect_species=[sp],
            dos=self.dos,
            volume=100,
            temperature=300,
            element_pools={"X": (target, [("P", 1.0)])},
        )
        concs = system._global_defect_concs(1.0)
        content = self.species_content(system, concs)
        self.assertAlmostEqual(content["P"] / target, 1.0, delta=1e-8)

    def test_coupled_pools_hit_targets_across_scales(self):
        for energy in (0.1, 0.5, 1.0):
            with self.subTest(energy=energy):
                system, target_x, target_y = self.coupled_system(energy=energy)
                concs = system._global_defect_concs(1.0)
                content = self.species_content(system, concs)
                self.assertAlmostEqual(
                    (content["P"] + content["Q"]) / target_x, 1.0, delta=1e-8
                )
                self.assertAlmostEqual(content["Q"] / target_y, 1.0, delta=1e-8)

    def test_coupled_dilute_weights_reach_order_one_targets(self):
        """O(1) targets with deeply dilute unconstrained populations
        (formation energies ~1 eV at 300 K): the solve must bridge ~16
        orders of magnitude between the Boltzmann populations and the
        targets."""
        sp_a = DefectSpecies(
            "A", nsites=20,
            charge_states=[DefectChargeState(charge=0, energy=1.0, degeneracy=1)],
        )
        sp_b = DefectSpecies(
            "B", nsites=20,
            charge_states=[DefectChargeState(charge=0, energy=1.0, degeneracy=1)],
        )
        sp_c = DefectSpecies(
            "C", nsites=20,
            charge_states=[DefectChargeState(charge=0, energy=1.0, degeneracy=1)],
        )
        system = DefectSystem(
            defect_species=[sp_a, sp_b, sp_c],
            dos=self.dos,
            volume=100,
            temperature=300,
            element_pools={
                "Mg": (8.0, [("A", 1.0), ("B", 1.0)]),
                "O": (5.0, [("A", 1.0), ("C", 1.0)]),
            },
        )
        concs = system._global_defect_concs(1.0)
        content = self.species_content(system, concs)
        self.assertAlmostEqual((content["A"] + content["B"]) / 8.0, 1.0, delta=1e-8)
        self.assertAlmostEqual((content["A"] + content["C"]) / 5.0, 1.0, delta=1e-8)

    def test_saturated_species_reaches_target_below_saturation(self):
        """A species with strongly negative formation energy saturates its
        sites at mu = 0; the solve must still reach a feasible target on
        the far side of the saturation plateau."""
        sp = DefectSpecies(
            "S", nsites=10,
            charge_states=[DefectChargeState(charge=0, energy=-2.0, degeneracy=1)],
        )
        system = DefectSystem(
            defect_species=[sp],
            dos=self.dos,
            volume=100,
            temperature=300,
            element_pools={"X": (0.5, [("S", 1.0)])},
        )
        concs = system._global_defect_concs(1.0)
        content = self.species_content(system, concs)
        self.assertAlmostEqual(content["S"] / 0.5, 1.0, delta=1e-8)

    def test_get_sc_fermi_with_charged_states_and_element_pool(self):
        """get_sc_fermi probes Fermi levels across the whole window, where
        charged states' Boltzmann weights swing through saturation; the
        pool solve must converge at every probe."""
        sp = DefectSpecies(
            "Mg", nsites=10,
            charge_states=[
                DefectChargeState(charge=0, energy=1.0, degeneracy=1),
                DefectChargeState(charge=1, energy=1.3, degeneracy=1),
            ],
        )
        target = 1e-8
        system = DefectSystem(
            defect_species=[sp],
            dos=self.dos,
            volume=100,
            temperature=300,
            element_pools={"Mg": (target, [("Mg", 1.0)])},
        )
        e_fermi, _ = system.get_sc_fermi()
        self.assertTrue(np.isfinite(e_fermi))
        concs = system._global_defect_concs(e_fermi)
        content = self.species_content(system, concs)
        self.assertAlmostEqual(content["Mg"] / target, 1.0, delta=1e-8)

    def test_get_sc_fermi_propagates_element_pool_diagnostics(self):
        """An infeasible pool fails identically at every probed Fermi
        level; get_sc_fermi must surface the pool diagnostic rather than
        reporting that no solution exists in the energy window."""
        sp = DefectSpecies(
            "S", nsites=2,
            charge_states=[DefectChargeState(charge=0, energy=1.0, degeneracy=1)],
        )
        system = DefectSystem(
            defect_species=[sp],
            dos=self.dos,
            volume=100,
            temperature=300,
            element_pools={"X": (5.0, [("S", 1.0)])},
        )
        with self.assertRaises(ElementPoolError) as ctx:
            system.get_sc_fermi()
        self.assertIn("Element pool", str(ctx.exception))
        self.assertNotIn("No solution found", str(ctx.exception))

    def test_underflowed_boltzmann_weights_reach_dilute_target(self):
        """E/kT large enough that the unconstrained populations underflow
        to exactly zero (2 eV at 30 K): the solve must still bridge to the
        target."""
        sp = DefectSpecies(
            "S", nsites=10,
            charge_states=[DefectChargeState(charge=0, energy=2.0, degeneracy=1)],
        )
        target = 1e-10
        system = DefectSystem(
            defect_species=[sp],
            dos=self.dos,
            volume=100,
            temperature=30,
            element_pools={"X": (target, [("S", 1.0)])},
        )
        concs = system._global_defect_concs(1.0)
        content = self.species_content(system, concs)
        self.assertAlmostEqual(content["S"] / target, 1.0, delta=1e-8)

    def test_coupled_shared_species_across_formation_energies_and_targets(self):
        """Coupled pools sharing a species, with 1.5-2.2 eV formation
        energies and targets from O(10) down to O(0.01): every combination
        is feasible and must be solved."""
        for energy in (1.5, 1.8, 2.2):
            for target_x, target_y in (
                (8.0, 5.0), (2.0, 0.5), (15.0, 1.0), (0.1, 0.05), (5.0, 5.0)
            ):
                with self.subTest(energy=energy, targets=(target_x, target_y)):
                    sp = {
                        name: DefectSpecies(
                            name, nsites=20,
                            charge_states=[
                                DefectChargeState(charge=0, energy=energy, degeneracy=1)
                            ],
                        )
                        for name in ("A", "B", "C")
                    }
                    system = DefectSystem(
                        defect_species=list(sp.values()),
                        dos=self.dos,
                        volume=100,
                        temperature=300,
                        element_pools={
                            "X": (target_x, [("A", 1.0), ("B", 1.0)]),
                            "Y": (target_y, [("A", 1.0), ("C", 1.0)]),
                        },
                    )
                    concs = system._global_defect_concs(1.0)
                    content = self.species_content(system, concs)
                    self.assertAlmostEqual(
                        (content["A"] + content["B"]) / target_x, 1.0, delta=1e-8
                    )
                    self.assertAlmostEqual(
                        (content["A"] + content["C"]) / target_y, 1.0, delta=1e-8
                    )

    def test_zero_remaining_target_zeroes_variable_states_and_solves_rest(self):
        """An element whose target is fully committed by fixed
        concentrations admits no further variable content: variable states
        of species with positive stoichiometry in that element get
        concentration 0, and the other elements solve without them."""
        target_y = 1e-9
        sp_a = DefectSpecies(
            "A", nsites=10, fixed_concentration=3.0,
            charge_states=[DefectChargeState(charge=0, energy=1.0, degeneracy=1)],
        )
        sp_b = DefectSpecies(
            "B", nsites=10,
            charge_states=[DefectChargeState(charge=0, energy=1.0, degeneracy=1)],
        )
        sp_c = DefectSpecies(
            "C", nsites=10,
            charge_states=[DefectChargeState(charge=0, energy=1.0, degeneracy=1)],
        )
        system = DefectSystem(
            defect_species=[sp_a, sp_b, sp_c],
            dos=self.dos,
            volume=100,
            temperature=300,
            element_pools={
                "X": (3.0, [("A", 1.0), ("B", 1.0)]),
                "Y": (target_y, [("B", 1.0), ("C", 1.0)]),
            },
        )
        concs = system._global_defect_concs(1.0)
        content = self.species_content(system, concs)
        self.assertEqual(content["B"], 0.0)
        self.assertEqual(content["A"], 3.0)
        self.assertAlmostEqual(content["C"] / target_y, 1.0, delta=1e-8)

    def test_three_coupled_pools_with_chained_sharing(self):
        """Three elements coupled in a chain (A supplies X and Y, C
        supplies Y and Z) at dilute scales."""
        u = np.exp(-1.0 / (kboltz * 300.0))
        targets = {"X": 3 * u, "Y": 4 * u, "Z": 2 * u}
        sp = {
            name: DefectSpecies(
                name, nsites=1,
                charge_states=[DefectChargeState(charge=0, energy=1.0, degeneracy=1)],
            )
            for name in ("A", "B", "C", "D")
        }
        system = DefectSystem(
            defect_species=list(sp.values()),
            dos=self.dos,
            volume=100,
            temperature=300,
            element_pools={
                "X": (targets["X"], [("A", 1.0), ("B", 1.0)]),
                "Y": (targets["Y"], [("A", 1.0), ("C", 1.0)]),
                "Z": (targets["Z"], [("C", 1.0), ("D", 1.0)]),
            },
        )
        concs = system._global_defect_concs(1.0)
        content = self.species_content(system, concs)
        self.assertAlmostEqual(
            (content["A"] + content["B"]) / targets["X"], 1.0, delta=1e-8
        )
        self.assertAlmostEqual(
            (content["A"] + content["C"]) / targets["Y"], 1.0, delta=1e-8
        )
        self.assertAlmostEqual(
            (content["C"] + content["D"]) / targets["Z"], 1.0, delta=1e-8
        )

    def test_dilute_targets_with_mixed_stoichiometries(self):
        """Stoichiometric coefficients of different magnitude (2.0 and
        0.5) set different characteristic scales per element."""
        u = np.exp(-1.0 / (kboltz * 300.0))
        target_x, target_y = 3 * u, 0.5 * u
        sp_p = DefectSpecies(
            "P", nsites=1,
            charge_states=[DefectChargeState(charge=0, energy=1.0, degeneracy=1)],
        )
        sp_q = DefectSpecies(
            "Q", nsites=1,
            charge_states=[DefectChargeState(charge=0, energy=1.0, degeneracy=1)],
        )
        system = DefectSystem(
            defect_species=[sp_p, sp_q],
            dos=self.dos,
            volume=100,
            temperature=300,
            element_pools={
                "X": (target_x, [("P", 2.0), ("Q", 1.0)]),
                "Y": (target_y, [("Q", 0.5)]),
            },
        )
        concs = system._global_defect_concs(1.0)
        content = self.species_content(system, concs)
        self.assertAlmostEqual(
            (2.0 * content["P"] + content["Q"]) / target_x, 1.0, delta=1e-8
        )
        self.assertAlmostEqual(0.5 * content["Q"] / target_y, 1.0, delta=1e-8)

    def test_feasible_grid_across_energies_stoichiometries_and_depths(self):
        """Systematic sweep of the axes that set solver difficulty --
        Boltzmann weight scale (formation energy), stoichiometry magnitude,
        and target depth as a fraction of capacity -- asserting the solve
        contract on every combination. All systems are jointly feasible by
        construction: each element can be supplied by its dedicated
        species alone."""
        for energy in (-1.0, 0.0, 1.0, 2.0):
            for stoich_b in (0.5, 1.0, 2.0):
                for fraction in (1e-12, 1e-4, 0.3, 0.9):
                    with self.subTest(
                        energy=energy, stoich=stoich_b, fraction=fraction
                    ):
                        sp = {
                            name: DefectSpecies(
                                name, nsites=20,
                                charge_states=[
                                    DefectChargeState(
                                        charge=0, energy=energy, degeneracy=1
                                    )
                                ],
                            )
                            for name in ("A", "B", "C")
                        }
                        target_x = fraction * 20 * stoich_b
                        target_y = fraction * 20
                        system = DefectSystem(
                            defect_species=list(sp.values()),
                            dos=self.dos,
                            volume=100,
                            temperature=300,
                            element_pools={
                                "X": (target_x, [("A", 1.0), ("B", stoich_b)]),
                                "Y": (target_y, [("A", 1.0), ("C", 1.0)]),
                            },
                        )
                        concs = system._global_defect_concs(1.0)
                        content = self.species_content(system, concs)
                        self.assertAlmostEqual(
                            (content["A"] + stoich_b * content["B"]) / target_x,
                            1.0,
                            delta=1e-6,
                        )
                        self.assertAlmostEqual(
                            (content["A"] + content["C"]) / target_y,
                            1.0,
                            delta=1e-6,
                        )

    def test_all_pools_exhausted_zeroes_pooled_and_frees_unpooled_species(self):
        """With every element fully committed by fixed concentrations there
        is nothing left to solve: pooled variable states get concentration
        0 and species outside any pool keep their Boltzmann populations."""
        sp_a = DefectSpecies(
            "A", nsites=10, fixed_concentration=3.0,
            charge_states=[DefectChargeState(charge=0, energy=1.0, degeneracy=1)],
        )
        sp_b = DefectSpecies(
            "B", nsites=10,
            charge_states=[DefectChargeState(charge=0, energy=1.0, degeneracy=1)],
        )
        sp_c = DefectSpecies(
            "C", nsites=10,
            charge_states=[DefectChargeState(charge=0, energy=1.0, degeneracy=1)],
        )
        system = DefectSystem(
            defect_species=[sp_a, sp_b, sp_c],
            dos=self.dos,
            volume=100,
            temperature=300,
            element_pools={"X": (3.0, [("A", 1.0), ("B", 1.0)])},
        )
        concs = system._global_defect_concs(1.0)
        content = self.species_content(system, concs)
        self.assertEqual(content["B"], 0.0)
        w = np.exp(-1.0 / (kboltz * 300.0))
        self.assertAlmostEqual(content["C"] / (10 * w / (1 + w)), 1.0, delta=1e-10)

    def test_exhausted_element_starving_another_pool_raises(self):
        """A species zeroed by an exhausted element cannot supply another
        pool: if it was that pool's only variable supplier, the now
        unreachable target must raise, naming the starved element."""
        sp_a = DefectSpecies(
            "A", nsites=10, fixed_concentration=3.0,
            charge_states=[DefectChargeState(charge=0, energy=1.0, degeneracy=1)],
        )
        sp_b = DefectSpecies(
            "B", nsites=10,
            charge_states=[DefectChargeState(charge=0, energy=1.0, degeneracy=1)],
        )
        system = DefectSystem(
            defect_species=[sp_a, sp_b],
            dos=self.dos,
            volume=100,
            temperature=300,
            element_pools={
                "X": (3.0, [("A", 1.0), ("B", 1.0)]),
                "Y": (1.0, [("B", 1.0)]),
            },
        )
        with self.assertRaises(ValueError) as ctx:
            system._global_defect_concs(1.0)
        self.assertIn("Y", str(ctx.exception))

    def test_exhausted_pool_tolerates_float_rounding_in_committed_total(self):
        """A pool target met exactly by fixed concentrations must be
        treated as exhausted even when float summation of the
        contributions does not land on the target exactly
        (0.1 + 0.2 != 0.3)."""
        sp_a = DefectSpecies(
            "A", nsites=10, fixed_concentration=0.1,
            charge_states=[DefectChargeState(charge=0, energy=1.0, degeneracy=1)],
        )
        sp_b = DefectSpecies(
            "B", nsites=10, fixed_concentration=0.2,
            charge_states=[DefectChargeState(charge=0, energy=1.0, degeneracy=1)],
        )
        sp_c = DefectSpecies(
            "C", nsites=10,
            charge_states=[DefectChargeState(charge=0, energy=1.0, degeneracy=1)],
        )
        system = DefectSystem(
            defect_species=[sp_a, sp_b, sp_c],
            dos=self.dos,
            volume=100,
            temperature=300,
            element_pools={"X": (0.3, [("A", 1.0), ("B", 1.0), ("C", 1.0)])},
        )
        concs = system._global_defect_concs(1.0)
        content = self.species_content(system, concs)
        self.assertEqual(content["C"], 0.0)

    def test_exhausted_pool_tolerates_rounding_with_negative_committed(self):
        """The rounding clamp must hold when the committed total is
        negative (fixed negative-stoichiometry states): a one-ulp negative
        remainder against an exactly-met negative target is exhaustion, not
        a constraint inconsistency."""
        sp_oi = DefectSpecies(
            "O_i", nsites=10,
            charge_states=[DefectChargeState(charge=0, energy=0.9, degeneracy=1)],
        )
        sp_vo = DefectSpecies(
            "V_O", nsites=10, fixed_concentration=0.3,
            charge_states=[DefectChargeState(charge=0, energy=0.9, degeneracy=1)],
        )
        system = DefectSystem(
            defect_species=[sp_oi, sp_vo],
            dos=self.dos,
            volume=100,
            temperature=300,
            element_pools={"dO": (-(0.1 + 0.2), [("O_i", 1.0), ("V_O", -1.0)])},
        )
        concs = system._global_defect_concs(1.0)
        content = self.species_content(system, concs)
        self.assertEqual(content["O_i"], 0.0)

    def test_small_net_target_at_half_saturation_converges(self):
        """A small net target over species near half-saturation needs the
        content resolved far finer than the line search's step floor; the
        solve must still reach the guard tolerance rather than stalling
        above it and reporting spurious infeasibility."""
        for target in (-1e-6, -1e-7, -1e-9):
            with self.subTest(target=target):
                sp_oi = DefectSpecies(
                    "O_i", nsites=10,
                    charge_states=[
                        DefectChargeState(charge=0, energy=0.0, degeneracy=1)
                    ],
                )
                sp_vo = DefectSpecies(
                    "V_O", nsites=10,
                    charge_states=[
                        DefectChargeState(charge=0, energy=0.0, degeneracy=1)
                    ],
                )
                system = DefectSystem(
                    defect_species=[sp_oi, sp_vo],
                    dos=self.dos,
                    volume=100,
                    temperature=300,
                    element_pools={"dO": (target, [("O_i", 1.0), ("V_O", -1.0)])},
                )
                concs = system._global_defect_concs(1.0)
                content = self.species_content(system, concs)
                self.assertAlmostEqual(
                    (content["O_i"] - content["V_O"]) / target, 1.0, delta=1e-6
                )

    def test_negative_net_content_target_balances_species(self):
        """A negative net-content target over mixed-sign stoichiometries
        (oxygen deficiency in an off-stoichiometry scan) is a balance
        tipped towards the negative-stoichiometry species."""
        target = -1e-7
        sp_oi = DefectSpecies(
            "O_i", nsites=10,
            charge_states=[DefectChargeState(charge=0, energy=0.9, degeneracy=1)],
        )
        sp_vo = DefectSpecies(
            "V_O", nsites=10,
            charge_states=[DefectChargeState(charge=0, energy=0.9, degeneracy=1)],
        )
        system = DefectSystem(
            defect_species=[sp_oi, sp_vo],
            dos=self.dos,
            volume=100,
            temperature=300,
            element_pools={"dO": (target, [("O_i", 1.0), ("V_O", -1.0)])},
        )
        concs = system._global_defect_concs(1.0)
        content = self.species_content(system, concs)
        self.assertGreater(content["O_i"], 0.0)
        self.assertAlmostEqual(
            (content["O_i"] - content["V_O"]) / target, 1.0, delta=1e-8
        )

    def test_negative_target_below_negative_capacity_raises(self):
        """A negative net-content target deeper than the
        negative-stoichiometry species can supply must raise the
        feasibility error, not a constraint-inconsistency error."""
        sp_oi = DefectSpecies(
            "O_i", nsites=10,
            charge_states=[DefectChargeState(charge=0, energy=0.9, degeneracy=1)],
        )
        sp_vo = DefectSpecies(
            "V_O", nsites=10,
            charge_states=[DefectChargeState(charge=0, energy=0.9, degeneracy=1)],
        )
        system = DefectSystem(
            defect_species=[sp_oi, sp_vo],
            dos=self.dos,
            volume=100,
            temperature=300,
            element_pools={"dO": (-20.0, [("O_i", 1.0), ("V_O", -1.0)])},
        )
        with self.assertRaises(ElementPoolError) as ctx:
            system._global_defect_concs(1.0)
        self.assertIn("minimum", str(ctx.exception))

    def test_mixed_zero_and_dilute_targets_in_one_solve(self):
        """A zero-target balance pool whose positive species also supplies
        a dilute dopant pool: both constraints must hold in one solve."""
        target_y = 1e-9
        sp_oi = DefectSpecies(
            "O_i", nsites=10,
            charge_states=[DefectChargeState(charge=0, energy=0.9, degeneracy=1)],
        )
        sp_vo = DefectSpecies(
            "V_O", nsites=10,
            charge_states=[DefectChargeState(charge=0, energy=0.9, degeneracy=1)],
        )
        sp_m = DefectSpecies(
            "M", nsites=10,
            charge_states=[DefectChargeState(charge=0, energy=0.9, degeneracy=1)],
        )
        system = DefectSystem(
            defect_species=[sp_oi, sp_vo, sp_m],
            dos=self.dos,
            volume=100,
            temperature=300,
            element_pools={
                "dO": (0.0, [("O_i", 1.0), ("V_O", -1.0)]),
                "Y": (target_y, [("O_i", 1.0), ("M", 1.0)]),
            },
        )
        concs = system._global_defect_concs(1.0)
        content = self.species_content(system, concs)
        gross = content["O_i"] + content["V_O"]
        self.assertLessEqual(abs(content["O_i"] - content["V_O"]), 1e-6 * gross)
        self.assertAlmostEqual(
            (content["O_i"] + content["M"]) / target_y, 1.0, delta=1e-8
        )

    def test_zero_target_with_negative_stoichiometry_balances_species(self):
        """A pool with target zero over species of opposite stoichiometry
        (the natural encoding of exact stoichiometry in an off-stoichiometry
        scan) must balance the species against each other, not zero one and
        leave the other unconstrained."""
        sp_oi = DefectSpecies(
            "O_i", nsites=10,
            charge_states=[DefectChargeState(charge=0, energy=0.5, degeneracy=1)],
        )
        sp_vo = DefectSpecies(
            "V_O", nsites=10,
            charge_states=[DefectChargeState(charge=0, energy=0.7, degeneracy=1)],
        )
        system = DefectSystem(
            defect_species=[sp_oi, sp_vo],
            dos=self.dos,
            volume=100,
            temperature=300,
            element_pools={"dO": (0.0, [("O_i", 1.0), ("V_O", -1.0)])},
        )
        concs = system._global_defect_concs(1.0)
        content = self.species_content(system, concs)
        self.assertGreater(content["O_i"], 0.0)
        self.assertAlmostEqual(content["O_i"] / content["V_O"], 1.0, delta=1e-8)

    def test_mixed_scale_pools_hit_both_targets(self):
        target_x, target_y = 1e-6, 1e-18
        sp_p = DefectSpecies(
            "P", nsites=1,
            charge_states=[DefectChargeState(charge=0, energy=0.5, degeneracy=1)],
        )
        sp_q = DefectSpecies(
            "Q", nsites=1,
            charge_states=[DefectChargeState(charge=0, energy=0.5, degeneracy=1)],
        )
        system = DefectSystem(
            defect_species=[sp_p, sp_q],
            dos=self.dos,
            volume=100,
            temperature=300,
            element_pools={
                "X": (target_x, [("P", 1.0), ("Q", 1.0)]),
                "Y": (target_y, [("Q", 1.0)]),
            },
        )
        concs = system._global_defect_concs(1.0)
        content = self.species_content(system, concs)
        self.assertAlmostEqual(
            (content["P"] + content["Q"]) / target_x, 1.0, delta=1e-8
        )
        self.assertAlmostEqual(content["Q"] / target_y, 1.0, delta=1e-8)


class TestDefectSystemBandEdgeCorrections(unittest.TestCase):
    def setUp(self):
        self.dos = DOS(
            dos=np.ones(101),
            edos=np.linspace(-5.0, 5.0, 101),
            bandgap=2.0,
            nelect=10,
        )
        self.cs0 = DefectChargeState(charge=0, energy=1.0, degeneracy=1)
        self.cs1 = DefectChargeState(charge=1, energy=1.5, degeneracy=2)
        self.species = DefectSpecies("V_O", nsites=1, charge_states=[self.cs0, self.cs1])

    def _make_system(self, **kwargs):
        return DefectSystem(
            defect_species=[self.species],
            dos=self.dos,
            volume=100,
            temperature=300,
            **kwargs,
        )

    def test_no_shifts_leaves_energies_and_bandgap_unchanged(self):
        system = self._make_system()
        self.assertEqual(system.defect_species[0].charge_states[0].energy, 1.0)
        self.assertEqual(system.defect_species[0].charge_states[1].energy, 1.5)
        self.assertIn("2 eV", repr(system))

    def test_rigid_shift_changes_displayed_bandgap_but_not_formation_energies(self):
        system = self._make_system(vbm_shift=0.05, cbm_shift=-0.02)
        self.assertEqual(system.defect_species[0].charge_states[0].energy, 1.0)
        self.assertEqual(system.defect_species[0].charge_states[1].energy, 1.5)
        self.assertIn("1.93 eV", repr(system))
        self.assertNotIn("→", repr(system))
        self.assertEqual(self.cs0.energy, 1.0)
        self.assertEqual(self.cs1.energy, 1.5)

    def test_non_rigid_shift_applies_charge_dependent_correction(self):
        d_vbm = 0.05
        system = self._make_system(vbm_shift=d_vbm, rigid_shift=False)
        self.assertAlmostEqual(
            system.defect_species[0].charge_states[0].energy,
            1.0 - self.cs0.charge * d_vbm,
        )
        self.assertAlmostEqual(
            system.defect_species[0].charge_states[1].energy,
            1.5 - self.cs1.charge * d_vbm,
        )
        self.assertEqual(self.cs0.energy, 1.0)
        self.assertEqual(self.cs1.energy, 1.5)

    def test_formation_energy_correction_distinguishes_metastable_states(self):
        cs_a = DefectChargeState(charge=1, energy=0.5, degeneracy=1)
        cs_b = DefectChargeState(charge=1, energy=0.9, degeneracy=1)
        species = DefectSpecies("X_i", nsites=1, charge_states=[cs_a, cs_b])
        system = DefectSystem(
            defect_species=[species],
            dos=self.dos,
            volume=100,
            temperature=300,
            formation_energy_corrections={cs_a: 0.1, cs_b: -0.05},
        )
        self.assertAlmostEqual(system.defect_species[0].charge_states[0].energy, 0.6)
        self.assertAlmostEqual(system.defect_species[0].charge_states[1].energy, 0.85)
        self.assertEqual(cs_a.energy, 0.5)
        self.assertEqual(cs_b.energy, 0.9)

    def test_formation_energy_correction_takes_priority_over_rigid_shift(self):
        system = self._make_system(
            vbm_shift=0.05,
            rigid_shift=False,
            formation_energy_corrections={self.cs1: 0.2},
        )
        self.assertAlmostEqual(
            system.defect_species[0].charge_states[1].energy, 1.5 + 0.2
        )
        self.assertAlmostEqual(system.defect_species[0].charge_states[0].energy, 1.0)
        self.assertEqual(self.cs0.energy, 1.0)
        self.assertEqual(self.cs1.energy, 1.5)

    def test_formation_energy_correction_with_unrecognised_charge_state_raises(self):
        unknown_cs = DefectChargeState(charge=2, energy=3.0, degeneracy=1)
        with self.assertRaises(ValueError):
            self._make_system(formation_energy_corrections={unknown_cs: 0.1})


class TestDefectSystemFactory(unittest.TestCase):
    def setUp(self):
        self.dos = DOS(
            dos=np.ones(101),
            edos=np.linspace(-5.0, 5.0, 101),
            bandgap=2.0,
            nelect=10,
        )
        self.cs0 = DefectChargeState(charge=0, energy=1.0, degeneracy=1)
        self.cs1 = DefectChargeState(charge=1, energy=1.5, degeneracy=2)
        self.species = DefectSpecies("V_O", nsites=1, charge_states=[self.cs0, self.cs1])

    def test_at_evaluates_shift_functions_at_given_temperature(self):
        factory = DefectSystemFactory(
            defect_species=[self.species],
            dos=self.dos,
            volume=100,
            vbm_shift_fn=lambda T: 0.001 * T,
            rigid_shift=False,
        )
        system = factory.at(300)
        self.assertAlmostEqual(system.vbm_shift, 0.3)
        self.assertAlmostEqual(
            system.defect_species[0].charge_states[1].energy,
            1.5 - self.cs1.charge * 0.3,
        )

    def test_at_sets_temperature(self):
        factory = DefectSystemFactory(defect_species=[self.species], dos=self.dos, volume=100)
        system = factory.at(300)
        self.assertEqual(system.temperature, 300)

    def test_at_supports_per_call_overrides(self):
        factory = DefectSystemFactory(defect_species=[self.species], dos=self.dos, volume=100)
        system = factory.at(300, convergence_tolerance=1e-12)
        self.assertEqual(system.convergence_tolerance, 1e-12)

    def test_at_with_formation_energy_correction_fns_per_charge_state(self):
        cs_a = DefectChargeState(charge=1, energy=0.5, degeneracy=1)
        cs_b = DefectChargeState(charge=1, energy=0.9, degeneracy=1)
        species = DefectSpecies("X_i", nsites=1, charge_states=[cs_a, cs_b])
        factory = DefectSystemFactory(
            defect_species=[species],
            dos=self.dos,
            volume=100,
            formation_energy_correction_fns={
                cs_a: lambda T: 0.001 * T,
                cs_b: lambda T: -0.0005 * T,
            },
        )
        system = factory.at(300)
        self.assertAlmostEqual(system.defect_species[0].charge_states[0].energy, 0.5 + 0.3)
        self.assertAlmostEqual(system.defect_species[0].charge_states[1].energy, 0.9 - 0.15)

    def test_repeated_at_calls_are_independent(self):
        factory = DefectSystemFactory(
            defect_species=[self.species],
            dos=self.dos,
            volume=100,
            vbm_shift_fn=lambda T: 0.001 * T,
            rigid_shift=False,
        )
        low = factory.at(200)
        high = factory.at(400)
        self.assertIsNot(
            low.defect_species[0].charge_states[1],
            high.defect_species[0].charge_states[1],
        )
        self.assertNotAlmostEqual(
            low.defect_species[0].charge_states[1].energy,
            high.defect_species[0].charge_states[1].energy,
        )


class TestDefectSystemReport(unittest.TestCase):
    def setUp(self):
        self.dos = DOS(
            dos=np.ones(101),
            edos=np.linspace(-5.0, 5.0, 101),
            bandgap=2.0,
            nelect=10,
        )
        self.species = DefectSpecies(
            "V_O",
            nsites=1,
            charge_states=[
                DefectChargeState(charge=0, energy=1.0, degeneracy=1),
                DefectChargeState(charge=1, energy=1.5, degeneracy=2),
            ],
        )

    def test_report_contains_expected_sections(self):
        system = DefectSystem(
            defect_species=[self.species], dos=self.dos, volume=100, temperature=300
        )
        with patch("sys.stdout", new=StringIO()):
            output = system.report()
        for snippet in (
            "DefectSystem",
            "SC Fermi energy:",
            "Carriers:",
            "Defect concentrations:",
            "V_O",
        ):
            self.assertIn(snippet, output)

    def test_report_shows_single_corrected_bandgap_value(self):
        system = DefectSystem(
            defect_species=[self.species],
            dos=self.dos,
            volume=100,
            temperature=300,
            vbm_shift=0.05,
            cbm_shift=-0.02,
        )
        with patch("sys.stdout", new=StringIO()):
            output = system.report()
        self.assertIn("1.93", output)
        self.assertNotIn("→", output)
        self.assertEqual(system.dos._bandgap, 2.0)


class TestDefectSystemPoolSerialisation(unittest.TestCase):
    def setUp(self):
        self.dos = DOS(
            dos=np.ones(101),
            edos=np.linspace(-5.0, 5.0, 101),
            bandgap=2.0,
            nelect=10,
        )
        self.species_a = DefectSpecies(
            "A",
            nsites=5,
            charge_states=[
                DefectChargeState(charge=0, energy=1.0, degeneracy=1),
                DefectChargeState(charge=1, energy=1.2, degeneracy=2),
            ],
        )
        self.species_b = DefectSpecies(
            "B",
            nsites=1,
            charge_states=[
                DefectChargeState(charge=0, energy=0.8, degeneracy=1),
                DefectChargeState(charge=-1, energy=1.5, degeneracy=2),
            ],
        )
        self.species_c = DefectSpecies(
            "C",
            nsites=3,
            charge_states=[
                DefectChargeState(charge=0, energy=0.5, degeneracy=1),
                DefectChargeState(charge=1, energy=0.9, degeneracy=2),
            ],
        )

    def _system(self):
        # site pool mixes object (A) and name ("B"); element pool references C
        # by object -- exercises both reference spellings through normalisation.
        return DefectSystem(
            defect_species=[self.species_a, self.species_b, self.species_c],
            dos=self.dos,
            volume=100,
            temperature=300,
            site_pools={"cation": (4.0, [self.species_a, "B"])},
            element_pools={"X": (0.3, [(self.species_c, 1.0)])},
        )

    @staticmethod
    def _concs_by_name(system, e_fermi):
        concs = system._global_defect_concs(e_fermi)
        return {
            ds.name: sum(concs.get(cs, 0.0) for cs in ds.charge_states)
            for ds in system.defect_species
        }

    def test_round_trip_through_json_preserves_pools_and_physics(self):
        system = self._system()
        as_dict = system.as_dict()
        self.assertEqual(
            as_dict["site_pools"],
            {"cation": {"n_sites": 4.0, "species": ["A", "B"]}},
        )
        self.assertEqual(
            as_dict["element_pools"],
            {"X": {"target": 0.3, "members": [{"species": "C", "stoichiometry": 1.0}]}},
        )
        reloaded = DefectSystem.from_dict(json.loads(json.dumps(as_dict)))
        self.assertEqual(reloaded.site_pools, {"cation": (4.0, ["A", "B"])})
        self.assertEqual(reloaded.element_pools, {"X": (0.3, [("C", 1.0)])})
        original = self._concs_by_name(system, 1.0)
        reloaded_totals = self._concs_by_name(reloaded, 1.0)
        self.assertEqual(set(original), set(reloaded_totals))
        for name, total in original.items():
            self.assertAlmostEqual(reloaded_totals[name], total, places=8)

    def test_round_trip_through_yaml_preserves_pools_and_physics(self):
        system = self._system()
        reloaded = DefectSystem.from_dict(
            yaml.safe_load(yaml.safe_dump(system.as_dict()))
        )
        self.assertEqual(reloaded.site_pools, {"cation": (4.0, ["A", "B"])})
        self.assertEqual(reloaded.element_pools, {"X": (0.3, [("C", 1.0)])})
        original = self._concs_by_name(system, 1.0)
        reloaded_totals = self._concs_by_name(reloaded, 1.0)
        for name, total in original.items():
            self.assertAlmostEqual(reloaded_totals[name], total, places=8)

    def test_system_without_pools_emits_neither_key(self):
        system = DefectSystem(
            defect_species=[self.species_a],
            dos=self.dos,
            volume=100,
            temperature=300,
        )
        as_dict = system.as_dict()
        self.assertNotIn("site_pools", as_dict)
        self.assertNotIn("element_pools", as_dict)

    def test_dict_without_pool_keys_still_loads(self):
        system = DefectSystem(
            defect_species=[self.species_a],
            dos=self.dos,
            volume=100,
            temperature=300,
        )
        as_dict = system.as_dict()
        as_dict.pop("site_pools", None)
        as_dict.pop("element_pools", None)
        reloaded = DefectSystem.from_dict(as_dict)
        self.assertEqual(reloaded.site_pools, {})
        self.assertEqual(reloaded.element_pools, {})

    def test_corrections_round_trip_via_baked_energies(self):
        donor_states = [
            DefectChargeState(charge=0, energy=1.0, degeneracy=1),
            DefectChargeState(charge=1, energy=0.6, degeneracy=2),
        ]
        acceptor_states = [
            DefectChargeState(charge=0, energy=1.0, degeneracy=1),
            DefectChargeState(charge=-1, energy=0.6, degeneracy=2),
        ]
        donor = DefectSpecies("D", nsites=2, charge_states=donor_states)
        acceptor = DefectSpecies("Acc", nsites=2, charge_states=acceptor_states)
        system = DefectSystem(
            defect_species=[donor, acceptor],
            dos=self.dos,
            volume=100,
            temperature=300,
            vbm_shift=0.1,
            cbm_shift=-0.05,
            formation_energy_corrections={donor_states[1]: 0.2},
            rigid_shift=False,
        )
        reloaded = DefectSystem.from_dict(json.loads(json.dumps(system.as_dict())))
        e_original = system.get_sc_fermi()[0]
        e_reloaded = reloaded.get_sc_fermi()[0]
        self.assertAlmostEqual(e_original, e_reloaded, places=6)
        original = self._concs_by_name(system, e_original)
        reloaded_totals = self._concs_by_name(reloaded, e_reloaded)
        for name, total in original.items():
            self.assertAlmostEqual(reloaded_totals[name], total, places=8)

    def test_multi_member_pools_round_trip_faithfully(self):
        # A multi-member element pool with asymmetric, opposite-sign
        # stoichiometries, and a species (A) that belongs to both a site pool
        # and an element pool: pins that reconstruction preserves member order,
        # the name-stoichiometry pairing, and the sign, with no cross-pool mix-up.
        system = DefectSystem(
            defect_species=[self.species_a, self.species_c],
            dos=self.dos,
            volume=100,
            temperature=300,
            site_pools={"shared": (6.0, [self.species_a, "C"])},
            element_pools={"dX": (0.0, [(self.species_a, 1.0), ("C", -2.0)])},
        )
        reloaded = DefectSystem.from_dict(json.loads(json.dumps(system.as_dict())))
        self.assertEqual(reloaded.site_pools, {"shared": (6.0, ["A", "C"])})
        self.assertEqual(
            reloaded.element_pools, {"dX": (0.0, [("A", 1.0), ("C", -2.0)])}
        )

    def test_numpy_convergence_tolerance_is_yaml_safe(self):
        system = DefectSystem(
            defect_species=[self.species_a],
            dos=self.dos,
            volume=100,
            temperature=300,
            convergence_tolerance=np.float64(1e-8),
        )
        as_dict = system.as_dict()
        self.assertIs(type(as_dict["convergence_tolerance"]), float)
        yaml.safe_dump(as_dict)

    def test_baked_numpy_correction_keeps_as_dict_yaml_safe(self):
        # The documented temperature-dependent-shift workflow: a numpy vbm_shift
        # with rigid_shift=False bakes -charge * vbm_shift into each energy,
        # promoting it to np.float64. The whole-system dict must stay YAML-safe.
        system = DefectSystem(
            defect_species=[self.species_a, self.species_b],
            dos=self.dos,
            volume=100,
            temperature=300,
            vbm_shift=np.float64(0.1),
            rigid_shift=False,
        )
        yaml.safe_dump(system.as_dict())

    def test_fully_numpy_system_as_dict_is_yaml_safe(self):
        # Contract guard at the boundary: every numeric field built from a numpy
        # scalar. Catches a YAML-safety leak in any field without relying on a
        # per-field test being remembered for it.
        species = DefectSpecies(
            "Z",
            nsites=np.int64(4),
            charge_states=[
                DefectChargeState(
                    charge=np.int64(0),
                    energy=np.float64(0.5),
                    degeneracy=np.float64(1),
                ),
                DefectChargeState(
                    charge=np.int64(1), fixed_concentration=np.float64(2.0)
                ),
            ],
        )
        system = DefectSystem(
            defect_species=[species],
            dos=self.dos,
            volume=np.float64(100.0),
            temperature=np.float64(300.0),
            convergence_tolerance=np.float64(1e-8),
            site_pools={"p": (np.float64(4.0), ["Z"])},
            element_pools={"X": (np.float64(0.5), [("Z", np.float64(1.0))])},
        )
        yaml.safe_dump(system.as_dict())


class TestDefectSystemFixedConcentrations(unittest.TestCase):
    def setUp(self):
        self.dos = DOS(
            dos=np.ones(101),
            edos=np.linspace(-5.0, 5.0, 101),
            bandgap=2.0,
            nelect=10,
        )

    def _donor(self, name="X", energy=0.5, nsites=1):
        return DefectSpecies(
            name=name,
            nsites=nsites,
            charge_states=[DefectChargeState(charge=1, energy=energy, degeneracy=1)],
        )

    def test_constructor_fixes_species_total_by_name(self):
        species = self._donor("X")
        system = DefectSystem(
            defect_species=[species],
            dos=self.dos,
            volume=100,
            temperature=300,
            fixed_concentrations={"X": 0.01},
        )
        self.assertAlmostEqual(
            system.concentration_dict(per_volume=False)["X"], 0.01
        )
        self.assertEqual(system.defect_species_by_name("X").fixed_concentration, 0.01)

    def test_factory_at_fixes_species_total_by_name(self):
        factory = DefectSystemFactory(
            defect_species=[self._donor("X")], dos=self.dos, volume=100
        )
        system = factory.at(300, fixed_concentrations={"X": 0.01})
        self.assertAlmostEqual(
            system.concentration_dict(per_volume=False)["X"], 0.01
        )
        self.assertEqual(system.defect_species_by_name("X").fixed_concentration, 0.01)

    def test_anneal_and_quench_freezes_some_species_and_re_equilibrates_rest(self):
        # A minority donor frozen at its high-temperature total, plus a major
        # donor/acceptor pair that re-equilibrates when the temperature drops.
        frozen = self._donor("X_frozen", energy=1.5)
        donor = self._donor("D", energy=0.5)
        acceptor = DefectSpecies(
            name="A",
            nsites=1,
            charge_states=[DefectChargeState(charge=-1, energy=0.5, degeneracy=1)],
        )
        factory = DefectSystemFactory(
            defect_species=[frozen, donor, acceptor], dos=self.dos, volume=100
        )
        T_high, T_low = 1000, 300

        high = factory.at(T_high).concentration_dict(per_volume=False)
        low = factory.at(
            T_low, fixed_concentrations={"X_frozen": high["X_frozen"]}
        ).concentration_dict(per_volume=False)

        # the frozen species keeps its high-temperature total
        self.assertAlmostEqual(
            low["X_frozen"], high["X_frozen"], delta=abs(high["X_frozen"]) * 1e-9
        )
        # a non-frozen species and the carriers re-equilibrate
        self.assertNotAlmostEqual(low["D"], high["D"])
        self.assertNotAlmostEqual(low["n0"], high["n0"])

    def test_fix_on_one_call_does_not_leak_to_another_or_to_the_factory(self):
        species = self._donor("X")
        factory = DefectSystemFactory(
            defect_species=[species], dos=self.dos, volume=100
        )

        fixed = factory.at(300, fixed_concentrations={"X": 0.01})
        free = factory.at(300)

        # the un-fixed call is unaffected by the earlier fixed call
        self.assertIsNone(free.defect_species_by_name("X").fixed_concentration)
        self.assertNotAlmostEqual(
            free.concentration_dict(per_volume=False)["X"], 0.01
        )
        self.assertAlmostEqual(
            fixed.concentration_dict(per_volume=False)["X"], 0.01
        )
        # the factory's own species object is never mutated
        self.assertIsNone(species.fixed_concentration)

    def test_fix_composes_with_formation_energy_corrections(self):
        cs_a = DefectChargeState(charge=1, energy=0.5, degeneracy=1)
        cs_b = DefectChargeState(charge=1, energy=0.9, degeneracy=1)
        corrected = DefectSpecies("X_i", nsites=1, charge_states=[cs_a, cs_b])
        acceptor = DefectSpecies(
            name="A",
            nsites=1,
            charge_states=[DefectChargeState(charge=-1, energy=0.5, degeneracy=1)],
        )
        factory = DefectSystemFactory(
            defect_species=[corrected, acceptor],
            dos=self.dos,
            volume=100,
            formation_energy_correction_fns={
                cs_a: lambda T: 0.1,
                cs_b: lambda T: -0.05,
            },
        )
        system = factory.at(300, fixed_concentrations={"X_i": 0.02})

        # the corrections (resolved by identity) survive applying the fix
        self.assertAlmostEqual(system.defect_species[0].charge_states[0].energy, 0.6)
        self.assertAlmostEqual(system.defect_species[0].charge_states[1].energy, 0.85)
        # and the fix (resolved by name) is in effect
        self.assertAlmostEqual(
            system.concentration_dict(per_volume=False)["X_i"], 0.02
        )

    def test_unknown_species_name_raises_value_error_naming_it(self):
        with self.assertRaises(ValueError) as ctx:
            DefectSystem(
                defect_species=[self._donor("X")],
                dos=self.dos,
                volume=100,
                temperature=300,
                fixed_concentrations={"NOPE": 0.01},
            )
        self.assertIn("NOPE", str(ctx.exception))

    def test_over_budget_fix_raises_value_error_at_construction(self):
        with self.assertRaises(ValueError):
            DefectSystem(
                defect_species=[self._donor("X", nsites=1)],
                dos=self.dos,
                volume=100,
                temperature=300,
                fixed_concentrations={"X": 5.0},
            )

    def test_nan_fixed_concentration_raises_at_construction(self):
        # a charge-neutral species never enters charge neutrality, so without
        # the guard a NaN survives the solve and surfaces silently as nan.
        neutral = DefectSpecies(
            "N", 1, [DefectChargeState(charge=0, energy=0.5, degeneracy=1)]
        )
        with self.assertRaises(ValueError) as ctx:
            DefectSystem(
                defect_species=[neutral],
                dos=self.dos,
                volume=100,
                temperature=300,
                fixed_concentrations={"N": float("nan")},
            )
        self.assertIn("N", str(ctx.exception))
        self.assertIn("finite", str(ctx.exception))

    def test_negative_fixed_concentration_raises_with_clear_message(self):
        with self.assertRaises(ValueError) as ctx:
            DefectSystem(
                defect_species=[self._donor("X")],
                dos=self.dos,
                volume=100,
                temperature=300,
                fixed_concentrations={"X": -0.01},
            )
        message = str(ctx.exception)
        self.assertIn("X", message)
        self.assertIn("non-negative", message)

    def test_overrides_a_constructed_species_level_fixed_concentration(self):
        species = DefectSpecies(
            "X",
            nsites=1,
            charge_states=[DefectChargeState(charge=1, energy=0.5, degeneracy=1)],
            fixed_concentration=0.05,
        )
        system = DefectSystem(
            defect_species=[species],
            dos=self.dos,
            volume=100,
            temperature=300,
            fixed_concentrations={"X": 0.01},
        )
        self.assertEqual(system.defect_species_by_name("X").fixed_concentration, 0.01)
        self.assertAlmostEqual(
            system.concentration_dict(per_volume=False)["X"], 0.01
        )

    def test_fixes_multiple_species_in_one_call(self):
        system = DefectSystem(
            defect_species=[self._donor("X", energy=0.5), self._donor("Y", energy=0.7)],
            dos=self.dos,
            volume=100,
            temperature=300,
            fixed_concentrations={"X": 0.01, "Y": 0.02},
        )
        cd = system.concentration_dict(per_volume=False)
        self.assertAlmostEqual(cd["X"], 0.01)
        self.assertAlmostEqual(cd["Y"], 0.02)

    def test_pooled_fix_exceeding_pool_budget_raises_naming_the_pool(self):
        a = DefectSpecies(
            "A", 5, [DefectChargeState(charge=0, energy=1.0, degeneracy=1)]
        )
        b = DefectSpecies(
            "B", 5, [DefectChargeState(charge=0, energy=1.0, degeneracy=1)]
        )
        with self.assertRaisesRegex(ValueError, "shared"):
            DefectSystem(
                defect_species=[a, b],
                dos=self.dos,
                volume=100,
                temperature=300,
                site_pools={"shared": (1.0, ["A", "B"])},
                fixed_concentrations={"A": 0.7, "B": 0.5},
            )


class TestDiluteLimitWarning(unittest.TestCase):
    def setUp(self):
        self.dos = DOS(
            dos=np.ones(101),
            edos=np.linspace(-5.0, 5.0, 101),
            bandgap=2.0,
            nelect=10,
        )

    def _make_system(self, defect_species, **kwargs):
        return DefectSystem(
            defect_species=defect_species,
            dos=self.dos,
            volume=100,
            temperature=300,
            **kwargs,
        )

    @staticmethod
    def _saturating_species(name="S", nsites=2, energy=-2.0):
        # A neutral state's occupancy w/(1+w) is e_fermi-independent; a low
        # energy drives it to ~100% of its sites.
        return DefectSpecies(
            name,
            nsites=nsites,
            charge_states=[DefectChargeState(charge=0, energy=energy, degeneracy=1)],
        )

    def test_threshold_defaults_to_one_percent(self):
        system = self._make_system([self._saturating_species()])
        self.assertEqual(system.occupancy_warning_threshold, 0.01)

    def test_threshold_is_stored(self):
        system = self._make_system(
            [self._saturating_species()], occupancy_warning_threshold=0.2
        )
        self.assertEqual(system.occupancy_warning_threshold, 0.2)

    def test_none_threshold_is_stored(self):
        system = self._make_system(
            [self._saturating_species()], occupancy_warning_threshold=None
        )
        self.assertIsNone(system.occupancy_warning_threshold)

    def test_out_of_range_threshold_raises(self):
        species = self._saturating_species()
        for bad in (0.0, 1.5, -0.1, float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                self._make_system([species], occupancy_warning_threshold=bad)

    def test_threshold_of_one_is_allowed(self):
        system = self._make_system(
            [self._saturating_species()], occupancy_warning_threshold=1.0
        )
        self.assertEqual(system.occupancy_warning_threshold, 1.0)

    def test_site_occupancy_fractions_match_site_percentages(self):
        # Pooled system: occupancy is measured against the pool size, and the
        # fractions are exactly site_percentages / 100.
        a = DefectSpecies(
            "A",
            nsites=5,
            charge_states=[DefectChargeState(charge=0, energy=0.057, degeneracy=1)],
        )
        b = DefectSpecies(
            "B",
            nsites=1,
            charge_states=[DefectChargeState(charge=0, energy=0.2, degeneracy=1)],
        )
        system = self._make_system([a, b], site_pools={"shared": (4.0, [a, b])})
        e_fermi = system.get_sc_fermi()[0]
        fractions = system._site_occupancy_fractions(e_fermi)
        percentages = system.site_percentages()
        for name in ("A", "B"):
            self.assertAlmostEqual(fractions[name] * 100, percentages[name], places=10)

    @staticmethod
    def _dilute_warnings(records):
        return [r for r in records if issubclass(r.category, DiluteLimitWarning)]

    def test_high_occupancy_species_emits_single_warning(self):
        system = self._make_system([self._saturating_species("S")])
        with warnings.catch_warnings(record=True) as records:
            warnings.simplefilter("always")
            system.get_sc_fermi()
        dilute = self._dilute_warnings(records)
        self.assertEqual(len(dilute), 1)
        self.assertIn("S", str(dilute[0].message))

    def test_dilute_species_emits_no_warning(self):
        species = DefectSpecies(
            "D",
            nsites=1,
            charge_states=[DefectChargeState(charge=0, energy=1.0, degeneracy=1)],
        )
        system = self._make_system([species])
        with warnings.catch_warnings(record=True) as records:
            warnings.simplefilter("always")
            system.get_sc_fermi()
        self.assertEqual(self._dilute_warnings(records), [])

    def test_none_threshold_suppresses_warning(self):
        system = self._make_system(
            [self._saturating_species("S")], occupancy_warning_threshold=None
        )
        with warnings.catch_warnings(record=True) as records:
            warnings.simplefilter("always")
            system.get_sc_fermi()
        self.assertEqual(self._dilute_warnings(records), [])

    def test_custom_threshold_changes_when_warning_fires(self):
        # Neutral state at ~0.057 eV -> ~10% site occupancy.
        def at_ten_percent():
            return DefectSpecies(
                "T",
                nsites=1,
                charge_states=[
                    DefectChargeState(charge=0, energy=0.057, degeneracy=1)
                ],
            )

        fires = self._make_system(
            [at_ten_percent()], occupancy_warning_threshold=0.01
        )
        silent = self._make_system(
            [at_ten_percent()], occupancy_warning_threshold=0.5
        )
        with warnings.catch_warnings(record=True) as fires_records:
            warnings.simplefilter("always")
            fires.get_sc_fermi()
        with warnings.catch_warnings(record=True) as silent_records:
            warnings.simplefilter("always")
            silent.get_sc_fermi()
        self.assertEqual(len(self._dilute_warnings(fires_records)), 1)
        self.assertEqual(self._dilute_warnings(silent_records), [])

    def test_multiple_high_occupancy_species_listed_in_one_warning(self):
        a = self._saturating_species("A", energy=-2.0)
        b = self._saturating_species("B", energy=-1.5)
        system = self._make_system([a, b])
        with warnings.catch_warnings(record=True) as records:
            warnings.simplefilter("always")
            system.get_sc_fermi()
        dilute = self._dilute_warnings(records)
        self.assertEqual(len(dilute), 1)
        message = str(dilute[0].message)
        self.assertIn("A", message)
        self.assertIn("B", message)
        self.assertIn(" and ", message)

    def test_warning_is_dilute_limit_category_with_expected_content(self):
        system = self._make_system([self._saturating_species("S")])
        with warnings.catch_warnings(record=True) as records:
            warnings.simplefilter("always")
            system.get_sc_fermi()
        dilute = self._dilute_warnings(records)
        self.assertEqual(len(dilute), 1)
        self.assertIs(dilute[0].category, DiluteLimitWarning)
        message = str(dilute[0].message)
        self.assertIn("S", message)
        self.assertIn("%", message)
        self.assertIn("dilute", message.lower())
        self.assertIn("non-physical", message.lower())

    def test_pooled_species_occupancy_measured_against_pool_size(self):
        # nsites=1 but the species shares a 10-site pool. At ~0.057 eV it holds
        # ~1 defect/cell -> ~10% of the POOL. If occupancy were (wrongly)
        # measured against nsites=1 it would read ~100%, so the reported figure
        # discriminates the denominator. site_percentages is the pool-based
        # reference the warning must agree with.
        species = DefectSpecies(
            "P",
            nsites=1,
            charge_states=[DefectChargeState(charge=0, energy=0.057, degeneracy=1)],
        )
        system = self._make_system(
            [species],
            site_pools={"shared": (10.0, [species])},
            occupancy_warning_threshold=0.01,
        )
        with warnings.catch_warnings(record=True) as records:
            warnings.simplefilter("always")
            system.get_sc_fermi()
        dilute = self._dilute_warnings(records)
        self.assertEqual(len(dilute), 1)
        pool_pct = system.site_percentages()["P"]
        self.assertLess(pool_pct, 100.0)
        self.assertIn(f"{pool_pct:.0f}%", str(dilute[0].message))


if __name__ == "__main__":
    unittest.main()
