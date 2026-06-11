import os
import unittest
from io import StringIO
from unittest.mock import Mock, patch

import numpy as np
from scipy.constants import physical_constants

from py_sc_fermi.defect_charge_state import DefectChargeState
from py_sc_fermi.defect_species import DefectSpecies
from py_sc_fermi.defect_system import DefectSystem, DefectSystemFactory
from py_sc_fermi.dos import DOS

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

    def test_site_percentages(self):
        self.defect_system.get_sc_fermi = Mock(return_value=[1, {}])
        self.defect_system.dos.carrier_concentrations = Mock(return_value=(1, 1))
        self.defect_system.defect_species[0].get_concentration = Mock(return_value=1)
        self.defect_system.defect_species[1].get_concentration = Mock(return_value=1)
        self.defect_system.defect_species[0].nsites = 1
        self.defect_system.defect_species[1].nsites = 1
        self.defect_system.defect_species[0].name = "v_O"
        self.defect_system.defect_species[1].name = "O_i"
        self.assertEqual(
            self.defect_system.site_percentages(), {"v_O": 100, "O_i": 100}
        )

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

    def test_pools_can_reference_species_by_name(self):
        n_pool = 10.0
        system = DefectSystem(
            defect_species=[self.species_a, self.species_b],
            dos=self.dos,
            volume=100,
            temperature=300,
            site_pools={"shared": (n_pool, ["A", "B"])},
        )
        concs = system._global_defect_concs(1.0)
        total_occupied = sum(
            concs[cs]
            for sp in system.defect_species
            for cs in sp.charge_states
        )
        self.assertLessEqual(total_occupied, n_pool)

    def test_pool_raises_when_fixed_concentrations_exceed_site_count(self):
        self.species_a.charge_states[0].fix_concentration(20.0)
        system = DefectSystem(
            defect_species=[self.species_a],
            dos=self.dos,
            volume=100,
            temperature=300,
            site_pools={"shared": (10.0, [self.species_a])},
        )
        with self.assertRaises(ValueError):
            system._global_defect_concs(1.0)

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
        system = DefectSystem(
            defect_species=[self.species_a],
            dos=self.dos,
            volume=100,
            temperature=300,
            site_pools={"shared": (10.0, [self.species_a])},
        )
        with self.assertRaises(ValueError):
            system._global_defect_concs(1.0)


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
        self.assertAlmostEqual(c_a, fixed_value, places=8)
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
        self.assertAlmostEqual(ratios[0] / ratios[1], 1.0, places=6)

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


class TestDefectSystemElementPoolConvergence(unittest.TestCase):
    """Convergence of the coupled element-pool chemical-potential solve
    across physically realistic concentration scales (~1e-2 to ~1e-18
    defects per unit cell)."""

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

    def test_coupled_dilute_pools_hit_targets_exactly(self):
        system, target_x, target_y = self.coupled_system(energy=1.0)
        concs = system._global_defect_concs(1.0)
        content = self.species_content(system, concs)
        self.assertAlmostEqual(
            (content["P"] + content["Q"]) / target_x, 1.0, delta=1e-8
        )
        self.assertAlmostEqual(content["Q"] / target_y, 1.0, delta=1e-8)

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


if __name__ == "__main__":
    unittest.main()
