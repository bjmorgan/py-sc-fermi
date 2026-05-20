import unittest
from unittest.mock import Mock, patch
from io import StringIO
import warnings

import numpy as np
import os
import textwrap

from py_sc_fermi.defect_species import DefectSpecies
from py_sc_fermi.dos import DOS
from py_sc_fermi.defect_system import DefectSystem
from py_sc_fermi.defect_charge_state import DefectChargeState


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
            convergence_tolerance=1e-6,
        )
        self.assertEqual(defect_system.volume, volume)
        self.assertEqual(defect_system.dos, dos)
        self.assertEqual(defect_system.temperature, temperature)
        self.assertEqual(defect_system.defect_species[0], mock_defect_species[0])
        self.assertEqual(defect_system.defect_species[1], mock_defect_species[1])
        self.assertEqual(defect_system.n_trial_steps, None)


class TestDefectSystem(unittest.TestCase):
    def setUp(self):
        volume = 100
        mock_defect_species = [Mock(spec=DefectSpecies), Mock(spec=DefectSpecies)]
        mock_defect_species[0].name = "v_O"
        mock_defect_species[1].name = "O_i"
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
            "n_trial_steps": 100,
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
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            defect_system = self.defect_system.from_dict(dictionary)
            self.defect_system.from_dict(dictionary)
        self.assertEqual(defect_system.volume, 100)
        self.assertEqual(defect_system.temperature, 100)
        self.assertEqual(defect_system.n_trial_steps, 100)
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
            charge_states={1: charge_state},
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
        self.defect_system.defect_species[0].charge_state_concentrations = Mock(return_value=[(cs_v_O, 1)])
        self.defect_system.defect_species[1].charge_state_concentrations = Mock(return_value=[(cs_O_i, 1)])
        self.defect_system.defect_species[0].charge_states = [cs_v_O]
        self.defect_system.defect_species[1].charge_states = [cs_O_i]
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
        
    def test_n_trial_steps_deprecation_warning(self):
        """Setting n_trial_steps should emit a DeprecationWarning."""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            DefectSystem(
                defect_species=[],
                dos=Mock(spec=DOS),
                volume=100,
                temperature=300,
                n_trial_steps=100,
            )
        
        deprecation_warnings = [w for w in caught if issubclass(w.category, DeprecationWarning)]
        self.assertEqual(len(deprecation_warnings), 1)
        self.assertIn("n_trial_steps", str(deprecation_warnings[0].message))


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
            nsites=1,
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

    def test_no_pools_matches_dilute_limit(self):
        system = DefectSystem(
            defect_species=[self.species_a, self.species_b],
            dos=self.dos,
            volume=100,
            temperature=300,
        )
        e_fermi = 1.0
        pooled_total = sum(system._global_defect_concs(e_fermi).values())
        dilute_total = sum(
            conc
            for sp in (self.species_a, self.species_b)
            for _, conc in sp.charge_state_concentrations(e_fermi, 300)
        )
        self.assertAlmostEqual(pooled_total, dilute_total, places=8)

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
            for sp in (self.species_a, self.species_b)
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
            for sp in (self.species_a, self.species_b)
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
            nsites=1,
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
        total_mg = sum(concs[cs] for cs in self.species.charge_states)
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
        total_mg = sum(concs[cs] for cs in self.species.charge_states)
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
        self.assertEqual(concs[self.species.charge_states[0]], fixed_value)
        total_mg = sum(concs[cs] for cs in self.species.charge_states)
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

    def test_no_shift_functions_is_a_no_op(self):
        system = self._make_system()
        with system._with_band_edge_corrections():
            self.assertEqual(self.cs0.energy, 1.0)
            self.assertEqual(self.cs1.energy, 1.5)
            self.assertEqual(system.dos._bandgap, 2.0)
        self.assertFalse(system._corrections_active)

    def test_rigid_shift_changes_bandgap_but_not_formation_energies(self):
        system = self._make_system(
            vbm_shift_fn=lambda T: 0.05,
            cbm_shift_fn=lambda T: -0.02,
            rigid_shift=True,
        )
        with system._with_band_edge_corrections():
            self.assertEqual(self.cs0.energy, 1.0)
            self.assertEqual(self.cs1.energy, 1.5)
            self.assertAlmostEqual(system.dos._bandgap, 2.0 + (-0.02 - 0.05))
        self.assertEqual(self.cs0.energy, 1.0)
        self.assertEqual(self.cs1.energy, 1.5)
        self.assertEqual(system.dos._bandgap, 2.0)
        self.assertFalse(system._corrections_active)

    def test_non_rigid_shift_applies_charge_dependent_correction(self):
        d_vbm = 0.05
        system = self._make_system(vbm_shift_fn=lambda T: d_vbm, rigid_shift=False)
        with system._with_band_edge_corrections():
            self.assertAlmostEqual(self.cs0.energy, 1.0 - self.cs0.charge * d_vbm)
            self.assertAlmostEqual(self.cs1.energy, 1.5 - self.cs1.charge * d_vbm)
        self.assertEqual(self.cs0.energy, 1.0)
        self.assertEqual(self.cs1.energy, 1.5)

    def test_formation_energy_correction_takes_priority_over_rigid_shift_flag(self):
        system = self._make_system(
            vbm_shift_fn=lambda T: 0.05,
            rigid_shift=False,
            formation_energy_corrections={("V_O", 1): lambda T: 0.2},
        )
        with system._with_band_edge_corrections():
            self.assertAlmostEqual(self.cs1.energy, 1.5 + 0.2)
            self.assertAlmostEqual(self.cs0.energy, 1.0)
        self.assertEqual(self.cs0.energy, 1.0)
        self.assertEqual(self.cs1.energy, 1.5)

    def test_corrections_are_not_re_applied_on_nested_entry(self):
        system = self._make_system(
            vbm_shift_fn=lambda T: 0.05, cbm_shift_fn=lambda T: -0.02
        )
        with system._with_band_edge_corrections():
            bandgap_once_applied = system.dos._bandgap
            with system._with_band_edge_corrections():
                self.assertTrue(system._corrections_active)
                self.assertEqual(system.dos._bandgap, bandgap_once_applied)
            self.assertEqual(system.dos._bandgap, bandgap_once_applied)
        self.assertEqual(system.dos._bandgap, 2.0)
        self.assertFalse(system._corrections_active)

    def test_exception_in_shift_function_resets_guard_and_state(self):
        def boom(temperature):
            raise RuntimeError("shift function failed")

        system = self._make_system(vbm_shift_fn=boom)
        with self.assertRaises(RuntimeError):
            with system._with_band_edge_corrections():
                pass
        self.assertFalse(system._corrections_active)
        self.assertEqual(system.dos._bandgap, 2.0)
        self.assertEqual(self.cs0.energy, 1.0)
        self.assertEqual(self.cs1.energy, 1.5)

    def test_exception_in_formation_energy_correction_rolls_back_partial_changes(self):
        def boom(temperature):
            raise RuntimeError("correction failed")

        system = self._make_system(
            vbm_shift_fn=lambda T: 0.05,
            formation_energy_corrections={("V_O", 1): boom},
        )
        with self.assertRaises(RuntimeError):
            with system._with_band_edge_corrections():
                pass
        self.assertFalse(system._corrections_active)
        self.assertEqual(system.dos._bandgap, 2.0)
        self.assertEqual(self.cs0.energy, 1.0)
        self.assertEqual(self.cs1.energy, 1.5)

    def test_exception_in_context_body_still_restores_state(self):
        system = self._make_system(
            vbm_shift_fn=lambda T: 0.05, cbm_shift_fn=lambda T: -0.02
        )
        with self.assertRaises(ValueError):
            with system._with_band_edge_corrections():
                raise ValueError("body failed")
        self.assertFalse(system._corrections_active)
        self.assertEqual(system.dos._bandgap, 2.0)
        self.assertEqual(self.cs0.energy, 1.0)
        self.assertEqual(self.cs1.energy, 1.5)


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

    def test_report_leaves_band_edge_corrections_restored(self):
        system = DefectSystem(
            defect_species=[self.species],
            dos=self.dos,
            volume=100,
            temperature=300,
            vbm_shift_fn=lambda T: 0.05,
            cbm_shift_fn=lambda T: -0.02,
        )
        with patch("sys.stdout", new=StringIO()):
            system.report()
        self.assertEqual(system.dos._bandgap, 2.0)
        self.assertFalse(system._corrections_active)


if __name__ == "__main__":
    unittest.main()
