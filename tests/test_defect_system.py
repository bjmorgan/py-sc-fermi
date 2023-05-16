import unittest
from unittest.mock import Mock, patch
from io import StringIO

import numpy as np
import os
import textwrap
from py_sc_fermi.defect_species import DefectSpecies
from py_sc_fermi.dos import DOS
from py_sc_fermi.defect_system import DefectSystem, CustomWarningManager
from py_sc_fermi.defect_charge_state import DefectChargeState


input_string = "1\n12\n0.1\n298\n1\nv_O 1 1\n 1 1 1\n1\nO_i 1e+22\n1\nO_i 1 1e+22\n"
input_string_spin = (
    "0\n12\n0.1\n298\n1\nv_O 1 1\n 1 1 1\n1\nO_i 1e+22\n1\nO_i 1 1e+22\n"
)
test_data_dir = "dummy_inputs/"
test_report_filename = os.path.join(
    os.path.dirname(__file__), test_data_dir, "report_string.txt"
)
test_yaml_filename = os.path.join(
    os.path.dirname(__file__), test_data_dir, "defect_system.yaml"
)
test_exception_yaml_filename = os.path.join(
    os.path.dirname(__file__), test_data_dir, "bad_yaml.yaml"
)
test_vasprun_filename = os.path.join(
    os.path.dirname(__file__), test_data_dir, "vasprun_nsp.xml"
)


class TestCustomWarningManager(unittest.TestCase):
    def setUp(self):
        self.warning_manager = CustomWarningManager()

    @patch('sys.stdout', new_callable=StringIO)
    def test_dos_overflow_warning(self, mock_stdout):
        self.warning_manager.custom_warning('overflow', RuntimeWarning, 'dos_file.py', 42)
        expected_warning = textwrap.dedent(
                        """DOSOverflowWarning: An overflow occurred during computation of
                        electron and hole concentrations. This is likely a natural result of the use of
                        a numerical solver for the Fermi energy search. This can likely be ignored
                        though you should always check the final results are reasonable.""")
        self.assertEqual(mock_stdout.getvalue().strip(), expected_warning.strip())

    @patch('sys.stdout', new_callable=StringIO)
    def test_defect_overflow_warning(self, mock_stdout):
        self.warning_manager.custom_warning('overflow', RuntimeWarning, 'defect_file.py', 42)
        expected_warning = textwrap.dedent(
                        """DefectOverflowWarning: An overflow occurred during computation of
                        defect concentrations. This is likely a natural result of the use of
                        a numerical solver for the Fermi energy search. This can likely be ignored
                        though you should always check the final results are reasonable.""")
        self.assertEqual(mock_stdout.getvalue().strip(), expected_warning.strip())

    @patch('sys.stdout', new_callable=StringIO)
    def test_other_warning(self, mock_stdout):
        self.warning_manager.custom_warning('other warning', RuntimeWarning, 'other_file.py', 42, None, None)
        expected_warning = "RuntimeWarning: other warning"
        self.assertEqual(mock_stdout.getvalue().strip(), expected_warning)


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
            n_trial_steps=100,
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
        self.defect_system.dos = DOS.from_vasprun(test_vasprun_filename, nelect=12)
        defect_dict = self.defect_system.as_dict()
        self.assertEqual(defect_dict["volume"], 100)
        self.assertEqual(defect_dict["temperature"], 298)
        self.assertEqual(defect_dict["n_trial_steps"], 1500)

    def test_from_dict(self):
        dictionary = {
            "volume": 100,
            "temperature": 100,
            "n_trial_steps": 100,
            "convergence_tolerance": 1,
            "defect_species": [{
                "name": "V_O",
                "nsites": 2,
                "charge_states": {1 : {"charge": 1, "energy": 0, "degeneracy": 1}},
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

    def test__get_report_string(self):
        self.defect_system.get_sc_fermi = Mock(return_value=[0.5, {}])
        self.defect_system.dos.carrier_concentrations = Mock(return_value=(100, 100))
        self.defect_system.defect_species[0].name = "v_O"
        self.defect_system.defect_species[0].charge_state_concentrations = Mock(
            return_value={+1: 1000}
        )
        self.defect_system.defect_species[0].charge_states = {
            1: DefectChargeState(charge=1, fixed_concentration=1000)
        }
        self.defect_system.defect_species[0].get_concentration = Mock(return_value=1000)
        self.defect_system.defect_species[1].get_concentration = Mock(return_value=1000)
        self.defect_system.defect_species[1].name = "O_i"
        self.defect_system.defect_species[1].charge_states = {
            -1: DefectChargeState(charge=-1, fixed_concentration=1000)
        }
        self.defect_system.defect_species[1].charge_state_concentrations = Mock(
            return_value={-1: 1000}
        )

        with open(test_report_filename, "r") as tst_string:
            test_string = tst_string.read()

        self.assertEqual(
            self.defect_system._get_report_string().strip(), test_string.strip()
        )

    def test_get_sc_fermi(self):
        self.defect_system.dos.emin = Mock(return_value=0)
        self.defect_system.dos.emax = Mock(return_value=1)
        self.defect_system.dos.carrier_concentrations = Mock(return_value=(1, 1))
        self.defect_system.q_tot = Mock(return_value=0)
        self.assertEqual(
            self.defect_system.get_sc_fermi(),
            (0.5, 0),
        )

    def test_get_sc_fermi_bottoms_out(self):
        self.defect_system.dos.emin = Mock(return_value=0)
        self.defect_system.dos.emax = Mock(return_value=1)
        self.defect_system.q_tot = Mock(return_value=(0.1))
        with self.assertRaises(RuntimeError):
            self.defect_system.get_sc_fermi()

    def test_get_sc_fermi_tops_out(self):
        self.defect_system.dos.emin = Mock(return_value=1)
        self.defect_system.dos.emax = Mock(return_value=0)
        self.defect_system.q_tot = Mock(return_value=(-0.1))
        with self.assertRaises(RuntimeError):
            self.defect_system.get_sc_fermi()

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
        self.defect_system.get_sc_fermi = Mock(return_value=[1, {}])
        self.defect_system.dos.carrier_concentrations = Mock(return_value=(1, 1))
        self.defect_system.defect_species[0].get_concentration = Mock(return_value=1)
        self.defect_system.defect_species[1].get_concentration = Mock(return_value=1)
        self.defect_system.defect_species[0].charge_state_concentrations = Mock(return_value={1: 1})
        self.defect_system.defect_species[1].charge_state_concentrations = Mock(return_value={-1: 1})
        self.defect_system.defect_species[0].charge_states = {1: 1}
        self.defect_system.defect_species[1].charge_states = {-1: 1}
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
            f"DefectSystem\n  nelect: 100 e\n  bandgap: 0.1 eV\n  volume: 100 A^3\n  temperature: 298 K\n\nContains defect species:\n".strip(),
        )


if __name__ == "__main__":
    unittest.main()
