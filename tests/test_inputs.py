import unittest
from unittest.mock import Mock
import numpy as np
from py_sc_fermi.defect_species import DefectSpecies
from py_sc_fermi.dos import DOS
import os

from py_sc_fermi.inputs import (
    volume_from_structure,
    read_dos_data,
    volume_from_unitcell,
    read_volume_from_structure_file,
    read_input_fermi,
    is_yaml,
    InputSet,
)
from pymatgen.core.structure import Structure

test_data_dir = "dummy_inputs/"
test_poscar_filename = os.path.join(os.path.dirname(__file__), test_data_dir, "POSCAR")
test_unitcell_filename = os.path.join(
    os.path.dirname(__file__), test_data_dir, "unitcell.dat"
)
test_sc_fermi_input_filename = os.path.join(
    os.path.dirname(__file__), test_data_dir, "input_fermi.dat"
)
test_frozen_sc_fermi_input_filename = os.path.join(
    os.path.dirname(__file__), test_data_dir, "frozen_charge_states.dat"
)
test_bad_frozen_sc_fermi_input_filename = os.path.join(
    os.path.dirname(__file__), test_data_dir, "bad_frozen_species.dat"
)
test_dos_filename = os.path.join(os.path.dirname(__file__), test_data_dir, "totdos.dat")
test_defect_system_yaml_filename = os.path.join(
    os.path.dirname(__file__), test_data_dir, "defect_system.yaml"
)
test_exception_yaml_filename = os.path.join(
    os.path.dirname(__file__), test_data_dir, "bad_yaml.yaml"
)

structure = Structure.from_file(test_poscar_filename)
volume = structure.volume


class TestInputsSetInit(unittest.TestCase):
    def test_input_set_is_initialised(self):
        dos = Mock(spec=DOS)
        volume = 100
        defect_species = [Mock(spec=DefectSpecies)]
        temperature = 298
        conv = 1
        n_trial = 100
        input_set = InputSet(dos, volume, defect_species, temperature, conv, n_trial)
        self.assertEqual(input_set.dos, dos)
        self.assertEqual(input_set.volume, volume)
        self.assertEqual(input_set.defect_species, defect_species)
        self.assertEqual(input_set.temperature, temperature)
        self.assertEqual(input_set.convergence_tolerance, conv)
        self.assertEqual(input_set.n_trial_steps, n_trial)


class TestInputSet(unittest.TestCase):
    def test_from_yaml(self):
        input_set = InputSet.from_yaml(test_defect_system_yaml_filename)
        self.assertEqual(input_set.volume, 59)
        self.assertEqual(input_set.dos.nelect, 18)
        self.assertEqual(input_set.dos.bandgap, 0.8084)
        self.assertEqual(input_set.temperature, 300)
        self.assertEqual(len(input_set.defect_species), 3)

    def test_from_yaml_raises(self):
        with self.assertRaises(ValueError):
            InputSet.from_yaml(test_exception_yaml_filename)

    def test_from_sc_fermi_inputs(self):
        input_set = InputSet.from_sc_fermi_inputs(
            test_sc_fermi_input_filename, test_unitcell_filename, test_dos_filename
        )
        self.assertEqual(input_set.volume, 544.7091796190017)
        self.assertEqual(input_set.dos.nelect, 18)
        self.assertEqual(input_set.dos.bandgap, 0.8084)
        self.assertEqual(input_set.temperature, 100)
        self.assertEqual(len(input_set.defect_species), 2)


class TestInputs(unittest.TestCase):
    def test_volume_from_structure(self):
        self.assertAlmostEqual(volume_from_structure(test_poscar_filename), volume)

    def test_volume_from_unitcell(self):
        self.assertAlmostEqual(volume_from_unitcell(test_unitcell_filename), volume)

    def test_is_yaml(self):
        self.assertTrue(is_yaml(test_defect_system_yaml_filename))
        self.assertFalse(is_yaml(test_sc_fermi_input_filename))

    def test_read_input_fermi(self):
        defect_data = read_input_fermi(test_sc_fermi_input_filename, volume=1)
        self.assertEqual(len(defect_data.defect_species), 2)
        self.assertEqual(defect_data.temperature, 100)
        self.assertEqual(defect_data.bandgap, 0.8084)
        self.assertEqual(defect_data.nelect, 18)

    def test_read_input_fermi_raises(self):
        with self.assertRaises(ValueError):
            read_input_fermi(test_sc_fermi_input_filename, volume=None, frozen=True)

    def test_read_input_fermi_frozen(self):
        read_input_fermi(test_frozen_sc_fermi_input_filename, volume=1, frozen=True)

    def test_bad_frozen_name_raises(self):
        with self.assertRaises(ValueError):
            read_input_fermi(
                test_bad_frozen_sc_fermi_input_filename, volume=1, frozen=True
            )

    def test_read_dos_data(self):
        dos_data = read_dos_data(filename=test_dos_filename, bandgap=1, nelect=1)
        self.assertEqual(type(dos_data), DOS)

    def test_read_volume_from_structure_file(self):
        self.assertAlmostEqual(
            read_volume_from_structure_file(test_poscar_filename), volume
        )
        self.assertAlmostEqual(
            read_volume_from_structure_file(test_unitcell_filename), volume
        )


if __name__ == "__main__":
    unittest.main()
