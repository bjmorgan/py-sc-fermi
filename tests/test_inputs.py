import unittest
import numpy as np
from py_sc_fermi.dos import DOS
import os

from py_sc_fermi.inputs import (
    volume_from_structure,
    read_unitcell_data,
    read_input_data,
    read_dos_data,
    inputs_from_files,
    dos_from_dict,
    defect_species_from_dict,
    defect_system_from_yaml
)
from pymatgen.core.structure import Structure

test_data_dir = "inputs/"
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
test_dos_filename = os.path.join(os.path.dirname(__file__), test_data_dir, "totdos.dat")
test_defect_system_yaml_filename = os.path.join(os.path.dirname(__file__), test_data_dir, "defect_system.yaml")

structure = Structure.from_file(test_poscar_filename)
volume = structure.volume


class TestInputs(unittest.TestCase):
    def test_volume_from_structure(self):
        self.assertAlmostEqual(volume_from_structure(test_poscar_filename), volume)

    def test_read_unitcell_data(self):
        self.assertAlmostEqual(read_unitcell_data(test_unitcell_filename), volume)

    def test_read_input_data(self):
        defect_data = read_input_data(test_sc_fermi_input_filename, volume=1)
        self.assertEqual(len(defect_data["defect_species"]), 2)
        self.assertEqual(defect_data["temperature"], 100)
        self.assertEqual(defect_data["bandgap"], 0.8084)
        self.assertEqual(defect_data["nelect"], 18)
        self.assertEqual(defect_data["spinpol"], 1)

    def test_read_dos_data(self):
        dos_data = read_dos_data(filename=test_dos_filename, bandgap=1, nelect=1)
        self.assertEqual(type(dos_data), DOS)

    def test_inputs_from_files(self):
        inputs = inputs_from_files(
            test_unitcell_filename,
            test_dos_filename,
            test_frozen_sc_fermi_input_filename,
            frozen=True,
        )
        self.assertEqual(len(inputs["defect_species"]), 4)
        self.assertEqual(inputs["temperature"], 300)
        self.assertEqual(inputs["bandgap"], 0.8084)
        self.assertEqual(inputs["nelect"], 18)
        self.assertEqual(inputs["spinpol"], 1)

    def test_dos_from_dict(self):
        dos_dict = {
            "dos": np.ones(101),
            "edos": np.linspace(-10.0, 10.0, 101),
            "bandgap": 1,
            "nelect": 10,
        }
        dos = dos_from_dict(dos_dict)
        self.assertEqual(dos_dict["bandgap"], dos.bandgap)
        self.assertEqual(dos_dict["nelect"], dos.nelect)
        self.assertEqual(False, dos.spin_polarised)
        np.testing.assert_equal(np.ones(101), dos.dos)
        np.testing.assert_equal(np.linspace(-10.0, 10.0, 101), dos.edos)

    def test_defect_species_from_dict(self):
         defect_dict = {"V_O" : {"nsites" :1, "charge_states" : {0 : {"degeneracy" : 1 , "formation_energy": 1}}}}
         defect_species = defect_species_from_dict(defect_dict, volume = 1)
         self.assertEqual(defect_species.name, "V_O")
         self.assertEqual(defect_species.nsites, 1)
         self.assertEqual(defect_species.charge_states[0].degeneracy, 1)
         self.assertEqual(defect_species.charge_states[0].energy, 1) 

    def test_defect_system_from_yaml(self):
        defect_system = defect_system_from_yaml(test_defect_system_yaml_filename)
        self.assertEqual(defect_system.dos.bandgap, 0.8084)
        self.assertEqual(defect_system.dos.nelect, 18)
        self.assertEqual(defect_system.volume, 59)
        self.assertEqual(defect_system.temperature, 300)
        self.assertEqual(len(defect_system.defect_species),3)


if __name__ == "__main__":
    unittest.main()
