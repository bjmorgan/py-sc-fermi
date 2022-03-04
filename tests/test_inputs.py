import unittest
from numpy.testing import assert_almost_equal, assert_equal
from py_sc_fermi.defect_charge_state import DefectChargeState
from py_sc_fermi.dos import DOS
from py_sc_fermi.defect_species import DefectSpecies
import os

from py_sc_fermi.inputs import (
    volume_from_structure,
    read_unitcell_data,
    read_input_data,
    read_dos_data,
    inputs_from_files,
    read_defect_species,
    update_frozen_defect_species,
    update_frozen_chgstates
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
    os.path.dirname(__file__), test_data_dir, "frozen_chg_states.dat"
)
test_dos_filename = os.path.join(os.path.dirname(__file__), test_data_dir, "totdos.dat")

structure = Structure.from_file(test_poscar_filename)
volume = structure.volume

class TestInputs(unittest.TestCase):
    def test_volume_from_structure(self):
        assert_almost_equal(volume_from_structure(test_poscar_filename), volume)

    def test_read_unitcell_data(self):
        assert_almost_equal(
            read_unitcell_data(test_unitcell_filename), volume
        )

    def test_read_input_data(self):
        defect_data = read_input_data(test_sc_fermi_input_filename, volume=1)
        assert_equal(len(defect_data["defect_species"]), 2)
        assert_equal(defect_data["temperature"], 100)
        assert_equal(defect_data["egap"], 0.8084)
        assert_equal(defect_data["nelect"], 18)
        assert_equal(defect_data["nspinpol"], 1)

    def test_read_dos_data(self):
        dos_data = read_dos_data(test_dos_filename, egap=1, nelect=1)
        assert type(dos_data) == DOS

    def test_inputs_from_files(self):
        inputs = inputs_from_files(
            test_unitcell_filename,
            test_dos_filename,
            test_frozen_sc_fermi_input_filename,
            frozen=True,
        )
        assert inputs

    def test_read_defect_species(self):
        string = "V_Ga 2 1\n0  2.4451 1\n-1  0.0265 1\nGa_Sb 2 1\n0  2.2649 1\n-1  2.0937 1"
        string = string.splitlines()
        defects = read_defect_species(string,2)
        assert len(defects) == 2
        assert defects[1].name == 'Ga_Sb'
        assert defects[1].nsites == 1

    def test_update_frozen_defect_species(self):
        string = "V_Ga  0.3285677364522E+20"
        string = string.splitlines()
        defects = [DefectSpecies('V_Ga', 1, [])]
        read_frozen_defect_species(string,defects,1,1)
        assert_almost_equal(defects[0].get_concentration(0.1, 298), 3.285677364522e-05)

    def test_update_frozen_chgstates(self):
        string = "V_Ga -1 0.19E+19\nGa_i 1 0.5E+20"
        string = string.splitlines()
        defects = [DefectSpecies('V_Ga', 1, [DefectChargeState(-1,1,1)]), DefectSpecies('Ga_i', 1, [DefectChargeState(1,1,1)])]
        read_frozen_chgstates(string,defects,1,2)
        assert_almost_equal(defects[0].get_concentration(0.1,300), 1.9e-06)
        assert_almost_equal(defects[1].get_concentration(0.1,300), 5e-05)


if __name__ == "__main__":
    unittest.main()
