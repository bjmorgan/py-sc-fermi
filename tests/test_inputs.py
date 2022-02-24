import unittest
from numpy.testing import assert_almost_equal
import os

from py_sc_fermi.inputs import volume_from_structure, read_unitcell_data
from pymatgen.core import Structure

test_data_dir = 'inputs'
test_poscar_filename = os.path.join( os.path.dirname( __file__ ), test_data_dir, 'POSCAR' )

structure = Structure.from_file(test_poscar_filename)
volume = Structure.volume

def test_volume_from_structure():
    assert volume == get_volume_from_structure(test_poscar_filename)

# def test_read_unitcell_data():
#     assert_almost_equal(volume,read_unitcell_data('inputs/unitcell.dat'))
    
if __name__ == '__main__':
    unittest.main()
