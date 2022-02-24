import unittest
from numpy.testing import assert_almost_equal
import os

from py_sc_fermi.inputs import volume_from_structure, read_unitcell_data
from pymatgen.core import Structure

test_data_dir = 'inputs'
test_poscar_filename = os.path.join( os.path.dirname( __file__ ), test_data_dir, 'POSCAR' )
test_unitcell_filename = os.path.join( os.path.dirname( __file__ ), test_data_dir, 'unitcell.dat' )

structure = Structure.from_file(test_poscar_filename)
volume = Structure.volume

def test_volume_from_structure():
    assert_almost_equal(volume,volume_from_structure(test_poscar_filename))

def test_read_unitcell_data():
     assert_almost_equal(volume,read_unitcell_data(test_unitcell_filename))
    
if __name__ == '__main__':
    unittest.main()
