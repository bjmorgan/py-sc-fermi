import unittest
from numpy.testing import assert_almost_equal

from py_sc_fermi.inputs import volume_from_structure, read_unitcell_data
from pymatgen.core import Structure

structure = Structure.from_file('inputs/POSCAR')
volume = Structure.volume

def test_volume_from_structure():
    assert volume == get_volume_from_structure('inputs/POSCAR')

def test_read_unitcell_data():
    assert_almost_equal(volume,read_unitcell_data('inputs/unitcell.dat'))
    
if __name__ == '__main__':
    unittest.main()
