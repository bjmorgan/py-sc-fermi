import unittest

from py_sc_fermi.inputs import volume_from_structure
from pymatgen.core import Structure

structure = Structure.from_file('/examples/v_F/POSCAR')
volume = Structure.volume

def test_volume_from_structure():
    assert volume == get_volume_from_structure('/examples/v_F/POSCAR')

if __name__ == '__main__':
    unittest.main()
