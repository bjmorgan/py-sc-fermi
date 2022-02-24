import unittest
from unittest.mock import Mock, patch
import numpy as np
from py_sc_fermi.dos import DOS

class TestDOSInit(unittest.TestCase):

    def test_DOS_is_initialised(self):
        dos_data = np.random.random(100)
        edos = np.linspace(-10.0, 10.0, 100)
        egap = 3.0
        nelect = 10
        with patch('py_sc_fermi.dos.DOS.normalise_dos', autospec=True) as mock_normalise_dos:
            dos = DOS(dos=dos_data,
                     edos=edos,
                     egap=egap,
                     nelect=nelect)
            self.assertEqual( mock_normalise_dos.call_count, 1 )
        np.testing.assert_equal( dos._dos, dos_data )
        np.testing.assert_equal( dos._edos, edos )
        self.assertEqual( dos._egap, egap )
        self.assertEqual( dos._nelect, nelect )

class TestDos(unittest.TestCase):

    def setUp(self):
        dos_data = np.ones(101)
        edos = np.linspace(-10.0, 10.0, 101)
        egap = 3.0
        nelect = 10       
        self.dos = DOS( dos=dos_data, edos=edos, egap=egap, nelect=nelect )

    def test_dos_property(self):
        np.testing.assert_equal( self.dos.dos, self.dos._dos )

    def test_edos_property(self):
        np.testing.assert_equal( self.dos.edos, self.dos._edos )

    def test_egap_property(self):
        self.assertEqual( self.dos.egap, self.dos._egap )

    def test_nelect_property(self):
        self.assertEqual( self.dos.nelect, self.dos._nelect )

    def test_normalise_dos(self):
        # The integral of rho(E) = 1, from -10.0 to 0, is equal to 10.0
        # With nelect = 10.0, the normaliss_dos() method should leave
        # all dos values = 1.0
        dos_data = np.ones(101)
        edos = np.linspace(-10.0, 10.0, 101)
        egap = 3.0
        nelect = 10
        dos = DOS(dos=dos_data,
                  edos=edos,
                  egap=egap,
                  nelect=nelect)
        np.testing.assert_equal(dos.dos, dos_data)
        
     def test_emin(self):
        dos_data = np.ones(101)
        edos = np.linspace(-10.0, 10.0, 101)
        egap = 3.0
        nelect = 10
        dos = DOS(dos=dos_data,
                  edos=edos,
                  egap=egap,
                  nelect=nelect)
        np.testing.assert_equal(dos.emin, -10)

if __name__ == '__main__':
    unittest.main()
