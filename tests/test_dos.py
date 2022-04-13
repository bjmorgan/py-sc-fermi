import unittest
from unittest import mock
from unittest.mock import Mock, patch
import numpy as np
import os
from py_sc_fermi.dos import DOS

test_data_dir = "inputs/"
test_vasprun_filename = os.path.join(
    os.path.dirname(__file__), test_data_dir, "vasprun_nsp.xml"
)


class TestDOSInit(unittest.TestCase):
    def test_DOS_is_initialised(self):
        dos_data = np.random.random(100)
        edos = np.linspace(-10.0, 10.0, 100)
        bandgap = 3.0
        nelect = 10
        with patch(
            "py_sc_fermi.dos.DOS.normalise_dos", autospec=True
        ) as mock_normalise_dos:
            dos = DOS(dos=dos_data, edos=edos, bandgap=bandgap, nelect=nelect)
            self.assertEqual(mock_normalise_dos.call_count, 1)
        np.testing.assert_equal(dos._dos, dos_data)
        np.testing.assert_equal(dos._edos, edos)
        self.assertEqual(dos._bandgap, bandgap)
        self.assertEqual(dos._nelect, nelect)


class TestDos(unittest.TestCase):
    def setUp(self):
        dos_data = np.ones(101)
        edos = np.linspace(-10.0, 10.0, 101)
        bandgap = 3.0
        nelect = 10
        self.dos = DOS(dos=dos_data, edos=edos, bandgap=bandgap, nelect=nelect)

    def test_dos_property(self):
        np.testing.assert_equal(self.dos.dos, self.dos._dos)

    def test_edos_property(self):
        np.testing.assert_equal(self.dos.edos, self.dos._edos)

    def test_bandgap_property(self):
        self.assertEqual(self.dos.bandgap, self.dos._bandgap)

    def test_nelect_property(self):
        self.assertEqual(self.dos.nelect, self.dos._nelect)

    def test_normalise_dos(self):
        # The integral of rho(E) = 1, from -10.0 to 0, is equal to 10.0
        # With nelect = 10.0, the normaliss_dos() method should leave
        # all dos values = 1.0
        dos_data = np.ones(101)
        np.testing.assert_equal(self.dos.dos, dos_data)

    def test_emin(self):
        self.assertEqual(self.dos.emin(), -10)

    def test_emax(self):
        self.assertEqual(self.dos.emax(), 10)

    def test_carrier_concentrations(self):
        np.testing.assert_almost_equal(
            self.dos.carrier_concentrations(1.5, 298),
            4.28893131265779e-27,
            1.7780649634855188e-30,
        )

    def test_from_vasprun(self):
        dos = self.dos.from_vasprun(test_vasprun_filename, nelect=320)
        self.assertEqual(dos.nelect, 320)
        self.assertEqual(dos.bandgap, 8.7342)


if __name__ == "__main__":
    unittest.main()
