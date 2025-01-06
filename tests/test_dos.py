import unittest
from unittest.mock import patch
import numpy as np
import os
from py_sc_fermi.dos import DOS

test_data_dir = "dummy_inputs/"
test_vasprun_filename = os.path.join(
    os.path.dirname(__file__), test_data_dir, "vasprun_nsp.xml"
)


class TestDOSInit(unittest.TestCase):
    def setUp(self):
        self.dos_data = np.random.random(100)
        self.edos = np.linspace(-10.0, 10.0, 100)
        self.bandgap = 3.0
        self.nelect = 10

    def test_initialisation_calls_normalise_dos(self):
        with patch(
            "py_sc_fermi.dos.DOS.normalise_dos", autospec=True
        ) as mock_normalise_dos:
            dos = DOS(dos=self.dos_data, edos=self.edos, bandgap=self.bandgap, nelect=self.nelect)
            self.assertEqual(mock_normalise_dos.call_count, 1)
        np.testing.assert_equal(dos._dos, self.dos_data)
        np.testing.assert_equal(dos._edos, self.edos)
        self.assertEqual(dos._bandgap, self.bandgap)
        self.assertEqual(dos._nelect, self.nelect)


class TestDos(unittest.TestCase):
    def setUp(self):
        dos_data = np.ones(101)
        edos = np.linspace(-10.0, 10.0, 101)
        bandgap = 3.0
        nelect = 10
        spin_polarised = False
        self.dos = DOS(
            dos=dos_data,
            edos=edos,
            bandgap=bandgap,
            nelect=nelect,
            spin_polarised=spin_polarised,
        )

    def tearDown(self):
        del self.dos

    def test_dos_property(self):
        np.testing.assert_equal(self.dos.dos, self.dos._dos)

    def test_edos_property(self):
        np.testing.assert_equal(self.dos.edos, self.dos._edos)

    def test_bandgap_property(self):
        self.assertEqual(self.dos.bandgap, self.dos._bandgap)

    def test_nelect_property(self):
        self.assertEqual(self.dos.nelect, self.dos._nelect)

    def test_spin_polarised_property(self):
        self.assertFalse(self.dos.spin_polarised)

    def test_normalise_dos(self):
        # The integral of rho(E) = 1, from -10.0 to 0, is equal to 10.0
        # With nelect = 10.0, the normaliss_dos() method should leave
        # all dos values = 1.0
        dos_data = np.ones(101)
        np.testing.assert_equal(self.dos.dos, dos_data)

    def test_sum_dos(self):
        dos_data = np.ones(101)
        dos = DOS(
            dos=dos_data, edos=np.linspace(-10.0, 10.0, 101), bandgap=1, nelect=10
        )
        dos.sum_dos()
        np.testing.assert_equal(dos.dos, dos_data)

    def test__p0_index(self):
        self.assertEqual(self.dos._p0_index(), 50)

    def test__n0_index(self):
        self.assertEqual(self.dos._n0_index(), 65)

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

    def test_from_dict(self):
        dos = self.dos.from_dict(
            {
                "dos": np.ones(101),
                "edos": np.linspace(-10.0, 10.0, 101),
                "bandgap": 3.0,
                "nelect": 10,
                "spin_pol": False
            }
        )
        self.assertEqual(dos.nelect, 10)
        self.assertEqual(dos.bandgap, 3.0)
        np.testing.assert_equal(dos.dos, np.ones(101))
        self.assertEqual(dos.spin_polarised, False)

    def test_from_dict_with_spin_polarised(self):
        dos = self.dos.from_dict(
            {
                "dos": np.array([np.ones(101), np.ones(101)]),
                "edos": np.linspace(-10.0, 10.0, 101),
                "bandgap": 3.0,
                "nelect": 10,
                "spin_pol": True
            }
        )
        self.assertEqual(dos.nelect, 10)
        self.assertEqual(dos.bandgap, 3.0)
        np.testing.assert_equal(
            dos.dos, np.sum([np.ones(101) / 2, np.ones(101) / 2], axis=0)
        )
        self.assertEqual(dos.spin_polarised, True)

    def test_as_dict(self):
        self.dos._dos = [1,2,3,4,5]
        self.dos._edos = [1,2,3,4,5]
        dictionary = self.dos.as_dict()
        self.assertEqual(dictionary["spin_pol"],False)
        self.assertEqual(dictionary["nelect"], 10)
        self.assertEqual(dictionary["bandgap"], 3)
        self.assertEqual(dictionary["edos"], [1,2,3,4,5])
        self.assertEqual(dictionary["dos"], [1,2,3,4,5])


if __name__ == "__main__":
    unittest.main()
