import unittest
from unittest.mock import patch, Mock, mock_open
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
        edos = np.linspace(-10.0, 10.0, 100)
        bandgap = 3.0
        # Semiconducting DOS: density of 1 in the valence band (E <= 0) and
        # conduction band (E >= bandgap), zero inside the gap:
        dos_data = np.where((edos <= 0) | (edos >= bandgap), 1.0, 0.0)
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
        expected = np.where(
            (self.dos.edos <= 0) | (self.dos.edos >= self.dos.bandgap), 1.0, 0.0
        )
        np.testing.assert_allclose(self.dos.dos, expected)

    def test_sum_dos(self):
        edos = np.linspace(-10.0, 10.0, 100)
        bandgap = 1.0
        dos_data = np.where((edos <= 0) | (edos >= bandgap), 1.0, 0.0)
        dos = DOS(dos=dos_data, edos=edos, bandgap=bandgap, nelect=10)
        dos.sum_dos()
        np.testing.assert_allclose(dos.dos, dos_data)

    def test__p0_idx(self):
        self.assertEqual(self.dos._p0_idx, 49)

    def test__n0_idx(self):
        self.assertEqual(self.dos._n0_idx, 64)

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

    def test_carrier_concentrations_robust_to_band_edge_index(self):
        # carrier concentrations should be insensitive to the exact, noisy,
        # auto-determined VBM/CBM indices. Here the physical DOS is held fixed
        # (band edges at 0 and 3 eV) while the supplied bandgap is perturbed,
        # which shifts the determined CBM index (_n0_idx) by two grid points,
        # however the result should be unchanged:
        edos = np.linspace(-10.0, 10.0, 100)
        dos_data = np.where((edos <= 0) | (edos >= 3.0), 1.0, 0.0)
        nelect = 10
        e_fermi, T = 1.5, 298

        result_exact = DOS(
            dos=dos_data, edos=edos, bandgap=3.0, nelect=nelect
        ).carrier_concentrations(e_fermi, T)
        result_noisy = DOS(
            dos=dos_data, edos=edos, bandgap=3.3, nelect=nelect
        ).carrier_concentrations(e_fermi, T)

        # sanity-check that the perturbed bandgap really does move the
        # auto-determined band-edge index that the integration bound is built on
        self.assertNotEqual(
            DOS(dos=dos_data, edos=edos, bandgap=3.0, nelect=nelect)._n0_idx,
            DOS(dos=dos_data, edos=edos, bandgap=3.3, nelect=nelect)._n0_idx,
        )
        np.testing.assert_allclose(result_exact, result_noisy, rtol=1e-3)

    def test_integration_indices(self):
        # default fixture (bandgap 3.0): p0_idx=49, n0_idx=64,
        # mid_gap=(49+64)//2=56, so integration runs to mid-gap.
        self.assertEqual(self.dos._p0_integration_idx, 56)
        self.assertEqual(self.dos._n0_integration_idx, 57)

    def test_integration_indices_narrow_bandgap_clamping(self):
        # For a bandgap narrower than one grid spacing, mid-gap collapses onto
        # the band-edge indices and the max(...)/min(...) clamps keep the hole
        # and electron integration ranges valid and non-overlapping.
        edos = np.linspace(-10.0, 10.0, 100)
        bandgap = 0.1  # narrower than the ~0.2 eV grid spacing
        dos_data = np.where((edos <= 0) | (edos >= bandgap), 1.0, 0.0)
        dos = DOS(dos=dos_data, edos=edos, bandgap=bandgap, nelect=10)
        self.assertEqual(dos._p0_idx, 49)
        self.assertEqual(dos._n0_idx, 50)
        # mid_gap=(49+50)//2 = 49 alone would make the hole range end before
        # p0_idx, the clamps push both integration bounds to 50:
        self.assertEqual(dos._p0_integration_idx, 50)  # (min, 50]
        self.assertEqual(dos._n0_integration_idx, 50)  # [50, max)

    def test_edos_must_bracket_zero(self):
        edos = np.linspace(1.0, 10.0, 100)  # does not bracket E=0
        dos_data = np.ones_like(edos)
        with self.assertRaises(ValueError) as context:
            DOS(dos=dos_data, edos=edos, bandgap=3.0, nelect=10)
        self.assertIn("DOS energy range must bracket zero", str(context.exception))

    def test_negative_bandgap_rejected(self):
        edos = np.linspace(-10.0, 10.0, 100)
        dos_data = np.ones_like(edos)
        with self.assertRaises(ValueError) as context:
            DOS(dos=dos_data, edos=edos, bandgap=-1.0, nelect=10)
        self.assertIn("bandgap must be non-negative", str(context.exception))

    def test_from_vasprun(self):
        dos = self.dos.from_vasprun(test_vasprun_filename, nelect=320)
        self.assertEqual(dos.nelect, 320)
        self.assertEqual(dos.bandgap, 8.7342)
        
    def test_from_vasprun_type_error(self):
        mock_vr = Mock()
        mock_vr.eigenvalue_band_properties = (1.0, 2.0, 3.0, 4.0)  # last value should be bool
        with patch(
            'py_sc_fermi.dos.Vasprun', autospec=True
        ) as mock_Vasprun:
            mock_Vasprun.return_value = mock_vr
            with self.assertRaises(TypeError) as context:
                DOS.from_vasprun("dummy_path.xml")
        self.assertIn("Expected tuple[float, float, float, bool]", str(context.exception))
        mock_Vasprun.assert_called_once_with(
            "dummy_path.xml",
            parse_potcar_file=False,
            separate_spins=False
        )
        
    def test_from_vasprun_type_error_spin_separated(self):
        mock_vr = Mock()
        mock_vr.eigenvalue_band_properties = (
            (1.0, 1.0),
            (2.0, 2.0),
            (3.0, 3.0),
            (True, True)
        )
        with patch(
            'py_sc_fermi.dos.Vasprun', autospec=True
        ) as mock_Vasprun:
            mock_Vasprun.return_value = mock_vr
            with self.assertRaises(TypeError) as context:
                DOS.from_vasprun("dummy_path.xml")
        self.assertIn("Expected tuple[float, float, float, bool]", str(context.exception))
        mock_Vasprun.assert_called_once_with(
            "dummy_path.xml",
            parse_potcar_file=False,
            separate_spins=False
        )

    def test_from_dict(self):
        edos = np.linspace(-10.0, 10.0, 100)
        bandgap = 3.0
        dos_data = np.where((edos <= 0) | (edos >= bandgap), 1.0, 0.0)
        dos = self.dos.from_dict(
            {
                "dos": dos_data,
                "edos": edos,
                "bandgap": bandgap,
                "nelect": 10,
                "spin_pol": False,
            }
        )
        self.assertEqual(dos.nelect, 10)
        self.assertEqual(dos.bandgap, bandgap)
        np.testing.assert_allclose(dos.dos, dos_data)
        self.assertEqual(dos.spin_polarised, False)

    def test_from_dict_with_spin_polarised(self):
        edos = np.linspace(-10.0, 10.0, 100)
        bandgap = 3.0
        dos_data = np.where((edos <= 0) | (edos >= bandgap), 1.0, 0.0)
        dos = self.dos.from_dict(
            {
                "dos": np.array([dos_data / 2, dos_data / 2]),
                "edos": edos,
                "bandgap": bandgap,
                "nelect": 10,
                "spin_pol": True,
            }
        )
        self.assertEqual(dos.nelect, 10)
        self.assertEqual(dos.bandgap, bandgap)
        np.testing.assert_allclose(dos.dos, dos_data)
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
