import unittest

from py_sc_fermi.defect_charge_state import DefectChargeState
from py_sc_fermi.defect_species import DefectSpecies
from py_sc_fermi.pools import SitePool


def _species(name: str) -> DefectSpecies:
    return DefectSpecies(
        name,
        nsites=4,
        charge_states=[DefectChargeState(charge=0, energy=1.0, degeneracy=1)],
    )


class TestSitePool(unittest.TestCase):
    def test_normalises_species_to_names(self):
        pool = SitePool(n_sites=10.0, species=[_species("A"), "B"])
        self.assertEqual(pool.n_sites, 10.0)
        self.assertEqual(pool.species, ["A", "B"])

    def test_species_property_returns_a_copy(self):
        pool = SitePool(n_sites=1.0, species=["A"])
        pool.species.append("X")
        self.assertEqual(pool.species, ["A"])

    def test_rejects_zero_n_sites(self):
        with self.assertRaisesRegex(ValueError, "n_sites must be > 0"):
            SitePool(n_sites=0.0, species=["A"])

    def test_rejects_negative_n_sites(self):
        with self.assertRaisesRegex(ValueError, "n_sites must be > 0"):
            SitePool(n_sites=-1.0, species=["A"])

    def test_equality_by_value(self):
        self.assertEqual(
            SitePool(n_sites=2.0, species=["A", "B"]),
            SitePool(n_sites=2.0, species=["A", "B"]),
        )
        self.assertNotEqual(
            SitePool(n_sites=2.0, species=["A"]),
            SitePool(n_sites=3.0, species=["A"]),
        )
