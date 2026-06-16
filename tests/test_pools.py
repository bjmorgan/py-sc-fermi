import math
import unittest

from py_sc_fermi.defect_charge_state import DefectChargeState
from py_sc_fermi.defect_species import DefectSpecies
from py_sc_fermi.pools import ElementPool, SitePool


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
        with self.assertRaisesRegex(ValueError, "n_sites must be finite and > 0"):
            SitePool(n_sites=0.0, species=["A"])

    def test_rejects_negative_n_sites(self):
        with self.assertRaisesRegex(ValueError, "n_sites must be finite and > 0"):
            SitePool(n_sites=-1.0, species=["A"])

    def test_rejects_non_finite_n_sites(self):
        for bad in (math.inf, math.nan):
            with self.assertRaisesRegex(ValueError, "n_sites must be finite and > 0"):
                SitePool(n_sites=bad, species=["A"])

    def test_rejects_empty_species(self):
        with self.assertRaisesRegex(ValueError, "at least one species"):
            SitePool(n_sites=1.0, species=[])

    def test_not_equal_to_foreign_type(self):
        self.assertNotEqual(SitePool(n_sites=1.0, species=["A"]), object())

    def test_equality_by_value(self):
        self.assertEqual(
            SitePool(n_sites=2.0, species=["A", "B"]),
            SitePool(n_sites=2.0, species=["A", "B"]),
        )
        self.assertNotEqual(
            SitePool(n_sites=2.0, species=["A"]),
            SitePool(n_sites=3.0, species=["A"]),
        )


class TestElementPool(unittest.TestCase):
    def test_normalises_members_to_names(self):
        pool = ElementPool(target=5.0, members={_species("A"): 1.0, "B": -2.0})
        self.assertEqual(pool.target, 5.0)
        self.assertEqual(pool.members, {"A": 1.0, "B": -2.0})

    def test_members_property_returns_a_copy(self):
        pool = ElementPool(target=1.0, members={"A": 1.0})
        pool.members["X"] = 9.0
        self.assertEqual(pool.members, {"A": 1.0})

    def test_accepts_negative_target(self):
        self.assertEqual(ElementPool(target=-3.0, members={"A": -1.0}).target, -3.0)

    def test_rejects_nan_target(self):
        with self.assertRaisesRegex(ValueError, "target must be finite"):
            ElementPool(target=math.nan, members={"A": 1.0})

    def test_rejects_inf_target(self):
        with self.assertRaisesRegex(ValueError, "target must be finite"):
            ElementPool(target=math.inf, members={"A": 1.0})

    def test_object_and_name_for_same_species_raise(self):
        with self.assertRaisesRegex(ValueError, "resolve to the same species name: A"):
            ElementPool(target=1.0, members={_species("A"): 1.0, "A": -1.0})

    def test_two_distinct_objects_sharing_a_name_raise(self):
        with self.assertRaisesRegex(ValueError, "resolve to the same species name: A"):
            ElementPool(target=1.0, members={_species("A"): 1.0, _species("A"): -1.0})

    def test_rejects_non_finite_stoichiometry(self):
        for bad in (math.nan, math.inf):
            with self.assertRaisesRegex(ValueError, "stoichiometries must be finite"):
                ElementPool(target=1.0, members={"A": bad})

    def test_rejects_empty_members(self):
        with self.assertRaisesRegex(ValueError, "at least one member"):
            ElementPool(target=1.0, members={})

    def test_not_equal_to_foreign_type(self):
        self.assertNotEqual(ElementPool(target=1.0, members={"A": 1.0}), object())

    def test_equality_by_value(self):
        self.assertEqual(
            ElementPool(target=1.0, members={"A": 1.0}),
            ElementPool(target=1.0, members={"A": 1.0}),
        )
        self.assertNotEqual(
            ElementPool(target=1.0, members={"A": 1.0}),
            ElementPool(target=2.0, members={"A": 1.0}),
        )


class TestPoolSerialisation(unittest.TestCase):
    def test_sitepool_round_trips_through_as_dict(self):
        pool = SitePool(n_sites=4.0, species=[_species("A"), "B"])
        self.assertEqual(pool.as_dict(), {"n_sites": 4.0, "species": ["A", "B"]})
        self.assertEqual(SitePool.from_dict(pool.as_dict()), pool)

    def test_elementpool_round_trips_with_members_mapping(self):
        pool = ElementPool(target=0.3, members={_species("C"): 1.0, "D": -2.0})
        self.assertEqual(
            pool.as_dict(), {"target": 0.3, "members": {"C": 1.0, "D": -2.0}}
        )
        self.assertEqual(ElementPool.from_dict(pool.as_dict()), pool)
