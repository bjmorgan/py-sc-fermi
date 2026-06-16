"""Pool data shapes, reference normalisation, and serialisation for
``DefectSystem`` site and element pools.

These describe how pools are accepted from the caller (each species given
as a ``DefectSpecies`` object or by name), how they are stored on a
``DefectSystem`` (references reduced to names), and how they serialise to
JSON/YAML. The element-pool chemical-potential solver is in
``py_sc_fermi.element_pools``.
"""

from __future__ import annotations

from typing import TypeAlias, TypedDict

from py_sc_fermi.defect_species import DefectSpecies

# Pools as accepted by the constructor: each species given as a DefectSpecies
# object or by its name.
SitePoolsInput: TypeAlias = dict[str, tuple[float, list[DefectSpecies | str]]]
ElementPoolsInput: TypeAlias = dict[str, tuple[float, list[tuple[DefectSpecies | str, float]]]]

# Pools as stored on a DefectSystem: references reduced to names.
SitePools: TypeAlias = dict[str, tuple[float, list[str]]]
ElementPools: TypeAlias = dict[str, tuple[float, list[tuple[str, float]]]]

# Element pools with names resolved back to roster DefectSpecies (solve time).
ElementPoolsResolved: TypeAlias = dict[str, tuple[float, list[tuple[DefectSpecies, float]]]]


class SerialisedSitePool(TypedDict):
    """JSON/YAML form of one site pool: a named site budget and its species."""

    n_sites: float
    species: list[str]


class SerialisedElementMember(TypedDict):
    """JSON/YAML form of one (species, stoichiometry) member of an element pool."""

    species: str
    stoichiometry: float


class SerialisedElementPool(TypedDict):
    """JSON/YAML form of one element pool: a target content and its members."""

    target: float
    members: list[SerialisedElementMember]


def _species_name(species: DefectSpecies | str) -> str:
    """Reduce a species reference (a ``DefectSpecies`` or its name) to the name."""
    return species.name if isinstance(species, DefectSpecies) else species


def _normalise_site_pools(site_pools: SitePoolsInput | None) -> SitePools:
    """Reduce every site-pool species reference to a species name."""
    if not site_pools:
        return {}
    return {
        pool_name: (n_sites, [_species_name(sp) for sp in species_list])
        for pool_name, (n_sites, species_list) in site_pools.items()
    }


def _normalise_element_pools(element_pools: ElementPoolsInput | None) -> ElementPools:
    """Reduce every element-pool species reference to a species name."""
    if not element_pools:
        return {}
    return {
        element: (target, [(_species_name(sp), stoich) for sp, stoich in pool_list])
        for element, (target, pool_list) in element_pools.items()
    }


def _site_pools_as_dict(site_pools: SitePools) -> dict[str, SerialisedSitePool]:
    """Serialise stored site pools to JSON/YAML-safe mappings."""
    return {
        pool_name: {"n_sites": float(n_sites), "species": list(species_names)}
        for pool_name, (n_sites, species_names) in site_pools.items()
    }


def _element_pools_as_dict(
    element_pools: ElementPools,
) -> dict[str, SerialisedElementPool]:
    """Serialise stored element pools to JSON/YAML-safe mappings."""
    return {
        element: {
            "target": float(target),
            "members": [
                {"species": name, "stoichiometry": float(stoich)}
                for name, stoich in pool_list
            ],
        }
        for element, (target, pool_list) in element_pools.items()
    }
