"""Pool classes for ``DefectSystem`` site and element constraints.

``SitePool`` and ``ElementPool`` are the input and stored form of a pool: each
accepts its species as ``DefectSpecies`` objects or by name and reduces them to
names, and each serialises to JSON/YAML via ``as_dict``/``from_dict``.
``ResolvedElementPool`` is the internal solve-time form, with names resolved
back to roster ``DefectSpecies``. The element-pool chemical-potential solver is
in ``py_sc_fermi.element_pools``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TypedDict

from py_sc_fermi.defect_species import DefectSpecies


class SerialisedSitePool(TypedDict):
    """JSON/YAML form of one site pool: a named site budget and its species."""

    n_sites: float
    species: list[str]


class SerialisedElementPool(TypedDict):
    """JSON/YAML form of one element pool: a target content and its members."""

    target: float
    members: dict[str, float]


def _species_name(species: DefectSpecies | str) -> str:
    """Reduce a species reference (a ``DefectSpecies`` or its name) to the name."""
    return species.name if isinstance(species, DefectSpecies) else species


class SitePool:
    """A named site budget shared by several defect species.

    Args:
        n_sites: total number of sites in this pool. Must be > 0.
        species: the defect species sharing the pool, each given as a
            ``DefectSpecies`` or by name. Stored reduced to names.
    """

    def __init__(self, n_sites: float, species: list[str | DefectSpecies]):
        self._init_canonical(n_sites, [_species_name(s) for s in species])

    def _init_canonical(self, n_sites: float, species: list[str]) -> None:
        """Validate and store a pool whose species are already names.

        Shared by ``__init__`` (after it reduces its references to names) and by
        ``from_dict`` (whose serialised species are already names), so both
        construct through the same validation.
        """
        if not (n_sites > 0):
            raise ValueError(f"SitePool n_sites must be > 0; got {n_sites}.")
        self._n_sites = n_sites
        self._species: list[str] = species

    @property
    def n_sites(self) -> float:
        """Total number of sites in the pool."""
        return self._n_sites

    @property
    def species(self) -> list[str]:
        """The pooled species, by name."""
        return list(self._species)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, SitePool):
            return NotImplemented
        return self._n_sites == other._n_sites and self._species == other._species

    def __repr__(self) -> str:
        return f"SitePool(n_sites={self._n_sites!r}, species={self._species!r})"

    def as_dict(self) -> SerialisedSitePool:
        """JSON/YAML-safe mapping: ``{"n_sites", "species"}`` (species by name)."""
        return {"n_sites": float(self._n_sites), "species": list(self._species)}

    @classmethod
    def from_dict(cls, d: SerialisedSitePool) -> SitePool:
        """Reconstruct a ``SitePool`` from ``as_dict`` output."""
        pool = cls.__new__(cls)
        pool._init_canonical(d["n_sites"], list(d["species"]))
        return pool


class ElementPool:
    """A net-content target for one element and the species supplying it.

    Args:
        target: the element content the pool is solved to (per unit cell). May
            be negative (intrinsic off-stoichiometry); only NaN/inf are rejected.
        members: mapping of each supplying species (a ``DefectSpecies`` or its
            name) to its stoichiometry for this element. Stored keyed by name.
    """

    def __init__(self, target: float, members: dict[str | DefectSpecies, float]):
        names = [_species_name(s) for s in members]
        duplicates = sorted({n for n in names if names.count(n) > 1})
        if duplicates:
            raise ValueError(
                "ElementPool members resolve to the same species name: "
                + ", ".join(duplicates)
            )
        self._init_canonical(target, dict(zip(names, members.values(), strict=True)))

    def _init_canonical(self, target: float, members: dict[str, float]) -> None:
        """Validate and store a pool whose members are already keyed by name.

        Shared by ``__init__`` (after it reduces its references to names) and by
        ``from_dict`` (whose serialised members are already names), so both
        construct through the same validation.
        """
        if not math.isfinite(target):
            raise ValueError(f"ElementPool target must be finite; got {target}.")
        self._target = target
        self._members: dict[str, float] = members

    @property
    def target(self) -> float:
        """The element content the pool is solved to."""
        return self._target

    @property
    def members(self) -> dict[str, float]:
        """Mapping of pooled species name to stoichiometry."""
        return dict(self._members)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ElementPool):
            return NotImplemented
        return self._target == other._target and self._members == other._members

    def __repr__(self) -> str:
        return f"ElementPool(target={self._target!r}, members={self._members!r})"

    def as_dict(self) -> SerialisedElementPool:
        """JSON/YAML-safe mapping: ``{"target", "members": {name: stoich}}``."""
        return {
            "target": float(self._target),
            "members": {name: float(stoich) for name, stoich in self._members.items()},
        }

    @classmethod
    def from_dict(cls, d: SerialisedElementPool) -> ElementPool:
        """Reconstruct an ``ElementPool`` from ``as_dict`` output."""
        pool = cls.__new__(cls)
        pool._init_canonical(d["target"], dict(d["members"]))
        return pool


@dataclass(frozen=True)
class ResolvedElementPool:
    """An element pool with names resolved to roster ``DefectSpecies``.

    Built internally by ``DefectSystem._resolve_element_pools`` from a validated
    ``ElementPool``; carries no validation of its own.
    """

    target: float
    members: dict[DefectSpecies, float]
