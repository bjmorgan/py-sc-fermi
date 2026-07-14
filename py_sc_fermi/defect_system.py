from __future__ import annotations

import copy
import math
import sys
import warnings
from collections import Counter
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Any, NamedTuple

import numpy as np
from scipy.constants import physical_constants as _physical_constants
from scipy.optimize import brentq

from py_sc_fermi import element_pools
from py_sc_fermi.defect_charge_state import DefectChargeState
from py_sc_fermi.defect_species import DefectSpecies
from py_sc_fermi.defect_system_result import DefectSystemResult
from py_sc_fermi.dos import DOS
from py_sc_fermi.element_pools import ElementPoolError
from py_sc_fermi.pools import (
    ElementPool,
    ResolvedElementPool,
    SitePool,
)
from py_sc_fermi.warnings import DiluteLimitWarning

_kboltz = _physical_constants["Boltzmann constant in eV/K"][0]

CorrectionKey = DefectChargeState | tuple[str, str]
"""A formation-energy-correction key: a ``DefectChargeState`` object, or a
``(species_name, charge_state_name)`` pair resolved against the system's
defect species."""

FixedConcentrationKey = str | tuple[str, str]
"""A fixed-concentration key: a species name (fixes the species' total
concentration), or a ``(species_name, charge_state_name)`` pair (fixes that
single charge state)."""


class _VariableState(NamedTuple):
    """A non-fixed-concentration charge state within an exclusion group,
    paired with its per-site Boltzmann weight.

    ``log_w = log(degeneracy) - Ef / (kB * T)`` is the *per-site* weight
    (independent of `nsites`/pool size), with `Ef` the charge state's
    formation energy at the current Fermi level.
    """

    charge_state: DefectChargeState
    species: DefectSpecies
    log_w: float


@dataclass
class _ExclusionGroup:
    """A set of charge states competing for a shared budget of `n_free`
    sites under Langmuir/site-exclusion statistics.

    `variable_states` holds one `_VariableState` per non-fixed-concentration
    charge state in the group.
    """

    n_free: float
    variable_states: list[_VariableState]


class DefectSystem:
    """This class is used to calculate the self consistent Fermi energy for
    a defective material, observing the condition of charge neutrality and
    therefore, point defect and carrier concentrations under equilibrium
    conditions.

    Args:
        defect_species (list[DefectSpecies]): List of ``DefectSpecies`` objects
          which are present in the ``DefectSystem``.
        volume (float): volume of the unit cell in Angstroms cubed
        dos (DOS): the ``DOS`` object associated with the unit cell
        temperature (float): temperature at which self-consistent Fermi energy
          will be solved for.
        convergence_tolerance (float, optional): Tolerance for the Fermi energy
          convergence in eV. If not specified, uses scipy's default.
        site_pools (dict[str, SitePool] | None, optional):
          Mapping of pool name to a `SitePool` (its total site count and the
          species sharing those sites, each given as a `DefectSpecies` object
          or by name and reduced to names by `SitePool`).
          By default (None), every DefectSpecies gets its own implicit
          exclusion group of `nsites` sites, so its charge states already
          compete with each other via Langmuir statistics; `site_pools` is
          only needed when several species must share one physical site
          budget.
        element_pools (dict[str, ElementPool] | None, optional):
          Mapping of element name to an `ElementPool` (its `target` content
          and `{species: stoichiometry}` members). Constrains the total amount
          of an element supplied by the member species: the chemical
          potentials are solved so that
          ``sum_i stoichiometry_i * concentration_i`` equals the target
          for every constrained element (a fixed-budget / closed-system
          constraint, as distinct from a fixed thermodynamic chemical
          potential). Each species may be given as a `DefectSpecies`
          object or by name, and one species may appear in several pools.
          The target is a content per unit cell on the same scale as the
          concentrations: fixed-concentration states count against it, and
          the remainder is distributed across the variable states. Mixed
          stoichiometry signs make the constraint a net balance, so a
          target of zero pins exact stoichiometry (e.g.
          ``{"dO": ElementPool(target=0.0, members={"O_i": 1.0, "V_O": -1.0})}``)
          and a negative target an off-stoichiometry deficiency; a scan can
          cross zero continuously. An `ElementPoolError` is raised if a target
          is unreachable (beyond the achievable content) or inconsistent with
          the fixed concentrations. Defaults to None (no element
          constraints). See `py_sc_fermi.element_pools` for the solver.
        vbm_shift (float, optional): a temperature-dependent shift of the
          valence-band maximum, in eV, evaluated by the caller for this
          system's `temperature`. The DOS is scissored at construction so the
          band gap (and hence the carrier concentrations from the Fermi-Dirac
          integration) changes by `cbm_shift - vbm_shift`; with the VBM pinned
          at E=0, `vbm_shift` narrows the gap. Separately, when
          `rigid_shift=False`, `-charge * vbm_shift` is added to every
          variable-concentration charge state's formation energy. Also sets the
          effective band gap shown by `__repr__`. Defaults to 0.0.
        cbm_shift (float, optional): a temperature-dependent shift of the
          conduction-band minimum, in eV, evaluated by the caller for this
          system's `temperature`. Enters only through the DOS scissor: the band
          gap changes by `cbm_shift - vbm_shift`, moving the conduction block
          (and the carrier concentrations) and setting the effective band gap
          shown by `__repr__`. Defaults to 0.0 (no shift).
        formation_energy_corrections (dict[CorrectionKey, float] | None, optional):
          per-charge-state formation-energy corrections (in eV), evaluated by
          the caller for this system's `temperature`. Each key identifies one
          charge state, either as the `DefectChargeState` object itself or as
          a `(species_name, charge_state_name)` pair; the two forms may be
          mixed. Keying per charge state (rather than per formal charge)
          allows different corrections for metastable states that share a
          formal charge. Every key must identify a charge state in
          `defect_species` that has a formation energy (a fixed-concentration
          state has nothing for a correction to act on), and no charge state
          may be referenced twice. Defaults to None (no per-state
          corrections).
        rigid_shift (bool, optional): gates only the formation-energy channel;
          the DOS scissor from `vbm_shift`/`cbm_shift` is applied regardless. If
          True (the default), the band structure and defect levels are assumed
          to move together as a rigid body, so every variable-concentration
          charge state not covered by `formation_energy_corrections` is left
          unchanged. If False, the defect levels are fixed in absolute energy
          while the band edges move, so such charge states have their formation
          energy shifted by `-charge * vbm_shift`.
        fixed_concentrations (dict[FixedConcentrationKey, float] | None, optional):
          mapping of species name -> fixed total concentration (per unit
          cell), or ``(species_name, charge_state_name)`` -> fixed
          concentration for that single charge state. A named species has its
          total concentration fixed at the given value, overriding any
          species-level `fixed_concentration` it was constructed with; a
          `(species_name, charge_state_name)` pair fixes that one charge
          state, leaving the rest of its species variable. Both forms are
          applied by name to this system's own copies of `defect_species`, so
          they compose with `formation_energy_corrections` (resolved against
          the passed-in `defect_species`) and never mutate the caller's
          species. A charge state fixed by either key form has its formation
          energy (including any correction applied to it) left unused by the
          solve. A non-finite or negative value, or one that cannot be
          hosted (above its site-exclusion group's site budget -- its own
          `nsites`, or its shared `site_pools` size -- or below its
          individually-fixed charge states), is rejected at construction by the
          same checks as a species-level fix. Defaults to None (no fixes).
        occupancy_warning_threshold (float | None, optional): the site-occupancy
          fraction (0.01 = 1%) above which a ``DiluteLimitWarning`` is emitted
          when the system is solved, naming any species whose solved occupancy
          exceeds it. py-sc-fermi assumes dilute, non-interacting defects, so a
          high occupancy flags a regime in which un-modelled defect-defect
          interactions may make the results non-physical. The default is a
          heuristic, provisional value that will be revisited; lower it where
          interactions matter sooner (e.g. charged defects in low-dielectric
          hosts), or set ``None`` to silence the warning. Must be ``None`` or a
          finite fraction in (0, 1]. This is a reporting preference, not part of
          the physical system, and is not serialised. Defaults to 0.01.
        label (str | None, optional): an optional human-readable tag for this
          system (e.g. a growth condition or sample name), useful for
          identifying it within a batch of related systems. Included in
          ``as_dict`` only when set. Defaults to None.

    Raises:
        ValueError: if two entries in `defect_species` share a name, a pool
          references a species not in `defect_species`, a pool lists a
          species more than once, a species appears in more than one site
          pool, `formation_energy_corrections` references a `DefectChargeState`
          that is not part of `defect_species`, a `fixed_concentrations` key
          names a species (or, for a pair key, a charge state) not in
          `defect_species`, or `occupancy_warning_threshold` is not ``None``
          or a finite fraction in (0, 1].

    Note:
        `DefectSystem` is a fixed-temperature snapshot whose public
        attributes are read-only once constructed: rebinding any of them raises
        `AttributeError`. `vbm_shift`, `cbm_shift`,
        `formation_energy_corrections`, `rigid_shift`, and
        `fixed_concentrations` are applied once at construction to copies of
        `defect_species` -- the objects passed in (including any
        `formation_energy_corrections` keys) are never modified, even
        temporarily. To build `DefectSystem`s at several temperatures from
        temperature-dependent shift/correction functions, see
        `DefectSystemFactory`.
    """

    def __init__(
        self,
        defect_species: list[DefectSpecies],
        dos: DOS,
        volume: float,
        temperature: float,
        convergence_tolerance: float | None = None,
        site_pools: dict[str, SitePool] | None = None,
        element_pools: dict[str, ElementPool] | None = None,
        vbm_shift: float = 0.0,
        cbm_shift: float = 0.0,
        formation_energy_corrections: dict[CorrectionKey, float] | None = None,
        rigid_shift: bool = True,
        fixed_concentrations: dict[FixedConcentrationKey, float] | None = None,
        occupancy_warning_threshold: float | None = 0.01,
        label: str | None = None,
    ):
        self._volume = volume
        self._dos = dos
        delta_gap = cbm_shift - vbm_shift
        if delta_gap != 0.0:
            self._dos = self._dos.scissored(delta_gap)
        self._temperature = temperature
        self._label = label
        self._result: DefectSystemResult | None = None
        self._convergence_tolerance = convergence_tolerance
        self._vbm_shift = vbm_shift
        self._cbm_shift = cbm_shift
        self._rigid_shift = rigid_shift
        self._occupancy_warning_threshold = self._validate_occupancy_warning_threshold(
            occupancy_warning_threshold
        )
        self._occupancy_warning_emitted = False

        self._defect_species = tuple(copy.deepcopy(defect_species))
        self._site_pools: dict[str, SitePool] = dict(site_pools or {})
        self._element_pools: dict[str, ElementPool] = dict(element_pools or {})
        self._validate_pools()

        self._apply_formation_energy_corrections(
            defect_species, formation_energy_corrections or {}
        )
        self._apply_fixed_concentrations(fixed_concentrations or {})
        self._validate_fixed_concentrations()

    @property
    def volume(self) -> float:
        return self._volume

    @property
    def dos(self) -> DOS:
        return self._dos

    @property
    def temperature(self) -> float:
        return self._temperature

    @property
    def label(self) -> str | None:
        """An optional human-readable tag for this system, copied into
        ``result`` and serialised when set."""
        return self._label

    @property
    def convergence_tolerance(self) -> float | None:
        return self._convergence_tolerance

    @property
    def vbm_shift(self) -> float:
        return self._vbm_shift

    @property
    def cbm_shift(self) -> float:
        return self._cbm_shift

    @property
    def rigid_shift(self) -> bool:
        return self._rigid_shift

    @property
    def defect_species(self) -> tuple[DefectSpecies, ...]:
        return self._defect_species

    @property
    def site_pools(self) -> Mapping[str, SitePool]:
        return MappingProxyType(self._site_pools)

    @property
    def element_pools(self) -> Mapping[str, ElementPool]:
        return MappingProxyType(self._element_pools)

    @property
    def occupancy_warning_threshold(self) -> float | None:
        return self._occupancy_warning_threshold

    @staticmethod
    def _validate_occupancy_warning_threshold(value: float | None) -> float | None:
        """Validate the site-occupancy warning threshold.

        Returns ``value`` unchanged if it is ``None`` (warning disabled) or a
        finite fraction in ``(0, 1]``.

        Raises:
            ValueError: if ``value`` is not ``None`` and not a finite number in
                ``(0, 1]``. Occupancy is a fraction that maxes at 1.0, so a
                threshold outside this range -- including NaN or infinity --
                could never sensibly fire.
        """
        if value is None:
            return None
        if not (0.0 < value <= 1.0):
            raise ValueError(
                "occupancy_warning_threshold must be a fraction in (0, 1] or "
                f"None; got {value}"
            )
        return value

    def _apply_fixed_concentrations(
        self, fixed_concentrations: dict[FixedConcentrationKey, float]
    ) -> None:
        """Fix concentrations on the deep copies, by name.

        A ``str`` key names a species and sets its species-level
        ``fixed_concentration``, overriding any it was constructed with. A
        ``(species_name, charge_state_name)`` pair fixes that single charge
        state. Both resolve against this system's own copies in
        ``self.defect_species`` (never the caller's originals). The
        concentration values are checked (finite and non-negative) by
        ``DefectChargeState._fix_concentration``/``DefectSpecies._fix_concentration``
        and, for species totals and their site-exclusion budgets, by
        ``_validate_fixed_concentrations``.

        Raises:
            ValueError: if a key is neither form, a name is not found (via
                the name lookups, listing the available names), or a value
                is not finite and non-negative.
        """
        for key, conc in fixed_concentrations.items():
            if isinstance(key, tuple) and len(key) == 2:
                species_name, cs_name = key
                species = self.defect_species_by_name(species_name)
                species.charge_state_by_name(cs_name)._fix_concentration(conc)
            elif isinstance(key, str):
                self.defect_species_by_name(key)._fix_concentration(conc)
            else:
                raise ValueError(
                    "fixed_concentrations keys must be a species name or "
                    f"(species_name, charge_state_name) pair; got {key!r}."
                )

    @staticmethod
    def _resolve_correction_keys(
        defect_species: list[DefectSpecies],
        corrections: dict[CorrectionKey, float],
    ) -> dict[DefectChargeState, float]:
        """Canonicalise correction keys to ``DefectChargeState`` objects.

        ``(species_name, charge_state_name)`` keys are resolved against
        ``defect_species``; ``DefectChargeState`` keys pass through unchanged.

        Raises:
            ValueError: if a key is neither form, a pair names an unknown
                species or charge state, the charge state has no formation
                energy for a correction to act on (fixed concentration), or
                two keys resolve to the same charge state.
        """
        species_by_name = {ds.name: ds for ds in defect_species}
        resolved: dict[DefectChargeState, float] = {}
        for key, value in corrections.items():
            if isinstance(key, tuple) and len(key) == 2:
                species_name, cs_name = key
                if species_name not in species_by_name:
                    available = ", ".join(species_by_name)
                    raise ValueError(
                        f"no defect species named '{species_name}'; "
                        f"available: {available}"
                    )
                cs = species_by_name[species_name].charge_state_by_name(cs_name)
            elif isinstance(key, DefectChargeState):
                cs = key
            else:
                raise ValueError(
                    "formation_energy_corrections keys must be "
                    "DefectChargeState objects or "
                    f"(species_name, charge_state_name) pairs; got {key!r}."
                )
            if cs.energy is None:
                owner = next(
                    (ds.name for ds in defect_species if cs in ds.charge_states),
                    "<unknown>",
                )
                raise ValueError(
                    f"formation_energy_corrections references charge state "
                    f"'{cs.name}' of species '{owner}', which has no formation "
                    "energy (fixed concentration) and cannot be corrected."
                )
            if cs in resolved:
                owner = next(
                    (ds.name for ds in defect_species if cs in ds.charge_states),
                    "<unknown>",
                )
                raise ValueError(
                    "formation_energy_corrections contains more than one key for "
                    f"charge state '{cs.name}' of species '{owner}'."
                )
            resolved[cs] = value
        return resolved

    def _apply_formation_energy_corrections(
        self,
        original_defect_species: list[DefectSpecies],
        formation_energy_corrections: dict[CorrectionKey, float],
    ) -> None:
        """Permanently shift each variable-concentration charge state's
        formation energy in `self.defect_species` (a copy of
        `original_defect_species`) by its entry in
        `formation_energy_corrections` if present, else
        `-cs.charge * self.vbm_shift` if `self.rigid_shift` is False, else 0.
        Keys are `DefectChargeState` objects or
        `(species_name, charge_state_name)` pairs, canonicalised to objects
        before use.
        """
        resolved_corrections = self._resolve_correction_keys(
            original_defect_species, formation_energy_corrections
        )
        original_states = [cs for ds in original_defect_species for cs in ds.charge_states]
        copied_states = [cs for ds in self.defect_species for cs in ds.charge_states]

        unrecognised = set(resolved_corrections) - set(original_states)
        if unrecognised:
            raise ValueError(
                f"formation_energy_corrections references "
                f"{len(unrecognised)} DefectChargeState object(s) that are "
                "not part of `defect_species`."
            )

        for original_cs, copied_cs in zip(original_states, copied_states, strict=True):
            if copied_cs._energy is None:
                continue
            if original_cs in resolved_corrections:
                delta = resolved_corrections[original_cs]
            elif not self.rigid_shift:
                delta = -copied_cs.charge * self.vbm_shift
            else:
                delta = 0.0
            copied_cs._energy += delta

    def _validate_pools(self) -> None:
        """Validate the species roster and pool references.

        Raises:
            ValueError: if two roster species share a name, a pool references
                a species not in the roster, a site pool lists a species more
                than once, or a species appears in more than one site pool.
        """
        name_counts = Counter(ds.name for ds in self.defect_species)
        duplicates = sorted(name for name, count in name_counts.items() if count > 1)
        if duplicates:
            raise ValueError(
                f"defect_species contains duplicate names: {', '.join(duplicates)}. "
                "Species names must be unique."
            )

        roster = set(name_counts)
        site_pool_of: dict[str, str] = {}
        for pool_name, pool in self.site_pools.items():
            self._check_unknown_species("site pool", pool_name, pool.species, roster)
            self._check_duplicate_species("site pool", pool_name, pool.species)
            for name in pool.species:
                first = site_pool_of.setdefault(name, pool_name)
                if first != pool_name:
                    raise ValueError(
                        f"species '{name}' appears in site pools '{first}' and "
                        f"'{pool_name}'; a species may belong to at most one "
                        "site pool"
                    )
        for element, element_pool in self.element_pools.items():
            self._check_unknown_species(
                "element pool", element, element_pool.members, roster
            )

    @staticmethod
    def _check_unknown_species(
        kind: str, pool_name: str, members: Iterable[str], roster: set[str]
    ) -> None:
        """Raise if any member name is not in ``roster``."""
        unknown = sorted(set(members) - roster)
        if unknown:
            raise ValueError(
                f"{kind} '{pool_name}' references species not in "
                f"defect_species: {', '.join(unknown)}"
            )

    @staticmethod
    def _check_duplicate_species(
        kind: str, pool_name: str, members: Iterable[str]
    ) -> None:
        """Raise if any member name appears more than once (site pools only)."""
        repeated = sorted(
            name for name, count in Counter(members).items() if count > 1
        )
        if repeated:
            raise ValueError(
                f"{kind} '{pool_name}' lists species more than once: "
                f"{', '.join(repeated)}"
            )

    def _validate_fixed_concentrations(self) -> None:
        """Reject fixed concentrations that are inconsistent or cannot be hosted.

        Per-object values are validated at their own boundary, by
        ``DefectSpecies``/``DefectChargeState`` (at construction, and by their
        internal ``_fix_concentration`` setters). The checks here are
        cross-cutting and independent of the Fermi level -- they depend only
        on the fixed concentrations and the site budgets -- so they are static
        properties of the system, caught here at construction rather than left
        to surface (wrapped as a Fermi-window error) from a later solve:

        * a species fixed at the total level must be consistent with its
          individually-fixed charge states: those cannot exceed the total, and
          if every charge state is fixed (none variable) they must sum to it,
          since there is then nothing to make up any shortfall; and
        * within a site-exclusion group, the members' total fixed
          concentration cannot exceed the group's site budget.

        Raises:
            ValueError: if a species' fixed charge states are inconsistent
                with its species-level fixed concentration, or a group's
                total fixed concentration exceeds its site budget (its
                ``site_pools`` size, or, for an unpooled species, its own
                ``nsites``).
        """
        for sp in self.defect_species:
            if sp.fixed_concentration is None:
                continue
            charge_state_total = self._charge_state_fixed_total(sp)
            if math.isclose(charge_state_total, sp.fixed_concentration):
                continue
            if charge_state_total > sp.fixed_concentration:
                raise ValueError(
                    f"'{sp.name}' is fixed at {sp.fixed_concentration} but its "
                    f"fixed charge states require {charge_state_total}"
                )
            if all(cs.fixed_concentration is not None for cs in sp.charge_states):
                raise ValueError(
                    f"'{sp.name}' is fixed at {sp.fixed_concentration} but its "
                    f"fixed charge states sum to {charge_state_total}, with no "
                    f"variable charge state to make up the difference"
                )

        for label, n_sites, species in self._exclusion_group_specs():
            fixed_total = sum((self._fixed_total(sp) for sp in species), 0.0)
            if fixed_total > n_sites:
                raise ValueError(
                    f"'{label}' has {n_sites} sites but fixed-concentration "
                    f"states occupy {fixed_total}"
                )

    @staticmethod
    def _charge_state_fixed_total(species: DefectSpecies) -> float:
        """Sum of the explicit fixed concentrations of ``species``' charge
        states (zero if none are individually fixed).
        """
        total = 0.0
        for cs in species.charge_states:
            if cs.fixed_concentration is not None:
                total += cs.fixed_concentration
        return total

    @staticmethod
    def _fixed_total(species: DefectSpecies) -> float:
        """Total fixed concentration ``species`` contributes to its group's
        occupancy: its species-level ``fixed_concentration`` if set, otherwise
        the sum of its charge states' fixed concentrations.
        """
        if species.fixed_concentration is not None:
            return species.fixed_concentration
        return DefectSystem._charge_state_fixed_total(species)

    def __repr__(self) -> str:
        # The scissor already baked (cbm_shift - vbm_shift) into self.dos at
        # construction; reading dos.bandgap directly avoids double-counting it.
        bandgap = self.dos.bandgap
        bandgap_str = f"{bandgap:.3g} eV"

        lines = [
            "DefectSystem",
            f"  bandgap:     {bandgap_str}    nelect: {self.dos.nelect}",
            f"  volume:      {self.volume:.4g} \u00c5\u00b3    temperature: {self.temperature} K",
            f"\n  {len(self.defect_species)} defect species:",
        ]
        for ds in self.defect_species:
            header = f"\n  {ds.name}  [nsites: {ds.nsites}]"
            if ds.fixed_concentration is not None:
                header += f"  [fixed conc: {ds.fixed_concentration:.3e}]"
            lines.append(header)
            for cs in ds.charge_states:
                if cs.fixed_concentration is not None:
                    lines.append(
                        f"    q = {cs.charge:+2d}   [c] = {cs.fixed_concentration:.3e}"
                        f"  (deg. {cs.degeneracy})"
                    )
                else:
                    # a variable (non-fixed) charge state always carries an energy
                    energy = cs.energy
                    assert energy is not None
                    lines.append(
                        f"    q = {cs.charge:+2d}   E = {energy:8.4f} eV  (deg. {cs.degeneracy})"
                    )
        if self.site_pools:
            lines.append("\n  site pools:")
            for pool_name, pool in self.site_pools.items():
                lines.append(
                    f"    {pool_name}: {pool.n_sites:.4g} sites  \u2192  "
                    f"[{', '.join(pool.species)}]"
                )
        if self.element_pools:
            lines.append("\n  element pools:")
            for elem, element_pool in self.element_pools.items():
                sp_names = ", ".join(
                    f"{name} \u00d7{stoich:g}"
                    for name, stoich in element_pool.members.items()
                )
                lines.append(
                    f"    {elem}: {element_pool.target:.4g} per cell  \u2192  [{sp_names}]"
                )
        return "\n".join(lines)

    @property
    def defect_species_names(self) -> list[str]:
        """list of the names of all ``DefectSpecies`` considered in the
        ``DefectSystem``.

        Returns:
            list[str]: list of names of ``DefectSpecies`` objects
        """
        return [ds.name for ds in self.defect_species]

    def _resolve_species(self, name: str) -> DefectSpecies:
        """Resolve a species name to its entry in the roster."""
        return self.defect_species_by_name(name)

    def _build_group(
        self,
        n_sites: float,
        species: list[DefectSpecies],
        e_fermi: float,
        fixed_concs: dict[DefectChargeState, float],
    ) -> _ExclusionGroup:
        """Build an exclusion group of `n_sites` sites shared by `species`.

        Species-level `fixed_concentration` and individually-fixed charge
        states are written into `fixed_concs`; every other charge state
        becomes a variable state of the group. The free-site budget is
        `n_sites` minus the group's total fixed concentration, which
        `_validate_fixed_concentrations` guarantees is non-negative at
        construction.
        """
        variable_states: list[_VariableState] = []
        for sp in species:
            if sp.fixed_concentration is not None:
                for cs, conc in sp.charge_state_concentrations(e_fermi, self.temperature):
                    fixed_concs[cs] = conc
                continue
            for cs in sp.charge_states:
                if cs.fixed_concentration is not None:
                    fixed_concs[cs] = cs.fixed_concentration
                else:
                    Ef = cs.get_formation_energy(e_fermi)
                    log_w = np.log(cs.degeneracy) - Ef / (_kboltz * self.temperature)
                    variable_states.append(_VariableState(cs, sp, log_w))

        n_free = n_sites - sum((self._fixed_total(sp) for sp in species), 0.0)
        return _ExclusionGroup(n_free=n_free, variable_states=variable_states)

    def _exclusion_group_specs(self) -> list[tuple[str, float, list[DefectSpecies]]]:
        """Enumerate the site-exclusion groups as
        ``(label, site budget, member species)``: one group per `site_pools`
        entry (its pool size and members), plus an implicit single-species
        group of `nsites` for every species not in a site pool.
        """
        specs: list[tuple[str, float, list[DefectSpecies]]] = []
        pooled_species: set[DefectSpecies] = set()
        for pool_name, pool in self.site_pools.items():
            sp_objs = [self._resolve_species(name) for name in pool.species]
            pooled_species.update(sp_objs)
            specs.append((pool_name, pool.n_sites, sp_objs))
        for sp in self.defect_species:
            if sp not in pooled_species:
                specs.append((sp.name, sp.nsites, [sp]))
        return specs

    def _build_exclusion_groups(
        self, e_fermi: float
    ) -> tuple[list[_ExclusionGroup], dict[DefectChargeState, float]]:
        """Partition all charge states into site-exclusion groups.

        Each `site_pools` entry forms a group sharing its pool size between
        its member species. Every other species forms an implicit
        single-species group of size `nsites` -- so an unpooled species and
        a `site_pools` entry of `(nsites, [species])` are equivalent.

        Returns the groups (covering every variable charge state) together
        with a dict of concentrations for every fixed and species-fixed
        charge state.
        """
        groups: list[_ExclusionGroup] = []
        fixed_concs: dict[DefectChargeState, float] = {}
        for _, n_sites, species in self._exclusion_group_specs():
            groups.append(self._build_group(n_sites, species, e_fermi, fixed_concs))
        return groups, fixed_concs

    def _resolve_element_pools(self) -> dict[str, ResolvedElementPool]:
        """Resolve string species references in `element_pools` to DefectSpecies."""
        return {
            elem: ResolvedElementPool(
                target=pool.target,
                members=MappingProxyType(
                    {
                        self._resolve_species(name): stoich
                        for name, stoich in pool.members.items()
                    }
                ),
            )
            for elem, pool in self.element_pools.items()
        }

    @staticmethod
    def _stoichiometry_lookup(
        pools: dict[str, ResolvedElementPool],
    ) -> dict[DefectSpecies, dict[str, float]]:
        """Map each species to {element: stoichiometry} for every element
        pool it participates in."""
        lookup: dict[DefectSpecies, dict[str, float]] = {}
        for elem, pool in pools.items():
            for sp, stoich in pool.members.items():
                lookup.setdefault(sp, {})[elem] = stoich
        return lookup

    @staticmethod
    def _remaining_element_targets(
        pools: dict[str, ResolvedElementPool],
    ) -> dict[str, float]:
        """For each element pool, the target content still to be supplied by
        variable-concentration states, after subtracting fixed contributions."""
        remaining: dict[str, float] = {}
        for elem, pool in pools.items():
            committed = 0.0
            for sp, stoich in pool.members.items():
                if sp.fixed_concentration is not None:
                    committed += stoich * sp.fixed_concentration
                else:
                    committed += stoich * sum(
                        cs.fixed_concentration
                        for cs in sp.charge_states
                        if cs.fixed_concentration is not None
                    )
            rem = pool.target - committed
            # Float summation of fixed contributions does not land exactly
            # on a target they are meant to meet (0.1 + 0.2 != 0.3): treat
            # a remainder within rounding distance of the committed total
            # as exactly zero. The committed total may be negative (fixed
            # negative-stoichiometry states), so scale by its magnitude.
            # Genuinely non-zero remainders are validated by the caller,
            # which knows whether the pool can shed content.
            if abs(rem) <= 1e-12 * abs(committed):
                rem = 0.0
            remaining[elem] = rem
        return remaining

    @staticmethod
    def _group_term_arrays(
        variable_states: list[_VariableState],
        elements: list[str],
        stoich: dict[DefectSpecies, dict[str, float]],
    ) -> tuple[np.ndarray, np.ndarray]:
        """Per-state log Boltzmann weights and stoichiometry vectors for a
        group's variable states, as arrays of shape (n,) and (n, len(elements))."""
        log_w = np.array([state.log_w for state in variable_states])
        s = np.array(
            [
                [stoich.get(state.species, {}).get(elem, 0.0) for elem in elements]
                for state in variable_states
            ]
        )
        return log_w, s

    def _solve_chemical_potentials(
        self,
        groups: list[_ExclusionGroup],
        elements: list[str],
        stoich: dict[DefectSpecies, dict[str, float]],
        remaining: dict[str, float],
        pools: dict[str, ResolvedElementPool],
    ) -> np.ndarray:
        """Solve for one chemical potential `mu_X` per element pool such
        that the total content of each pooled element, summed over every
        group, equals `remaining[X]`. Assembles the per-group state arrays
        and delegates the numerical solve to
        `element_pools.solve_chemical_potentials`."""
        group_data = [
            (group.n_free, *self._group_term_arrays(group.variable_states, elements, stoich))
            for group in groups
            if group.n_free > 0 and group.variable_states
        ]
        return element_pools.solve_chemical_potentials(
            group_data,
            elements,
            np.array([remaining[e] for e in elements]),
            [pools[e].target for e in elements],
        )

    def _solve_element_pools(
        self,
        groups: list[_ExclusionGroup],
        pools: dict[str, ResolvedElementPool],
        elements: list[str],
        stoich: dict[DefectSpecies, dict[str, float]],
    ) -> tuple[np.ndarray, list[str], set[DefectSpecies]]:
        """Dispose of the element pools and solve for their chemical
        potentials.

        Args:
            groups: the site-exclusion groups for this Fermi level.
            pools: resolved element pools (`_resolve_element_pools`).
            elements: the constrained element names, ``list(pools)``.
            stoich: per-species ``{element: stoichiometry}`` lookup.

        Returns:
            ``(mu, elements, forced_zero)``. `mu` holds the chemical
            potential of each element in the returned `elements`, a copy
            of the input narrowed to those that remain in the solve (empty
            when there are no pools or all are exhausted). `forced_zero` is
            the set of species whose variable states the caller must
            assign concentration zero.

        A pool with a negative-stoichiometry variable state can shed
        content, so zero and negative remaining targets are balance
        conditions at finite ``lambda_X = exp(mu_X)`` and the element stays
        in the solve. Without one, variable states can only add content: a
        negative remainder is inconsistent with the fixed concentrations
        (raising `ElementPoolError`), and a zero remainder is the
        ``lambda_X -> 0`` limit, in which every variable state with
        positive stoichiometry in X has concentration 0 (collected into
        `forced_zero`) and X drops out of the chemical-potential solve.
        """
        forced_zero: set[DefectSpecies] = set()
        if not elements:
            return np.array([]), elements, forced_zero

        remaining = self._remaining_element_targets(pools)
        variable_species = {
            state.species for group in groups for state in group.variable_states
        }
        balance_capable = {
            e
            for e in elements
            if any(
                s_by_elem.get(e, 0.0) < 0.0 and sp in variable_species
                for sp, s_by_elem in stoich.items()
            )
        }
        for e in elements:
            if remaining[e] < 0.0 and e not in balance_capable:
                target = pools[e].target
                raise ElementPoolError(
                    f"Element pool '{e}': fixed-concentration states "
                    f"contribute {target - remaining[e]:.3e}, exceeding "
                    f"the target {target:.3e} by {-remaining[e]:.3e}, and "
                    "no variable state can remove content. Your "
                    "constraints are mutually inconsistent."
                )
        exhausted = {
            e for e in elements if remaining[e] == 0.0 and e not in balance_capable
        }
        solve_groups = groups
        if exhausted:
            forced_zero = {
                sp
                for sp, s_by_elem in stoich.items()
                if any(s_by_elem.get(e, 0.0) > 0.0 for e in exhausted)
            }
            elements = [e for e in elements if e not in exhausted]
            solve_groups = [
                _ExclusionGroup(
                    group.n_free,
                    [
                        state
                        for state in group.variable_states
                        if state.species not in forced_zero
                    ],
                )
                for group in groups
            ]
        if not elements:
            return np.array([]), elements, forced_zero
        mu = self._solve_chemical_potentials(
            solve_groups, elements, stoich, remaining, pools
        )
        return mu, elements, forced_zero

    def _global_defect_concs(self, e_fermi: float) -> dict[DefectChargeState, float]:
        """
        Returns a dict mapping each DefectChargeState -> concentration per cell.

        Every species is assigned to a site-exclusion group (a `site_pools`
        entry, or its own implicit single-species group of `nsites` sites)
        and its variable charge states are given Langmuir/site-exclusion
        statistics, ``c_i = n_free * w_i / (1 + sum_j w_j)``. If
        `element_pools` are defined, the element chemical potentials are
        solved for first (`_solve_element_pools`) and folded into each
        state's weight as ``w_i * exp(s_i . mu)``.
        """
        groups, all_concs = self._build_exclusion_groups(e_fermi)

        pools = self._resolve_element_pools()
        elements = list(pools.keys())
        stoich = self._stoichiometry_lookup(pools)

        mu, elements, forced_zero = self._solve_element_pools(
            groups, pools, elements, stoich
        )

        for group in groups:
            if not group.variable_states:
                continue
            if group.n_free == 0:
                for cs, _, _ in group.variable_states:
                    all_concs[cs] = 0.0
                continue
            kept = []
            for state in group.variable_states:
                if state.species in forced_zero:
                    all_concs[state.charge_state] = 0.0
                else:
                    kept.append(state)
            if not kept:
                continue
            log_w, s = self._group_term_arrays(kept, elements, stoich)
            for state, c_i in zip(
                kept,
                element_pools.group_concs(group.n_free, log_w, s, mu),
                strict=True,
            ):
                all_concs[state.charge_state] = c_i

        return all_concs

    @classmethod
    def from_dict(cls, dictionary: dict) -> DefectSystem:
        """Generate a DefectSystem from a dictionary.

        Reads ``site_pools`` and ``element_pools`` when present; dictionaries
        without those keys (e.g. older serialisations) load unchanged.

        Args:
            dictionary (dict): Dictionary containing the DefectSystem data, as
              produced by ``DefectSystem.as_dict``.

        Returns:
            DefectSystem: DefectSystem corresponding to the provided dictionary.
        """
        site_pools = {
            name: SitePool.from_dict(pool)
            for name, pool in dictionary.get("site_pools", {}).items()
        }
        element_pools = {
            element: ElementPool.from_dict(pool)
            for element, pool in dictionary.get("element_pools", {}).items()
        }
        return cls(
            dos=DOS.from_dict(dictionary["dos"]),
            volume=dictionary["volume"],
            temperature=dictionary["temperature"],
            convergence_tolerance=dictionary.get("convergence_tolerance"),
            defect_species=[
                DefectSpecies.from_dict(defect_species)
                for defect_species in dictionary["defect_species"]
            ],
            site_pools=site_pools,
            element_pools=element_pools,
            label=dictionary.get("label"),
        )

    def defect_species_by_name(self, name: str) -> DefectSpecies:
        """return a ``DefectSpecies`` contained within the ``DefectSystem``
        via its name.

        Args:
            name (str): name of the ``DefectSpecies`` to return

        Returns:
            DefectSpecies: ``DefectSpecies`` where ``DefectSpecies.name == name``

        Raises:
            ValueError: if no ``DefectSpecies`` in the system has that name.
        """
        for ds in self.defect_species:
            if ds.name == name:
                return ds
        available = ", ".join(self.defect_species_names)
        raise ValueError(f"no defect species named '{name}'; available: {available}")
    
    @property
    def result(self) -> DefectSystemResult:
        """The solved state of this system as a ``DefectSystemResult``.

        Solves for the self-consistent Fermi energy on first access and caches
        the result; subsequent accesses return the same object without
        re-solving. Safe to cache unconditionally because the system is an
        immutable snapshot.

        Returns:
            DefectSystemResult: the Fermi level, carrier and defect
            concentrations at the self-consistent Fermi energy.
        """
        # Deliberately a manual cache, not functools.cached_property: the latter
        # inserts a functools frame that DiluteLimitWarning's attribution walk
        # (_warning_stacklevel) stops at, misattributing the warning to
        # functools.py instead of the caller. See
        # test_result_warning_is_attributed_to_caller.
        if self._result is None:
            e_fermi, _ = self.get_sc_fermi()
            p0, n0 = self.dos.carrier_concentrations(e_fermi, self.temperature)
            per_cell = self._per_charge_state_concs(e_fermi)
            self._result = DefectSystemResult(
                temperature=float(self.temperature),
                fermi_energy=float(e_fermi),
                volume=float(self.volume),
                label=self.label,
                p0_per_cell=float(p0),
                n0_per_cell=float(n0),
                charge_state_concentrations_per_cell=MappingProxyType(
                    {
                        species: MappingProxyType(states)
                        for species, states in per_cell.items()
                    }
                ),
            )
        return self._result

    def get_sc_fermi(self) -> tuple[float, float]:
        """Calculate the self-consistent Fermi energy.
        
        Finds the Fermi energy at which charge neutrality is achieved,
        using Brent's method for root finding.
        
        Returns:
            tuple[float, float]: The self-consistent Fermi energy and the
            absolute residual charge density at that energy.
        
        Raises:
            RuntimeError: If no solution is found within the DOS energy range.
            ElementPoolError: If element-pool constraints cannot be satisfied
              at a probed Fermi level.
        """
        emin = self.dos.emin()
        emax = self.dos.emax()

        try:
            kwargs = {}
            if self.convergence_tolerance is not None:
                kwargs['xtol'] = self.convergence_tolerance

            e_fermi = brentq(
                self.q_tot,
                emin,
                emax,
                **kwargs,
            )  # type: ignore[call-overload]
        except ElementPoolError:
            # An element-pool failure diagnoses the constraints, not the
            # energy window; surface it undisguised.
            raise
        except ValueError as err:
            raise RuntimeError(
                f"No solution found between {emin} and {emax}: {err}"
            ) from err
        
        # q_tot is continuous wherever finite, so a sign change across the DOS
        # window brackets a genuine root. The returned residual is a diagnostic
        # on that solution.
        residual = abs(self.q_tot(e_fermi))
        self._warn_if_high_occupancy(e_fermi)
        return e_fermi, residual

    def _warn_if_high_occupancy(self, e_fermi: float) -> None:
        """Emit a ``DiluteLimitWarning`` if any species' solved site occupancy
        exceeds ``occupancy_warning_threshold``.

        Does nothing when the threshold is ``None`` or when this system has
        already warned. Otherwise computes the per-species occupancy fractions
        at the converged ``e_fermi`` and, if any species is above the threshold,
        emits one warning naming every offending species and its occupancy,
        worst first. The warning fires at most once per ``DefectSystem``: the
        occupancy verdict is deterministic for an immutable snapshot, so a user
        inspecting one solved system through several routes (a direct
        ``get_sc_fermi``, or the deeper ``result``, ``site_percentages`` or
        ``element_chemical_potential_shifts``) sees a single warning, attributed
        to their own call rather than to an internal solve method.

        Args:
            e_fermi (float): the self-consistent Fermi energy from
              ``get_sc_fermi``.
        """
        threshold = self.occupancy_warning_threshold
        if self._occupancy_warning_emitted or threshold is None:
            return
        fractions = self._site_occupancy_fractions(e_fermi)
        offenders = sorted(
            (
                (name, fraction)
                for name, fraction in fractions.items()
                if fraction > threshold
            ),
            key=lambda item: item[1],
            reverse=True,
        )
        if not offenders:
            return
        warnings.warn(
            self._high_occupancy_message(offenders),
            DiluteLimitWarning,
            stacklevel=self._warning_stacklevel(),
        )
        self._occupancy_warning_emitted = True

    @staticmethod
    def _warning_stacklevel() -> int:
        """``stacklevel`` for the high-occupancy ``warnings.warn`` call that
        points at the first frame outside this module.

        The warning is reached at different stack depths -- directly via
        ``get_sc_fermi``, or more deeply through the ``result`` property
        (accessed as ``system.result``, or via ``site_percentages`` or
        ``element_chemical_potential_shifts``, which read
        ``result.fermi_energy``) -- so no fixed ``stacklevel`` blames the
        caller on every path. Counting frames up to the first one outside
        ``py_sc_fermi.defect_system`` attributes the warning to the user's
        call. ``warnings.warn`` is invoked one frame above this helper, where
        the count begins.
        """
        frame: Any = sys._getframe(1)
        level = 1
        while frame is not None and frame.f_globals.get("__name__") == __name__:
            frame = frame.f_back
            level += 1
        return level

    @staticmethod
    def _high_occupancy_message(offenders: list[tuple[str, float]]) -> str:
        """Build the ``DiluteLimitWarning`` message for the offending species.

        Args:
            offenders: ``(name, occupancy fraction)`` pairs, ordered worst first.

        Returns:
            str: a two-sentence message reporting the occupancies, naming the
            dilute assumption, and warning the results may be non-physical.
        """
        caveat = (
            "py-sc-fermi assumes dilute, non-interacting defects, so results at "
            "this occupancy may be non-physical."
        )
        if len(offenders) == 1:
            name, fraction = offenders[0]
            return (
                f"'{name}' reaches {fraction * 100:.3g}% site occupancy at the "
                f"self-consistent Fermi level. {caveat}"
            )
        listed = [f"'{name}' ({fraction * 100:.3g}%)" for name, fraction in offenders]
        joined = f"{', '.join(listed[:-1])} and {listed[-1]}"
        return (
            f"{joined} reach high site occupancy at the self-consistent Fermi "
            f"level. {caveat}"
        )

    def total_defect_charge_contributions(self, e_fermi: float) -> tuple[float,float]:
        lhs = rhs = 0.0
        for cs, conc in self._global_defect_concs(e_fermi).items():
            if cs.charge > 0:
                lhs += conc * cs.charge
            elif cs.charge < 0:
                rhs += conc * abs(cs.charge)
        return lhs, rhs
    
    def q_tot(self, e_fermi: float) -> float:
        """for a given Fermi energy, calculate the net charge density of the
        ``DefectSystem`` as the difference between charge contributions from all
        positive species (including holes) and all negative species (including
        electrons).

        Args:
            e_fermi (float): Fermi energy

        Returns:
            float: net charge density of the ``DefectSystem`` at ``e_fermi``
        """
        p0, n0 = self.dos.carrier_concentrations(e_fermi, self.temperature)
        lhs_def, rhs_def = self.total_defect_charge_contributions(e_fermi)
        lhs = p0 + lhs_def
        rhs = n0 + rhs_def
        diff = rhs - lhs
        return diff

    def get_transition_levels(self) -> dict[str, list[list]]:
        """Return transition_levels transition levels profiles of all ``DefectSpecies``
        all defects as dictionary of ``{DefectSpecies.name : [e_fermi, e_formation]}``
        over the whole density of states energy range.

        Returns:
            dict[str, list[list]]: Dictionary giving per-defect transition-level
            profiles.
        """
        transition_levels = {}
        for defect_species in self.defect_species_names:
            transition_level = self.defect_species_by_name(defect_species).tl_profile(
                self.dos.emin(), self.dos.emax()
            )
            x = [[x_value][0][0] for x_value in transition_level]
            y = [[y_value][0][1] for y_value in transition_level]
            transition_levels.update({defect_species: [x, y]})
        return transition_levels

    def _per_charge_state_concs(self, e_fermi: float) -> dict[str, dict[str, float]]:
        """Per-species, per-charge-state concentrations per unit cell at
        ``e_fermi``, keyed by species name then charge-state name.
        """
        cs_concs = self._global_defect_concs(e_fermi)
        return {
            ds.name: {cs.name: float(cs_concs[cs]) for cs in ds.charge_states}
            for ds in self.defect_species
        }

    def _site_occupancy_fractions(self, e_fermi: float) -> dict[str, float]:
        """Return each ``DefectSpecies``' site occupancy as a fraction of the
        sites available to it, at an already-solved Fermi level.

        The concentrations are the same pool- and exclusion-aware values
        underlying ``result`` (``_global_defect_concs`` at ``e_fermi``). The
        denominator is the sites available to the species:
        the pool's total site count for a species in a ``site_pools`` entry,
        otherwise its own ``nsites``. This is the fractional form of
        ``site_percentages`` (which scales it by 100); both draw on a single
        solved ``e_fermi`` rather than re-solving, so the high-occupancy warning
        and ``site_percentages`` share one source of truth.

        Args:
            e_fermi (float): a self-consistent Fermi energy, as returned by
              ``get_sc_fermi``.

        Returns:
            dict[str, float]: mapping of ``DefectSpecies`` name to its site
            occupancy as a fraction (0.0-1.0).
        """
        concs = self._global_defect_concs(e_fermi)
        pool_sites = {
            name: pool.n_sites
            for pool in self.site_pools.values()
            for name in pool.species
        }
        return {
            str(ds.name): float(
                sum(concs[cs] for cs in ds.charge_states)
                / pool_sites.get(ds.name, ds.nsites)
            )
            for ds in self.defect_species
        }

    def site_percentages(self) -> dict[str, float]:
        """Return each ``DefectSpecies``' solved site occupancy as a
        percentage of the sites available to it.

        The occupancies are drawn from the same solved, pool- and
        exclusion-aware concentrations as ``result`` (``_global_defect_concs``
        at ``result.fermi_energy``). The denominator
        is the sites available to the species -- its own ``nsites``, or, for a
        species in a ``site_pools`` entry, the pool's total site count (so a
        pool's members together occupy at most 100% of it). A species'
        concentration cannot exceed its available sites, so every occupancy is
        at most 100%.

        Reads the cached ``result`` rather than re-solving.

        Returns:
            dict[str, float]: mapping of ``DefectSpecies`` name to its site
            occupancy as a percentage.
        """
        e_fermi = self.result.fermi_energy
        return {
            name: fraction * 100
            for name, fraction in self._site_occupancy_fractions(e_fermi).items()
        }

    def element_chemical_potential_shifts(self) -> dict[str, float]:
        """Solved chemical-potential shifts for each constrained element.

        Returns the chemical-potential shift ``delta_mu`` (in eV) of every
        element constrained by ``element_pools``, evaluated at
        ``result.fermi_energy``. ``delta_mu`` is measured relative to the
        chemical-potential reference at which the formation energies were
        defined: a target above the element's unconstrained content gives
        ``delta_mu > 0``, a target below it ``delta_mu < 0``, and a target
        equal to it ``delta_mu ~ 0``. These are shifts, not absolute chemical
        potentials.

        py-sc-fermi reports the shift but cannot judge whether it is
        physically accessible, having no competing-phase data. The value is
        intended for an external stability-region check: a shift inside the
        host's stability region corresponds to an accessible concentration,
        one outside to a supersaturated or otherwise unphysical state.

        The shift is re-derived from the same deterministic element-pool solve
        used for the concentrations at this Fermi level, so it is consistent
        with ``result``.

        Reads the cached ``result`` rather than re-solving.

        Returns:
            dict[str, float]: each constrained element mapped to its
            ``delta_mu`` in eV, in the order the pools were declared. Empty
            when no ``element_pools`` are defined. An element whose target is
            met with zero variable content and that has no content-removing
            (negative-stoichiometry) state is fully excluded; its shift is the
            ``delta_mu -> -inf`` limit, reported as ``-math.inf``.

        Raises:
            ElementPoolError: if the element-pool constraints cannot be
                satisfied -- the targets are mutually inconsistent or exceed
                the available sites.
        """
        if not self.element_pools:
            return {}
        e_fermi = self.result.fermi_energy
        groups, _ = self._build_exclusion_groups(e_fermi)
        pools = self._resolve_element_pools()
        elements = list(pools)
        stoich = self._stoichiometry_lookup(pools)
        mu, solved_elements, _ = self._solve_element_pools(
            groups, pools, elements, stoich
        )
        shifts = {
            elem: float(mu[i] * _kboltz * self.temperature)
            for i, elem in enumerate(solved_elements)
        }
        return {elem: shifts.get(elem, -math.inf) for elem in pools}

    def as_dict(self) -> dict:
        """Return a dictionary representation of the ``DefectSystem``.

        The returned dictionary contains only JSON/YAML-safe types and round-
        trips through ``DefectSystem.from_dict``. ``site_pools`` and
        ``element_pools`` are included only when non-empty; each pool is
        serialised as a mapping with named fields, species referenced by name
        (a site pool as ``{"n_sites", "species"}``; an element pool as
        ``{"target", "members": {species: stoichiometry}}``).

        The serialised DOS and formation energies already include any
        corrections applied at construction: ``vbm_shift``/``cbm_shift`` are
        baked into the DOS as a band-gap scissor, and (under
        ``rigid_shift=False``) ``vbm_shift`` plus any
        ``formation_energy_corrections`` are baked into the formation energies.
        Those constructor parameters are deliberately not round-tripped, since
        re-applying a baked correction on load would double-count it.

        Returns:
            dict: dictionary representation of the ``DefectSystem``.
        """

        defect_system_dict = dict(
            volume=float(self.volume),
            temperature=float(self.temperature),
            defect_species=[
                defect_species.as_dict() for defect_species in self.defect_species
            ],
            dos=self.dos.as_dict(),
        )
        if self.convergence_tolerance is not None:
            defect_system_dict["convergence_tolerance"] = float(self.convergence_tolerance)
        if self.label is not None:
            defect_system_dict["label"] = self.label
        if self.site_pools:
            defect_system_dict["site_pools"] = {
                name: pool.as_dict() for name, pool in self.site_pools.items()
            }
        if self.element_pools:
            defect_system_dict["element_pools"] = {
                element: pool.as_dict()
                for element, pool in self.element_pools.items()
            }
        return defect_system_dict


class DefectSystemFactory:
    """Builds `DefectSystem` snapshots at different temperatures from
    temperature-dependent band-edge shift and formation-energy correction
    functions.

    `DefectSystem` itself takes pre-evaluated `vbm_shift`/`cbm_shift`/
    `formation_energy_corrections` for a single temperature. This factory
    stores the temperature-dependent functions instead, and `at(temperature)`
    evaluates them and constructs a fresh, independent `DefectSystem` for
    that temperature.

    Args:
        defect_species (list[DefectSpecies]): List of ``DefectSpecies`` objects
          which are present in the ``DefectSystem``.
        dos (DOS): the ``DOS`` object associated with the unit cell.
        volume (float): volume of the unit cell in Angstroms cubed.
        convergence_tolerance (float, optional): Tolerance for the Fermi energy
          convergence in eV, passed to every `DefectSystem` built by `at()`.
        site_pools (dict[str, SitePool] | None, optional):
          a mapping of pool name to a `SitePool`, passed to every
          `DefectSystem` built by `at()`.
        element_pools (dict[str, ElementPool] | None, optional):
          a mapping of element name to an `ElementPool`, passed to every
          `DefectSystem` built by `at()`. Defaults to None.
        vbm_shift_fn (Callable[[float], float] | None, optional): a function
          of temperature returning the valence-band-maximum shift (in eV) at
          that temperature. Defaults to None (no shift).
        cbm_shift_fn (Callable[[float], float] | None, optional): a function
          of temperature returning the conduction-band-minimum shift (in eV)
          at that temperature. Defaults to None (no shift).
        formation_energy_correction_fns (dict | None, optional): per-charge-state
          temperature-dependent formation-energy corrections, as a mapping
          ``dict[CorrectionKey, Callable[[float], float]]``. Each callable
          takes a temperature in K and returns a correction in eV to add to
          that charge state's formation energy at each snapshot. Keys
          identify a charge state either as the `DefectChargeState` object or
          as a `(species_name, charge_state_name)` pair; they are resolved
          when `at()` builds each `DefectSystem`. Defaults to None.
        rigid_shift (bool, optional): passed to every `DefectSystem` built by
          `at()`. Defaults to True.
        occupancy_warning_threshold (float | None, optional): passed to every
          `DefectSystem` built by `at()`, so a temperature sweep warns
          consistently. See `DefectSystem` for its meaning and validation.
          Defaults to 0.01.
        label (str | None, optional): passed to every `DefectSystem` built by
          `at()`, unless overridden per call. Defaults to None.
    """

    def __init__(
        self,
        defect_species: list[DefectSpecies],
        dos: DOS,
        volume: float,
        convergence_tolerance: float | None = None,
        site_pools: dict[str, SitePool] | None = None,
        element_pools: dict[str, ElementPool] | None = None,
        vbm_shift_fn: Callable[[float], float] | None = None,
        cbm_shift_fn: Callable[[float], float] | None = None,
        formation_energy_correction_fns: (
            dict[CorrectionKey, Callable[[float], float]] | None
        ) = None,
        rigid_shift: bool = True,
        occupancy_warning_threshold: float | None = 0.01,
        label: str | None = None,
    ):
        self.defect_species = defect_species
        self.dos = dos
        self.volume = volume
        self.convergence_tolerance = convergence_tolerance
        self.site_pools = site_pools
        self.element_pools = element_pools
        self.vbm_shift_fn = vbm_shift_fn
        self.cbm_shift_fn = cbm_shift_fn
        self.formation_energy_correction_fns = formation_energy_correction_fns or {}
        self.rigid_shift = rigid_shift
        self.occupancy_warning_threshold = occupancy_warning_threshold
        self.label = label

    def at(self, temperature: float, **overrides: Any) -> DefectSystem:
        """Build a `DefectSystem` snapshot at `temperature`.

        Evaluates `vbm_shift_fn`, `cbm_shift_fn`, and
        `formation_energy_correction_fns` at `temperature`, then constructs a
        `DefectSystem` from those values together with this factory's other
        attributes. Each call deep-copies `defect_species` independently (via
        `DefectSystem.__init__`), so repeated calls (e.g. for an
        anneal-and-quench sequence) do not affect one another.

        Args:
            temperature (float): temperature (in K) at which to evaluate the
              shift/correction functions and solve the resulting
              `DefectSystem`.
            **overrides: keyword arguments forwarded to `DefectSystem.__init__`,
              overriding any of this factory's corresponding attributes (e.g.
              `temperature` itself, or `fixed_concentrations`, a
              `dict[FixedConcentrationKey, float]` mapping a species name (or a
              ``(species_name, charge_state_name)`` pair, fixing that single
              charge state) to a fixed concentration per unit cell). Passing
              `fixed_concentrations` freezes those species' totals (or
              per-charge-state values) at the given values, as in the
              anneal-and-quench workflow: take the totals (or per-charge-state
              values) solved at a high temperature and re-solve at a lower one
              with them held fixed.

        Returns:
            DefectSystem: a new, independent `DefectSystem` for `temperature`.
        """
        vbm_shift = self.vbm_shift_fn(temperature) if self.vbm_shift_fn else 0.0
        cbm_shift = self.cbm_shift_fn(temperature) if self.cbm_shift_fn else 0.0
        formation_energy_corrections = {
            cs: fn(temperature)
            for cs, fn in self.formation_energy_correction_fns.items()
        }
        kwargs: dict[str, Any] = dict(
            defect_species=self.defect_species,
            dos=self.dos,
            volume=self.volume,
            temperature=temperature,
            convergence_tolerance=self.convergence_tolerance,
            site_pools=self.site_pools,
            element_pools=self.element_pools,
            vbm_shift=vbm_shift,
            cbm_shift=cbm_shift,
            formation_energy_corrections=formation_energy_corrections,
            rigid_shift=self.rigid_shift,
            occupancy_warning_threshold=self.occupancy_warning_threshold,
            label=self.label,
        )
        kwargs.update(overrides)
        return DefectSystem(**kwargs)
