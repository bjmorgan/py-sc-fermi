from __future__ import annotations

import copy
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.constants import physical_constants as _physical_constants
from scipy.optimize import brentq, root
from scipy.special import logsumexp

from py_sc_fermi.defect_charge_state import DefectChargeState
from py_sc_fermi.defect_species import DefectSpecies
from py_sc_fermi.dos import DOS

_kboltz = _physical_constants["Boltzmann constant in eV/K"][0]

# Maximum acceptable relative deviation |content / target - 1| for a
# successful element-pool solve; anything looser is raised as an error
# rather than returned as silently unconstrained concentrations.
_element_pool_tolerance = 1e-6


@dataclass
class _ExclusionGroup:
    """A set of charge states competing for a shared budget of `n_free`
    sites under Langmuir/site-exclusion statistics.

    `variable_states` holds ``(charge_state, species, log_w)`` for every
    non-fixed-concentration charge state in the group, where
    ``log_w = log(degeneracy) - Ef / (kB * T)`` is the *per-site* Boltzmann
    weight (independent of `nsites`/pool size).
    """

    n_free: float
    variable_states: list[tuple[DefectChargeState, DefectSpecies, float]]


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
        site_pools (dict[str, tuple[float, list[DefectSpecies]]] | None, optional):
          Mapping of pool name -> (total sites in that pool, list of
          DefectSpecies sharing those sites). By default (None), every
          DefectSpecies gets its own implicit exclusion group of `nsites`
          sites, so its charge states already compete with each other via
          Langmuir statistics; `site_pools` is only needed when several
          species must share one physical site budget.
        vbm_shift (float, optional): a temperature-dependent shift of the
          valence-band maximum, in eV, evaluated by the caller for this
          system's `temperature`. Used (with `cbm_shift`) to compute the
          effective band gap shown by `__repr__`/`report`, and (when
          `rigid_shift=True`, the default) added to every variable-concentration
          charge state's formation energy as `-charge * vbm_shift`. Defaults
          to 0.0 (no shift).
        cbm_shift (float, optional): a temperature-dependent shift of the
          conduction-band minimum, in eV, evaluated by the caller for this
          system's `temperature`. Used with `vbm_shift` to compute the
          effective band gap shown by `__repr__`/`report`. Defaults to 0.0
          (no shift).
        formation_energy_corrections (dict[DefectChargeState, float] | None, optional):
          per-charge-state formation-energy corrections (in eV), evaluated by
          the caller for this system's `temperature`. Keying by the
          `DefectChargeState` object itself (rather than e.g.
          `(species_name, charge)`) allows different corrections for
          metastable `DefectChargeState`s that share a formal charge. Every
          key must be one of the `DefectChargeState`s in `defect_species`.
          Defaults to None (no per-state corrections).
        rigid_shift (bool, optional): if True (the default), the band
          structure and defect levels are assumed to move together as a
          rigid body, so `vbm_shift`/`cbm_shift` only affect the displayed
          band gap and every variable-concentration charge state not covered
          by `formation_energy_corrections` is left unchanged. If False, the
          defect levels are fixed in absolute energy while the band edges
          move, so such charge states have their formation energy shifted by
          `-charge * vbm_shift`.

    Note:
        `DefectSystem` is an immutable, fixed-temperature snapshot:
        `vbm_shift`, `cbm_shift`, `formation_energy_corrections`, and
        `rigid_shift` are evaluated once at construction and applied to
        copies of `defect_species` -- the objects passed in (including any
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
        site_pools: dict[str, tuple[float, list[DefectSpecies]]] | None = None,
        element_pools: dict[str, tuple[float, list[tuple[Any, float]]]] | None = None,
        vbm_shift: float = 0.0,
        cbm_shift: float = 0.0,
        formation_energy_corrections: dict[DefectChargeState, float] | None = None,
        rigid_shift: bool = True,
    ):
        self.volume = volume
        self.dos = dos
        self.temperature = temperature
        self.convergence_tolerance = convergence_tolerance
        self.vbm_shift = vbm_shift
        self.cbm_shift = cbm_shift
        self.rigid_shift = rigid_shift

        memo: dict[int, Any] = {}
        self.defect_species = copy.deepcopy(defect_species, memo)
        self._apply_formation_energy_corrections(
            defect_species, formation_energy_corrections or {}
        )

        # share `memo` so any DefectSpecies objects referenced by both
        # `defect_species` and `site_pools`/`element_pools` resolve to the
        # same copies as `self.defect_species`.
        self.site_pools = copy.deepcopy(site_pools, memo) if site_pools else {}
        self.element_pools = copy.deepcopy(element_pools, memo) if element_pools else {}

    def _apply_formation_energy_corrections(
        self,
        original_defect_species: list[DefectSpecies],
        formation_energy_corrections: dict[DefectChargeState, float],
    ) -> None:
        """Permanently shift each variable-concentration charge state's
        formation energy in `self.defect_species` (a copy of
        `original_defect_species`) by `formation_energy_corrections[cs]` if
        present, else `-cs.charge * self.vbm_shift` if `self.rigid_shift` is
        False, else 0.
        """
        original_states = [cs for ds in original_defect_species for cs in ds.charge_states]
        copied_states = [cs for ds in self.defect_species for cs in ds.charge_states]

        unrecognised = set(formation_energy_corrections) - set(original_states)
        if unrecognised:
            raise ValueError(
                f"formation_energy_corrections references "
                f"{len(unrecognised)} DefectChargeState object(s) that are "
                "not part of `defect_species`."
            )

        for original_cs, copied_cs in zip(original_states, copied_states, strict=True):
            if copied_cs._energy is None:
                continue
            if original_cs in formation_energy_corrections:
                delta = formation_energy_corrections[original_cs]
            elif not self.rigid_shift:
                delta = -copied_cs.charge * self.vbm_shift
            else:
                delta = 0.0
            copied_cs._energy += delta

    def __repr__(self) -> str:
        bandgap = self.dos.bandgap + (self.cbm_shift - self.vbm_shift)
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
            for pool_name, (n, species) in self.site_pools.items():
                sp_names = ", ".join(
                    sp.name if isinstance(sp, DefectSpecies) else sp
                    for sp in species
                )
                lines.append(f"    {pool_name}: {n:.4g} sites  \u2192  [{sp_names}]")
        if self.element_pools:
            lines.append("\n  element pools:")
            for elem, (n, pool_list) in self.element_pools.items():
                sp_names = ", ".join(
                    (sp if isinstance(sp, str) else sp.name) + f" \u00d7{stoich:g}"
                    for sp, stoich in pool_list
                )
                lines.append(f"    {elem}: {n:.4g} per cell  \u2192  [{sp_names}]")
        return "\n".join(lines)

    @property
    def defect_species_names(self) -> list[str]:
        """list of the names of all ``DefectSpecies`` considered in the
        ``DefectSystem``.

        Returns:
            list[str]: list of names of ``DefectSpecies`` objects
        """
        return [ds.name for ds in self.defect_species]

    def report(self) -> str:
        """Solve for the self-consistent Fermi energy and print a summary of
        the system: SC Fermi energy, carrier concentrations, and per-species
        and per-charge-state defect concentrations.

        Returns:
            str: the formatted report (also printed to stdout)
        """
        scale = 1e24 / self.volume
        e_fermi, _ = self.get_sc_fermi()
        p0, n0 = self.dos.carrier_concentrations(e_fermi, self.temperature)

        # pool-aware per-charge-state concentrations (per unit cell)
        cs_concs = self._global_defect_concs(e_fermi)

        # group concentrations back to their parent DefectSpecies
        sp_cs: dict[str, list[tuple[DefectChargeState, float]]] = {
            ds.name: [] for ds in self.defect_species
        }
        for ds in self.defect_species:
            for cs in ds.charge_states:
                if cs in cs_concs:
                    sp_cs[ds.name].append((cs, cs_concs[cs]))

        lines = [
            repr(self),
            "",
            f"SC Fermi energy:  {e_fermi:.6f} eV",
            "",
            "Carriers:",
            f"  p0 (holes):     {p0 * scale:.3e} cm\u207b\u00b3",
            f"  n0 (electrons): {n0 * scale:.3e} cm\u207b\u00b3",
            "",
            "Defect concentrations:",
        ]
        for ds in self.defect_species:
            cs_list = sp_cs[ds.name]
            sp_total = sum(c for _, c in cs_list) * scale
            header = f"  {ds.name:<16s}  {sp_total:.3e} cm\u207b\u00b3"
            if ds.fixed_concentration is not None:
                header += "  [fixed]"
            lines.append(header)
            total_per_cell = sum(c for _, c in cs_list)
            for cs, conc in sorted(cs_list, key=lambda x: x[0].charge):
                pct = 100 * conc / total_per_cell if total_per_cell > 0 else 0.0
                flag = "  [fixed]" if cs.fixed_concentration is not None else ""
                lines.append(
                    f"    q = {cs.charge:+2d}   {conc * scale:.3e} cm\u207b\u00b3"
                    f"  ({pct:6.2f}%){flag}"
                )

        output = "\n".join(lines)
        print(output)
        return output

    def _resolve_species(self, sp: DefectSpecies | str) -> DefectSpecies:
        """Resolve a DefectSpecies or its name to a DefectSpecies instance."""
        return self.defect_species_by_name(sp) if isinstance(sp, str) else sp

    def _build_group(
        self,
        label: str,
        n_sites: float,
        species: list[DefectSpecies],
        e_fermi: float,
        fixed_concs: dict[DefectChargeState, float],
    ) -> _ExclusionGroup:
        """Build an exclusion group of `n_sites` sites shared by `species`.

        Species-level `fixed_concentration` and individually-fixed charge
        states are written into `fixed_concs` and subtracted from `n_sites`
        to give the group's free-site budget; every other charge state
        becomes a variable state of the group.
        """
        fixed_total = 0.0
        variable_states: list[tuple[DefectChargeState, DefectSpecies, float]] = []
        for sp in species:
            if sp.fixed_concentration is not None:
                fixed_total += sp.fixed_concentration
                for cs, conc in sp.charge_state_concentrations(e_fermi, self.temperature):
                    fixed_concs[cs] = conc
                continue
            for cs in sp.charge_states:
                if cs.fixed_concentration is not None:
                    fixed_concs[cs] = cs.fixed_concentration
                    fixed_total += cs.fixed_concentration
                else:
                    Ef = cs.get_formation_energy(e_fermi)
                    log_w = np.log(cs.degeneracy) - Ef / (_kboltz * self.temperature)
                    variable_states.append((cs, sp, log_w))

        n_free = n_sites - fixed_total
        if n_free < 0:
            raise ValueError(
                f"'{label}' has {n_sites} sites but fixed-concentration "
                f"states occupy {fixed_total}"
            )
        return _ExclusionGroup(n_free=n_free, variable_states=variable_states)

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
        pooled_species: set[DefectSpecies] = set()

        for pool_name, (n_pool, species_list) in self.site_pools.items():
            sp_objs = [self._resolve_species(sp) for sp in species_list]
            pooled_species.update(sp_objs)
            groups.append(self._build_group(pool_name, n_pool, sp_objs, e_fermi, fixed_concs))

        for sp in self.defect_species:
            if sp not in pooled_species:
                groups.append(self._build_group(sp.name, sp.nsites, [sp], e_fermi, fixed_concs))

        return groups, fixed_concs

    def _resolve_element_pools(
        self,
    ) -> dict[str, tuple[float, list[tuple[DefectSpecies, float]]]]:
        """Resolve string species references in `element_pools` to DefectSpecies."""
        return {
            elem: (target, [(self._resolve_species(sp), stoich) for sp, stoich in pool_list])
            for elem, (target, pool_list) in self.element_pools.items()
        }

    @staticmethod
    def _stoichiometry_lookup(
        pools: dict[str, tuple[float, list[tuple[DefectSpecies, float]]]],
    ) -> dict[DefectSpecies, dict[str, float]]:
        """Map each species to {element: stoichiometry} for every element
        pool it participates in."""
        lookup: dict[DefectSpecies, dict[str, float]] = {}
        for elem, (_, pool_list) in pools.items():
            for sp, stoich in pool_list:
                lookup.setdefault(sp, {})[elem] = stoich
        return lookup

    @staticmethod
    def _remaining_element_targets(
        pools: dict[str, tuple[float, list[tuple[DefectSpecies, float]]]],
    ) -> dict[str, float]:
        """For each element pool, the target content still to be supplied by
        variable-concentration states, after subtracting fixed contributions."""
        remaining: dict[str, float] = {}
        for elem, (target, pool_list) in pools.items():
            committed = 0.0
            for sp, stoich in pool_list:
                if sp.fixed_concentration is not None:
                    committed += stoich * sp.fixed_concentration
                else:
                    committed += stoich * sum(
                        cs.fixed_concentration
                        for cs in sp.charge_states
                        if cs.fixed_concentration is not None
                    )
            rem = target - committed
            if rem < 0:
                raise ValueError(
                    f"Element pool '{elem}': fixed-concentration states "
                    f"already contribute {committed:.3e} which exceeds the "
                    f"target {target:.3e}. Your constraints are mutually "
                    "inconsistent."
                )
            remaining[elem] = rem
        return remaining

    @staticmethod
    def _group_term_arrays(
        variable_states: list[tuple[DefectChargeState, DefectSpecies, float]],
        elements: list[str],
        stoich: dict[DefectSpecies, dict[str, float]],
    ) -> tuple[np.ndarray, np.ndarray]:
        """Per-state log Boltzmann weights and stoichiometry vectors for a
        group's variable states, as arrays of shape (n,) and (n, len(elements))."""
        log_w = np.array([lw for _, _, lw in variable_states])
        s = np.array(
            [
                [stoich.get(sp, {}).get(elem, 0.0) for elem in elements]
                for _, sp, _ in variable_states
            ]
        )
        return log_w, s

    @staticmethod
    def _group_concs(
        n_free: float, log_w: np.ndarray, s: np.ndarray, mu: np.ndarray
    ) -> np.ndarray:
        """c_i = n_free * w_i * exp(s_i . mu) / (1 + sum_j w_j * exp(s_j . mu))
        for every variable state in a group, given the element chemical
        potentials `mu`."""
        exponents = log_w + s @ mu
        log_z = logsumexp(np.concatenate(([0.0], exponents)))
        return n_free * np.exp(exponents - log_z)

    def _solve_chemical_potentials(
        self,
        groups: list[_ExclusionGroup],
        elements: list[str],
        stoich: dict[DefectSpecies, dict[str, float]],
        remaining: dict[str, float],
        pools: dict[str, tuple[float, list[tuple[DefectSpecies, float]]]],
    ) -> np.ndarray:
        """Solve for one chemical potential `mu_X` per element pool such
        that, summing `_group_concs` over every group, the total content of
        each pooled element equals `remaining[X]`.

        This is the stationarity condition of the convex grand potential
        ``F(mu) = sum_g n_free_g * log(Z_g(mu)) - mu . remaining``, whose
        gradient is the element-content vector and whose Hessian is, group
        by group, the covariance matrix of the states' stoichiometry vectors
        weighted by their occupancy -- symmetric positive semi-definite, so
        `F` has a single minimum whenever `remaining` is achievable. The
        stationarity condition is solved as the per-element-scaled root
        problem ``content_X(mu) / remaining[X] - 1 = 0``, with the scaled
        Hessian as its analytic Jacobian.
        """
        K = len(elements)
        remaining_vec = np.array([remaining[e] for e in elements])

        group_data = [
            (group.n_free, *self._group_term_arrays(group.variable_states, elements, stoich))
            for group in groups
            if group.n_free > 0 and group.variable_states
        ]

        # feasibility: the most of element X any group can supply is
        # n_free * (largest stoichiometry among its variable states).
        for k, elem in enumerate(elements):
            max_content = sum(
                n_free * max(0.0, s[:, k].max()) for n_free, _, s in group_data
            )
            if remaining_vec[k] > max_content:
                target, _ = pools[elem]
                committed = target - remaining_vec[k]
                raise ValueError(
                    f"Element pool '{elem}': target {target:.3e} exceeds "
                    f"the maximum achievable content "
                    f"{committed + max_content:.3e} (committed "
                    f"{committed:.3e} plus up to {max_content:.3e} from "
                    "variable states)."
                )

        def content_and_hessian(mu: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
            content = np.zeros(K)
            hessian = np.zeros((K, K))
            for n_free, log_w, s in group_data:
                c = self._group_concs(n_free, log_w, s, mu)
                content += s.T @ c
                mean_s = (s * c[:, None]).sum(axis=0) / n_free
                ds = s - mean_s
                hessian += (ds * c[:, None]).T @ ds
            return content, hessian

        def residual_and_jacobian(mu: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
            content, hessian = content_and_hessian(mu)
            return content / remaining_vec - 1.0, hessian / remaining_vec[:, None]

        # Defect concentrations natively span ~1e-30..1 per cell, and
        # scipy's default convergence criteria assume O(1) problems:
        # bracketing methods are scale-free, but gradient-norm and
        # absolute-residual thresholds are not, and would accept mu = 0
        # unchanged for any dilute target. Scaling each equation by its own
        # target keeps every residual O(1), however dilute the targets and
        # however many orders of magnitude separate them.
        #
        # The initial guess is one diagonal Newton step of the log-content
        # system ``log(content_X(mu)) = log(remaining[X])`` from mu = 0,
        # which inverts the dilute limit ``content_X ~ content_X(0) *
        # exp(sbar_X * mu_X)`` with characteristic stoichiometry
        # ``sbar_X = H_XX / content_X``. Started from mu = 0 itself, the
        # solver can sit on an underflow plateau where the Jacobian
        # vanishes.
        c0, h0 = content_and_hessian(np.zeros(K))
        with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
            x0 = np.log(remaining_vec / c0) * c0 / np.diag(h0)
        x0 = np.clip(
            np.nan_to_num(x0, nan=700.0, posinf=700.0, neginf=-700.0), -700.0, 700.0
        )
        result = root(residual_and_jacobian, x0=x0, jac=True, method="hybr")
        if not result.success:
            targets = ", ".join(f"{elem}={remaining[elem]:.3e}" for elem in elements)
            raise ValueError(
                f"Element pool targets ({targets}) could not be satisfied: "
                f"the chemical-potential solve did not converge "
                f"({result.message}). The targets may be jointly infeasible."
            )

        # Re-check the solution independently of the solver's own verdict:
        # an optimiser returning success without moving from its starting
        # point is exactly the failure mode that motivated the
        # dimensionless formulation, and any recurrence must be loud.
        achieved, _ = content_and_hessian(result.x)
        rel_err = np.abs(achieved / remaining_vec - 1.0)
        worst = int(np.argmax(rel_err))
        if rel_err[worst] > _element_pool_tolerance:
            elem = elements[worst]
            target, _ = pools[elem]
            committed = target - remaining_vec[worst]
            raise ValueError(
                f"Element pool '{elem}': the chemical-potential solve "
                f"reported success but the target was not met (target "
                f"{target:.6e}, achieved {committed + achieved[worst]:.6e})."
            )
        return result.x

    def _global_defect_concs(self, e_fermi: float) -> dict[DefectChargeState, float]:
        """
        Returns a dict mapping each DefectChargeState -> concentration per cell.

        Every species is assigned to a site-exclusion group (a `site_pools`
        entry, or its own implicit single-species group of `nsites` sites)
        and its variable charge states are given Langmuir/site-exclusion
        statistics, ``c_i = n_free * w_i / (1 + sum_j w_j)``. If
        `element_pools` are defined, the element chemical potentials are
        solved for first (`_solve_chemical_potentials`) and folded into each
        state's weight as ``w_i * exp(s_i . mu)``.
        """
        groups, all_concs = self._build_exclusion_groups(e_fermi)

        pools = self._resolve_element_pools()
        elements = list(pools.keys())
        stoich = self._stoichiometry_lookup(pools)

        forced_zero: set[DefectSpecies] = set()
        mu = np.array([])
        if elements:
            remaining = self._remaining_element_targets(pools)
            # An element whose target is already fully committed by fixed
            # concentrations admits no further variable content. In the
            # lambda_X -> 0 limit, every variable state of a species with
            # positive stoichiometry in X has concentration 0, and X drops
            # out of the chemical-potential solve.
            exhausted = {e for e in elements if remaining[e] == 0.0}
            solve_groups = groups
            if exhausted:
                forced_zero = {
                    sp
                    for sp, s_by_elem in stoich.items()
                    if any(s_by_elem.get(e, 0.0) > 0.0 for e in exhausted)
                }
                elements = [e for e in elements if e not in exhausted]
                remaining = {e: remaining[e] for e in elements}
                solve_groups = [
                    _ExclusionGroup(
                        group.n_free,
                        [
                            state
                            for state in group.variable_states
                            if state[1] not in forced_zero
                        ],
                    )
                    for group in groups
                ]
            if elements:
                mu = self._solve_chemical_potentials(
                    solve_groups, elements, stoich, remaining, pools
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
                if state[1] in forced_zero:
                    all_concs[state[0]] = 0.0
                else:
                    kept.append(state)
            if not kept:
                continue
            log_w, s = self._group_term_arrays(kept, elements, stoich)
            for (cs, _, _), c_i in zip(
                kept,
                self._group_concs(group.n_free, log_w, s, mu),
                strict=True,
            ):
                all_concs[cs] = c_i

        return all_concs

    @classmethod
    def from_dict(cls, dictionary: dict) -> DefectSystem:
        """Generate a DefectSystem from a dictionary.
    
        Args:
            dictionary (dict): Dictionary containing the DefectSystem data.
    
        Returns:
            DefectSystem: DefectSystem corresponding to the provided dictionary.
        """
        return cls(
            dos=DOS.from_dict(dictionary["dos"]),
            volume=dictionary["volume"],
            temperature=dictionary["temperature"],
            convergence_tolerance=dictionary.get("convergence_tolerance"),
            defect_species=[
                DefectSpecies.from_dict(defect_species)
                for defect_species in dictionary["defect_species"]
            ],
        )
        
    def defect_species_by_name(self, name: str) -> DefectSpecies:
        """return a ``DefectSpecies`` contained within the ``DefectSystem``
        via its name.

        Args:
            name (str): name of the ``DefectSpecies`` to return

        Returns:
            DefectSpecies: ``DefectSpecies`` where ``DefectSpecies.name == name``
        """
        return [ds for ds in self.defect_species if ds.name == name][0]
    
    def get_sc_fermi(self) -> tuple[float, float]:
        """Calculate the self-consistent Fermi energy.
        
        Finds the Fermi energy at which charge neutrality is achieved,
        using Brent's method for root finding.
        
        Returns:
            tuple[float, float]: The self-consistent Fermi energy and the
            absolute residual charge density at that energy.
        
        Raises:
            RuntimeError: If no solution is found within the DOS energy range.
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
        except ValueError as err:
            raise RuntimeError(
                f"No solution found between {emin} and {emax}: {err}"
            ) from err
        
        # q_tot is continuous wherever finite, so a sign change across the DOS
        # window brackets a genuine root. The returned residual is a diagnostic
        # on that solution.
        residual = abs(self.q_tot(e_fermi))
        return e_fermi, residual

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

    def concentration_dict(
        self,
        decomposed: bool = False,
        per_volume: bool = True,
    ) -> dict[str, Any]:
        """Returns a dictionary of the properties of the ``DefectSystem`` object
        after solving for the self-consistent Fermi energy.

        Args:
            decomposed (bool, optional): if True, return a dictionary in which the
              concentration of each ``DefectChargeState`` is given explicitly,
              rather than as a sum over all ``DefectChargeState`` objects in the
              each ``DefectSpecies``. Defaults to False.
            per_volume (bool, optional): if True, return concentrations in units
              of cm^-3, else returns concentration per unit cell. Defaults to True.

        Returns:
            dict[str, Any]: dictionary specifying the Fermi Energy,
            hole concentration (``"p0"``), electron concentration
            (``"n0"``), temperature, and the defect concentrations.
        """
        if per_volume:
            scale = 1e24 / self.volume
        else:
            scale = 1

        e_fermi = self.get_sc_fermi()[0]
        p0, n0 = self.dos.carrier_concentrations(e_fermi, self.temperature)
        run_stats = {
            "Fermi Energy": float(e_fermi),
            "p0": float(p0 * scale),
            "n0": float(n0 * scale),
        }
        cs_concs = self._global_defect_concs(e_fermi)
        if not decomposed:
            sum_concs = {}
            for ds in self.defect_species:
                total = sum(cs_concs.get(cs, 0.0) for cs in ds.charge_states)
                sum_concs[str(ds.name)] = float(total * scale)
            return {**run_stats, **sum_concs}
        else:
            decomp_concs: dict[str, dict[int, float]] = {}
            for ds in self.defect_species:
                by_charge: dict[int, float] = {}
                for cs in ds.charge_states:
                    if cs in cs_concs:
                        by_charge[cs.charge] = (
                            by_charge.get(cs.charge, 0.0) + float(cs_concs[cs] * scale)
                        )
                decomp_concs[str(ds.name)] = by_charge
            return {**run_stats, **decomp_concs}

    def site_percentages(
        self, 
    ) -> dict[str, float]:
        """Returns a dictionary of the DefectSpecies in the DefectSystem which
        giving the percentage of the sites in the structure that will host that 
        defect.

        Returns:
            dict[str, Any]: dictionary specifying the per-DefectSpecies site
            concentrations.
        """
        e_fermi = self.get_sc_fermi()[0]

        sum_concs = {
                str(ds.name): float(
                    (ds.get_concentration(e_fermi, self.temperature) / ds.nsites) * 100
                )
                for ds in self.defect_species
            }
        return sum_concs

    def as_dict(self) -> dict:
        """Return a dictionary representation of the ``DefectSystem``.

        Returns:
            dict: dictionary representation of the ``DefectSystem``, suitable for
            round-tripping through ``DefectSystem.from_dict``.
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
            defect_system_dict["convergence_tolerance"] = self.convergence_tolerance
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
        site_pools (dict[str, tuple[float, list[DefectSpecies]]] | None, optional):
          passed to every `DefectSystem` built by `at()`.
        element_pools (dict[str, tuple[float, list[tuple[Any, float]]]] | None, optional):
          passed to every `DefectSystem` built by `at()`.
        vbm_shift_fn (Callable[[float], float] | None, optional): a function
          of temperature returning the valence-band-maximum shift (in eV) at
          that temperature. Defaults to None (no shift).
        cbm_shift_fn (Callable[[float], float] | None, optional): a function
          of temperature returning the conduction-band-minimum shift (in eV)
          at that temperature. Defaults to None (no shift).
        formation_energy_correction_fns (dict[DefectChargeState, Callable] | None, optional):
          per-charge-state functions of temperature returning a
          formation-energy correction (in eV) at that temperature. Keys must
          be `DefectChargeState`s in `defect_species`. Defaults to None.
        rigid_shift (bool, optional): passed to every `DefectSystem` built by
          `at()`. Defaults to True.
    """

    def __init__(
        self,
        defect_species: list[DefectSpecies],
        dos: DOS,
        volume: float,
        convergence_tolerance: float | None = None,
        site_pools: dict[str, tuple[float, list[DefectSpecies]]] | None = None,
        element_pools: dict[str, tuple[float, list[tuple[Any, float]]]] | None = None,
        vbm_shift_fn: Callable[[float], float] | None = None,
        cbm_shift_fn: Callable[[float], float] | None = None,
        formation_energy_correction_fns: (
            dict[DefectChargeState, Callable[[float], float]] | None
        ) = None,
        rigid_shift: bool = True,
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
              `fixed_concentrations` on individual `DefectSpecies`, or
              `temperature` itself).

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
        )
        kwargs.update(overrides)
        return DefectSystem(**kwargs)
