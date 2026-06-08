import warnings
from collections.abc import Callable
from contextlib import contextmanager
from typing import Any

import numpy as np
from scipy.constants import physical_constants as _physical_constants
from scipy.optimize import brentq  # type: ignore[call-overload]
from scipy.special import logsumexp

from py_sc_fermi.defect_charge_state import DefectChargeState
from py_sc_fermi.defect_species import DefectSpecies
from py_sc_fermi.dos import DOS

_kboltz = _physical_constants["Boltzmann constant in eV/K"][0]


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
        n_trial_steps (int, optional): Deprecated. Previously set the maximum
          number of solver iterations. The solver now uses Brent's method
          which converges reliably without this parameter.
        site_pools (dict[str, tuple[float, list[DefectSpecies]]] | None, optional):
                Mapping of pool name → (total sites in that pool, list of
                DefectSpecies sharing those sites). If None, no site competition
                is applied and each species is treated independently.
    """

    def __init__(
        self,
        defect_species: list[DefectSpecies],
        dos: DOS,
        volume: float,
        temperature: float,
        convergence_tolerance: float | None = None,
        n_trial_steps: int | None = None, # deprecated
        site_pools: dict[str, tuple[float, list[DefectSpecies]]] | None = None,
        element_pools: dict[str, tuple[float, list[tuple[Any, float]]]] | None = None,
        vbm_shift_fn: Callable[[float], float] | None = None,
        cbm_shift_fn: Callable[[float], float] | None = None,
        formation_energy_corrections: dict[tuple[str, int], Callable[[float], float]] | None = None,
        rigid_shift: bool = True,
    ):
        self.defect_species = defect_species
        self.volume = volume
        self.dos = dos
        self.temperature = temperature
        self.convergence_tolerance = convergence_tolerance
        self.vbm_shift_fn = vbm_shift_fn
        self.cbm_shift_fn = cbm_shift_fn
        self.formation_energy_corrections = formation_energy_corrections or {}
        self.rigid_shift = rigid_shift
        self._corrections_active = False
        if n_trial_steps is not None:
            warnings.warn(
                "n_trial_steps is deprecated and will be removed in a future version. "
                "The solver now uses Brent's method which converges reliably without "
                "this parameter.",
                DeprecationWarning,
                stacklevel=2,
            )
        self.n_trial_steps = n_trial_steps

        self.site_pools = site_pools or {}
        self.element_pools = element_pools or {}

    def __repr__(self):
        # Compute corrections at current temperature (if any)
        has_corrections = self.vbm_shift_fn is not None or self.cbm_shift_fn is not None
        if has_corrections:
            d_vbm = self.vbm_shift_fn(self.temperature) if self.vbm_shift_fn else 0.0
            d_cbm = self.cbm_shift_fn(self.temperature) if self.cbm_shift_fn else 0.0
            corrected_bandgap = self.dos.bandgap + (d_cbm - d_vbm)
        else:
            d_vbm = 0.0
            corrected_bandgap = self.dos.bandgap
        has_fe_corrections = bool(self.formation_energy_corrections) or not self.rigid_shift

        if has_corrections:
            bandgap_str = (
                f"{self.dos.bandgap:.3g} eV  \u2192  {corrected_bandgap:.3g} eV"
                f" at {self.temperature} K"
            )
        else:
            bandgap_str = f"{self.dos.bandgap:.3g} eV"

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
                    e_stored = cs.energy
                    if has_fe_corrections:
                        key = (ds.name, cs.charge)
                        if key in self.formation_energy_corrections:
                            delta = self.formation_energy_corrections[key](self.temperature)
                        elif not self.rigid_shift:
                            delta = -cs.charge * d_vbm
                        else:
                            delta = 0.0
                        e_corr = e_stored + delta
                        energy_str = (
                            f"E = {e_stored:8.4f} eV  \u2192  {e_corr:8.4f} eV"
                            f" at {self.temperature} K"
                        )
                    else:
                        energy_str = f"E = {e_stored:8.4f} eV"
                    lines.append(
                        f"    q = {cs.charge:+2d}   {energy_str}  (deg. {cs.degeneracy})"
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

    @contextmanager
    def _with_band_edge_corrections(self):
        """Temporarily apply VBM/CBM shift corrections at the current temperature.

        Re-entrant: nested calls (e.g. report → get_sc_fermi) are no-ops once
        the corrections are already active, so corrections are never applied twice.
        """
        if self._corrections_active or (
            self.vbm_shift_fn is None and self.cbm_shift_fn is None
        ):
            yield
            return

        self._corrections_active = True
        applied: dict[tuple[str, int], float] = {}
        bandgap_delta = 0.0
        try:
            d_vbm = self.vbm_shift_fn(self.temperature) if self.vbm_shift_fn else 0.0
            d_cbm = self.cbm_shift_fn(self.temperature) if self.cbm_shift_fn else 0.0

            # Apply formation energy corrections per charge state.
            # By default (rigid_shift=True) only the band gap changes; formation
            # energies are unchanged unless explicit per-charge-state corrections
            # are provided. With rigid_shift=False the legacy q·d_vbm correction is
            # used as a fallback.
            for ds in self.defect_species:
                for cs in ds.charge_states:
                    if cs._energy is None:
                        continue
                    key = (ds.name, cs.charge)
                    if key in self.formation_energy_corrections:
                        delta = self.formation_energy_corrections[key](self.temperature)
                    elif not self.rigid_shift:
                        delta = -cs.charge * d_vbm
                    else:
                        delta = 0.0
                    cs._energy += delta
                    applied[key] = delta

            bandgap_delta = d_cbm - d_vbm
            self.dos._bandgap += bandgap_delta

            yield
        finally:
            for ds in self.defect_species:
                for cs in ds.charge_states:
                    if cs._energy is not None:
                        cs._energy -= applied.get((ds.name, cs.charge), 0.0)
            self.dos._bandgap -= bandgap_delta
            self._corrections_active = False

    def report(self) -> str:
        """Solve for the self-consistent Fermi energy and print a summary of
        the system: SC Fermi energy, carrier concentrations, and per-species
        and per-charge-state defect concentrations.

        Returns:
            str: the formatted report (also printed to stdout)
        """
        with self._with_band_edge_corrections():
            return self._report_inner()

    def _report_inner(self) -> str:
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

    def _apply_element_constraints(
        self,
        concs: dict[DefectChargeState, float],
        e_fermi: float,
    ) -> None:
        """Rescale variable-concentration charge state concentrations so that each
        element pool’s total defect content matches its fixed target.

        Site pool competition (if any) is applied before this method is called, so
        site pools take priority. This method then scales all variable-concentration
        states of the element’s species together to satisfy the elemental content
        constraint. Fixed-concentration charge states are left unchanged.

        Args:
            concs: mapping from every DefectChargeState → concentration per cell.
                   Modified in place.
        """
        # build set of species already governed by a site pool
        site_pooled_species: set = {
            self.defect_species_by_name(sp) if isinstance(sp, str) else sp
            for _, (_, sps) in self.site_pools.items()
            for sp in sps
        }

        for elem, (fixed_total, pool_list) in self.element_pools.items():
            # resolve string species names
            resolved: list[tuple[DefectSpecies, float]] = [
                (self.defect_species_by_name(sp) if isinstance(sp, str) else sp, stoich)
                for sp, stoich in pool_list
            ]

            # split species into site-pool-governed (fixed for our purposes) and free
            site_governed = [(sp, stoich) for sp, stoich in resolved if sp in site_pooled_species]
            free = [(sp, stoich) for sp, stoich in resolved if sp not in site_pooled_species]

            # 1) elemental content already committed by site-pool species and fixed-cs species
            committed = sum(
                sum(concs.get(cs, 0.0) for cs in sp.charge_states) * stoich
                for sp, stoich in site_governed
            )
            committed += sum(
                sum(
                    concs.get(cs, 0.0)
                    for cs in sp.charge_states
                    if cs.fixed_concentration is not None
                )
                * stoich
                for sp, stoich in free
            )

            remaining = fixed_total - committed
            if remaining < 0:
                raise ValueError(
                    f"Element pool ‘{elem}’: site pool and fixed-concentration states "
                    f"already contribute {committed:.3e} which exceeds the target "
                    f"{fixed_total:.3e}. Your constraints are mutually inconsistent."
                )

            # 2) compute log-weights for all free variable states (log-space to avoid overflow)
            log_contributions: list[tuple] = []  # (cs, stoich, log_w)
            for sp, stoich in free:
                for cs in sp.charge_states:
                    if cs in concs and cs.fixed_concentration is None and cs._energy is not None:
                        Ef = cs.get_formation_energy(e_fermi)
                        log_w = (
                            np.log(sp.nsites)
                            + np.log(cs.degeneracy)
                            - Ef / (_kboltz * self.temperature)
                            + np.log(stoich)
                        )
                        log_contributions.append((cs, stoich, log_w))

            if not log_contributions:
                continue

            # log of total elemental content from free variable states
            log_current_var = logsumexp([lw for _, _, lw in log_contributions])

            # 3) assign concentrations: remaining * (unnorm_w / total_elem_w)
            #    where unnorm_w = nsites * deg * exp(-Ef/kT), total_elem_w sums stoich * unnorm_w
            log_remaining = np.log(remaining)
            for cs, stoich, log_w in log_contributions:
                concs[cs] = np.exp(log_remaining + (log_w - np.log(stoich)) - log_current_var)

    def _global_defect_concs(self, e_fermi: float):
        """
        Returns a dict mapping each DefectChargeState → concentration per cell,
        applying site-competition for any pools defined in self.site_pools.
        Pool entries may be either DefectSpecies instances or their name strings.
        """
        all_concs= {}

        # 1) Handle each pool
        for pool_name, (N_pool, species_list) in self.site_pools.items():
            # normalize list to DefectSpecies objects
            sp_objs = []
            for sp in species_list:
                if isinstance(sp, str):
                    sp_objs.append(self.defect_species_by_name(sp))
                else:
                    sp_objs.append(sp)

            # a) fixed occupancy per species
            fixed_per_sp = {
                sp: sum(
                    conc
                    for cs, conc in sp.charge_state_concentrations(
                        e_fermi, self.temperature
                    )
                    if cs.fixed_concentration is not None
                )
                for sp in sp_objs
            }
            total_fixed = sum(fixed_per_sp.values())
            free_sites = N_pool - total_fixed
            if free_sites < 0:
                raise ValueError(
                    f"Pool '{pool_name}' has {N_pool} sites but fixed states "
                    f"occupy {total_fixed}"
                )

            # b) log Boltzmann weight per species and per variable state
            #    (log-space to avoid overflow at extreme Fermi energies)
            sp_log_ws: dict = {}  # sp -> list of (cs, log_w)
            for sp in sp_objs:
                sp_log_ws[sp] = []
                for cs in sp.variable_conc_charge_states():
                    Ef = cs.get_formation_energy(e_fermi)
                    log_w = (
                        np.log(sp.nsites)
                        + np.log(cs.degeneracy)
                        - Ef / (_kboltz * self.temperature)
                    )
                    sp_log_ws[sp].append((cs, log_w))

            # log of total weight per species (-inf if no variable states)
            sp_log_total = {
                sp: logsumexp([lw for _, lw in pairs]) if pairs else -np.inf
                for sp, pairs in sp_log_ws.items()
            }

            # c) partition function: log(1 + sum_sp w_sp)
            log_Z = logsumexp([0.0] + list(sp_log_total.values()))

            # d) assign each species its share
            for sp in sp_objs:
                # fixed states pass through unchanged
                for cs, conc in sp.charge_state_concentrations(
                    e_fermi, self.temperature
                ):
                    if cs.fixed_concentration is not None:
                        all_concs[cs] = conc
                # variable states: share proportional to Boltzmann weight
                if sp_log_ws[sp]:
                    log_share = np.log(free_sites) + sp_log_total[sp] - log_Z
                    for cs, log_w_i in sp_log_ws[sp]:
                        all_concs[cs] = np.exp(log_share + log_w_i - sp_log_total[sp])

        # 2) Species not in any pool: old dilute‐limit
        pooled = {
            self.defect_species_by_name(sp) if isinstance(sp, str) else sp
            for _, (_, sps) in self.site_pools.items() for sp in sps
        }
        for sp in self.defect_species:
            if sp not in pooled:
                for cs, conc in sp.charge_state_concentrations(e_fermi, self.temperature):
                    all_concs[cs] = conc

        if self.element_pools:
            self._apply_element_constraints(all_concs, e_fermi)

        return all_concs

    @classmethod
    def from_dict(cls, dictionary: dict) -> "DefectSystem":
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
            n_trial_steps=dictionary.get("n_trial_steps"),
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
        with self._with_band_edge_corrections():
            return self._sc_fermi_solve()

    def _sc_fermi_solve(self) -> tuple[float, float]:
        emin = self.dos.emin()
        emax = self.dos.emax()
        
        try:
            kwargs = {}
            if self.convergence_tolerance is not None:
                kwargs['xtol'] = self.convergence_tolerance
            if self.n_trial_steps is not None:
                kwargs['maxiter'] = self.n_trial_steps
            
            e_fermi = brentq(
                self.q_tot,
                emin,
                emax,
                **kwargs,
            )  # type: ignore[call-overload]
        except ValueError as err:
            raise RuntimeError(f"No solution found between {emin} and {emax}") from err
        
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
        with self._with_band_edge_corrections():
            return self._transition_levels_inner()

    def _transition_levels_inner(self) -> dict[str, list[list]]:
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
        with self._with_band_edge_corrections():
            return self._concentration_dict_inner(decomposed=decomposed, per_volume=per_volume)

    def _concentration_dict_inner(
        self, decomposed: bool = False, per_volume: bool = True
    ) -> dict[str, Any]:
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
        with self._with_band_edge_corrections():
            return self._site_percentages_inner()

    def _site_percentages_inner(self) -> dict[str, float]:
        e_fermi = self.get_sc_fermi()[0]

        sum_concs = {
                str(ds.name): float(
                    (ds.get_concentration(e_fermi, self.temperature) / ds.nsites) * 100
                )
                for ds in self.defect_species
            }
        return sum_concs

    def as_dict(self) -> dict:
        """

        Returns:
            dict: _description_
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
        if self.n_trial_steps is not None:
            defect_system_dict["n_trial_steps"] = self.n_trial_steps
        return defect_system_dict
