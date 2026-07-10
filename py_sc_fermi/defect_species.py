from __future__ import annotations

import math
from collections import Counter

import numpy as np
from scipy.constants import physical_constants
from scipy.special import logsumexp

from py_sc_fermi.defect_charge_state import DefectChargeState

kboltz = physical_constants["Boltzmann constant in eV/K"][0]

class DefectSpecies:
    """Class for individual defect species.

    Args:
        name (str): A unique identifying string for this defect species,
           e.g. ``"V_O"`` might be used for an oxygen vacancy.
        nsites (int): Number of sites energetically degenerate sites where this
         defect can form in the unit cell (the site degeneracy).
        charge_states (list[DefectChargeState]): A list of
           ``DefectChargeState`` objects belonging to this defect species.
           Multiple charge states may share the same formal charge, to
           represent metastable defect configurations; such states must be
           given explicit, distinct names (charge-state names, including the
           charge-derived defaults, must be unique within a species).

    """

    def __init__(
        self,
        name: str,
        nsites: int,
        charge_states: list[DefectChargeState],
        fixed_concentration: float | None = None,
    ):
        """Instantiate a DefectSpecies object."""

        if not charge_states:
            raise ValueError(
                f"DefectSpecies '{name}' must have at least one charge state."
            )
        if nsites <= 0:
            raise ValueError(
                f"DefectSpecies '{name}' must have nsites > 0; got {nsites}."
            )
        name_counts = Counter(cs.name for cs in charge_states)
        duplicates = sorted(n for n, count in name_counts.items() if count > 1)
        if duplicates:
            raise ValueError(
                f"DefectSpecies '{name}' has duplicate charge-state names: "
                f"{', '.join(duplicates)}. Charge-state names must be unique "
                "within a species; metastable states sharing a formal charge "
                "must be given explicit names."
            )
        self._name = name
        self._nsites = nsites
        self._charge_states = list(charge_states)
        self._fixed_concentration = fixed_concentration

    def fix_concentration(self, concentration: float) -> None:
        """fix the concentration of this ``DefectSpecies``

        Args:
            concentration (float): concentration per unit cell
        """
        self._fixed_concentration = concentration

    @property
    def name(self) -> str:
        """identifying string for this ``DefectSpecies``

        Returns:
            str: label for ``DefectSpecies``
        """
        return self._name

    @property
    def nsites(self) -> int:
        """site degeneracy of this ``DefectSpecies`` in the unit cell.

        Returns:
            int: site degeneracy for ``DefectSpecies``
        """
        return self._nsites

    @property
    def charge_states(
        self,
    ) -> tuple[DefectChargeState, ...]:
        """

        Returns:
            tuple[DefectChargeState, ...]: the ``DefectChargeState`` objects
            that comprise this ``DefectSpecies``
        """
        return tuple(self._charge_states)

    def charge_state_by_name(self, name: str) -> DefectChargeState:
        """Return the ``DefectChargeState`` in this species with the given name.

        Args:
            name (str): name of the ``DefectChargeState`` to return.

        Returns:
            DefectChargeState: the charge state where ``cs.name == name``.

        Raises:
            ValueError: if no charge state in this species has that name.
        """
        for cs in self._charge_states:
            if cs.name == name:
                return cs
        available = ", ".join(cs.name for cs in self._charge_states)
        raise ValueError(
            f"DefectSpecies '{self.name}' has no charge state named '{name}'; "
            f"available: {available}"
        )

    @property
    def charges(self) -> list[int]:
        """list of all the charges of the ``DefectChargeState`` objects that
        comprise this ``DefectSpecies``

        Returns:
            list[int]: list of charge states of this ``DefectSpecies``
        """
        return [cs.charge for cs in self._charge_states]
    
    def site_weights(
        self, e_fermi: float, temperature: float
    ) -> list[tuple[float, DefectChargeState, float]]:
        """
        For *variable* charge‐states only, return
            (nsites, DefectChargeState, weight)
        where
            weight = g * exp( - E_f(e_fermi) / (k_B * T) )
        """
        weights: list[tuple[float, DefectChargeState, float]] = []
        for cs in self.variable_conc_charge_states():
            Ef = cs.get_formation_energy(e_fermi)
            w  = cs.degeneracy * np.exp(-Ef / (kboltz * temperature))
            weights.append((self.nsites, cs, w))
        return weights


    @property
    def fixed_concentration(self) -> float | None:
        """fixed concentration of the ``DefectSpecies``. ``None`` if the
        concentration of this defect is variable.

        Returns:
            float | None: fixed concentration per unit cell of the
            ``DefectSpecies``
        """
        return self._fixed_concentration

    def __repr__(self) -> str:
        to_return = f"\n{self.name}, nsites={self.nsites}"
        if self.fixed_concentration is not None:
            to_return += f"\nfixed [c] = {self.fixed_concentration}"
        to_return += "\n" + "".join(
            [f"  {cs.__repr__()}\n" for cs in self.charge_states]
        )
        return to_return

    @classmethod
    def from_dict(cls, d: dict) -> DefectSpecies:
        """return a ``DefectSpecies`` object from a dictionary containing the defect
        species data, as produced by ``as_dict``.

        Args:
            defect_species_dict (dict): dictionary containing the defect species
               data.

        Raises:
            ValueError: if the dictionary specifies no charge states, or if any
               of the specified ``DefectChargeState`` objects have no fixed
               concentration and no formation energy

        Returns:
            DefectSpecies: as specified by the provided dictionary
        """
        states = [
            DefectChargeState.from_dict(cs_dict) for cs_dict in d["charge_states"]
        ]
        return cls(
            name=d["name"],
            nsites=d["nsites"],
            charge_states=states,
            fixed_concentration=d.get("fixed_concentration", None),
        )

    def charge_states_by_formation_energy(
        self, e_fermi: float
    ) -> list[DefectChargeState]:
        """
        Return all *variable* DefectChargeState objects sorted by formation
        energy at a given Fermi energy.

        Args:
            e_fermi (float): Fermi energy at which to evaluate formation energies.

        Returns:
            list[DefectChargeState]: all variable charge‐states of this species,
            sorted from lowest to highest formation energy at e_fermi.
        """
        return sorted(
            self.variable_conc_charge_states(),
            key=lambda st: st.get_formation_energy(e_fermi),
        )

    def as_dict(self) -> dict:
        """get representation of ``DefectSpecies`` as a dictionary

        Returns:
            dict: dictionary representation of ``DefectSpecies``
        """

        defect_dict = {
            "name": str(self.name),
            "nsites": int(self.nsites),
            "charge_states": [cs.as_dict() for cs in self._charge_states]
        }
        if self.fixed_concentration is not None:
            defect_dict.update({"fixed_concentration": float(self.fixed_concentration)})

        return defect_dict

    def min_energy_charge_state(self, e_fermi: float) -> DefectChargeState:
        """Returns the defect charge state with the minimum energy at a given
        Fermi energy.

        Args:
            e_fermi (float): Fermi Energy

        Returns:
            DefectChargeState: the ``DefectChargeState`` of this ``DefectSpecies``
            with the lowest energy at ``e_fermi``.

        Raises:
            ValueError: if this ``DefectSpecies`` has no variable-concentration
                charge states, as the minimum-energy charge state is then
                undefined.
        """
        charge_states = self.charge_states_by_formation_energy(e_fermi)
        if not charge_states:
            raise ValueError(
                f"DefectSpecies '{self.name}' has no variable-concentration "
                "charge states, so a minimum-energy charge state is undefined."
            )
        return charge_states[0]

    def effective_formation_energy(
        self, e_fermi: float, temperature: float = 0.0
    ) -> float:
        """Effective formation energy of this ``DefectSpecies``, summed over
        all of its variable-concentration charge states and metastable forms,
        at a given Fermi energy and temperature.

        Unlike ``get_formation_energies``, which groups charge states by
        formal charge, this combines *every* variable-concentration charge
        state into a single value:
        ``F_d(E_F) = -kT * ln(sum_i g_i * exp(-E_i(E_F) / kT))``.

        At ``temperature=0`` (the default), this is the minimum formation
        energy over all charge states at ``e_fermi`` -- i.e. the standard
        "lower envelope" formation-energy-vs-Fermi-energy curve, evaluable at
        any ``e_fermi`` rather than only at the kinks returned by
        ``tl_profile``. At ``temperature>0``, this is the smooth curve
        ``F_d = -kT * ln(c_total / nsites)`` implied by the total
        concentration of this species.

        Args:
            e_fermi (float): Fermi energy at which to evaluate the effective
                formation energy.
            temperature (float, optional): temperature for the
                Boltzmann-weighted sum. Defaults to 0.

        Returns:
            float: the effective formation energy of this ``DefectSpecies``
            at ``e_fermi``.

        Raises:
            ValueError: if this ``DefectSpecies`` has no variable-concentration
                charge states, so an effective formation energy is undefined.
        """
        states = self.variable_conc_charge_states()
        if not states:
            raise ValueError(
                f"DefectSpecies '{self.name}' has no variable-concentration "
                "charge states, so an effective formation energy is undefined."
            )
        return self._effective_energy(states, e_fermi, temperature)

    def _effective_energy(
        self, states: list[DefectChargeState], e_fermi: float, temperature: float
    ) -> float:
        """Boltzmann-weighted effective formation energy of a group of
        ``DefectChargeState`` objects, at a given Fermi energy and temperature.

        Args:
            states (list[DefectChargeState]): the states to combine.
            e_fermi (float): Fermi energy at which to evaluate formation
                energies.
            temperature (float): if 0, the minimum formation energy among
                ``states`` is returned (the T -> 0 limit, in which only the
                lowest-energy state is occupied). Otherwise, the
                Boltzmann-weighted effective formation energy
                ``-kT * ln(sum_i g_i * exp(-E_i / kT))`` is returned.

        Returns:
            float: the effective formation energy of ``states`` at ``e_fermi``.
        """
        energies = np.array([cs.get_formation_energy(e_fermi) for cs in states])
        if temperature == 0.0:
            return float(np.min(energies))
        kT = kboltz * temperature
        log_weights = np.array([np.log(cs.degeneracy) for cs in states]) - energies / kT
        return float(-kT * logsumexp(log_weights))

    def get_formation_energies(
        self, e_fermi: float, temperature: float = 0.0
    ) -> dict[int, float]:
        """
        Return the effective formation energy of each formal charge at a given
        Fermi energy, considering only *variable*-concentration charge states.

        Args:
            e_fermi (float): Fermi energy at which to calculate formation energies.
            temperature (float, optional): if 0 (the default), the formation
                energy of a charge is that of its lowest-energy
                ``DefectChargeState``. If > 0, charges with multiple
                ``DefectChargeState``s (metastable forms) instead use the
                Boltzmann-weighted effective formation energy of those forms.
                Defaults to 0.

        Returns:
            dict[int, float]: mapping from `DefectChargeState.charge` to its
            (effective) formation energy at e_fermi.
        """
        grouped: dict[int, list[DefectChargeState]] = {}
        for cs in self.variable_conc_charge_states():
            grouped.setdefault(cs.charge, []).append(cs)
        return {
            charge: self._effective_energy(states, e_fermi, temperature)
            for charge, states in grouped.items()
        }

    def tl_profile(
        self, efermi_min: float, efermi_max: float, temperature: float = 0.0
    ) -> np.ndarray:
        """get transition level profile for this ``DefectSpecies`` between a
        minimum and maximum Fermi energy.

        Args:
            efermi_min (float): minimum Fermi energy
            efermi_max (float): maximum Fermi energy
            temperature (float, optional): temperature used to combine
                metastable ``DefectChargeState``s sharing a formal charge, via
                ``get_formation_energies``/``get_transition_level_and_energy``.
                Defaults to 0, in which case each charge is represented by its
                lowest-energy state.

        Returns:
            np.ndarray: transition level profile between efermi_min
            and efermi_max.

        Raises:
            ValueError: if this ``DefectSpecies`` has no variable-concentration
                charge states, so a transition-level profile is undefined.
        """
        form_eners = self.get_formation_energies(efermi_min, temperature=temperature)
        if not form_eners:
            raise ValueError(
                f"DefectSpecies '{self.name}' has no variable-concentration "
                "charge states, so a transition-level profile is undefined."
            )
        q1 = min(form_eners, key=form_eners.__getitem__)
        points = [(efermi_min, form_eners[q1])]
        while q1 != min(self.charges):
            qlist = [q for q in self.charges if q < q1]
            nextp, nextq = min(
                (
                    (
                        self.get_transition_level_and_energy(
                            q1, q2, temperature=temperature
                        ),
                        q2,
                    )
                    for q2 in qlist
                ),
                key=lambda p: p[0][0],
            )
            if nextp[0] < efermi_max:
                points.append(nextp)
                q1 = nextq
            else:
                break
        form_eners_max = self.get_formation_energies(efermi_max, temperature=temperature)
        q_end = min(form_eners_max, key=form_eners_max.__getitem__)
        points.append((efermi_max, form_eners_max[q_end]))
        return np.array(points)

    def get_transition_level_and_energy(
        self, q1: int, q2: int, temperature: float = 0.0
    ) -> tuple[float, float]:
        """Calculates the Fermi energy and formation
        energy for the transition level between charge states q1 and q2.

        Args:
            q1 (int): charge on first ``DefectChargeState`` of interest
            q2 (int): charge on second ``DefectChargeState`` of interest
            temperature (float, optional): temperature used to combine
                metastable ``DefectChargeState``s sharing a formal charge, via
                ``get_formation_energies``. Defaults to 0.

        Returns:
            tuple[float, float]: Fermi energy and formation energy of the
            transition level between ``DefectChargeState`` objects with charges
            q1 and q2.
        """

        form_eners = self.get_formation_energies(0, temperature=temperature)
        trans_level = (form_eners[q2] - form_eners[q1]) / (q1 - q2)
        energy = q1 * trans_level + form_eners[q1]
        return (trans_level, energy)

    def get_concentration(self, e_fermi: float, temperature: float) -> float:
        """Determine the net concentration for this ``DefectSpecies`` at a
        specific Fermi energy and temperature.

        Args:
            e_fermi (float): fermi energy
            temperature (float): temperature

        Returns:
            float: concentration per calculation cell of this ``DefectSpecies``

        Note:
            If this ``DefectSpecies`` has a set fixed concentration, then this
            will be returned.
        """
        # Honor a forced total immediately
        if self.fixed_concentration is not None:
            return self.fixed_concentration

        # Sum over every (state, concentration) pair
        total = 0.0
        for _cs, conc in self.charge_state_concentrations(e_fermi, temperature):
            total += conc
        return total

    def fixed_conc_charge_states(self) -> list[DefectChargeState]:
        """get ``DefectChargeState`` objects of this ``DefectSpecies`` with fixed
        concentration
        (i.e those for which ``DefectChargeState.fixed_concentration != None``)

        Returns:
            list[DefectChargeState]: list of ``DefectChargeState`` objects within
            this ``DefectSpecies`` with fixed concentration
        """
        return [cs for cs in self._charge_states if cs.fixed_concentration is not None]

    def variable_conc_charge_states(self) -> list[DefectChargeState]:
        """get ``DefectChargeState`` objects in this ``DefectSpecies`` with variable
        concentration (i.e those with ``DefectChargeState.fixed_concentration == None``)

        Returns:
            list[DefectChargeState]: list of ``DefectChargeState`` objects within
            this ``DefectSpecies`` with variable concentration
        """
        return [cs for cs in self._charge_states if cs.fixed_concentration is None]

    def charge_state_concentrations(
        self, e_fermi: float, temperature: float
    ) -> list[tuple[DefectChargeState, float]]:
        """
        Compute per-cell concentrations for every charge state in this defect species.

        This method returns a list of (DefectChargeState, concentration) pairs,
        where each DefectChargeState is treated individually—so you can have
        multiple states with the same formal charge.  Concentrations are:

        - For states with `fixed_concentration` set: that value is used directly.
        - For variable states: c_cell = c_site * nsites, where
            c_site = exp(–E_f/EkBT)·degeneracy per site.

        If the species itself has `fixed_concentration` set, all variable-state
        concentrations are rescaled as a group so that

            sum(all state concentrations) == species.fixed_concentration.

        Args:
            e_fermi (float):
                The Fermi energy (relative to VBM) in electron-volts at which to
                evaluate formation energies and Boltzmann factors.
            temperature (float):
                The absolute temperature in Kelvin to use in the Boltzmann
                exponent (kB·T).

        Returns:
            list[tuple[DefectChargeState, float]]:
                A list of tuples.  Each tuple contains one DefectChargeState
                instance and its computed concentration per unit cell (float).
                All states in `self._charge_states` appear exactly once.

        Raises:
            ValueError:
                If `self.fixed_concentration` is set but the sum of all fixed-charge
                concentrations alone already exceeds it, or if there are no variable
                states left to satisfy the remaining occupancy.

        Example:
            >>> ds = DefectSpecies("V_O", nsites=12, charge_states=[...])
            >>> conc_list = ds.charge_state_concentrations(e_fermi=1.2, temperature=800)
            >>> for state, c in conc_list:
            ...     print(state.charge, c)
        """
        results: list[tuple[DefectChargeState, float]] = []
        for cs in self._charge_states:
            if cs.fixed_concentration is not None:
                c_cell = cs.fixed_concentration
            else:
                c_cell = cs.get_concentration(e_fermi, temperature) * self.nsites
            results.append((cs, c_cell))

        if self.fixed_concentration is not None:
            var_states = [
                (i, cs) for i, (cs, _) in enumerate(results) if cs.fixed_concentration is None
            ]
            fixed_conc_total = sum(c for (cs, c) in results if cs.fixed_concentration is not None)
            constrained_conc = self.fixed_concentration - fixed_conc_total
            if math.isclose(fixed_conc_total, self.fixed_concentration):
                # Fixed charge states already account for the whole species
                # total (up to floating-point noise); leave no remainder, so
                # the variable-state rescale below cannot produce a tiny
                # negative concentration.
                constrained_conc = 0.0
            elif constrained_conc < 0:
                raise ValueError(
                    f"Fixed charge state concentrations ({fixed_conc_total}) exceed "
                    f"total species concentration ({self.fixed_concentration})"
                )
            if var_states:
                # Use logsumexp for numerically stable Boltzmann proportions
                log_weights = np.array([
                    np.log(cs.degeneracy)
                    - cs.get_formation_energy(e_fermi) / (kboltz * temperature)
                    for _, cs in var_states
                ])
                log_total = logsumexp(log_weights)
                for (idx, cs), log_w in zip(var_states, log_weights, strict=True):
                    results[idx] = (cs, np.exp(log_w - log_total) * constrained_conc)
            elif constrained_conc > 0:
                raise ValueError(
                    f"Fixed charge state concentrations ({fixed_conc_total}) are "
                    f"below total species concentration ({self.fixed_concentration}) "
                    "with no variable charge state to make up the difference"
                )

        return results

    def defect_charge_contributions(
        self, e_fermi: float, temperature: float
    ) -> tuple[float, float]:
        """
        Calculate the defect charge contributions to the total charge of this
        ``DefectSpecies`` at a given Fermi energy and temperature.

        Args:
            e_fermi (float): Fermi energy.
            temperature (float): temperature

        Returns:
            tuple[float, float]: charge contributions of the
            ``DefectChargeState`` objects that comprise this ``DefectSpecies``
            at the given Fermi energy and temperature.
        """

        lhs = rhs = 0.0
        # charge_state_concentrations now returns list[tuple[DefectChargeState, float]]
        for cs, conc in self.charge_state_concentrations(e_fermi, temperature):
            q = cs.charge
            if q > 0:
                lhs += conc * q
            elif q < 0:
                rhs += conc * abs(q)
        return lhs, rhs
