import numpy as np
from typing import List, Dict, Tuple, Optional
from py_sc_fermi.defect_charge_state import DefectChargeState


class DefectSpecies(object):
    """Class for individual defect species.

    Args:
        name (str): A unique identifying string for this defect species,
           e.g. ``"V_O"`` might be used for an oxygen vacancy.
        nsites (int): Number of sites energetically degenerate sites where this
         defect can form in the unit cell (the site degeneracy).
        charge_states (Dict[int, DefectChargeState]): A dictionary of
           ``DefectChargeState`` with their charge as the key, i.e.
           {charge : ``DefectChargeState``}

    """

    def __init__(
        self,
        name: str,
        nsites: int,
        charge_states: List[DefectChargeState],
        fixed_concentration: Optional[float] = None,
    ):
        """Instantiate a DefectSpecies object."""

        self._name = name
        self._nsites = nsites
        self._charge_states = charge_states

        if isinstance(charge_states, dict):
             self._charge_states = list(charge_states.values())
        elif isinstance(charge_states, list):
             self._charge_states = charge_states
        else:
            raise TypeError("charge_states must be a dict or a list")
        
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
            int: site degeneracy fot ``DefectSpecies``
        """
        return self._nsites

    @property
    def charge_states(
        self,
    ) -> List[DefectChargeState]:
        """

        Returns:
            List[DefectChargeState]: list of ``DefectChargeState`` objects that
            comprise this ``DefectSpecies``
        """
        return self._charge_states

    @property
    def charges(self) -> List[int]:
        """list of all the charges of the ``DefectChargeState`` objects that
        comprise this ``DefectSpecies``

        Returns:
            List[int]: list of charge states of this ``DefectSpecies``
        """
        return [cs.charge for cs in self._charge_states]

    @property
    def fixed_concentration(self) -> Optional[float]:
        """fixed concentration of the ``DefectSpecies``. ``None`` if the
        concentration of this defect is variable.

        Returns:
            Optional[float]: fixed concentration per unit cell of the
            ``DefectSpecies``
        """
        return self._fixed_concentration

    def __repr__(self):
        to_return = f"\n{self.name}, nsites={self.nsites}"
        if self.fixed_concentration is not None:
            to_return += f"\nfixed [c] = {self.fixed_concentration}"
        to_return += "\n" + "".join(
            [f"  {cs.__repr__()}\n" for cs in self.charge_states.values()]
        )
        return to_return

    @classmethod
    def from_dict(cls, d: dict) -> "DefectSpecies":
        """return a ``DefectSpecies`` object from a dictionary containing the defect
        species data. Primarily for use defining a full ``DefectSystem`` from a
        .yaml file.

        Args:
            defect_species_dict (dict): dictionary containing the defect species
               data.

        Raises:
            ValueError: if any of the ``DefectChargeState`` objects specified have no
               fixed concentration and no formation energy

        Returns:
            DefectChargeState: as specified by the provided dictionary
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
    ) -> List[DefectChargeState]:
        """
        Return all *variable* DefectChargeState objects sorted by formation
        energy at a given Fermi energy.

        Args:
            e_fermi (float): Fermi energy at which to evaluate formation energies.

        Returns:
            List[DefectChargeState]: all variable charge‐states of this species,
            sorted from lowest to highest formation energy at e_fermi.
        """
        # variable_conc_charge_states() now returns List[DefectChargeState]
        return sorted(
            self.variable_conc_charge_states(),
            key=lambda st: st.get_formation_energy(e_fermi),
        )

    def as_dict(self) -> dict:
        """get representation of ``DefectSpecies`` as a dictionary

        Returns:
            dict: dictionary representation of ``DefectChargeState``
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
        """
        return self.charge_states_by_formation_energy(e_fermi)[0]

    def get_formation_energies(self, e_fermi: float) -> Dict[int, float]:
        """
        Return the formation energy of each *variable* charge‐state at a given
        Fermi energy.

        Args:
            e_fermi (float): Fermi energy at which to calculate formation energies.

        Returns:
            Dict[int, float]: mapping from `DefectChargeState.charge` to its
            formation energy at e_fermi.
        """
        # variable_conc_charge_states() now returns List[DefectChargeState]
        return {
            cs.charge: cs.get_formation_energy(e_fermi)
            for cs in self.variable_conc_charge_states()
        }

    def tl_profile(self, efermi_min: float, efermi_max: float) -> np.ndarray:
        """get transition level profile for this ``DefectSpecies`` between a
        minimum and maximum Fermi energy.

        Args:
            efermi_min (float): minimum Fermi energy
            efermi_max (float): maximum Fermi energy

        Returns:
            np.ndarray: transition level profile between efermi_min
            and efermi_max.
        """
        cs = self.min_energy_charge_state(efermi_min)
        q1 = cs.charge
        form_e = cs.get_formation_energy(efermi_min)
        points = [(efermi_min, form_e)]
        while q1 != min(self.charges):
            qlist = [q for q in self.charges if q < q1]
            nextp, nextq = min(
                ((self.get_transition_level_and_energy(q1, q2), q2) for q2 in qlist),
                key=lambda p: p[0][0],
            )
            if nextp[0] < efermi_max:
                points.append(nextp)
                q1 = nextq
            else:
                break
        cs = self.min_energy_charge_state(efermi_max)
        form_e = cs.get_formation_energy(efermi_max)
        points.append((efermi_max, form_e))
        return np.array(points)

    def get_transition_level_and_energy(self, q1: int, q2: int) -> Tuple[float, float]:
        """Calculates the Fermi energy and formation
        energy for the transition level between charge states q1 and q2.

        Args:
            q1 (int): charge on first ``DefectChargeState`` of interest
            q2 (int): charge on second ``DefectChargeState`` of interest

        Returns:
            Tuple[float, float]: Fermi energy and formation energy of the
            transition level between ``DefectChargeState`` objects with charges
            q1 and q2.
        """

        form_eners = self.get_formation_energies(0)
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
        for cs, conc in self.charge_state_concentrations(e_fermi, temperature):
            total += conc
        return total

    def fixed_conc_charge_states(self) -> List[DefectChargeState]:
        """get ``DefectChargeState`` objects of this ``DefectSpecies`` with fixed
        concentration
        (i.e those for which ``DefectChargeState.fixed_concentration != None``)

        Returns:
            List[DefectChargeState]: list of ``DefectChargeState`` objects within
            this ``DefectSpecies`` with fixed concentration
        """
        return [cs for cs in self._charge_states if cs.fixed_concentration is not None]

    def variable_conc_charge_states(self) -> List[DefectChargeState]:
        """get ``DefectChargeState`` objects in this ``DefectSpecies`` with variable
        concentration (i.e those with ``DefectChargeState.fixed_concentration == None``)

        Returns:
            List[DefectChargeState]: list of ``DefectChargeState`` objects within
            this ``DefectSpecies`` with variable concentration
        """
        return [cs for cs in self._charge_states if cs.fixed_concentration is None]

    def charge_state_concentrations(
        self, e_fermi: float, temperature: float
    ) -> List[Tuple[DefectChargeState, float]]:
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
            List[Tuple[DefectChargeState, float]]:
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
        results: List[Tuple[DefectChargeState, float]] = []
        for cs in self._charge_states:
            c_site = cs.get_concentration(e_fermi, temperature)
            if cs.fixed_concentration is not None:
                c_cell = cs.fixed_concentration
            else:
                c_cell = c_site * self.nsites
            results.append((cs, c_cell))

        if self.fixed_concentration is not None:
            fixed_sum = sum(c for (cs, c) in results if cs.fixed_concentration is not None)
            var_sum   = sum(c for (cs, c) in results if cs.fixed_concentration is None)
            to_alloc  = self.fixed_concentration - fixed_sum

            if var_sum > 0:
                scale = to_alloc / var_sum
                for idx, (cs, c) in enumerate(results):
                    if cs.fixed_concentration is None:
                        results[idx] = (cs, c * scale)
            else:
                if abs(to_alloc) > 0:
                    raise ValueError(
                        f"{self.name}: fixed_concentration {self.fixed_concentration} "
                        f"cannot be satisfied (fixed states sum to {fixed_sum})."
                    )

        return results

    def defect_charge_contributions(
        self, e_fermi: float, temperature: float
    ) -> Tuple[float, float]:
        """
        Calculate the defect charge contributions to the total charge of this
        ``DefectSpecies`` at a given Fermi energy and temperature.

        Args:
            e_fermi (float): Fermi energy.
            temperature (float): temperature

        Returns:
            Tuple[float, float]: charge contributions of the
            ``DefectChargeState`` objects that comprise this ``DefectSpecies``
            at the given Fermi energy and temperature.
        """

        lhs = rhs = 0.0
        # charge_state_concentrations now returns List[Tuple[DefectChargeState, float]]
        for cs, conc in self.charge_state_concentrations(e_fermi, temperature):
            q = cs.charge
            if q > 0:
                lhs += conc * q
            elif q < 0:
                rhs += conc * abs(q)
        return lhs, rhs
