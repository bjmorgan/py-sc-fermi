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
        charge_states: Dict[int, DefectChargeState],
        fixed_concentration: float = None,
    ):
        """Instantiate a DefectSpecies object."""

        self._name = name
        self._nsites = nsites
        self._charge_states = charge_states
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
    ) -> Dict[int, DefectChargeState]:
        """

        Returns:
            Dict[int, DefectChargeState]: The charge states of this defect species as s dictionary of
            ``{charge (int): DefectChargeState}`` key-value pairs"""
        return self._charge_states

    @property
    def charges(self) -> List[int]:
        """list of all the charges of the ``DefectChargeState`` objects that comprise
        this ``DefectSpecies``

        Returns:
            List[int]: list of charge states of this ``DefectSpecies``
        """
        return list(self.charge_states.keys())

    @property
    def fixed_concentration(self) -> Optional[float]:
        """fixed concentration of the ``DefectSpecies``. ``None`` if the
        concentration of this defect is variable.

        Returns:
            Optional[float]: fixed concentration per unit cell of the ``DefectSpecies``
        """
        return self._fixed_concentration

    def __repr__(self):
        to_return = f"\n{self.name}, nsites={self.nsites}"
        if self.fixed_concentration != None:
            to_return += f"\nfixed [c] = {self.fixed_concentration}"
        to_return += "\n" + "".join(
            [f"  {cs.__repr__()}\n" for cs in self.charge_states.values()]
        )
        return to_return

    @classmethod
    def from_dict(cls, defect_species_dict: dict, volume: Optional[float] = None):
        """return a ``DefectSpecies`` object from a dictionary containing the defect
        species data. Primarily for use defining a full ``DefectSystem`` from a
        .yaml file.

        Args:
            defect_species_dict (dict): dictionary containing the defect species
               data.
            volume (Optional[float], optional): volume of the unit cell.
               Defaults to ``None``.

        Raises:
            ValueError: if any of the ``DefectChargeState`` objects specified have no
               fixed concentration and no formation energy

        Returns:
            DefectChargeState: as specified by the provided dictionary
        """
        charge_states = []
        name = list(defect_species_dict.keys())[0]
        for n, c in defect_species_dict[name]["charge_states"].items():
            if "fixed_concentration" not in list(c.keys()):
                fixed_concentration = None
            elif volume is not None:
                fixed_concentration = float(c["fixed_concentration"]) / 1e24 * volume
            if "formation_energy" not in list(c.keys()):
                formation_energy = None
            else:
                formation_energy = float(c["formation_energy"])
            if formation_energy == None and fixed_concentration == None:
                raise ValueError(
                    f"{name, n} must have one or both fixed concentration or formation energy"
                )
            charge_state = DefectChargeState(
                charge=n,
                energy=formation_energy,
                degeneracy=c["degeneracy"],
                fixed_concentration=fixed_concentration,
            )
            charge_states.append(charge_state)

        if (
            "fixed_concentration" in defect_species_dict[name].keys()
            and volume is not None
        ):
            fixed_concentration = (
                float(defect_species_dict[name]["fixed_concentration"]) / 1e24 * volume
            )
            return cls(
                name,
                defect_species_dict[name]["nsites"],
                charge_states={cs.charge: cs for cs in charge_states},
                fixed_concentration=fixed_concentration,
            )
        else:
            return cls(
                name,
                defect_species_dict[name]["nsites"],
                charge_states={cs.charge: cs for cs in charge_states},
            )

    @classmethod
    def _from_list_of_strings(cls, defect_string: List[str]):
        """generate a ``DefectSpecies`` object from a string containing the defect
        species data. Only intended for use reading defect species from a
        SC-Fermi input file.

        Args:
            defect_string (List[str]): list of strings describing the
            ``DefectSpecies``

        Returns:
            DefectSpecies: returns a ``DefectSpecies`` object as defined by
            the input list of strings
        """
        defect_species = defect_string.pop(0).split()
        name = defect_species[0]
        n_charge_states = int(defect_species[1])
        nsites = int(defect_species[2])
        charge_states = []
        for i in range(n_charge_states):
            string = defect_string.pop(0)
            charge_state = DefectChargeState.from_string(string)
            charge_states.append(charge_state)
        return cls(name, nsites, {cs.charge: cs for cs in charge_states})

    def charge_states_by_formation_energy(
        self, e_fermi: float
    ) -> List[DefectChargeState]:
        """Returns a list of charge states sorted by formation energy at a given
        Fermi energy.

        Args:
            e_fermi (float): Fermi energy

        Returns:
            List[DefectChargeState]: list of ``DefectChargeState`` objects in this
            ``DefectSpecies`` sorted by formation energy at ``e_fermi``.

        Note:
            ``DefectChargeState`` objects with fixed-concentration are not
            included in the returned list, even if their formation energy
            is specified.
        """
        return sorted(
            self.variable_conc_charge_states().values(),
            key=lambda x: x.get_formation_energy(e_fermi),
        )

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
        """Returns a dictionary of formation energies for all
        ``DefectChargeState`` objects in this ``DefectSpecies``.
        Formation energies are calculated at E[Fermi] relative to E[VBM].

        Args:
            e_fermi (float): fermi energy

        Returns:
            Dict[int, float]: key-value pairs of charge on ``DefectChargeState``
            objects that comprise this ``DefectSpecies`` and their formation energy,
            i.e ``{DefectChargeState.charge : formation_energy}``

        Note:
            Fixed-concentration ``DefectChargeState`` objects are not included
            in the returned list, even if their formation energy is specified.
        """
        form_eners = {}
        for q, cs in self.variable_conc_charge_states().items():
            form_eners[q] = cs.get_formation_energy(e_fermi)
        return form_eners

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
        if self.fixed_concentration:
            return self.fixed_concentration
        else:
            return sum(self.charge_state_concentrations(e_fermi, temperature).values())

    def fixed_conc_charge_states(
        self,
    ) -> Dict[int, DefectChargeState]:
        """get ``DefectChargeState`` objects of this ``DefectSpecies`` with fixed
        concentration (i.e those for which ``DefectChargeState.fixed_concentration != None``)

        Returns:
            Dict[int, DefectChargeState]: key-value pairs of charge on fixed
            concentration ``DefectChargeState`` objects, and the charge state
            which is variable, i.e. ``{DefectChargeState.charge : DefectChargeState}``
        """
        return {
            q: cs
            for q, cs in self.charge_states.items()
            if cs.fixed_concentration is not None
        }

    def variable_conc_charge_states(self) -> Dict[int, DefectChargeState]:
        """get ``DefectChargeState`` objects in this ``DefectSpecies`` with variable
        concentration (i.e those with ``DefectChargeState.fixed_concentration == None``)

        Returns:
            Dict[int, DefectChargeState]: key-value pairs of charge on variable
            concentration ``DefectChargeState`` objects within this ``DefectSpecies``,
            and the charge state which is variable, i.e.
            ``{DefectChargeState.charge : DefectChargeState}``
        """
        return {
            q: cs
            for q, cs in self.charge_states.items()
            if cs.fixed_concentration is None
        }

    def charge_state_concentrations(
        self, e_fermi: float, temperature: float
    ) -> Dict[int, float]:
        """at a given Fermi energy and temperature, calculate the concentrations
        of the different ``DefectChargeStates`` of this ``DefectSpecies``.

        Args:
            e_fermi (float): Fermi energy
            temperature (float): temperature

        Returns:
            Dict[int, float]: key-value pairs of charge of each
            ``DefectChargeState`` and the concentration of the
            ``DefectChargeState`` with that charge, i.e.
            {``DefectChargeState.charge``: concentration}
        """

        var_concs = self.variable_conc_charge_states()
        fixed_concs = self.fixed_conc_charge_states()

        cs_concentrations = {
            cs.charge: cs.get_concentration(e_fermi, temperature) * self.nsites
            for cs in var_concs.values()
        }
        for q, cs in fixed_concs.items():
            cs_concentrations[q] = cs.get_concentration(e_fermi, temperature)

        if self.fixed_concentration is not None:
            fixed_conc_chg_states = sum(
                [
                    c
                    for q, c in cs_concentrations.items()
                    if q in self.fixed_conc_charge_states()
                ]
            )
            variable_conc_chg_states = sum(
                [
                    c
                    for q, c in cs_concentrations.items()
                    if q in self.variable_conc_charge_states()
                ]
            )
            constrained_conc = self.fixed_concentration - fixed_conc_chg_states
            scaling = constrained_conc / variable_conc_chg_states
            for q in cs_concentrations:
                if q in self.variable_conc_charge_states():
                    cs_concentrations[q] *= scaling
        return cs_concentrations

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

        lhs = 0.0
        rhs = 0.0
        for q, concd in self.charge_state_concentrations(e_fermi, temperature).items():
            if q < 0:
                rhs += concd * abs(q)
            if q > 0:
                lhs += concd * abs(q)
        return lhs, rhs
