import numpy as np
from typing import List, Dict, Tuple, Optional
from py_sc_fermi.defect_charge_state import DefectChargeState


class DefectSpecies(object):
    """Class for individual defect species.


    :param str name: A unique identifying string for this defect species,
        e.g. `V_O` might be an oxygen vacancy.
    :param int nsites: Number of sites energetically degenerate sites in the
        calculation cell where this defect can form.
    :param dict charge_states: A dictionary of :py:class:`DefectChargeState`
        for this defect. Each key-value pair takes the form
        ``{charge (int): :py:class:`DefectChargeState``
    :param float fixed_concentration: If set, this fixes the total concentration
        of this defect species to this value (per calculation cell).
        The concetration of the charge states will be able to change,
        but not the total concentration of the defect species.
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
        """fixed the net concentration (per calculation cell) of this defect
        species to the specified value.

        :param float concentration: Fixed concentration to be enforced
        (per calculation cell)
        """
        self._fixed_concentration = concentration

    @property
    def name(self) -> str:
        """:return: The identifying string for this defect species."""
        return self._name

    @property
    def nsites(self) -> int:
        """:return: The number of sites energetically degenerate sites for this
        defect species."""
        return self._nsites

    @property
    def charge_states(
        self,
    ) -> Dict[int, DefectChargeState]:
        """:return: The charge states of this defect species as s dictionary of
        ``{charge (int): :py:class:`DefectChargeState`}`` key-value pairs"""
        return self._charge_states

    @property
    def charges(self) -> List[int]:
        """
        :return: A list of all the charges of this defect species.
        """
        return list(self.charge_states.keys())

    @property
    def fixed_concentration(self) -> Optional[float]:
        """:return: The fixed concentration (per unit cell) of this defect
        species, or `None` if the defect concentrations are free to change"""
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
    def from_dict(
        cls, defect_species_dict: dict, volume: Optional[float] = None):
        """
        return a DefectSpecies object from a dictionary containing the defect
        species data.

        :param dict defect_species_dict: dictionary containing the defect species
            data.
        :param float volume: volume of the defect system in Angstrom^3.
        :return: :py:class:`DefectSpecies`
        :rtype: py_sc_fermi.defect_species.DefectSpecies
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

        if "fixed_concentration" in defect_species_dict[name].keys() and volume is not None:
            fixed_concentration = (
                float(defect_species_dict[name]["fixed_concentration"]) / 1e24 * volume
            )
            return cls(
                name,
                defect_species_dict[name]["nsites"],
                charge_states={cs.charge : cs for cs in charge_states},
                fixed_concentration=fixed_concentration,
            )
        else:
            return cls(
                name,
                defect_species_dict[name]["nsites"],
                charge_states={cs.charge : cs for cs in charge_states},
            )

    @classmethod
    def _from_list_of_strings(
        cls, defect_string: List[str]
    ):
        """ 
        generate a DefectSpecies object from a string containing the defect
        species data. Only intended for use reading defect species from a
        SC-Fermi input file.

        :param str defect_string: string containing the defect species data.
        :return: :py:class:`DefectSpecies`
        :rtype: py_sc_fermi.defect_species.DefectSpecies
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
        return cls(name, nsites, {cs.charge : cs for cs in charge_states})

    def charge_states_by_formation_energy(
        self, e_fermi: float
    ) -> List[DefectChargeState]:
        """
        Returns a list of charge states sorted by formation energy at a given
        Fermi energy.

        :param float e_fermi: Fermi energy relative to E(VBM) (in eV).
        :return: A list of :py:class:`DefectChargeState` sorted by formation energy.
        :rtype: list[:py:class:`DefectChargeState`]

        .. note::
            Fixed-concentration defect charge states are not included in the
            returned list, even if their formation energy is specified.

        """
        return sorted(
            self.variable_conc_charge_states().values(),
            key=lambda x: x.get_formation_energy(e_fermi),
        )

    def min_energy_charge_state(
        self, e_fermi: float
    ) -> DefectChargeState:
        """Returns the defect charge state with the minimum energy at a given
        Fermi energy.

        :param float e_fermi: Fermi energy relative to E(VBM) (in eV).
        :return: The :py:class:`DefectChargeState` with the minimum energy at
            e_fermi.
        :rtype: :py:class:`DefectChargeState`
        """
        return self.charge_states_by_formation_energy(e_fermi)[0]

    def get_formation_energies(self, e_fermi: float) -> Dict[int, float]:
        """Returns a dictionary of formation energies for all charge states.
        Formation energies are calculated at E_Fermi relative to E_VBM.

        :param float e_fermi: Fermi energy relative to E(VBM) (in eV).
        :return: A dictionary of formation energies for all charge states.
        :rtype: dict[int, float]
            e_fermi (float): Fermi energy relative to E(VBM) (in eV).

        .. note::
            Fixed-concentration defect charge states are not included in the
            returned list, even if their formation energy is specified.
        """
        form_eners = {}
        for q, cs in self.variable_conc_charge_states().items():
            form_eners[q] = cs.get_formation_energy(e_fermi)
        return form_eners

    def tl_profile(self, efermi_min: float, efermi_max: float) -> np.ndarray:
        """
        get transition level profile between a minimum and maximum Fermi energy.

        :param float efermi_min: Minimum Fermi energy (in eV).
        :param float efermi_max: Maximum Fermi energy (in eV).
        :return: Transition level profile between efermi_min and efermi_max.
        :rtype: np.array
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
        """
        Calculates the Fermi energy (relative to the host VBM) and formation
        energy for the transition level between charge states q1 and q2.

        :param int q1: first charge state of interest.
        :param int q2: second charge state of interest.
        :return: Fermi energy (in eV) and formation energy (in eV) of the
            transition level between charge states q1 and q2.
        :rtype: Tuple[float, float]
        """

        form_eners = self.get_formation_energies(0)
        trans_level = (form_eners[q2] - form_eners[q1]) / (q1 - q2)
        energy = q1 * trans_level + form_eners[q1]
        return (trans_level, energy)

    def get_concentration(self, e_fermi: float, temperature: float) -> float:
        """
        Determine the net concentration (per calculation cell for this defect
        species at a specific Fermi energy (in eV relative to E(VBM)) and
        temperature.

        :param float e_fermi: Fermi energy relative to E(VBM) (in eV).
        :param float temperature: Temperature (in K).
        :return: concentration
        :rtype: float

        .. note:
            If this ``DefectSpecies`` has a fixed concentration, then this
            will be returned.

        """
        if self.fixed_concentration:
            return self.fixed_concentration
        else:
            return sum(self.charge_state_concentrations(e_fermi, temperature).values())

    def fixed_conc_charge_states(
        self,
    ) -> Dict[int, DefectChargeState]:
        """:return: A dictionary of fixed-concentration defect charge states
        of this ``DefectSpecies.``. Each key-value pair is of the form
        ``{charge_state: DefectChargeState}``."""
        return {
            q: cs
            for q, cs in self.charge_states.items()
            if cs.fixed_concentration is not None
        }

    def variable_conc_charge_states(
        self,
    ) -> Dict[int, DefectChargeState]:
        """:return: A dictionary of variable-concentration defect charge states
        of this ``DefectSpecies.``. Each key-value pair is of the form
        ``{charge_state: DefectChargeState}``."""
        return {
            q: cs
            for q, cs in self.charge_states.items()
            if cs.fixed_concentration is None
        }

    def charge_state_concentrations(
        self, e_fermi: float, temperature: float
    ) -> Dict[int, float]:
        """
        At a given Fermi energy and temperature, calculate the concentrations
        of the different charge states of ``DefectSpecies``.

        :param float e_fermi: Fermi energy relative to E(VBM) (in eV).
        :param float temperature: Temperature (in K).
        :return: A dictionary of concentrations of the different charge states
            of ``DefectSpecies``. Each key-value pair is of the form
            ``{charge state (int): concentration (float)}``.
        :rtype: Dict[int, float]
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
        Calculate the defect charge contributions to the total charge of the
        defect at a given Fermi energy (in eV) and temperature.

        :param float e_fermi: Fermi energy relative to E(VBM) (in eV).
        :param float temperature: Temperature (in K).
        :return: The total charge of the defect (in e) and the charge
            contribution of each charge state (in e).
        :rtype: Tuple[float, float]
        """
        lhs = 0.0
        rhs = 0.0
        for q, concd in self.charge_state_concentrations(e_fermi, temperature).items():
            if q < 0:
                rhs += concd * abs(q)
            if q > 0:
                lhs += concd * abs(q)
        return lhs, rhs
