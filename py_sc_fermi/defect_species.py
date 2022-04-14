import numpy as np


class DefectSpecies(object):
    """Class for individual defect species.

    Attribute:
        name (str): An identifying string for this defect species, e.g. `V_O`.
        nsites (int): The number of sites in the unit cell where this defect may reside.
        charge_states (dict): A dictionary of charge states. Each key-value pair
            is charge (int): charge-state (`DefectChargeState` or
            `FrozenDefectChargeState`).
        fixed_concentration (optional, :obj:`float`): If set, this fixes the net
            concentration of this defect species per unit cell. If not set, the
            concentration of this defect species will depend on the system conditions.

    """

    def __init__(
        self,
        name: str,
        nsites: int,
        charge_states: list,
        fixed_concentration: float = None,
    ):
        """Instantiate a DefectSpecies object.

        Args:
            name (str): An identifying string for this defect species, e.g. `V_O`.
            nsites (int): Number of sites in the unit cell where this defect may exist.
            charge_states (list of :obj:`DefectChargeState`
                                or :obj:`FrozenDefectChargeState` objects): A list of
                          defect charge state objects (these may be fixed-concentration
                          charge states)

        Returns:
            None

        """

        self._name = name
        self._nsites = nsites
        self._charge_states = {cs.charge: cs for cs in charge_states}
        self._fixed_concentration = fixed_concentration

    def fix_concentration(self, concentration):
        """fixed the net concentration (per unit cell) of this defect species"""
        self._fixed_concentration = concentration

    @property
    def name(self) -> str:
        """The identifying name for this defect species."""
        return self._name

    @property
    def nsites(self) -> int:
        """The number of sites per unit cell that this defect species may occupy."""
        return self._nsites

    @property
    def charge_states(
        self,
    ) -> dict[int, "py_sc_fermi.defect_charge_state.DefectChargeState"]:
        """The charge states for this defect species as s dictionary of
        charge: defect_state key-value pairs"""
        return self._charge_states

    @property
    def fixed_concentration(self) -> float:
        """The fixed net concentration (per unit cell) of this defect species,
        or `None` if the defect concentrations are free to change"""
        return self._fixed_concentration

    def __repr__(self):
        to_return = f"\n{self.name}, nsites={self.nsites}"
        if self.fixed_concentration != None:
            to_return += f"\nfixed [c] = {self.fixed_concentration}"
        to_return += "\n" + "".join(
            [f"  {cs.__repr__()}\n" for cs in self.charge_states.values()]
        )
        return to_return

    def charge_states_by_formation_energy(
        self, e_fermi: float
    ) -> list["py_sc_fermi.defect_charge_state.DefectChargeState"]:
        """Returns a list of defect charge states, sorted by increasing formation energy
        at E_Fermi (relative to E(VBM)).

        Args:
            e_fermi (float): Fermi energy relative to E(VBM) (in eV).

        Returns:
            (list(DefectChargeState)): Ordered list of all defect charge states.

        Note:
            Defect charge states with fixed concentrations do not have a defined
            formation energy, and are not returned in the sorted list.

        """
        return sorted(
            self.variable_conc_charge_states().values(),
            key=lambda x: x.get_formation_energy(e_fermi),
        )

    def min_energy_charge_state(
        self, e_fermi: float
    ) -> "py_sc_fermi.defect_charge_state.DefectChargeState":
        """Returns the defect charge state with the minimum energy at E_Fermi
        (relative to E(VBM)).

        Args:
            e_fermi (float): Fermi energy relative to e(VBM) (in eV).

        Returns:
            DefectChargeState: The defect charge state with the lowest formation energy.

        """
        return self.charge_states_by_formation_energy(e_fermi)[0]

    def get_formation_energies(self, e_fermi: float) -> dict[int, float]:
        """Returns a dictionary of formation energies for all charge states. Formation
        energies are calculated at E_Fermi relative to E_VBM.

        Args:
            e_fermi (float): Fermi energy relative to E(VBM) (in eV).

        Returns:
            (dict(int:float)): Dictionary of {charge: formation_energy} pairs.

        Notes:
            Fixed-concentration defect charge states do not have a defined formation
            energy and are not returned.

        """
        form_eners = {}
        for q, cs in self.variable_conc_charge_states().items():
            form_eners[q] = cs.get_formation_energy(e_fermi)
        return form_eners

    def tl_profile(self, efermi_min: float, efermi_max: float) -> "np.array":
        """
        get transition level profile between a range of specified
        Fermi energies

        Args:
            efermi_min (float): minimum Fermi energy for range of interest
            efermi_max (float): maximum Fermi energy for range of interest

        Returns:
            points (numpy.array): transtion levels
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

    def get_transition_level_and_energy(self, q1: int, q2: int) -> tuple[float, float]:
        """
        Get the transition level between two charge states.

        Calculates the Fermi energy (relative to the host VBM) and formation
        energy for the transition level between charge states q1 and q2.

        Args:
            q1 (int): Charge of charge state 1.
            q2 (int): Charge of charge state 2.

        Returns:
            (tuple): Fermi energy (relative to the host VBM) and
            formation energy of the transition level as a tuple of floats
            ``(trans_level, energy)``.
        """

        form_eners = self.get_formation_energies(0)
        trans_level = (form_eners[q2] - form_eners[q1]) / (q1 - q2)
        energy = q1 * trans_level + form_eners[q1]
        return (trans_level, energy)

    @property
    def charges(self) -> list[int]:
        """
        Return a list of charges for the charge states of this defect.

        Args:
            None

        Returns:
            (list(int)): A list of integer charges.

        """
        return self.charge_states.keys()

    def get_concentration(self, e_fermi: float, temperature: float) -> float:
        """
        Get the net concentration per unit cell for this defect species at a specific
        Fermi energy (relative to E(VBM)) and temperature.

        Args:
            e_fermi (float): Fermi energy, relative to E(VBM) (in eV).
            temperature (float): Temperature (in K).

        Returns:
            (float): The concentration (per unit cell), summed over all charge states.

        Note:
            If this DefectSpecies has a fixed concentration, then this
            will be returned.

        """
        if self.fixed_concentration:
            return self.fixed_concentration
        else:
            return sum(self.charge_state_concentrations(e_fermi, temperature).values())

    def fixed_conc_charge_states(
        self,
    ) -> dict[int, "py_sc_fermi.defect_charge_state.DefectChargeState"]:
        """return dict of fixed concentration charge states for this DefectSpecies"""
        return {
            q: cs
            for q, cs in self.charge_states.items()
            if cs.fixed_concentration != None
        }

    def variable_conc_charge_states(
        self,
    ) -> dict[int, "py_sc_fermi.defect_charge_state.DefectChargeState"]:
        """return dict of charge states with no fixed concentration for this defect species"""
        return {
            q: cs
            for q, cs in self.charge_states.items()
            if cs.fixed_concentration == None
        }

    def charge_state_concentrations(
        self, e_fermi: float, temperature: float
    ) -> dict[int, float]:
        """
        At a given Fermi energy and temperature, calculate the concentrations
        of the different charge states of this defect species

        Args:
            e_fermi (float): the Fermi energy (relative to the VBM) (in eV)
            temperature (flot): Temperature (in K)

        Returns:
            cs_concentrations (Dict): a dictionary of charge state (int)
            concentration (float) pairs.
        """
        
        var_concs = self.variable_conc_charge_states()
        fixed_concs = self.fixed_conc_charge_states()

        cs_concentrations = {cs.charge : cs.get_concentration(e_fermi, temperature) * self.nsites for cs in var_concs.values()}
        for q, cs in fixed_concs.items():
            cs_concentrations[q] = cs.get_concentration(e_fermi, temperature)
        
        if self.fixed_concentration != None:
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
    ) -> tuple[float, float]:
        """get total charge contribution for this DefectSpecies in all charge states"""
        lhs = 0.0
        rhs = 0.0
        for q, concd in self.charge_state_concentrations(e_fermi, temperature).items():
            if q < 0:
                rhs += concd * abs(q)
            if q > 0:
                lhs += concd * abs(q)
        return lhs, rhs
