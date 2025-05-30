from typing import Dict, List, Tuple, Any, Optional
from py_sc_fermi.dos import DOS
from py_sc_fermi.defect_species import DefectSpecies
from py_sc_fermi.defect_charge_state import DefectChargeState
from py_sc_fermi.inputs import InputSet
import warnings



class CustomWarningManager:
    def __init__(self):
        self.dos_overflow_warning_issued = False
        self.defect_overflow_warning_issued = False

    def custom_warning(self, message, category, filename, lineno, file=None, line=None):
        if category == RuntimeWarning:
            if "dos" in str(filename) and "overflow" in str(message):
                if not self.dos_overflow_warning_issued:
                    print(
                        """DOSOverflowWarning: An overflow occurred during computation of
                        electron and hole concentrations. This is likely a natural result of the use of
                        a numerical solver for the Fermi energy search. This can likely be ignored
                        though you should always check the final results are reasonable."""
                    )
                    self.dos_overflow_warning_issued = True
            elif "defect" in str(filename) and "overflow" in str(message):
                if not self.defect_overflow_warning_issued:
                    print(
                        """DefectOverflowWarning: An overflow occurred during computation of
                        defect concentrations. This is likely a natural result of the use of
                        a numerical solver for the Fermi energy search. This can likely be ignored
                        though you should always check the final results are reasonable."""
                    )
                    self.defect_overflow_warning_issued = True
            else:
                print(f"RuntimeWarning: {message}")


# Create a CustomWarningManager and set the custom_warning method as the warning handler
warning_manager = CustomWarningManager()
warnings.showwarning = warning_manager.custom_warning


class DefectSystem(object):
    """
    Compute self-consistent Fermi energy, carrier and defect concentrations
    in a material, with optional site- and element-level occupancy constraints.

    Args:
        defect_species (List[DefectSpecies]): 
            List of species in the system.
        dos (DOS): 
            Density of states object for the host.
        volume (float): 
            Unit cell volume (Å³).
        temperature (float): 
            Absolute temperature (K).
        site_pools (Dict[str,(float,List[Union[str,DefectSpecies]]]], optional):
            Maps pool names to (number of sites per cell, list of species sharing
            those sites).  If None, no site competition is applied.
        element_pools (Dict[str,(float,List[Union[str,DefectSpecies]]]], optional):
            Maps element symbols to (number of dopant atoms per cell, 
            list of (species, stoich) pairs).  If None, no elemental constraint
            is applied.
        convergence_tolerance (float, optional): 
            Charge-neutrality tolerance for the solver.
        n_trial_steps (int, optional): 
            Maximum iterations for the solver.
    """

    def __init__(
        self,
        defect_species: List[DefectSpecies],
        dos: DOS,
        volume: float,
        temperature: float,
        site_pools: Optional[Dict[str, Tuple[float, List[DefectSpecies]]]] = None,
        element_pools: Optional[Dict[str, Tuple[float, List[DefectSpecies]]]] = None,
        convergence_tolerance: float = 1e-18,
        n_trial_steps: int = 1500,
    ):
        self.defect_species = defect_species
        self.volume = volume
        self.dos = dos
        self.temperature = temperature
        self.element_pools = element_pools or {}
        self.site_pools = site_pools or {}
        self.convergence_tolerance = convergence_tolerance
        self.n_trial_steps = n_trial_steps

    def __repr__(self):
        to_return = [
            f"DefectSystem\n",
            f"  nelect: {self.dos.nelect} e\n",
            f"  bandgap: {self.dos.bandgap} eV\n",
            f"  volume: {self.volume} A^3\n",
            f"  temperature: {self.temperature} K\n",
            f"\nContains defect species:\n",
        ]
        for ds in self.defect_species:
            to_return.append(str(ds))
        return "".join(to_return)

    @property
    def defect_species_names(self) -> List[str]:
        """list of the names of all ``DefectSpecies`` considered in the
        ``DefectSystem``.

        Returns:
            List[str]: list of names of ``DefectSpecies`` objects
        """
        return [ds.name for ds in self.defect_species]

    @classmethod
    def from_input_set(cls, input_set: InputSet) -> "DefectSystem":
        """generate ``DefectSystem`` from ``InputSet``

        Args:
            input_set (InputSet): ``InputSet`` defining input parameters

        Returns:
            DefectSystem: ``DefectSystem`` corresponding to provided ``InputSet``
        """

        return cls(
            defect_species=input_set.defect_species,
            dos=input_set.dos,
            volume=input_set.volume,
            temperature=input_set.temperature,
            convergence_tolerance=input_set.convergence_tolerance,
            n_trial_steps=input_set.n_trial_steps,
        )

    @classmethod
    def from_yaml(cls, filename: str, structure_file="", dos_file="") -> "DefectSystem":
        """generate ``DefectSystem`` via a yaml file.

        Args:
            filename (str): path to yaml file containing the ``DefectSystem``
              data
            structure_file (str): path to file containing volume information.
              Defaults to an empty string.
            dos_file (str): path to file containing dos information. Defaults
              to an empty string.

        Returns:
            DefectSystem: ``DefectSystem`` corresponding to provided yaml file
        """

        input_set = InputSet.from_yaml(filename, structure_file, dos_file)
        return cls(
            defect_species=input_set.defect_species,
            dos=input_set.dos,
            volume=input_set.volume,
            temperature=input_set.temperature,
            convergence_tolerance=input_set.convergence_tolerance,
            n_trial_steps=input_set.n_trial_steps,
        )
    
    def _apply_element_constraints(
        self,
        concs: Dict[DefectChargeState, float],
        e_fermi: float
    ) -> None:
        """
        Rescale or seed every DefectChargeState in each element‐pool so that
        the *sum* of their concentrations (per cell) matches fixed_total.

        If the current dilute‐limit sum is zero, we *seed* them from their
        Boltzmann (site) weights at the given Fermi level, then scale up.
        """
        for elem, (fixed_total, pool_list) in self.element_pools.items():

            # 1) normalize pool_list → List[(DefectSpecies, stoich)]
            norm: List[Tuple[DefectSpecies, float]] = []
            for sp_id, stoich in pool_list:
                sp = (
                    self.defect_species_by_name(sp_id)
                    if isinstance(sp_id, str)
                    else sp_id
                )
                norm.append((sp, float(stoich)))

            # 2) current total dopant‐atoms in those states
            current = 0.0
            for sp, stoich in norm:
                tot_sp = sum(concs.get(cs, 0.0) for cs in sp.charge_states)
                current += tot_sp * stoich

            # 3a) if any already present, just scale them
            if current > 0:
                scale = fixed_total / current
                for sp, _ in norm:
                    for cs in sp.charge_states:
                        if cs in concs:
                            concs[cs] *= scale
                continue

            # 3b) otherwise seed from Boltzmann/site‐weights at e_fermi
            weights: List[Tuple[DefectChargeState, float]] = []
            for sp, stoich in norm:
                # sp.site_weights returns (pool_name, cs, weight)
                for _, cs, w in sp.site_weights(e_fermi, self.temperature):
                    weights.append((cs, w * stoich))

            total_w = sum(w for _, w in weights)
            if total_w <= 0:
                continue

            scale = fixed_total / total_w
            for cs, w in weights:
                concs[cs] = w * scale


    def _global_defect_concs(self, e_fermi: float) -> Dict[DefectChargeState, float]:
        """
        Build per-cell concentrations for every DefectChargeState,
        **first** seeding+rescaling dopant states via element-pools at the
        current e_fermi, then applying the site-competition denominator.
        """
        # 1) start from dilute limit + any per‐state fixed_concentration
        all_concs: Dict[DefectChargeState, float] = {}
        for sp in self.defect_species:
            for cs, c_cell in sp.charge_state_concentrations(e_fermi, self.temperature):
                all_concs[cs] = c_cell

        # 2) enforce element-pools *with* the current e_fermi
        self._apply_element_constraints(all_concs, e_fermi)

        # 3) site‐competition (Boltzmann‐denominator approach)
        for pool_name, (N_pool, species_list) in self.site_pools.items():
            # resolve species_list → DefectSpecies objects
            sp_objs = [
                self.defect_species_by_name(x) if isinstance(x, str) else x
                for x in species_list
            ]

            # collect only the variable (non-fixed) states and their weights
            weights: List[Tuple[DefectChargeState, float]] = []
            for sp in sp_objs:
                for _, cs, w in sp.site_weights(e_fermi, self.temperature):
                    if cs.fixed_concentration is None:
                        weights.append((cs, w))

            if not weights:
                continue

            Z = 1.0 + sum(w for _, w in weights)  # partition function

            # reassign each variable state → N_pool * (w / Z)
            for cs, w in weights:
                all_concs[cs] = N_pool * (w / Z)

            # any cs.fixed_concentration remain untouched

        return all_concs


    @classmethod
    def from_dict(cls, dictionary: dict) -> "DefectSystem":
        """generate ``DefectSystem`` from a dictionary

        Args:
            filename (str): path to yaml file containing the ``DefectSystem``
              data
            structure_file (str): path to file containing volume information.
              Defaults to an empty string.
            dos_file (str): path to file containing dos information. Defaults
              to an empty string.

        Returns:
            DefectSystem: ``DefectSystem`` corresponding to provided yaml file
        """
        return cls(
            dos=DOS.from_dict(dictionary["dos"]),
            volume=dictionary["volume"],
            temperature=dictionary["temperature"],
            convergence_tolerance=dictionary["convergence_tolerance"],
            n_trial_steps=dictionary["n_trial_steps"],
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

    def get_sc_fermi(self) -> Tuple[float, float]:
        """
        Solve to find Fermi energy in for which the ``DefectSystem`` is charge neutral

        Returns:
           Tuple[float, float]: Fermi energy, residual

        Raises:
          RuntimeError: if the solver fails does not find a valid solution within
            ``self.dos.emin`` and ``self.dos.emax``

        Note:
            The solver will return the Fermi energy either when
            ``self.convergence_tolerance`` is satisfied or when the solver has
            attempted ``self.n_trial_steps``.
            The residual is the the absolute charge density of
            the solver at the end of the last step. Please ensure the residual
            is satisfactorily low if convergence is not reached. It may be
            prudent to investigate the convergence of the solver with respect to
            ``self.n_trial_steps`` and ``self.convergence_tolerance``.
        """
        # initial guess
        emin = self.dos.emin()
        emax = self.dos.emax()
        direction = +1.0
        e_fermi = (emin + emax) / 2.0
        step = 1.0
        reached_e_min = False
        reached_e_max = False

        # loop until convergence or max number of steps reached
        with warnings.catch_warnings():
            warnings.filterwarnings("once")
            for i in range(self.n_trial_steps):
                q_tot = self.q_tot(e_fermi=e_fermi)
                if e_fermi > emax:
                    if reached_e_min or reached_e_max:
                        raise RuntimeError(
                            f"No solution found between {emin} and {emax}"
                        )
                    reached_e_max = True
                    direction = -1.0
                if e_fermi < emin:
                    if reached_e_max or reached_e_min:
                        raise RuntimeError(
                            f"No solution found between {emin} and {emax}"
                        )
                    reached_e_min = True
                    direction = +1.0
                if abs(q_tot) < self.convergence_tolerance:
                    break
                if q_tot > 0.0:
                    if direction == +1.0:
                        step *= 0.25
                        direction = -1.0
                elif q_tot < 0.0:
                    if direction == -1.0:
                        step *= 0.25
                        direction = +1.0
                e_fermi += step * direction

        # return results
        residual = abs(q_tot)
        return e_fermi, residual

    def total_defect_charge_contributions(self, e_fermi: float) -> Tuple[float,float]:
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

    def get_transition_levels(self) -> Dict[str, List[List]]:
        """Return transition_levels transition levels profiles of all ``DefectSpecies``
        all defects as dictionary of ``{DefectSpecies.name : [e_fermi, e_formation]}``
        over the whole density of states energy range.

        Returns:
            Dict[str, List[List]]: Dictionary giving per-defect transition-level
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
    ) -> Dict[str, Any]:
        """
        Solve for the self-consistent Fermi energy and return carrier + defect
        concentrations.

        This method now pulls its defect concentrations from `_global_defect_concs`,
        so both site-competition and element-pools are applied before anything is
        summed or returned.

        Args:
            decomposed (bool): If False, returns one total per DefectSpecies.
                               If True, returns a dict of charge → conc for each species.
            per_volume (bool): If True, scales concentrations to cm⁻³ (1e24/volume).
                               If False, returns raw per-cell counts.

        Returns:
            Dict[str, Any]: A dict containing
                - "Fermi Energy": float (eV)
                - "p0": hole conc (cm⁻³ or per-cell)
                - "n0": electron conc (cm⁻³ or per-cell)
                - one entry per species, either float or dict[int,float].
        """
        # 1) find the self-consistent Fermi level
        e_fermi = self.get_sc_fermi()[0]
        # 2) get carrier concentrations
        p0, n0 = self.dos.carrier_concentrations(e_fermi, self.temperature)
        # 3) get all defect concentrations per cell (with pools applied)
        all_concs = self._global_defect_concs(e_fermi)

        # 4) decide scale for per-volume vs per-cell
        scale = 1e24 / self.volume if per_volume else 1.0

        # 5) assemble the output
        result: Dict[str, Any] = {
            "Fermi Energy": float(e_fermi),
            "p0": float(p0 * scale),
            "n0": float(n0 * scale),
        }

        if not decomposed:
            # sum each species’ total concentration
            for sp in self.defect_species:
                total_sp = sum(
                    conc
                    for cs, conc in all_concs.items()
                    if cs in sp.charge_states
                )
                result[sp.name] = float(total_sp * scale)
        else:
            # break down by charge-state
            for sp in self.defect_species:
                breakdown: Dict[int, float] = {}
                for cs, conc in all_concs.items():
                    if cs in sp.charge_states:
                        breakdown[cs.charge] = float(conc * scale)
                result[sp.name] = breakdown

        return result

    def site_percentages(
        self, 
    ) -> Dict[str, float]:
        """Returns a dictionary of the DefectSpecies in the DefectSystem which
        giving the percentage of the sites in the structure that will host that 
        defect.

        Returns:
            Dict[str, Any]: dictionary specifying the per-DefectSpecies site
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
        """

        Returns:
            dict: _description_
        """

        defect_system_dict = dict(
            volume=float(self.volume),
            temperature=float(self.temperature),
            n_trial_steps=int(self.n_trial_steps),
            defect_species=[
                defect_species.as_dict() for defect_species in self.defect_species
            ],
            convergence_tolerance=float(self.convergence_tolerance),
            dos=self.dos.as_dict(),
        )
        return defect_system_dict
