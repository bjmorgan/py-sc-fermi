from typing import Dict, List, Tuple, Any, Optional, Union, Tuple
import warnings

from py_sc_fermi.dos import DOS
from py_sc_fermi.defect_species import DefectSpecies
from py_sc_fermi.defect_charge_state import DefectChargeState
from py_sc_fermi.inputs import InputSet


class CustomWarningManager:
    def __init__(self):
        self.dos_overflow_warning_issued = False
        self.defect_overflow_warning_issued = False

    def custom_warning(self, message, category, filename, lineno, file=None, line=None):
        if category == RuntimeWarning:
            if "dos" in str(filename).lower() and "overflow" in str(message).lower():
                if not self.dos_overflow_warning_issued:
                    print(
                        "DOSOverflowWarning: An overflow occurred during computation of "
                        "electron and hole concentrations. This is likely a numerical "
                        "artifact of the Fermi‐level search. You can usually ignore it, "
                        "but please verify your results."
                    )
                    self.dos_overflow_warning_issued = True
            elif "defect" in str(filename).lower() and "overflow" in str(message).lower():
                if not self.defect_overflow_warning_issued:
                    print(
                        "DefectOverflowWarning: An overflow occurred during computation of "
                        "defect concentrations. This is likely a numerical artifact of the "
                        "Fermi‐level search. You can usually ignore it, but please verify "
                        "your results."
                    )
                    self.defect_overflow_warning_issued = True
            else:
                print(f"RuntimeWarning: {message}")


# Install our custom warning handler
warning_manager = CustomWarningManager()
warnings.showwarning = warning_manager.custom_warning


class DefectSystem:
    """
    Compute self‐consistent Fermi energy, carrier, and defect concentrations
    in a material, with optional site‐competition and element‐pool constraints.

    Args:
        defect_species (List[DefectSpecies]):
            List of DefectSpecies objects in this system.
        dos (DOS):
            The DOS object for the host material.
        volume (float):
            Unit cell volume in Å³.
        temperature (float):
            Absolute temperature in Kelvin.
        site_pools (Optional[Dict[str, Tuple[float, List[Union[str,DefectSpecies]]]]]):
            Maps a "pool name" → (N_pool, list of species that compete for these sites).
            If None or empty, no site competition is applied.
        element_pools (Optional[Dict[str, Tuple[float, List[Tuple[Union[str,DefectSpecies], float]]]]]):
            Maps an element symbol → (N_fixed, list of (species, stoichiometry) tuples),
            enforcing that the sum of those species’ concentrations (× stoich) per cell
            equals exactly N_fixed.  If None or empty, no element constraints are applied.
        convergence_tolerance (float, optional):
            Convergence tolerance for the charge‐neutrality solver. Default 1e‐18.
        n_trial_steps (int, optional):
            Maximum number of iterations for the Fermi‐level solver. Default 1500.
    """

    def __init__(
        self,
        defect_species: List[DefectSpecies],
        dos: DOS,
        volume: float,
        temperature: float,
        site_pools: Optional[Dict[str, Tuple[float, List[Union[str,DefectSpecies]]]]] = None,
        element_pools: Optional[Dict[str, Tuple[float, List[Tuple[Union[str,DefectSpecies], float]]]]] = None,
        convergence_tolerance: float = 1e-18,
        n_trial_steps: int = 1500,
    ):
        self.defect_species = defect_species
        self.dos = dos
        self.volume = volume
        self.temperature = temperature
        self.site_pools = site_pools or {}
        self.element_pools = element_pools or {}
        self.convergence_tolerance = convergence_tolerance
        self.n_trial_steps = n_trial_steps

    def __repr__(self) -> str:
        lines = [
            "DefectSystem\n",
            f"  nelect: {self.dos.nelect} e\n",
            f"  bandgap: {self.dos.bandgap:.3f} eV\n",
            f"  volume: {self.volume:.3f} Å³\n",
            f"  temperature: {self.temperature:.1f} K\n",
            "\nContains defect species:\n",
        ]
        for ds in self.defect_species:
            lines.append(str(ds))
        return "".join(lines)

    @property
    def defect_species_names(self) -> List[str]:
        """List of the names of all DefectSpecies in this system."""
        return [ds.name for ds in self.defect_species]

    @classmethod
    def from_input_set(cls, input_set: InputSet) -> "DefectSystem":
        """Construct a DefectSystem from an InputSet."""
        return cls(
            defect_species=input_set.defect_species,
            dos=input_set.dos,
            volume=input_set.volume,
            temperature=input_set.temperature,
            site_pools=input_set.site_pools if hasattr(input_set, "site_pools") else {},
            element_pools=input_set.element_pools if hasattr(input_set, "element_pools") else {},
            convergence_tolerance=input_set.convergence_tolerance,
            n_trial_steps=input_set.n_trial_steps,
        )

    @classmethod
    def from_yaml(cls, filename: str, structure_file: str = "", dos_file: str = "") -> "DefectSystem":
        """
        Construct a DefectSystem from a YAML file via InputSet.from_yaml.
        (This assumes your YAML has keys for defect_species, dos, etc.)
        """
        input_set = InputSet.from_yaml(filename, structure_file, dos_file)
        return cls.from_input_set(input_set)

    def defect_species_by_name(self, name: str) -> DefectSpecies:
        """Return the DefectSpecies object whose name matches `name`."""
        for ds in self.defect_species:
            if ds.name == name:
                return ds
        raise KeyError(f"DefectSpecies '{name}' not found in DefectSystem.")

    def _apply_element_constraints(
        self,
        concs: Dict[DefectChargeState, float],
        e_fermi: float
    ) -> None:
        """
        Seed or rescale *only* the dopant states specified in element_pools so that
        their total per‐cell concentration equals exactly `fixed_total`.

        For each (elem → (fixed_total, pool_list)):
          1) If the dilute‐limit sum over those species is > 0, simply scale uniformly
             so that sum_i [c_i * stoich_i] == fixed_total.
          2) If that sum is zero, "seed" them from their Boltzmann (site) weights
             at the given e_fermi, then scale up so that sum == fixed_total.
        """
        for elem, (fixed_total, pool_list) in self.element_pools.items():

            # Normalize pool_list → List[(DefectSpecies, stoich)]
            norm: List[Tuple[DefectSpecies, float]] = []
            for sp_id, stoich in pool_list:
                sp_obj = self.defect_species_by_name(sp_id) if isinstance(sp_id, str) else sp_id
                norm.append((sp_obj, float(stoich)))

            # Compute current total dopant‐atoms in those species (per‐cell)
            current = 0.0
            for sp_obj, stoich in norm:
                sum_sp = sum(concs.get(cs, 0.0) for cs in sp_obj.charge_states)
                current += sum_sp * stoich

            # (2a) If any dopant is already nonzero, scale all dopant states by fixed_total/current
            if current > 0.0:
                scale = fixed_total / current
                for sp_obj, _ in norm:
                    for cs in sp_obj.charge_states:
                        if cs in concs:
                            concs[cs] *= scale
                continue

            # (2b) Otherwise, seed dopant states from their Boltzmann/site‐weights at e_fermi
            weights: List[Tuple[DefectChargeState, float]] = []
            for sp_obj, stoich in norm:
                for _, cs, w in sp_obj.site_weights(e_fermi, self.temperature):
                    weights.append((cs, w * stoich))

            total_w = sum(w for _, w in weights)
            if total_w <= 0.0:
                # No valid Boltzmann weights → leave as zero
                continue

            scale = fixed_total / total_w
            for cs, w in weights:
                concs[cs] = w * scale

    def _apply_element_constraints(
        self,
        concs: Dict[DefectChargeState, float],
        e_fermi: float
    ) -> None:
        """
        Explicitly seed *all* dopant states in the pool from their
        Boltzmann/site‐weights at the current Fermi level, then rescale
        so that the sum matches fixed_total.  This ensures that even states
        with zero dilute‐limit conc will be populated.
        """
        for elem, (fixed_total, pool_list) in self.element_pools.items():
            # 1) normalize each entry to a DefectSpecies + stoich
            norm: List[Tuple[DefectSpecies, float]] = []
            for sp_id, stoich in pool_list:
                sp = (
                    self.defect_species_by_name(sp_id)
                    if isinstance(sp_id, str)
                    else sp_id
                )
                norm.append((sp, float(stoich)))

            # 2) gather Boltzmann/site‐weights for *all* charge states in pool
            weights: List[Tuple[DefectChargeState, float]] = []
            for sp, stoich in norm:
                # site_weights returns (pool_name, cs, weight)
                for _, cs, w in sp.site_weights(e_fermi, self.temperature):
                    weights.append((cs, w * stoich))

            total_w = sum(w for _, w in weights)
            if total_w <= 0:
                # nothing to seed
                continue

            # 3) seed + scale into your fixed_total
            scale = fixed_total / total_w
            for cs, w in weights:
                concs[cs] = w * scale

    def _global_defect_concs(self, e_fermi: float) -> Dict[DefectChargeState, float]:
        """
        Build per‐cell concentrations for every DefectChargeState:
          1) Start from dilute‐limit + any per‐state fixed_concentration
          2) Call _apply_element_constraints(...) to seed & cap dopant states
          3) Perform site‐competition (Fermi‐Dirac denominator) ONLY over native
             (non‐dopant, non‐fixed) charge states.
          4) Return the final per‐cell dictionary all_concs.

        After this, each dopant species (e.g. Mg_Li, Mg_Ni, Mg_i) will be EXACTLY
        the number of atoms determined by element_pools, divided among its own
        charge‐states.  Those dopant charge states are never changed by the 
        site-competition step.
        """
        from warnings import warn

        # 1) Get dilute‐limit + per‐state fixed concentrations
        all_concs: Dict[DefectChargeState, float] = {}
        for sp in self.defect_species:
            for cs, c_cell in sp.charge_state_concentrations(e_fermi, self.temperature):
                all_concs[cs] = c_cell

        # 2) Seed/scale all dopant charge‐states so that 
        #    sum_{dopant_cs}(conc[cs]·stoich) == fixed_total, obeying priority.
        self._apply_element_constraints(all_concs, e_fermi)

        # Build set of all dopant charge‐states, so we never touch them in site‐competition:
        dopant_cs = set()
        for _, (_, pool_list) in self.element_pools.items():
            for sp_id, _ in pool_list:
                sp_obj = (
                    self.defect_species_by_name(sp_id)
                    if isinstance(sp_id, str)
                    else sp_id
                )
                dopant_cs.update(sp_obj.charge_states)

        # 3) Site‐competition (Fermi–Dirac denominator) over only the NAT IVE states 
        #    i.e. (cs.fixed_concentration is None) AND (cs not in dopant_cs).
        for pool_name, (N_pool, species_list) in self.site_pools.items():
            # Normalize species_list → DefectSpecies objects
            sp_objs: List[DefectSpecies] = []
            for entry in species_list:
                sp_objs.append(
                    self.defect_species_by_name(entry)
                    if isinstance(entry, str)
                    else entry
                )

            # (a) Compute “fixed_sum” = any truly fixed_concentration states + all dopant_cs
            fixed_sum = 0.0
            for sp in sp_objs:
                for cs in sp.charge_states:
                    if (cs.fixed_concentration is not None) or (cs in dopant_cs):
                        fixed_sum += all_concs.get(cs, 0.0)

            free_sites = N_pool - fixed_sum
            if free_sites < 0.0:
                # Dopants (or truly fixed) have already overfilled this pool.
                warn(
                    f"Pool '{pool_name}' is overfilled by {-free_sites:.3g} sites.\n"
                    "→ Zeroing out all native (non‐dopant, non‐fixed) states in this pool."
                )
                for sp in sp_objs:
                    for cs in sp.charge_states:
                        if (cs.fixed_concentration is None) and (cs not in dopant_cs):
                            all_concs[cs] = 0.0
                # Skip redistributing “free_sites”
                continue

            # (b) Gather Boltzmann weights for only the native (variable, non‐dopant) states
            weights: List[Tuple[DefectChargeState, float]] = []
            for sp in sp_objs:
                for _, cs, w in sp.site_weights(e_fermi, self.temperature):
                    if (cs.fixed_concentration is None) and (cs not in dopant_cs):
                        weights.append((cs, w))

            if not weights:
                # No variable states to distribute
                continue

            # (c) Partition function: Z = 1 + Σ w_native
            Z = 1.0 + sum(w for _, w in weights)

            # (d) Redistribute “free_sites” among those native states:
            for cs, w in weights:
                all_concs[cs] = free_sites * (w / Z)
            # (Any dopant_cs or truly fixed_concentration states remain untouched.)

        # 4) At this point, all dopant charge‐states are exactly what we placed in step 2,
        #    and no site‐competition step will have moved them.  We can simply return all_concs.
        return all_concs


    def total_defect_charge_contributions(self, e_fermi: float) -> Tuple[float, float]:
        """
        Return (lhs_defect, rhs_defect) = total positive defect charge,
        total negative defect charge, at the given e_fermi.
        """
        lhs = 0.0
        rhs = 0.0
        all_concs = self._global_defect_concs(e_fermi)
        for cs, conc_cell in all_concs.items():
            if cs.charge > 0:
                lhs += conc_cell * cs.charge
            elif cs.charge < 0:
                rhs += conc_cell * abs(cs.charge)
        return lhs, rhs

    def q_tot(self, e_fermi: float) -> float:
        """
        Net charge density = (electrons + negative defect charge) - (holes + positive defect charge).
        We return (n0 + rhs_def) - (p0 + lhs_def) so that zero crossing → neutrality.
        """
        p0, n0 = self.dos.carrier_concentrations(e_fermi, self.temperature)
        lhs_def, rhs_def = self.total_defect_charge_contributions(e_fermi)
        # lhs = p0 + lhs_def, rhs = n0 + rhs_def.  We want q_tot = rhs - lhs = (n0+rhs_def) - (p0+lhs_def).
        return (n0 + rhs_def) - (p0 + lhs_def)

    def get_sc_fermi(self) -> Tuple[float, float]:
        """
        Solve self‐consistently for the Fermi energy such that q_tot(e_fermi) ≈ 0.
        Returns (e_fermi, residual).  Raises RuntimeError if no solution found.
        """
        emin = self.dos.emin()
        emax = self.dos.emax()
        e_fermi = 0.5 * (emin + emax)
        step = (emax - emin) / 2.0
        direction = +1.0
        reached_e_min = False
        reached_e_max = False

        with warnings.catch_warnings():
            warnings.filterwarnings("once")
            for _ in range(self.n_trial_steps):
                qt = self.q_tot(e_fermi)
                if abs(qt) < self.convergence_tolerance:
                    break
                if qt > 0:  # net negative charge → lower e_fermi
                    if direction > 0:
                        step *= 0.25
                        direction = -1.0
                    e_fermi += direction * step
                else:      # net positive charge → raise e_fermi
                    if direction < 0:
                        step *= 0.25
                        direction = +1.0
                    e_fermi += direction * step

                if e_fermi < emin:
                    if reached_e_min or reached_e_max:
                        raise RuntimeError(f"No Fermi solution in [{emin}, {emax}]")
                    reached_e_min = True
                    direction = +1.0
                    e_fermi = emin + 0.5 * step
                if e_fermi > emax:
                    if reached_e_min or reached_e_max:
                        raise RuntimeError(f"No Fermi solution in [{emin}, {emax}]")
                    reached_e_max = True
                    direction = -1.0
                    e_fermi = emax - 0.5 * step

        return e_fermi, abs(qt)

    def get_transition_levels(self) -> Dict[str, List[List[float]]]:
        """
        Return a dict {species_name: [list_of_e_fermi_points, list_of_formation_energies]}
        representing the full transition‐level profile for each defect.
        """
        tl_dict: Dict[str, List[List[float]]] = {}
        for name in self.defect_species_names:
            sp_obj = self.defect_species_by_name(name)
            tl_profile = sp_obj.tl_profile(self.dos.emin(), self.dos.emax())
            e_fermi_vals = [pt[0][0] for pt in tl_profile]
            e_form_vals = [pt[0][1] for pt in tl_profile]
            tl_dict[name] = [e_fermi_vals, e_form_vals]
        return tl_dict

    def concentration_dict(
        self,
        decomposed: bool = False,
        per_volume: bool = True
    ) -> Dict[str, Any]:
        """
        Solve for the self‐consistent Fermi energy and return carrier + defect concentrations.

        Args:
          decomposed (bool): If False, returns one total concentration per DefectSpecies
                              (summed over charge states).  If True, returns a dict
                              of {charge → concentration} for each species.
          per_volume (bool): If True, scale concentrations to cm⁻³ via (1e24/volume);
                             If False, return raw per‐cell counts.

        Returns:
          Dict with keys:
            "Fermi Energy" → float (eV)
            "p0" → hole concentration (cm⁻³ or per‐cell)
            "n0" → electron concentration (cm⁻³ or per‐cell)
            one entry per species_name → either float (total) or Dict[int,float] (decomposed).
        """
        e_fermi, _ = self.get_sc_fermi()
        p0, n0 = self.dos.carrier_concentrations(e_fermi, self.temperature)
        all_concs = self._global_defect_concs(e_fermi)

        scale = (1e24 / self.volume) if per_volume else 1.0
        result: Dict[str, Any] = {
            "Fermi Energy": float(e_fermi),
            "p0": float(p0 * scale),
            "n0": float(n0 * scale),
        }

        if not decomposed:
            for sp_obj in self.defect_species:
                total_sp = sum(
                    conc_cell
                    for cs, conc_cell in all_concs.items()
                    if cs in sp_obj.charge_states
                )
                result[sp_obj.name] = float(total_sp * scale)
        else:
            for sp_obj in self.defect_species:
                breakdown: Dict[int, float] = {}
                for cs, conc_cell in all_concs.items():
                    if cs in sp_obj.charge_states:
                        breakdown[cs.charge] = float(conc_cell * scale)
                result[sp_obj.name] = breakdown

        return result

    def site_percentages(self) -> Dict[str, float]:
        """
        Return the percent occupancy of each DefectSpecies, i.e.
          (c_site / nsites) * 100 for each species’ lowest‐energy charge state,
          evaluated at the self‐consistent Fermi level.
        """
        e_fermi, _ = self.get_sc_fermi()
        return {
            ds.name: float(
                (ds.get_concentration(e_fermi, self.temperature) / ds.nsites) * 100.0
            )
            for ds in self.defect_species
        }

    def as_dict(self) -> dict:
        """
        Serialize this DefectSystem to a dictionary (for JSON/YAML output).
        """
        return {
            "volume": float(self.volume),
            "temperature": float(self.temperature),
            "n_trial_steps": int(self.n_trial_steps),
            "convergence_tolerance": float(self.convergence_tolerance),
            "dos": self.dos.as_dict(),
            "defect_species": [ds.as_dict() for ds in self.defect_species],
            "site_pools": self.site_pools,
            "element_pools": self.element_pools,
        }
