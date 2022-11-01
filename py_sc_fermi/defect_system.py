from typing import Dict, List, Tuple, Any
from py_sc_fermi.dos import DOS
from py_sc_fermi.defect_species import DefectSpecies
from py_sc_fermi.inputs import InputSet
import numpy as np


class DefectSystem(object):
    """This class is used to calculate the self consistent Fermi energy for
    a defective material, observing the condition of charge neutrality and
    therefore, point defect and carrier concentrations under equilibrium
    conditions.

    Args:
        defect_species (List[DefectSpecies]): List of ``DefectSpecies`` objects
          which are present in the ``DefectSystem``.
        volume (float): volume of the unit cell in Angstroms cubed
        dos (DOS): the ``DOS`` object associated with the unit cell
        temperature (float): temperature at which self-consentient Fermi energy
          will be solved for.
        convergence_tolerance (float): the charge neutrality tolerance for the
          self-consistent Fermi energy solver. Defaults to ``1e-18``.
        n_trial_steps (int): the maximum number of steps to take in the
          self-consistent Fermi energy solver. Defaults to 1500.
    """

    def __init__(
        self,
        defect_species: List[DefectSpecies],
        dos: DOS,
        volume: float,
        temperature: float,
        convergence_tolerance: float = 1e-18,
        n_trial_steps: int = 1500,
    ):

        self.defect_species = defect_species
        self.volume = volume
        self.dos = dos
        self.temperature = temperature
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
    def from_yaml(cls, filename: str) -> "DefectSystem":
        """generate ``DefectSystem`` via a yaml file.

        Args:
            filename (str): path to yaml file containing the ``DefectSystem``
              data

        Returns:
            DefectSystem: ``DefectSystem`` corresponding to provided yaml file
        """

        input_set = InputSet.from_yaml(filename)
        return cls(
            defect_species=input_set.defect_species,
            dos=input_set.dos,
            volume=input_set.volume,
            temperature=input_set.temperature,
            convergence_tolerance=input_set.convergence_tolerance,
            n_trial_steps=input_set.n_trial_steps,
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
        converged = False
        reached_e_min = False
        reached_e_max = False

        # loop until convergence or max number of steps reached
        for i in range(self.n_trial_steps):
            q_tot = self.q_tot(e_fermi=e_fermi)
            if e_fermi > emax:
                if reached_e_min or reached_e_max:
                    raise RuntimeError(f"No solution found between {emin} and {emax}")
                reached_e_max = True
                direction = -1.0
            if e_fermi < emin:
                if reached_e_max or reached_e_min:
                    raise RuntimeError(f"No solution found between {emin} and {emax}")
                reached_e_min = True
                direction = +1.0
            if abs(q_tot) < self.convergence_tolerance:
                converged = True
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

    def report(self) -> None:
        """print a report in the style of `SC-Fermi <https://github.com/jbuckeridge/sc-fermi>`_
        which summarises key properties of the defect system."""
        print(self._get_report_string())

    def _get_report_string(self) -> str:
        """generate string to facilitate self.report()"""
        string = ""
        e_fermi = self.get_sc_fermi()[0]
        string += f"SC Fermi level :      {e_fermi}  (eV)\n"
        p0, n0 = self.dos.carrier_concentrations(e_fermi, self.temperature)
        string += "Concentrations:\n"
        string += f"n (electrons)  : {n0 * 1e24 / self.volume} cm^-3\n"
        string += f"p (holes)      : {p0 * 1e24 / self.volume} cm^-3\n"
        for ds in self.defect_species:
            concall = ds.get_concentration(e_fermi, self.temperature)
            if ds.fixed_concentration == None:
                string += f"{ds.name:9}      : {concall * 1e24 / self.volume} cm^-3\n"
            else:
                string += (
                    f"{ds.name:9}      : {concall * 1e24 / self.volume} cm^-3 [fixed]\n"
                )
        string += "\nBreakdown of concentrations for each defect charge state:\n"
        for ds in self.defect_species:
            concall = ds.get_concentration(e_fermi, self.temperature)
            string += "---------------------------------------------------------\n"
            if concall == 0.0:
                string += f"{ds.name:11}: Zero total - cannot give breakdown\n"
                continue
            string += f"{ds.name:11}: Charge Concentration(cm^-3) Total\n"
            for q, conc in ds.charge_state_concentrations(
                e_fermi, self.temperature
            ).items():
                if ds.charge_states[q].fixed_concentration:
                    fix_str = " [fixed]"
                else:
                    fix_str = ""

                string += f"           : {q: 1}  {conc * 1e24 / self.volume:5e}          {(conc * 100 / concall):.2f} {fix_str}\n"
        return string

    def total_defect_charge_contributions(self, e_fermi: float) -> Tuple[float, float]:
        """
        Calculate the charge contributions from each ``DefectSpecies`` in all charge
        states to the total charge density

        Args:
            e_fermi (float): Fermi energy

        Returns:
            Tuple[float, float]: charge contributions of positive (lhs) and
            negative (rhs) charge states of all defects
        """
        contrib = np.array(
            [
                ds.defect_charge_contributions(e_fermi, self.temperature)
                for ds in self.defect_species
            ]
        )
        lhs = np.sum(contrib[:, 0])
        rhs = np.sum(contrib[:, 1])
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

    def as_dict(
        self,
        decomposed: bool = False,
        per_volume: bool = True,
    ) -> Dict[str, Any]:
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
            Dict[str, Any]: dictionary specifying the Fermi Energy,
            hole concentration (``"p0"``), electron concentration
            (``"n0"``), temperature, and the defect concentrations.
        """
        if per_volume == True:
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
        if decomposed == False:
            sum_concs = {
                str(ds.name): float(
                    ds.get_concentration(e_fermi, self.temperature) * scale
                )
                for ds in self.defect_species
            }
            return {**run_stats, **sum_concs}
        else:
            decomp_concs = {
                str(ds.name): {
                    int(q): float(
                        ds.charge_state_concentrations(e_fermi, self.temperature)[q]
                        * scale
                    )
                    for q in ds.charge_states
                }
                for ds in self.defect_species
            }
            return {**run_stats, **decomp_concs}
