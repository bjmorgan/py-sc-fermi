import numpy as np
from py_sc_fermi.constants import kboltz


class DefectChargeState:
    """Class for individual defect charge states"""

    def __init__(
        self,
        charge: int,
        energy: float = None,
        fixed_concentration: float = None,
        degeneracy: int = 1,
    ):
        """Instantiate a DefectChargeState.

        Args:
            charge (int): Charge of this charge state.
            energy (float): Formation energy of this charge state when E_Fermi = E(VBM)
            degeneracy (int): Degeneracy of this charge state (e.g. spin/intrinsic degeneracy).
            fixed_concentration (float): the fixed concentration of the defect charge state (default = None)

        Returns:
            None

        Raises:
            ValueError: If `energy` and `fixed concentration` == None.

        """
        if energy == None and fixed_concentration == None:
            raise ValueError(
                "You must specify either a fixed concentration or energy for this defect! \n Note, if you specify both, the concentration will treated as fixed"
            )
        self._charge = charge
        self._energy = energy
        self._degeneracy = degeneracy
        self._fixed_concentration = fixed_concentration

    def fix_concentration(self, concentration: float) -> None:
        """fix the net concentration (per unit cell) of this defect species"""
        self._fixed_concentration = concentration

    @property
    def energy(self) -> float:
        """Get the energy of this charge state at E_Fermi = E(VBM) (0)."""
        return self._energy

    @property
    def charge(self) -> int:
        """Get the charge of this charge state."""
        return self._charge

    @property
    def degeneracy(self) -> int:
        """Get the degeneracy of this charge state."""
        return self._degeneracy

    @property
    def fixed_concentration(self) -> float:
        """The fixed net concentration (per unit cell) of this defect species,
           or `None` if the defect concentrations are free to change"""
        return self._fixed_concentration


    def get_formation_energy(self, e_fermi: float) -> float:
        """Calculate the formation energy of this charge state at a
        specified Fermi energy.

        Args:
            e_fermi (float): Position of the Fermi energy, relative to the VBM (in eV).

        Returns:
            (float): The charge state energy (in eV).

        """
        return self.energy + self.charge * e_fermi

    def get_concentration(self, e_fermi: float, temperature: float) -> float:
        """Calculate the concentration of this charge state at a
        specified Fermi energy and temperature, per site in the unit
        cell..

        Args:
            e_fermi (float): Fermi energy relative to the VBM (in eV).
            temperature (float): Temperature (in K).

        Returns:
            concentration (float): Per calculation-cell concentration of this charge state.

        """
        if self.fixed_concentration == None:
            expfac = -self.get_formation_energy(e_fermi) / (kboltz * temperature)
            concentration = self.degeneracy * np.exp(expfac)
        else:
            concentration = self.fixed_concentration
        return concentration

    def __repr__(self):
        if self.fixed_concentration == None:
            return f"q={self.charge:+2}, e={self.energy}, deg={self.degeneracy}"
        else:
            return f"q={self.charge:+2}, [c]={self.fixed_concentration}, deg={self.degeneracy}"
