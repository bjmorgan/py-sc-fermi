import numpy as np
from scipy.constants import physical_constants  # type: ignore

kboltz = physical_constants["Boltzmann constant in eV/K"][0]


class DefectChargeState:
    """Class for individual defect charge states

    :param int charge: Charge of this charge state.
    :param float energy: Formation energy of this charge state when
        E_Fermi = E(VBM)
    :param int degeneracy: Degeneracy of this charge state 
        (e.g. spin/intrinsic degeneracy).
    :param float fixed_concentration: the fixed concentration of the defect
        charge state (default = None)

    :raises ValueError: If `energy` and `fixed concentration` == None.

    .. note:: If both a formation energy and fixed_concentration are specified,
        the concentration of the ``DefectChargeState`` will treated as fixed.

    """

    def __init__(
        self,
        charge: int,
        energy: float = None,
        fixed_concentration: float = None,
        degeneracy: int = 1,
    ):
        """Initialise a ``DefectChargeState`` instance."""

        if energy == None and fixed_concentration == None:
            raise ValueError(
                "You must specify either a fixed concentration or energy for this defect! \n Note, if you specify both, the concentration will treated as fixed"
            )
        self._charge = charge
        self._energy = energy
        self._degeneracy = degeneracy
        self._fixed_concentration = fixed_concentration

    def fix_concentration(self, concentration: float) -> None:
        """fix the net concentration (per calculation cell) of this defect
           species
           :param float concentration: the fixed concentration of this defect
           """
        self._fixed_concentration = concentration

    @property
    def energy(self) -> float:
        """:return: Formation energy of this charge state when E_Fermi = E(VBM)
        """
        return self._energy

    @property
    def charge(self) -> int:
        """:return: Charge of this charge state."""
        return self._charge

    @property
    def degeneracy(self) -> int:
        """:return: The number of energetically degenerate states for this
           charge state."""
        return self._degeneracy

    @property
    def fixed_concentration(self) -> float:
        """:return: the fixed concentration of this defect charge state, or None
           if the concentration is free to vary."""
        return self._fixed_concentration

    def get_formation_energy(self, e_fermi: float) -> float:
        """Calculate the formation energy of this charge state at a
        specified Fermi energy.

        :param float e_fermi: Fermi energy relative to the VBM (in eV).
        :return: Formation energy of this charge state when E_Fermi = E(VBM)
        :rtype: float
        """
        return self.energy + self.charge * e_fermi

    def get_concentration(self, e_fermi: float, temperature: float) -> float:
        """Calculate the concentration of this charge state at a
        specified Fermi energy and temperature, per site in the unit
        cell.

        :param float e_fermi: Fermi energy relative to the VBM (in eV).
        :param float temperature: Temperature (in K).
        :return concentration: Concentration of this charge state at the specified
            Fermi energy and temperature.
        :rtype: float

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
