import numpy as np  # type: ignore
from scipy.constants import physical_constants  # type: ignore
from typing import Optional

kboltz = physical_constants["Boltzmann constant in eV/K"][0]


class DefectChargeState:
    """Class describing an individual defect charge state

    :param int charge: integer charge of the defect
    :param float energy: Formation energy of this charge state at E(VBM)
    :param int degeneracy: Degeneracy of this charge state
        (e.g. spin/intrinsic degeneracy).
    :param float fixed_concentration: the fixed concentration of the defect
        charge state per unit cell (Default ``None``)

    :raises ValueError: If ``energy == None`` and ``fixed concentration == None``.
     ``DefectChargeState`` may possess either a formation energy 
     or a fixed concentration, or both.

    .. note:: If both a formation energy and fixed_concentration are specified,
        the concentration of the ``DefectChargeState`` will treated as fixed.

    """

    def __init__(
        self,
        charge: int,
        degeneracy: int = 1,
        energy: float = None,
        fixed_concentration: float = None,
    ):
        """Initialise a ``DefectChargeState`` instance."""

        if energy == None and fixed_concentration == None:
            raise ValueError(
                """You must specify either a fixed concentration or energy for 
                this defect! \n Note, if you specify both, the concentration
                will treated as fixed"""
            )
        self._charge = charge
        self._degeneracy = degeneracy
        self._energy = energy
        self._fixed_concentration = fixed_concentration

    @property
    def energy(self) -> Optional[float]:
        """:return: Formation energy of this charge state at E[VBM] (E_Fermi = 0)"""
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
    def fixed_concentration(self) -> Optional[float]:
        """:return: the fixed concentration of this defect charge state, 
        or ``None`` if the concentration is free to vary."""
        return self._fixed_concentration

    @classmethod
    def from_string(
        cls, string: str, volume: Optional[float] = None, frozen: bool = False
    ):
        """
        Create a ``DefectChargeState`` from a given string. This method was
        envisaged for use as a way to read in defect charge states from an input
        file for SC-Fermi. 
        
        If a user does wish to specify a defect charge state using this
        functionality, the string should be in the form:
             **charge formation_energy degeneracy**
        i.e. a defect with charge 2, formation energy of 0.1 eV and degeneracy
        of 2 would be specified as:
            ``"2 0.1 2"``
        if the charge state has a fixed concentration, the string should be in
        the form:
            **charge concentration**
        i.e. a defect with charge 2, concentration of 1e21 per cm-3
        would be specified as:
            ``"2 1e21"``

        :param str string: String representation of a defect charge state.
        :param float volume: Volume of the unit cell. Only required if the
            charge state has a fixed concentration.
        :param bool frozen: If ``True``, the parser will expect to read the string as a
            fixed concentration ``DefectChargeState``.
        """
        string = string.strip()
        stripped_string = string.split()
        if frozen is False:
            return cls(
                charge=int(stripped_string[0]),
                energy=float(stripped_string[1]),
                degeneracy=int(stripped_string[2]),
            )
        else:
            if volume is None:
                raise ValueError(
                    "You must specify a real, positive cell volume if passing a frozen concentration!"
                )
            else:
                return cls(
                    charge=int(stripped_string[1]),
                    fixed_concentration=float(stripped_string[2]) / 1e24 * volume,
                )

    def fix_concentration(self, concentration: float) -> None:
        """fix the net concentration (per unit cell) of this defect charge state

        :param ``float`` concentration: the fixed concentration of this defect
        """
        self._fixed_concentration = concentration

    def get_formation_energy(self, e_fermi: float) -> float:
        """Calculate the formation energy of this charge state at a
        specified Fermi energy.

        :param float e_fermi: Fermi energy relative to the VBM (in eV).
        :return: Formation energy of this charge state when E_Fermi = E(VBM)
        :rtype: float

        :raises ValueError: If ``self.energy == None``.
        """
        if self.energy is not None:
            return self.energy + self.charge * e_fermi
        else:
            raise ValueError(
                "Cannot calculate formation energy as a function of `e_fermi` without a defined formation energy!"
            )

    def get_concentration(self, e_fermi: float, temperature: float) -> float:
        """Calculate the concentration of this ``DefectChargeState`` at a
        specified Fermi energy and temperature, per site in the unit
        cell.

        :param float e_fermi: Fermi energy.
        :param float temperature: Temperature.
        :return concentration: Concentration at the specified
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
