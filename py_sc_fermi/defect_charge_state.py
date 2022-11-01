import numpy as np  # type: ignore
from scipy.constants import physical_constants  # type: ignore
from typing import Optional

kboltz = physical_constants["Boltzmann constant in eV/K"][0]


class DefectChargeState:
    """Class describing a single charge state (of a ``DefectSpecies``).

    Args:
         charge (int): charge of this ``DefectChargeState``
         degeneracy (int): degeneracy per unit cell
         energy (float): formation energy at E[Fermi] = 0
         fixed_concentration (float): fixed concentration per unit cell
    """

    def __init__(
        self,
        charge: int,
        degeneracy: int = 1,
        energy: float = None,
        fixed_concentration: float = None,
    ):
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
        """formation energy of the ``DefectChargeState`` at E[vbm] (E[Fermi] = 0)

        Returns:
            Optional[float]: formation energy
        """
        return self._energy

    @property
    def charge(self) -> int:
        """charge of the ``DefectChargeState``

        Returns:
            int: charge
        """
        return self._charge

    @property
    def degeneracy(self) -> int:
        """The degeneracy of the ``DefectChargeState`` (e.g. spin degeneracy)

        Returns:
            int: degeneracy per unit cell
        """
        return self._degeneracy

    @property
    def fixed_concentration(self) -> Optional[float]:
        """fixed concentration of this ``DefectChargeState`` or ``None`` if the
        concentration is free to vary.

        Returns:
            Optional[float]: fixed concentration per unit cell
        """
        return self._fixed_concentration

    @classmethod
    def from_string(
        cls, string: str, volume: Optional[float] = None, frozen: bool = False
    ) -> "DefectChargeState":
        """
        Create a ``DefectChargeState`` from a given string. This method was
        envisaged for use as a way to read in defect charge states from an input
        file for `SC-Fermi <https://github.com/jbuckeridge/sc-fermi>`_.

        If a user does wish to specify a defect charge state using this
        functionality, the string should be in the form:

        `charge  formation_energy  degeneracy`

        i.e. a defect with charge 2, formation energy of 0.1 eV and degeneracy
        of 2 would be specified as:

        ``"2 0.1 2"``

        if the charge state has a fixed concentration, the string should be in
        the form:

        `charge  concentration`

        i.e. a defect with charge 2, concentration of 1e21 per cm-3
        would be specified as:

        ``"2 1e21"``

        Args:
            string (str): string representation of the ``DefectChargeState``
            volume (Optional[float], optional): volume of the unit cell, only
                if ``frozen == True``. Defaults to ``None``.
            frozen (bool, optional): if ``True`` the concentration of this
                ``DefectChargeState`` cannot change when solving for a self
                consistent Fermi energy. Defaults to ``False``.

        Raises:
            ValueError: if defect concentration is fixed, but ``volume == None``

        Returns:
            ``DefectChargeState``: relevant ``DefectChargeState`` object
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
        """Fixes the concentration (per unit cell) of this ``DefectChargeState``

        Args:
            concentration (float): ``DefectChargeState`` concentration per unit cell
        """
        self._fixed_concentration = concentration

    def get_formation_energy(self, e_fermi: float) -> float:
        """get the formation energy of this ``DefectChargeState`` at a given Fermi
        energy

        Args:
            e_fermi (float): Fermi energy at which to calculate the formation energy

        Raises:
            ValueError: if ``DefectChargeState.energy == None``

        Returns:
            float: formation energy of ``DefectChargeState`` at ``e_fermi``
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

        Args:
            e_fermi (float): Fermi energy.
            temperature (float): Temperature.

        Returns:
            float: Concentration at the specified Fermi energy and temperature.
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
