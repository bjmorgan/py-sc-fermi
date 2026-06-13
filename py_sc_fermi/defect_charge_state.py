from __future__ import annotations

import warnings

import numpy as np
from scipy.constants import physical_constants

from py_sc_fermi.warnings import suppresses_numpy_overflow

kboltz = physical_constants["Boltzmann constant in eV/K"][0]


class DefectChargeState:
    """Class describing a single charge state (of a ``DefectSpecies``).

    Args:
         charge (int): charge of this ``DefectChargeState``
         degeneracy (float): degeneracy of this charge state
         energy (float): formation energy at E[Fermi] = 0
         fixed_concentration (float): fixed concentration per unit cell
    """

    def __init__(
        self,
        charge: int,
        degeneracy: float = 1,
        energy: float | None = None,
        fixed_concentration: float | None = None,
    ):
        if energy is None and fixed_concentration is None:
            raise ValueError(
                "You must specify either a fixed concentration or an energy for "
                "this defect. If you specify both, the concentration is treated as fixed."
            )
        if degeneracy <= 0:
            raise ValueError("degeneracy must be a positive number")
        self._charge = charge
        self._degeneracy = degeneracy
        self._energy = energy
        self._fixed_concentration = fixed_concentration

    @property
    def energy(self) -> float | None:
        """formation energy of the ``DefectChargeState`` at E[vbm] (E[Fermi] = 0)

        Returns:
            float | None: formation energy
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
    def degeneracy(self) -> float:
        """The degeneracy of the ``DefectChargeState`` (e.g. spin degeneracy)

        Returns:
            float: degeneracy of this charge state
        """
        return self._degeneracy

    @property
    def fixed_concentration(self) -> float | None:
        """fixed concentration of this ``DefectChargeState`` or ``None`` if the
        concentration is free to vary.

        Returns:
            float | None: fixed concentration per unit cell
        """
        return self._fixed_concentration

    @classmethod
    def from_dict(cls, dictionary: dict) -> DefectChargeState:
        """generate a ``DefectChargeState`` object from a dictionary

        Args:
            dictionary (dict): dictionary defining ``DefectChargeState``. Any
              fixed concentration given should be provided per-unit cell

        Returns:
            DefectChargeState: object described by `dictionary`
        """

        valid_keys = ["degeneracy", "energy", "charge", "fixed_concentration"]
        unrecognised_keys = set(dictionary.keys()) - set(valid_keys)
        if unrecognised_keys:
            warnings.warn(
                f"Ignoring unrecognised keys: {', '.join(unrecognised_keys)}",
                stacklevel=2,
            )

        if "fixed_concentration" in dictionary.keys():
            return DefectChargeState(
                degeneracy=dictionary["degeneracy"],
                charge=dictionary["charge"],
                fixed_concentration=dictionary["fixed_concentration"],
            )
        else:
            return DefectChargeState(
                degeneracy=dictionary["degeneracy"],
                energy=dictionary["energy"],
                charge=dictionary["charge"],
            )

    def as_dict(self) -> dict:
        """generate a dictionary representation of the ``DefectChargeState``

        Returns:
            dict: dictionary representation of the ``DefectChargeState``
        """

        defect_dict = {
            "degeneracy": float(self.degeneracy),
            "energy": float(self.energy) if self.energy is not None else None,
            "charge": int(self.charge),
        }
        if self.fixed_concentration is not None:
            defect_dict.update(
                {"fixed_concentration": float(self.fixed_concentration)}
            )

        return defect_dict

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
                "Cannot calculate formation energy as a function of `e_fermi` without "
                "a defined formation energy!"
            )

    @suppresses_numpy_overflow
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
        if self.fixed_concentration is None:
            expfac = -self.get_formation_energy(e_fermi) / (kboltz * temperature)
            concentration = self.degeneracy * np.exp(expfac)
        else:
            concentration = self.fixed_concentration
        return concentration

    def __repr__(self) -> str:
        if self.fixed_concentration is None:
            return f"q={self.charge:+2}, e={self.energy}, deg={self.degeneracy}"
        else:
            return f"q={self.charge:+2}, [c]={self.fixed_concentration}, deg={self.degeneracy}"
