from __future__ import annotations

import math
import warnings
from typing import Any

import numpy as np
from scipy.constants import physical_constants

from py_sc_fermi.warnings import UnrecognisedKeyWarning, suppresses_numpy_overflow

kboltz = physical_constants["Boltzmann constant in eV/K"][0]


class DefectChargeState:
    """Class describing a single charge state (of a ``DefectSpecies``).

    Args:
         charge (int): charge of this ``DefectChargeState``
         degeneracy (float, optional): degeneracy of this charge state.
           Defaults to 1.
         energy (float | None, optional): formation energy at E[Fermi] = 0.
           Defaults to None; either ``energy`` or ``fixed_concentration`` must
           be given.
         fixed_concentration (float | None, optional): fixed concentration per
           unit cell. Must be finite and non-negative. Defaults to None.
         name (str | None): identifying label for this charge state, used as
           the key in ``DefectSystem.result.charge_state_concentrations``
           (keyed by species name then charge-state name). Defaults to the
           charge string (e.g. ``"q+2"``). Metastable configurations sharing a
           formal charge must be given explicit names (e.g. ``"V_O_2+_tet"``
           vs ``"V_O_2+_oct"``), since names must be unique within a
           ``DefectSpecies``.
    """

    def __init__(
        self,
        charge: int,
        degeneracy: float = 1,
        energy: float | None = None,
        fixed_concentration: float | None = None,
        name: str | None = None,
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
        self._name = name
        self._fixed_concentration = (
            self._validated_fixed_concentration(fixed_concentration)
            if fixed_concentration is not None
            else None
        )

    def _validated_fixed_concentration(self, concentration: float) -> float:
        """Return ``concentration`` if it is a valid fixed concentration.

        Raises:
            ValueError: if ``concentration`` is not finite and non-negative.
        """
        if not math.isfinite(concentration) or concentration < 0:
            raise ValueError(
                f"DefectChargeState '{self.name}' has an invalid fixed "
                f"concentration {concentration}; it must be finite and "
                "non-negative"
            )
        return concentration

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
    def name(self) -> str:
        """Identifying label for this charge state.

        Returns:
            str: the explicit name if one was set, otherwise the
            charge-derived default (e.g. ``"q+2"``).
        """
        if self._name is not None:
            return self._name
        return f"q{self._charge:+d}"

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

        valid_keys = ["degeneracy", "energy", "charge", "fixed_concentration", "name"]
        unrecognised_keys = set(dictionary.keys()) - set(valid_keys)
        if unrecognised_keys:
            warnings.warn(
                f"Ignoring unrecognised keys: {', '.join(unrecognised_keys)}",
                UnrecognisedKeyWarning,
                stacklevel=2,
            )

        name = dictionary.get("name", None)
        if "fixed_concentration" in dictionary.keys():
            return DefectChargeState(
                degeneracy=dictionary["degeneracy"],
                charge=dictionary["charge"],
                energy=dictionary.get("energy"),
                fixed_concentration=dictionary["fixed_concentration"],
                name=name,
            )
        else:
            return DefectChargeState(
                degeneracy=dictionary["degeneracy"],
                energy=dictionary["energy"],
                charge=dictionary["charge"],
                name=name,
            )

    def as_dict(self) -> dict:
        """generate a dictionary representation of the ``DefectChargeState``

        The ``name`` key is included only when an explicit name was set;
        charge-derived default names are omitted and regenerated on load.

        Returns:
            dict: dictionary representation of the ``DefectChargeState``
        """

        defect_dict: dict[str, Any] = {
            "degeneracy": float(self.degeneracy),
            "energy": float(self.energy) if self.energy is not None else None,
            "charge": int(self.charge),
        }
        if self.fixed_concentration is not None:
            defect_dict.update(
                {"fixed_concentration": float(self.fixed_concentration)}
            )
        if self._name is not None:
            defect_dict["name"] = self._name

        return defect_dict

    def _fix_concentration(self, concentration: float) -> None:
        """Set the fixed concentration (per unit cell); internal setter used
        by construction-time fixing.

        Args:
            concentration (float): ``DefectChargeState`` concentration per unit cell

        Raises:
            ValueError: if ``concentration`` is not finite and non-negative.
        """
        self._fixed_concentration = self._validated_fixed_concentration(concentration)

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
    def _dilute_site_concentration(self, e_fermi: float, temperature: float) -> float:
        """Calculate the concentration of this ``DefectChargeState`` at a
        specified Fermi energy and temperature. A variable-concentration
        state returns a per-site concentration; a fixed-concentration state
        returns its fixed per-cell concentration unchanged.

        The variable-state value is the dilute-limit (Boltzmann) expression
        ``degeneracy * exp(-E_formation / kT)``, with ``E_formation`` the
        formation energy at ``e_fermi``. It applies no site exclusion, so it
        is unbounded and can exceed one per site; a solved ``DefectSystem``
        applies site-exclusion statistics instead, which agree with this
        expression only at low occupancy.

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
        name_part = f", name={self._name}" if self._name is not None else ""
        if self.fixed_concentration is None:
            return f"q={self.charge:+2}, e={self.energy}, deg={self.degeneracy}{name_part}"
        return (
            f"q={self.charge:+2}, [c]={self.fixed_concentration},"
            f" deg={self.degeneracy}{name_part}"
        )
