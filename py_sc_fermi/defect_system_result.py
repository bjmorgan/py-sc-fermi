from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class DefectSystemResult:
    """The solved state of a ``DefectSystem``: Fermi level, carrier and defect
    concentrations at the self-consistent Fermi energy.

    Constructed by ``DefectSystem.result``; users normally read it rather than
    build it. Canonical concentrations are stored per unit cell; the cm^-3 views
    (the common default) are derived using ``volume``. Per-species totals are
    derived from the per-charge-state breakdown, so the two cannot disagree.

    Supports value equality, but is unhashable: ``charge_state_concentrations_per_cell``
    nests plain ``dict``s, so no field-derived hash can be computed.

    Attributes:
        temperature: temperature of the solve, in K.
        fermi_energy: the self-consistent Fermi energy, in eV.
        volume: volume of the unit cell in Angstroms cubed, used for the cm^-3
            views.
        label: an optional human tag copied from the ``DefectSystem``.
        p0_per_cell: hole concentration per unit cell.
        n0_per_cell: electron concentration per unit cell.
        charge_state_concentrations_per_cell: per-unit-cell concentrations keyed
            by species name then charge-state name.
    """

    temperature: float
    fermi_energy: float
    volume: float
    label: str | None
    p0_per_cell: float
    n0_per_cell: float
    charge_state_concentrations_per_cell: dict[str, dict[str, float]]

    # Left unhashable rather than inheriting a dataclass-synthesised __hash__
    # that would raise at call time over the nested dict field.
    __hash__ = None  # type: ignore[assignment]

    @property
    def _scale(self) -> float:
        """Multiplier converting a per-unit-cell concentration to cm^-3."""
        return 1e24 / self.volume

    @property
    def p0(self) -> float:
        """Hole concentration in cm^-3."""
        return self.p0_per_cell * self._scale

    @property
    def n0(self) -> float:
        """Electron concentration in cm^-3."""
        return self.n0_per_cell * self._scale

    @property
    def charge_state_concentrations(self) -> dict[str, dict[str, float]]:
        """Per-charge-state concentrations in cm^-3, keyed by species name then
        charge-state name."""
        scale = self._scale
        return {
            species: {name: conc * scale for name, conc in states.items()}
            for species, states in self.charge_state_concentrations_per_cell.items()
        }

    @property
    def concentrations_per_cell(self) -> dict[str, float]:
        """Per-species total concentrations per unit cell (summed over charge
        states)."""
        return {
            species: sum(states.values())
            for species, states in self.charge_state_concentrations_per_cell.items()
        }

    @property
    def concentrations(self) -> dict[str, float]:
        """Per-species total concentrations in cm^-3 (summed over charge
        states)."""
        return {
            species: sum(states.values())
            for species, states in self.charge_state_concentrations.items()
        }

    def as_dict(self) -> dict:
        """Return a structured, JSON-safe record of the solved state.

        Concentrations are in cm^-3; ``volume`` is included so per-cell values
        are recoverable. Keys are snake_case. This is a one-way export: there is
        no ``from_dict`` (a result is derived by solving a ``DefectSystem``, not
        reconstructed).

        Returns:
            dict: ``fermi_energy`` (eV), ``p0`` and ``n0`` (cm^-3),
            ``temperature`` (K), ``volume`` (cubic Angstrom), ``label``,
            ``concentrations`` and ``charge_state_concentrations`` (per-species
            totals and per-charge-state, in cm^-3).
        """
        return {
            "fermi_energy": self.fermi_energy,
            "p0": self.p0,
            "n0": self.n0,
            "temperature": self.temperature,
            "volume": self.volume,
            "label": self.label,
            "concentrations": self.concentrations,
            "charge_state_concentrations": self.charge_state_concentrations,
        }

    def __repr__(self) -> str:
        label = f", label={self.label!r}" if self.label is not None else ""
        return (
            f"DefectSystemResult(fermi_energy={self.fermi_energy:.4f} eV, "
            f"temperature={self.temperature:g} K{label})"
        )
