from __future__ import annotations

from collections.abc import Mapping
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
    nests mappings, so no field-derived hash can be computed.

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
        high_occupancy_species: species whose solved site occupancy exceeds the
            originating system's ``occupancy_warning_threshold``, mapped to that
            occupancy fraction (0.0-1.0) and ordered worst (highest occupancy)
            first. Empty when no species exceeds the threshold, or when the
            threshold is ``None`` (the dilute-limit check is disabled). This is
            the same verdict the ``DiluteLimitWarning`` reports, exposed as data
            so it is always inspectable regardless of the warnings filter -- a
            solve run under suppressed warnings still records it here.
    """

    temperature: float
    fermi_energy: float
    volume: float
    label: str | None
    p0_per_cell: float
    n0_per_cell: float
    charge_state_concentrations_per_cell: Mapping[str, Mapping[str, float]]
    high_occupancy_species: Mapping[str, float]

    # @dataclass synthesises __hash__ onto this class; setting __hash__ = None
    # keeps instances unhashable, as they hold dict/mapping fields that a
    # synthesised hash would raise over at call time.
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

    def __str__(self) -> str:
        """Return a human-readable report of the solved state: Fermi energy,
        carrier concentrations, and per-species and per-charge-state defect
        concentrations (cm^-3) with each charge state's share of its species
        total. Describes the solution only; ``repr()`` of the originating
        ``DefectSystem`` describes the system.
        """
        cm3 = "cm\u207b\u00b3"  # cm to the minus three
        lines = [
            f"SC Fermi energy:  {self.fermi_energy:.6f} eV",
            "",
            "Carriers:",
            f"  p0 (holes):     {self.p0:.3e} {cm3}",
            f"  n0 (electrons): {self.n0:.3e} {cm3}",
            "",
            "Defect concentrations:",
        ]
        charge_state_concs = self.charge_state_concentrations
        for species, states in charge_state_concs.items():
            species_total = sum(states.values())
            lines.append(f"  {species:<16s}  {species_total:.3e} {cm3}")
            for name, conc in states.items():
                pct = 100 * conc / species_total if species_total > 0 else 0.0
                lines.append(f"    {name:<6s}   {conc:.3e} {cm3}  ({pct:6.2f}%)")
        return "\n".join(lines)
