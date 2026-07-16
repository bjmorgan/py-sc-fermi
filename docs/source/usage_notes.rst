Usage notes
===========

Conventions and units. For the model behind the calculation (the statistics,
the charge-neutrality solve, and the pool constraints), see :doc:`theory`.

- "Unit cell" throughout means the cell for which the density-of-states data
  was calculated. Volumes, degeneracies, and electron counts must all be
  consistent with that same cell, or the reported Fermi energy and
  concentrations will be wrong.
- The reported Fermi energy and transition levels are referenced to 0 eV at the
  valence-band maximum; the input density of states must be zeroed there.
- All user input must use these units:

  - energy: electron volts
  - temperature: Kelvin
  - volume: Angstroms\ :superscript:`3`

- Concentrations are reported in cm\ :superscript:`-3`. Internally py-sc-fermi
  works in per-unit-cell counts, and every concentration on
  ``DefectSystem.result`` has a ``_per_cell`` counterpart. Each function's
  documentation states which it expects or returns; if one does not, that is a
  bug. Please report it.
