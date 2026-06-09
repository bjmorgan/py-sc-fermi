# Change log

## V2.2.2

### Bug Fixes

- `DefectSpecies` now rejects an empty set of charge states at construction,
  raising a `ValueError` naming the species, rather than constructing
  successfully and failing later with a confusing `IndexError` (e.g. in
  `min_energy_charge_state`).

### Documentation

- Corrected the reversed `DefectChargeState.from_dict` docstring summary,
  which described building a dictionary from an object rather than the
  reverse.
- Tidied the `DefectChargeState.__init__` "no concentration or energy"
  `ValueError` message, removing the embedded indentation and stray newline
  and fixing the grammar.
- Removed the misleading "per unit cell" qualifier from the `degeneracy`
  documentation; degeneracy is a dimensionless per-charge-state factor.
- Corrected the `DefectSpecies.as_dict` and `DefectSpecies.from_dict`
  return-type descriptions, which mislabelled the result as a
  `DefectChargeState` representation.
- Filled in the placeholder `DefectSystem.as_dict` docstring.

## V2.2.1

### Improvements

- When `DefectSystem.get_sc_fermi` cannot bracket a self-consistent Fermi
  level, the raised `RuntimeError` now includes the underlying
  `scipy.optimize.brentq` message and chains it as the cause, rather than
  replacing it with a generic "no solution found" string, so the reason a
  solve failed is visible (#57).

## V2.2.0

### Improvements

- DOS integrations are now more robust to noise / band-edge effects (#55).
  Carrier-density integrations terminate near mid-gap rather than at
  auto-determined VBM/CBM indices, making the result largely independent
  of band-edge index location. Band-edge indices are now determined by
  closest grid point to E=0 (VBM) and E=bandgap (CBM), eliminating
  sensitivity to small noise in the energy grid. Aligns with the
  equivalent change in pymatgen's FermiDos.

### Bug Fixes

- DOS construction now validates that the energy range brackets zero
  (the VBM convention) and that bandgap ≤ max(edos), raising a
  `ValueError` on invalid input (#55).
- DOS construction now rejects negative bandgaps, which previously
  produced silently incorrect carrier concentrations due to overlapping
  hole and electron integration ranges (#56).

### Build & Packaging

- Minimum Python version is now 3.11. Python 3.8, 3.9, and 3.10 are
  no longer supported; users on older Python versions can pin to
  `py-sc-fermi <2.2`.
- Switched from `pymatgen` to `pymatgen-core`. All py-sc-fermi imports
  resolve unchanged; user code does not need updating.
- Dependencies are now declared inline in `pyproject.toml`;
  `requirements.txt` is removed.
- Added `[docs]` and `[dev]` install extras: `pip install ".[docs]"`
  for documentation building, `pip install ".[dev]"` for development.
- Added a PEP 561 `py.typed` marker so downstream type checkers
  pick up py-sc-fermi's type annotations.

### Notes

- All internal type hints have been updated to modern syntax
  (PEP 585 / PEP 604).

## V2.1.0

### Bug Fixes

- Fixed bug where importing py-sc-fermi globally replaced `warnings.showwarning`, silencing all non-RuntimeWarning warnings including user warnings and deprecation warnings. The `CustomWarningManager` class has been removed and numpy overflow is now suppressed locally in functions where it is expected during Fermi energy solving.
- Fixed numerical instability in charge state concentration scaling when using fixed species concentrations. The calculation now uses logsumexp for numerically stable proportions.
- Added validation for invalid degeneracy values in `DefectChargeState`.
- Added validation when fixed charge state concentrations exceed total species concentration.

### Improvements

- Refactored `get_sc_fermi` to use `scipy.optimize.brentq`, providing ~40x faster performance.
- `n_trial_steps` parameter is now deprecated and will be removed in a future version.
- `convergence_tolerance` is now optional and controls Fermi energy precision (passed to scipy's `xtol`).

### Notes

- Results may differ from previous versions at the ~1e-12 eV level due to the new solver. 

## V2.0.0

- Every object now has `as_dict/from_dict` methods to help save work as, e.g., 
  json files.
  - This introduces a breaking change! The functionality previously accessible 
    via `DefectSystem.as_dict()` in versions < 2.0.0 is now accessible via
    `DefectSystem.concentration_dict()`
  - This has also resulted in some changes to the CLI input file format. See
    the docs for new formatting requirements.

## V1.1.0
- New warning handling for overflow errors

## V1.0.0 full as-reviewed release
