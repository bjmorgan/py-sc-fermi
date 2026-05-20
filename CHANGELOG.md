# Change log

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
