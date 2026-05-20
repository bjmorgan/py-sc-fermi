# Change log

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
