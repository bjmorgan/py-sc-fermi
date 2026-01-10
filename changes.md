# Change log

## V2.0.5

- Fixed global warning suppression caused by importing py-sc-fermi (#50)
- Added custom warning classes (`DOSOverflowWarning`, `DefectOverflowWarning`) 
  that integrate properly with Python's warning system
- Fixed type hints for `DOS._p_func` and `DOS._n_func`

## V2.0.0

- Every object now has `as_dict/from_dict` methods to help save work as e.g. 
  json files.
  - This introduces a breaking change! The functionality previously accessible 
    via `DefectSystem.as_dict()` in versions < 2.0.0 is now accessible via
    `DefectSystem.concentration_dict()`
  - This has also resulted in some changes to the CLI input file format. See
    the docs for new formatting requirements.

## V1.1.0
- New warning handling for overflow errors

## V1.0.0 full as-reviewed release
