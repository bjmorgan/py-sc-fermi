# Change log

## V2.0.0

- every object now has `as_dict/from_dict` methods to help save work as e.g. 
  json files.
  - this introduces a breaking change! The functionality previously accessible 
    via `DefectSystem.as_dict()` in versions < 2.0.0 is now accessible via
    `DefectSystem.concentration_dict()`
  - this has also resulted in some changes to the CLI input file format. See
    the docs for new formatting requirements.

## V1.1.0
- new warning handling for overflow errors

## V1.0.0 full as-reviewed release
