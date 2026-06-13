# Change log

## V3.0.0

### Breaking Changes

- `DefectSpecies.charge_states` is now a `list[DefectChargeState]` rather than
  a `dict[int, DefectChargeState]`, and `charge_state_concentrations` now
  returns `list[tuple[DefectChargeState, float]]` rather than
  `dict[int, float]`. This allows a `DefectSpecies` to contain multiple
  `DefectChargeState` objects with the same formal charge, to support
  metastable defects.
- Removed `py_sc_fermi.inputs` (`InputSet`) and the `sc_fermi_solve` CLI tool,
  dropping support for reading the legacy SC-Fermi text input-file format.
  `DefectSystem`/`DefectSpecies`/`DefectChargeState` should now be constructed
  directly, via `from_dict`, or from a `.yaml` file.
- Removed `DefectChargeState.from_string` and
  `DefectSpecies._from_list_of_strings`, which existed solely to support the
  removed text input-file format.
- Removed the `n_trial_steps` parameter from `DefectSystem`, deprecated since
  v2.1.0. The self-consistent Fermi energy solver uses `scipy.optimize.brentq`
  and no longer takes a maximum-iterations argument; pass
  `convergence_tolerance` to control solver precision instead.
- `DefectSpecies` names within a `DefectSystem` must now be unique;
  constructing a `DefectSystem` with two species sharing a name raises a
  `ValueError`. Names key `concentration_dict` and `defect_species_by_name`, so
  duplicates were already ambiguous -- this now fails loudly at construction
  rather than silently.

### Improvements

- `DefectSpecies.get_formation_energies`, `get_transition_level_and_energy`,
  and `tl_profile` now correctly support multiple `DefectChargeState`s sharing
  a formal charge (metastable defects, see above): each charge is represented
  by its lowest-energy state by default. These methods also gained an optional
  `temperature` argument; when set, charges with multiple forms instead use
  the Boltzmann-weighted effective formation energy
  `-kT * ln(sum_i g_i * exp(-E_i / kT))` across those forms.
- Added `DefectSpecies.effective_formation_energy`, the Boltzmann-weighted
  formation energy summed over *all* of a species' charge states and
  metastable forms, `-kT * ln(sum_i g_i * exp(-(E_i + q_i*E_F)/kT))`. At
  `temperature=0` this reproduces the standard "lower envelope"
  formation-energy-vs-Fermi-energy curve at any Fermi energy (not just the
  kinks returned by `tl_profile`); at `temperature>0` it gives the smooth
  curve `-kT * ln(c_total/nsites)` implied by the species' total
  concentration.
- Site-exclusion (Langmuir) statistics, `c_i = N_free * w_i / (1 + sum_j w_j)`,
  are now the default for every `DefectSpecies`, so every charge state's
  concentration is capped by its species' `nsites` (the dilute Boltzmann
  formula is recovered automatically when `sum_j w_j << 1`). Added site pools
  and element pools (`DefectSystem(site_pools=..., element_pools=...)`):
  `site_pools` merges several `DefectSpecies` into one shared budget of sites
  (`N_free` becomes the pool size rather than each species' own `nsites`), and
  `element_pools` constrains the total content of a dopant element across one
  or more species to a fixed target, by solving for the element's chemical
  potential. Multiple `element_pools` are solved jointly via
  `scipy.optimize.minimize`, so targets that couple through a shared species
  are satisfied simultaneously. Note: within a shared `site_pools` group, only
  each charge state's `degeneracy` (not its species' `nsites`) sets its
  relative weight.
- `DefectSystem.as_dict`/`from_dict` now serialise `site_pools` and
  `element_pools`, so a system with site or element constraints round-trips
  through JSON and YAML (previously the pools were silently dropped). Each pool
  is written as a self-describing mapping with species referenced by name.
  `DOS.as_dict` now emits native Python floats, so a `DefectSystem` dictionary
  is YAML-safe.
- `site_pools` and `element_pools` references are validated at construction: a
  `ValueError` is raised for a reference to a species not in `defect_species`,
  a species listed more than once within a pool, or a species placed in more
  than one site pool. References may be given as `DefectSpecies` objects or
  names interchangeably (both are normalised to names internally).
  `DefectSystem.defect_species_by_name` now raises a `ValueError` listing the
  available names, rather than an `IndexError`, when given an unknown name.
- Added band-edge corrections (`vbm_shift`, `cbm_shift`,
  `formation_energy_corrections`, `rigid_shift`) to `DefectSystem`. `vbm_shift`
  and `cbm_shift` are pre-evaluated shifts (in eV) used to compute the
  effective band gap shown by `__repr__`/`report`; `formation_energy_corrections`
  is a `dict[DefectChargeState, float]` of per-charge-state formation-energy
  corrections, keyed by the `DefectChargeState` objects themselves so that
  metastable states sharing a formal charge can be corrected independently. If
  `rigid_shift` is True (the default), the band structure and defect levels
  are assumed to move together as a rigid body, so `vbm_shift`/`cbm_shift`
  only affect the displayed band gap and any variable-concentration charge
  state not covered by `formation_energy_corrections` is unchanged; if False,
  such charge states have their formation energy shifted by
  `-charge * vbm_shift`. `DefectSystem` is an immutable, fixed-temperature snapshot:
  corrections are applied once at construction to copies of `defect_species`,
  the caller's objects are never modified, and `report()`/`as_dict()`/
  `from_dict()` always agree.
- Added `DefectSystemFactory`, for building `DefectSystem` snapshots at a
  series of temperatures from temperature-dependent `vbm_shift_fn`,
  `cbm_shift_fn` and `formation_energy_correction_fns` (each a function of
  temperature, the latter keyed by `DefectChargeState`). `factory.at(T,
  **overrides)` evaluates these functions at `T` and returns an independent
  `DefectSystem`, e.g. `{T: factory.at(T).concentration_dict() for T in
  temperatures}`.

## V2.2.2

### Bug Fixes

- `DefectSpecies` now rejects an empty set of charge states at construction,
  raising a `ValueError` naming the species, rather than constructing
  successfully and failing later with a confusing `IndexError` (e.g. in
  `min_energy_charge_state`).
- `DefectSpecies.min_energy_charge_state` now raises a `ValueError` naming the
  species when there are no variable-concentration charge states (e.g. every
  charge state has a fixed concentration), rather than failing with a bare
  `IndexError`. This also covers the `tl_profile` analysis path, which calls
  it.
- `InputSet.from_yaml` now validates the `fixed_conc_units` argument and raises
  a `ValueError` for any value other than the two allowed units, `"cm^-3"`
  (convert from cm^-3 to per-unit-cell) and `"per_unit_cell"` (no conversion).
  Previously any string other than the exact `"cm^-3"` silently skipped the
  conversion, so a mistyped unit such as `"cm-3"` gave fixed concentrations
  wrong by a factor of roughly 1e24 / volume with no error or warning. Callers
  that relied on an arbitrary string to get the no-conversion path should now
  pass `"per_unit_cell"`. The parameter is also now documented in the method
  docstring.
- `InputSet.from_yaml` now raises a clear `ValueError` when `dos_file` is given
  with an unsupported extension (anything other than `.dat` or `.xml`).
  Previously such a path left the local `dos` variable unbound and failed later
  with a cryptic `UnboundLocalError`.

### Improvements

- Completed static type-hint coverage of the package. Added the missing
  return-type and parameter annotations on the `__repr__` methods,
  `InputSet.from_yaml`, `DefectSystem.from_yaml`, and the CLI entry points,
  and gave the `suppresses_numpy_overflow` decorator a signature-preserving
  `ParamSpec`/`TypeVar` typing so it no longer erases the types of the public
  methods it wraps (`DefectChargeState.get_concentration` and
  `DOS.carrier_concentrations`). These are annotation-only changes with no
  runtime effect.
- Adopted `from __future__ import annotations` in the modules that define
  classmethod constructors, so their return-type forward references (e.g.
  `-> "DefectSpecies"`) are now written unquoted. Annotation-only, with no
  runtime effect.

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
- Added a `DefectSpecies` return-type annotation to `DefectSpecies.from_dict`
  (and the private `_from_list_of_strings`), so the return contract is checked
  by mypy and surfaced to downstream type checkers rather than only described
  in the docstring.
- Filled in the placeholder `DefectSystem.as_dict` docstring.

### Build & Packaging

- Enabled `disallow_untyped_defs` and `disallow_incomplete_defs` in the mypy
  configuration, so missing or incomplete annotations are caught in CI.

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
