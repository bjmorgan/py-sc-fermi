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
  directly or via `from_dict`.
- Removed `DefectChargeState.from_string` and
  `DefectSpecies._from_list_of_strings`, which existed solely to support the
  removed text input-file format.
- Removed the `n_trial_steps` parameter from `DefectSystem`, deprecated since
  v2.1.0. The self-consistent Fermi energy solver uses `scipy.optimize.brentq`
  and no longer takes a maximum-iterations argument; pass
  `convergence_tolerance` to control solver precision instead.
- `DefectSpecies` names within a `DefectSystem` must now be unique;
  constructing a `DefectSystem` with two species sharing a name raises a
  `ValueError`. Names key `defect_species_by_name` and the per-species entries
  in `DefectSystem.result`, so duplicates were already ambiguous -- this now
  fails loudly at construction rather than silently.
- `DefectSystem.result.charge_state_concentrations` keys the inner dict by
  charge-state name rather than charge integer, with one entry per
  `DefectChargeState`. `DefectChargeState.name` defaults to the charge string
  (`"q+2"`, `"q-1"`, `"q+0"`); metastable states sharing a formal charge must
  be given explicit names, and duplicate names within a species raise
  `ValueError` at construction. Previously the per-charge-state read-out was
  keyed by charge integer and summed states sharing a formal charge into a
  single entry. Code that indexed the old dict by charge integer (e.g.
  `d["V_O"][2]`) must use the name (e.g.
  `system.result.charge_state_concentrations["V_O"]["q+2"]`, or
  `["V_O"]["V_O_2+"]` for a named state).
- `DefectSystem`'s physical public attributes (`volume`, `dos`, `temperature`,
  `convergence_tolerance`, `vbm_shift`, `cbm_shift`, `rigid_shift`,
  `defect_species`, `site_pools`, `element_pools`) are now read-only; rebinding
  any of them after construction raises `AttributeError`, enforcing the
  documented immutable-snapshot contract. Construct a new `DefectSystem` (or use
  `DefectSystemFactory.at(...)`) for a different temperature, DOS, or correction
  set. `occupancy_warning_threshold`, a non-physical reporting preference, is
  likewise set only at construction, validated there and read-only thereafter.
- `DefectChargeState.fix_concentration` and `DefectSpecies.fix_concentration`
  are no longer public. A concentration is fixed only at construction: via the
  `fixed_concentration` argument of `DefectChargeState`/`DefectSpecies`, or via
  `DefectSystem(fixed_concentrations=...)` / `DefectSystemFactory.at(...,
  fixed_concentrations=...)`. This removes the last route that could mutate a
  constructed `DefectSystem` in place.
- `DefectSystem.report()`, `concentration_dict()` and
  `charge_state_concentration_dict()` have been removed. A solved system is now
  read through a single cached `DefectSystem.result` (a `DefectSystemResult`):
  `result.concentrations` and `result.charge_state_concentrations` give the
  per-species and per-charge-state concentrations in cm^-3, with
  `concentrations_per_cell` / `charge_state_concentrations_per_cell`
  counterparts; `result.fermi_energy`, `result.p0` and `result.n0` are the
  solved scalars; `print(result)` gives the human-readable report and
  `result.as_dict()` a JSON-safe record. The `decomposed` and `per_volume`
  flags are gone, and `as_dict` uses snake_case keys (`fermi_energy`, not
  `"Fermi Energy"`).

### Improvements

- `DefectChargeState` has a `name` attribute, used as the key in
  `DefectSystem.result.charge_state_concentrations`. It defaults to the charge
  string (`"q+2"`); explicit names are required for metastable states sharing a
  formal charge, and must be unique within a species. Explicit names appear in
  `__repr__` and round-trip through `as_dict` / `from_dict`.
- Added `DefectSpecies.charge_state_by_name`, mirroring
  `DefectSystem.defect_species_by_name`; raises `ValueError` listing the
  available names when no charge state matches.
- `formation_energy_corrections` (`DefectSystem`) and
  `formation_energy_correction_fns` (`DefectSystemFactory`) accept
  `(species_name, charge_state_name)` pairs as keys, alongside
  `DefectChargeState` objects. Referencing the same charge state through two
  keys, keying a fixed-concentration state (which has no formation energy to
  correct), or passing a key that is neither form raises `ValueError`.
- `fixed_concentrations` (`DefectSystem`, and as a `factory.at()` override)
  accepts `(species_name, charge_state_name)` pairs as keys, fixing that
  single charge state at construction, alongside species-name keys for
  species totals. Charge-state-level quenching no longer requires mutating
  a constructed system.
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
- The consistency of fixed concentrations is checked at construction. A
  `ValueError` is raised if a species' individually-fixed charge states are
  inconsistent with its species-level `fixed_concentration` -- they exceed it,
  or, when every charge state is fixed so none can absorb a shortfall, fall
  below it -- or if the total fixed concentration in a site-exclusion group
  exceeds its available sites (an unpooled species' own `nsites`, or a
  `site_pools` entry's shared size). These conditions are all independent of the
  Fermi level; they previously surfaced only at solve time (wrapped as a
  `RuntimeError` about the Fermi-energy search window) or, for a shortfall, were
  silently under-reported. They now fail fast at construction with a clear
  message.
- `DefectChargeState` and `DefectSpecies` validate their `fixed_concentration`
  at construction: a value that is not finite and non-negative raises
  `ValueError`, rather than passing silently into budget sums (where a NaN
  defeats every comparison and a negative silently reduces the total).
- Added band-edge corrections (`vbm_shift`, `cbm_shift`,
  `formation_energy_corrections`, `rigid_shift`) to `DefectSystem`. `vbm_shift`
  and `cbm_shift` are pre-evaluated shifts (in eV) that scissor the DOS at
  construction: the conduction block is rigidly moved so the band gap changes
  by `cbm_shift - vbm_shift` (the VBM is pinned at E=0), which changes the
  carrier concentrations from the Fermi-Dirac integration and the
  self-consistent Fermi level, as well as the effective band gap shown by
  `__repr__`. The new `DOS.scissored(delta_gap)` performs this shift and
  returns a new `DOS`. `formation_energy_corrections` is a
  `dict` of per-charge-state formation-energy corrections, keyed by the
  `DefectChargeState` object or a `(species_name, charge_state_name)` pair,
  so that metastable states sharing a formal charge can be corrected
  independently. If `rigid_shift` is True (the default), the band structure
  and defect levels are assumed to move together as a rigid body, so any
  variable-concentration charge state not covered by
  `formation_energy_corrections` is unchanged; if False, such charge states
  have their formation energy shifted by `-charge * vbm_shift`. The DOS
  scissor and the formation-energy channel are independent. `DefectSystem`
  is an immutable, fixed-temperature snapshot:
  corrections are applied once at construction to copies of `defect_species`
  (and to a private scissored DOS), the caller's objects are never modified,
  and `as_dict()`/`from_dict()` and the solved `result` always agree.
- Added `DefectSystemFactory`, for building `DefectSystem` snapshots at a
  series of temperatures from temperature-dependent `vbm_shift_fn`,
  `cbm_shift_fn` and `formation_energy_correction_fns` (each a function of
  temperature; the latter a `dict` of temperature-dependent formation-energy
  corrections keyed per charge state, e.g. vibrational free-energy
  contributions). `factory.at(T, **overrides)` evaluates these functions at `T`
  and returns an independent `DefectSystem`, e.g. `{T:
  factory.at(T).result.concentrations for T in temperatures}`.
- `DefectSystem` gained a `fixed_concentrations` argument: a mapping of species
  name -- or `(species_name, charge_state_name)` pair, fixing that single
  charge state -- to a fixed concentration per unit cell, applied by
  name to the constructor's own copies of `defect_species` (overriding any
  species-level `fixed_concentration` they were constructed with, and composing
  with `formation_energy_corrections`, keyed by object or by name pair). This
  makes the anneal-and-quench override documented on `DefectSystemFactory.at`
  work: freeze some species' totals at their high-temperature values and
  re-solve at a lower temperature, e.g.
  `factory.at(T_low, fixed_concentrations={n: high[n] for n in frozen})`.
- `DefectSystem` and `DefectSystemFactory` now emit a `DiluteLimitWarning`
  (a dedicated warning subclass), at most once per system, when a defect
  species' solved site occupancy exceeds `occupancy_warning_threshold`
  (default 0.01, i.e. 1%; pass `None` to disable, or a finite fraction in
  (0, 1] to retune). py-sc-fermi models dilute, non-interacting defects, so a
  high occupancy flags a regime in which un-modelled defect-defect interactions
  may make the results non-physical. The threshold is a reporting preference
  and is not serialised.
- Added a `PyScFermiWarning` base class so all py-sc-fermi warnings can be
  filtered as a group. `DiluteLimitWarning` now subclasses it, and
  `DefectChargeState.from_dict` emits the new `UnrecognisedKeyWarning` for
  ignored keys rather than a bare `UserWarning`. All remain `UserWarning`
  subclasses, so existing filters are unaffected.
- Added `DefectSystem.element_chemical_potential_shifts`, which reports the
  solved chemical-potential shift (`delta_mu`, in eV) of each element
  constrained by `element_pools`, evaluated at the self-consistent Fermi
  level and relative to the reference at which the formation energies were
  defined. A target above an element's unconstrained content gives a positive
  shift, below it a negative shift, and equal to it a near-zero shift. The
  shift is re-derived from the same element-pool solve used for the
  concentrations, so it is consistent with the concentrations in `result`; an
  element driven to the complete-exclusion limit is reported as `-inf`.
  Intended for external stability-region checks (py-sc-fermi has no
  competing-phase data).
- Added `DefectSystem.result`, a cached `DefectSystemResult` holding the solved
  self-consistent Fermi level, carrier concentrations, and per-species and
  per-charge-state defect concentrations. Concentrations are exposed in cm^-3
  (the default) and per unit cell, with per-species totals derived from the
  per-charge-state breakdown so the two cannot disagree; `str(result)` is the
  human-readable report and `result.as_dict()` a JSON-safe record. The system
  solves once on first access and caches the result (it is an immutable
  snapshot). `DefectSystem` and `DefectSystemFactory` also gained an optional
  `label`, a human-readable tag copied into `result` and serialised when set.
- `DefectSystem` and the pool classes now return read-only views from their
  container accessors, so the immutable-snapshot contract holds at the
  container level (previously some accessors leaked live-mutable internals
  and others returned defensive copies). `DefectSystem.defect_species` and
  `SitePool.species` return tuples; `DefectSystem.site_pools` /
  `element_pools` and `ElementPool.members` return read-only mappings. All
  documented use is read-only, so this is non-breaking; to obtain a mutable
  copy, wrap explicitly (`list(...)`, `dict(...)`).
- `DefectSystem` now validates `volume` and `temperature` at construction,
  raising `ValueError` (naming the offending parameter and its value) when
  either is not finite and > 0. `DefectSystemFactory` applies the same check to
  `volume` at build time and to the `temperature` passed to `at()`, before any
  temperature-dependent shift function is evaluated. Previously these bad
  inputs failed silently or opaquely, far from the construction call. A negative
  `temperature` inverted the Boltzmann factors, giving a converged-looking but
  wrong Fermi level and defect concentrations; a negative `volume` left the
  Fermi level and per-cell concentrations correct but flipped the sign of the
  reported cm^-3 concentrations. A zero `temperature` divided by zero in the
  solve, and a zero `volume` divided by zero when a cm^-3 concentration was
  read. These now fail immediately at the boundary, matching the existing
  checks on the DOS, `DefectChargeState`, and `SitePool`.

### Bug Fixes

- `DefectSystem.site_percentages` is now computed from the solved,
  site-exclusion- and pool-aware concentrations (`_global_defect_concs` at
  the self-consistent Fermi level) rather than the unbounded dilute
  Boltzmann expression. A near-saturation or pooled species previously reported
  occupancies far above 100%, contradicting the concentrations in `result`;
  each species' occupancy is now divided by the sites available to it (its own
  `nsites`, or the shared `site_pools` size for a pooled species), so every
  reported occupancy is at most 100%.
- `DefectSpecies.tl_profile` and `DefectSystem.get_transition_levels` now ignore
  a species' fixed-concentration charge states, which carry no formation energy
  and so define no transition level; previously such a state could raise
  `KeyError`. Calling `DefectSpecies.get_transition_level_and_energy` directly
  with a non-variable charge now raises a descriptive `ValueError` rather than a
  bare `KeyError`.

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
