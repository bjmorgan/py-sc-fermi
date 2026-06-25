Migrating to v3
===============

py-sc-fermi v3 is a major release. It adds site and element pools, a
fixed-temperature snapshot model with ``DefectSystemFactory``, band-edge
corrections, and support for metastable charge states, and it removes the
legacy FORTRAN SC-Fermi text-input and command-line compatibility layer. This
page covers the changes that affect existing code and explains why solved
concentrations can differ from v2. The complete list of changes is in the
project changelog.


Why are my concentrations different from v2?
--------------------------------------------

The default statistical model has changed. v2 and the original SC-Fermi code
use dilute Boltzmann statistics; v3 uses site-exclusion (Langmuir) statistics
by default. Each defect species has a finite number of sites per unit cell
(``nsites``), and its charge states now share that budget::

    c_i = N_free * w_i / (1 + sum_j w_j)

where ``w_i`` is the Boltzmann weight of charge state ``i`` and ``N_free`` is
the species' available sites. In the dilute limit, where the total weight
``sum_j w_j`` is much smaller than one, this reduces to the familiar Boltzmann
expression ``c_i = N_free * w_i``, so concentrations for dilute defects are
unchanged.

Results therefore differ from v2 only when a defect approaches its site budget,
at low formation energies, high temperatures, or otherwise high occupancy.
There, site-exclusion caps the concentration at the available sites, giving a
lower and more physical population than the unbounded dilute expression. If a
v2 calculation placed a defect at a substantial fraction of its sites, expect
v3 to report a lower concentration for it.

The difference grows smoothly with occupancy, with no threshold below which the
two models agree. py-sc-fermi raises a ``DiluteLimitWarning`` once a solved
occupancy exceeds ``DefectSystem.occupancy_warning_threshold`` (1% by default,
and adjustable): a prompt to review the result and check whether the
dilute / lattice-gas defect model holds for your system.


Breaking API changes
---------------------

Charge states are lists
~~~~~~~~~~~~~~~~~~~~~~~~~

``DefectSpecies.charge_states`` is now a ``list[DefectChargeState]`` rather than
a ``dict`` keyed by charge, and ``charge_state_concentrations`` returns a list
of ``(DefectChargeState, concentration)`` pairs rather than a
``{charge: concentration}`` dict. Construct a species with a list, and iterate
that list directly:

.. code-block:: python

    species = DefectSpecies(
        name="V_O",
        nsites=1,
        charge_states=[
            DefectChargeState(charge=0, energy=0.0, degeneracy=1),
            DefectChargeState(charge=2, energy=1.2, degeneracy=1),
        ],
    )

    for charge_state in species.charge_states:   # was species.charge_states.values()
        ...

Because the charge states are a list, a single species may now hold several
states with the same formal charge, which is how v3 represents metastable
defects.

The legacy SC-Fermi input format has been removed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``py_sc_fermi.inputs`` module (``InputSet``), the ``sc_fermi_solve``
command-line tool, and ``DefectChargeState.from_string`` /
``DefectSpecies._from_list_of_strings`` have been removed, together with support
for reading the FORTRAN SC-Fermi text input files. Build the objects directly in
Python, from a dictionary with ``from_dict``, or from a ``.yaml`` file; the
tutorial shows each route.

``n_trial_steps`` has been removed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The self-consistent Fermi-energy solver uses ``scipy.optimize.brentq`` and no
longer takes a maximum-iterations argument. Remove ``n_trial_steps`` from
``DefectSystem`` construction, and pass ``convergence_tolerance`` to control
solver precision instead.

Defect species names must be unique
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A ``DefectSystem`` raises ``ValueError`` if two ``DefectSpecies`` share a name.
Names key ``concentration_dict`` and ``defect_species_by_name``, so duplicate
names were already ambiguous; this now fails at construction rather than
silently. Give each species a distinct name.

``DefectSystem`` attributes are read-only
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A ``DefectSystem`` is an immutable snapshot at a fixed temperature. Its physical
attributes (``temperature``, ``dos``, ``vbm_shift``, ``defect_species``, the
pools, and the rest) can no longer be reassigned after construction; doing so
raises ``AttributeError``. To compute at a different temperature, DOS, or
correction set, construct a new ``DefectSystem``, or use
``DefectSystemFactory.at(...)``, which builds one snapshot per temperature.
``occupancy_warning_threshold``, a reporting preference rather than part of the
physical system, remains settable.

.. code-block:: python

    # v2: mutate the temperature in place
    system.temperature = 300

    # v3: build a new snapshot, e.g. with the factory
    system = factory.at(300)


New in v3
---------

A few of the new capabilities, with the patterns you are most likely to reach
for when porting existing scripts. See the tutorial for the rest.

Temperature sweeps
~~~~~~~~~~~~~~~~~~~~

``DefectSystemFactory`` builds an independent ``DefectSystem`` snapshot per
temperature, so a sweep is a single comprehension. Temperature-dependent
band-edge shifts and formation-energy corrections are given as functions of
temperature, evaluated at each snapshot:

.. code-block:: python

    from py_sc_fermi.defect_system import DefectSystemFactory

    factory = DefectSystemFactory(
        defect_species=defect_species,
        dos=dos,
        volume=volume,
        vbm_shift_fn=lambda T: ...,   # optional band-edge shift (eV) at temperature T
        cbm_shift_fn=lambda T: ...,
    )

    results = {T: factory.at(T).concentration_dict() for T in temperatures}

Frozen defects (anneal and quench)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To model defect populations set at a high annealing temperature and frozen in
as the material cools, solve at the annealing temperature, take the totals to
hold fixed, and re-solve at the quenched temperature with those totals fixed.
``fixed_concentrations`` holds each named species' total (per unit cell) while
its charge states re-equilibrate:

.. code-block:: python

    factory = DefectSystemFactory(defect_species=defect_species, dos=dos, volume=volume)

    # equilibrium at the annealing temperature
    annealed = factory.at(annealing_temperature)
    totals = annealed.concentration_dict(per_volume=False)   # per unit cell

    # quench: freeze chosen species' totals, re-solve at the lower temperature
    frozen = {name: totals[name] for name in frozen_species}
    quenched = factory.at(quenched_temperature, fixed_concentrations=frozen)

Temperature-dependent band edges need no special handling here: set
``vbm_shift_fn`` and ``cbm_shift_fn`` on the factory (as above) and they are
evaluated at every snapshot, so the quench gets its own band-edge correction
automatically. Pass ``dos=`` to ``at()`` only when the quench needs a genuinely
different density of states, for example a custom electronic structure rather
than a rigid band-edge shift.

Temperature-dependent formation-energy corrections
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``DefectSystemFactory`` accepts ``formation_energy_correction_fns``, a mapping
from a ``DefectChargeState`` to an arbitrary function of temperature. The
function is evaluated at each snapshot and the result is added to that charge
state's formation energy. This is the natural place for corrections whose
magnitude varies with temperature — for example, vibrational free-energy
corrections from phonon calculations:

.. code-block:: python

    # Keep a reference to the charge state so it can be used as a key below.
    # Matching is by object identity, so this must be the same object passed
    # into DefectSpecies. Giving it a name also makes it appear clearly in
    # repr output and concentration results.
    v_o_2plus = DefectChargeState(charge=2, energy=1.2, degeneracy=1, name="V_O_2+")

    v_O = DefectSpecies(
        name="V_O",
        nsites=1,
        charge_states=[
            DefectChargeState(charge=0, energy=0.0, degeneracy=1, name="V_O_0"),
            v_o_2plus,
        ],
    )

    factory = DefectSystemFactory(
        defect_species=[v_O],
        dos=dos,
        volume=volume,
        formation_energy_correction_fns={
            v_o_2plus: lambda T: vibrational_free_energy(T),
        },
    )

    results = {T: factory.at(T).concentration_dict() for T in temperatures}

At each temperature the correction is re-evaluated and baked into the formation
energy before the self-consistent Fermi level is solved. ``DefectSystem``
itself (for a single temperature) accepts the pre-evaluated version as
``formation_energy_corrections: dict[DefectChargeState, float]`` if you do not
need a sweep.

Named charge states also make the decomposed concentration output unambiguous.
For the ``v_O`` species above, ``concentration_dict(decomposed=True)`` returns
``{"V_O_0": ..., "V_O_2+": ...}`` rather than the generated ``"q+0"`` /
``"q+2"`` fallback keys. The same keys appear in
``charge_state_concentration_dict()``.

Also new
~~~~~~~~

- Band-edge corrections (``vbm_shift`` / ``cbm_shift``) that scissor the density
  of states, shifting the carrier gap and the self-consistent Fermi level.
- Metastable charge states, and Boltzmann-weighted effective formation energies
  for charges with more than one form.
- ``DefectSystem.charge_state_concentration_dict(per_volume=True)`` returns
  concentrations as ``{species_name: [(DefectChargeState, conc), ...]}``,
  preserving every charge state as a separate entry. Unlike
  ``concentration_dict(decomposed=True)``, which groups by formal charge and
  sums metastable states that share one, this method lets you read the
  individual concentration of each metastable configuration.
- Site pools (a shared site budget across several species) and element pools (a
  target content for an element, solved through its chemical potential).
