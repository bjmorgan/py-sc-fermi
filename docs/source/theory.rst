The defect-equilibrium model
============================

py-sc-fermi computes the equilibrium concentrations of point defects and free
carriers in a semiconductor at a fixed temperature. Given the formation
energies of a set of defect charge states and the host density of states, it
finds the Fermi level at which the crystal is charge-neutral, and reports the
defect, electron, and hole concentrations there.

This page sets out the model behind that calculation: how a defect's
concentration follows from its formation energy, how charge neutrality fixes
the Fermi level, how site and element budgets enter as constraints, and how
metastable states and temperature-dependent inputs are handled. Throughout, the
formation energies and the density of states are taken as given; computing them
is outside py-sc-fermi's scope.


Defect concentrations and site exclusion
----------------------------------------

A point defect forms on one of a finite set of equivalent sites in the unit
cell. A species declares how many such sites the cell provides through
``nsites``, and each of its charge states carries a formation energy that
varies linearly with the Fermi level::

    E_i(E_F) = E_i(0) + q_i * E_F

with ``q_i`` the charge and ``E_i(0)`` the formation energy at the
valence-band maximum (``E_F = 0``). The tendency of a site to adopt charge
state ``i`` is set by its Boltzmann weight::

    w_i = g_i * exp(-E_i(E_F) / kT)

where ``g_i`` is the state's degeneracy and ``kT`` the thermal energy.

Because the charge states of a species draw on the same finite set of sites,
they compete for occupancy. Each site is either empty or in one of the
competing states, so its occupation follows a single-site partition function
``1 + sum_j w_j``: the ``1`` for the empty site, one ``w_j`` for each state
that can occupy it. The concentration of state ``i`` is then::

    c_i = N_free * w_i / (1 + sum_j w_j)

with ``N_free`` the number of available sites. This is the lattice-gas
(Langmuir) result, and it is py-sc-fermi's default. It caps the occupancy of
the site set at ``N_free``: no defect can exceed the sites that host it.


The dilute limit
~~~~~~~~~~~~~~~~~

When the sites are mostly empty the total weight ``sum_j w_j`` is much smaller
than one, the denominator approaches ``1``, and the concentration reduces to
the familiar dilute expression::

    c_i = N_free * w_i = N_free * g_i * exp(-E_i(E_F) / kT)

This is the defect concentration of the standard dilute model (`Zhang and
Northrup, 1991 <https://doi.org/10.1103/PhysRevLett.67.2339>`_), and the
statistics used by py-sc-fermi v2 and the original SC-Fermi code. Site
exclusion and the dilute limit therefore agree wherever defects are dilute, and
part ways only as a species approaches its site budget (at low formation
energies, high temperatures, or otherwise high occupancy), where site exclusion
holds the concentration below the unbounded dilute value. The difference grows
smoothly with occupancy, vanishing only in the dilute limit.


When site exclusion is not enough
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Site exclusion corrects the *counting* of configurations: it stops a defect
from occupying more sites than exist. It leaves the *energetics* untouched —
every formation energy is still that of an isolated defect, independent of how
many others are present. Once a defect occupies an appreciable fraction of its
sites, that independence is the assumption that fails: neighbouring defects
begin to perturb one another's formation energies, and concentrations computed
from isolated-defect energies become unreliable.

py-sc-fermi cannot model those interactions, but it can flag when they are
likely to matter. When a solve leaves a species above
``occupancy_warning_threshold`` (one per cent of its sites by default), it
raises a ``DiluteLimitWarning`` naming the species and its occupancy. The same
verdict is available without the warning as ``result.high_occupancy_species``,
and the solved occupancy of every species as ``system.site_percentages()``. The
threshold is a prompt to check the model against your system, not a boundary
between valid and invalid results.


Shared site pools
~~~~~~~~~~~~~~~~~~

By default each species has its own budget of ``nsites`` sites, and only its
own charge states compete for them. When several species genuinely share one
set of physical sites, as substitutional species sharing a sublattice do, group
them in a ``site_pools`` entry. Site exclusion is then applied across the group:
the members draw on one shared ``N_free`` and together occupy at most all of it.


The self-consistent Fermi level
-------------------------------

The Fermi level is not an input to the calculation: it is the one unknown that
py-sc-fermi solves for. Every defect concentration depends on it through the
formation energies ``E_i(E_F)``, the free-carrier concentrations depend on it
too, and the value that makes the crystal charge-neutral is what fixes them
all.

The free carriers come from the density of states. Integrating the electronic
DOS against the Fermi-Dirac distribution gives the hole concentration ``p0`` in
the valence band and the electron concentration ``n0`` in the conduction band,
each a function of the Fermi level and temperature. The valence-band maximum
sits at ``E_F = 0`` and the conduction-band minimum at the band gap, so raising
the Fermi level fills electrons and depletes holes.

Charge neutrality balances every positive charge against every negative one::

    p0 + sum(positive defect charge) = n0 + sum(negative defect charge)

the left side collecting the holes and the positively-charged defects, the
right the electrons and the negatively-charged defects, each defect weighted by
its concentration and the magnitude of its charge. Written as a single
residual::

    q_tot(E_F) = [n0 + negative] - [p0 + positive]

this is zero exactly at the self-consistent Fermi level. As the Fermi level
sweeps the DOS energy window the balance tips from hole-rich at the low-energy
end to electron-rich at the high-energy end, so ``q_tot`` changes sign once
between them. py-sc-fermi brackets the window and locates that crossing with
Brent's method. ``convergence_tolerance`` sets the tolerance on the returned
Fermi level; ``DefectSystem.result`` reports it as ``fermi_energy``, referenced
to the valence-band maximum, alongside ``p0``, ``n0``, and the defect
concentrations evaluated there.


Element constraints as chemical potentials
------------------------------------------

By default a charge state's formation energy already fixes everything its
concentration needs: ``E_i(0)`` was defined at some chosen set of atomic
chemical potentials, and the concentration follows directly. This is an
open-system picture — each element is drawn from an implicit reservoir held at a
fixed chemical potential.

Sometimes the constraint runs the other way: not the chemical potential of an
element, but its total amount. You might hold the total oxygen content of the
cell fixed, or pin an exact stoichiometry, and ask what chemical potential that
implies. An ``element_pools`` entry expresses this. It names an element, a
target content per cell, and the member species that carry the element, each
with a signed stoichiometry — positive where the species adds atoms of the
element (an interstitial), negative where it removes them (a vacancy). The
constraint is::

    sum_i stoichiometry_i * c_i = target

py-sc-fermi enforces it by shifting the element's chemical potential by
``delta_mu`` and folding that into each member's weight,
``w_i -> w_i * exp(s_i * delta_mu / kT)``, then solving for the ``delta_mu`` at
which the content meets the target. This is mass action: raising ``delta_mu``
makes the element's interstitials cheaper and its vacancies more costly, feeding
atoms into the cell; lowering it does the reverse. A target of zero with
mixed-sign stoichiometry pins exact stoichiometry, with interstitials balancing
vacancies; a negative target sets an off-stoichiometric deficiency; and a scan
can cross zero continuously.


What the model can and cannot judge
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``element_chemical_potential_shifts`` returns the solved ``delta_mu``, measured
relative to the reference at which the formation energies were defined: a target
above the element's unconstrained content gives ``delta_mu > 0``, one below it
``delta_mu < 0``.

Two limits bound the target, and only one of them is the model's to enforce.
Site exclusion is: a target needing more atoms than there are sites to host
them is unreachable, and py-sc-fermi rejects it. Thermodynamic accessibility is
not: whether the element's chemical potential can actually reach ``delta_mu``
without precipitating a competing phase depends on the phase diagram, which
py-sc-fermi does not have. The reported shift is meant to be checked against an
external stability region — a shift inside it corresponds to an accessible
concentration, one outside it to a supersaturated or otherwise unphysical state
that the equilibrium solve will still report if asked.


Metastable states and transition levels
----------------------------------------

A single formal charge can correspond to more than one atomic configuration,
such as two geometries of a 2+ oxygen vacancy. py-sc-fermi represents these
metastable states as separate ``DefectChargeState`` objects that share a formal
charge and are given distinct names.

In the concentration solve, no special treatment is needed: each metastable
state is an independent competitor in its site-exclusion group, with its own
weight ``w_i``, so its population emerges from the same Langmuir expression as
any other state. The Boltzmann weighting between configurations is not a
separate step; it falls out of the competition.

The combination becomes explicit only when a *single* formation energy is
wanted for a formal charge — to draw a transition-level diagram, or to report
an effective formation energy. There the states sharing a charge are collapsed
into one value::

    F(E_F) = -kT * ln( sum_i g_i * exp(-E_i(E_F) / kT) )

The ``temperature`` argument of ``effective_formation_energy``,
``get_formation_energies``, ``tl_profile`` and
``get_transition_level_and_energy`` controls this collapse. At
``temperature = 0`` (the default) it reduces to the lowest-energy state at each
Fermi level — the usual lower-envelope formation-energy curve. Above zero it is
the smooth Boltzmann-weighted combination, in which a low-lying metastable state
contributes before it becomes the ground state. This argument sets only how
metastable states are combined for reporting; it is independent of the
temperature at which the system is solved.

A transition level is the Fermi energy at which two charges are equal in
formation energy — where the lower envelope switches from one charge to the
next. Between charges ``q1`` and ``q2`` it lies at::

    E_F = (E_q2 - E_q1) / (q1 - q2)

with the formation energies evaluated at the band edge.
``get_transition_levels`` returns these profiles for every species over the DOS
energy range.


The fixed-temperature snapshot model
------------------------------------

A ``DefectSystem`` is a snapshot at one temperature. Everything that defines the
equilibrium — the temperature, the density of states, any band-edge shifts,
formation-energy corrections, and fixed concentrations — is supplied at
construction and cannot be changed afterwards; the public attributes are
read-only. A snapshot therefore has a single well-defined solved state, computed
on first access to ``result`` and cached.

A temperature series is a series of such snapshots. ``DefectSystemFactory``
holds the temperature-dependent inputs (band-edge shifts and formation-energy
corrections as functions of temperature), and ``at(T)`` evaluates them and
builds a fresh, independent ``DefectSystem`` for each temperature. The snapshots
share no state, so a sweep is a straightforward series of independent solves.

Band-edge shifts enter a snapshot through the density of states. ``vbm_shift``
and ``cbm_shift`` scissor the DOS, changing the band gap by
``cbm_shift - vbm_shift`` and so moving the carrier concentrations and the Fermi
level. Optionally, ``vbm_shift`` also shifts the defect formation energies, for
the case where the defect levels are taken to stay fixed in absolute energy as
the bands move rather than riding rigidly with them. Supplied as functions on
the factory, the shifts are re-evaluated at every temperature in a sweep.

Anneal and quench uses the snapshot model directly. Defect populations are often
set at a high growth or annealing temperature and then frozen in as the material
cools faster than the defects can re-equilibrate. To model this, solve one
snapshot at the annealing temperature, take the total concentrations it
predicts, and solve a second snapshot at the lower temperature with those totals
held fixed through ``fixed_concentrations``; each frozen total is held while its
species' charge states re-equilibrate to the low-temperature Fermi level. The
migration guide and the tutorial give the worked idiom.
