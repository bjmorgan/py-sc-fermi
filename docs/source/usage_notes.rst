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


What ``energy`` and ``degeneracy`` mean physically
----------------------------------------------------

py-sc-fermi solves the dilute-limit (or site-exclusion) statistics for a set
of charge states you supply; it does not itself compute a formation energy or
a degeneracy. It is worth being precise about what those two numbers are
taken to mean, since supplying the wrong physical quantity gives a
self-consistent-looking but wrong answer.

In the notation of Mosquera-Lois *et al.* (2023) [MosqueraLois2023]_, the
equilibrium concentration of a defect follows a Boltzmann factor built from a
free energy of formation, ``g_f = h_f - T*s_f``, which separates into an
enthalpy and a set of entropy contributions (configurational is handled
separately by the site-exclusion statistics; see :doc:`theory`). Rewriting
the Boltzmann factor as a prefactor times a purely enthalpic exponential::

    c_eq = exp(s_f / k_B) * exp(-h_f / (k_B T)) = (Z_d / Z_b) * exp(-h_f / (k_B T))

gives exactly the two numbers py-sc-fermi asks for:

- ``energy`` is ``h_f``, a formation enthalpy referenced to the
  valence-band maximum (their Eq. 15). ``DefectChargeState.energy`` itself is
  typically the *athermal, static* 0 K value — but it need not stay that way
  through a solve: ``DefectSystem`` accepts a per-charge-state
  ``formation_energy_corrections`` (in eV), and ``DefectSystemFactory.at(T)``
  evaluates temperature-dependent correction functions (plus ``vbm_shift``/
  ``cbm_shift`` for band-edge renormalisation) at each snapshot's
  temperature. That is the intended place to add finite-temperature physics
  you have computed elsewhere — e.g. the vibrational free-energy term
  ``a_vib_f,P(T)`` or a temperature-dependent band gap — so that ``h_f`` at
  each solved temperature is the full (or as-complete-as-you've-modelled)
  formation enthalpy at that temperature, not a fixed 0 K number reused
  everywhere.
- ``degeneracy`` is the pre-exponential ratio ``Z_d / Z_b``, i.e. every
  entropic contribution folded into a single number, entering the
  Boltzmann weight as ``w_i = g_i * exp(-E_i(E_F) / kT)`` (see
  :doc:`theory`). Unlike ``energy``, it has no per-temperature hook: a
  ``DefectChargeState``'s ``degeneracy`` is fixed once and reused at every
  temperature a ``DefectSystem`` is solved at.

py-sc-fermi does not compute any of these corrections from a structure or a
phonon calculation; both the base ``energy``/``degeneracy`` and any
temperature-dependent corrections are inputs you must have already evaluated
(e.g. from DFT total energies, symmetry analysis, or a phonon calculation)
before constructing a ``DefectChargeState`` or ``DefectSystemFactory``. For
``degeneracy`` specifically, `doped <https://github.com/SMTG-Bham/doped>`_
computes the spin and orientational (symmetry) contributions to ``Z_d / Z_b``
automatically from the pristine and relaxed defect structures, alongside its
defect-generation and finite-size-correction workflow, and can be used to
obtain the values to pass in here.

.. [MosqueraLois2023] I. Mosquera-Lois, S. R. Kavanagh, J. Klarbring, K. Tolborg
   and A. Walsh, *Imperfections are not 0 K: free energy of point defects in
   crystals*, Chemical Society Reviews, 52(17), 5812-5826 (2023).
   `doi:10.1039/D3CS00432E <https://doi.org/10.1039/D3CS00432E>`_
