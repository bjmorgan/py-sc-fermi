Usage notes
-------------------------------------
- At different points in the documentation, a "unit cell" is referred to, this is the cell for which the density
  of states data was calculated, and any volumes, degeneracies and numbers of electrons provided should be 
  consitent with this structure for the reported defect concentrations to be accurate. In most cases, inconsistent
  specification of these values will lead to incorrect Fermi energies.
- The reported Fermi energy and transition levels reported are referenced to 0 eV; the code expects that the input
  density of states data is zeroed on the valence band maximum.
- The code operates using the following units,
  and all user-input must be consistent:
  
  - energy: electron volts. 
  
  - temperature: Kelvin.  
  
  - volume: Angstroms :superscript:`3`. 
  
- Concentrations are a special case, internally, the :py:mod:`py-sc-fermi` operates in the concentration of sites in the unit cell 
  which are defective, but will report in cm :superscript:`-3`. The documentation should always specify what kind of concentration
  is expected by a particular function, if it does not, this is a bug! Please report it. 
 
- :py:mod:`py-sc-fermi` calculates defect concentrations from the formation energies it is given, following
  the model of `Zhang and Northrup`_https://doi.org/10.1103/PhysRevLett.67.2339. Two independent mechanisms are
  available for bringing temperature-dependence into those formation energies:

  - :py:class:`~py_sc_fermi.defect_system.DefectSystem` accepts ``vbm_shift``, ``cbm_shift``,
    ``formation_energy_corrections`` and ``rigid_shift`` arguments, for applying temperature-dependent
    band-edge or per-charge-state formation-energy corrections (e.g. from a phonon-renormalisation
    calculation) at a given temperature. :py:class:`~py_sc_fermi.defect_system.DefectSystemFactory` builds a
    series of such ``DefectSystem`` snapshots from temperature-dependent correction functions, for running
    a temperature sweep. See the tutorial for worked examples.

  - :py:meth:`~py_sc_fermi.defect_species.DefectSpecies.effective_formation_energy`,
    :py:meth:`~py_sc_fermi.defect_species.DefectSpecies.get_formation_energies`,
    :py:meth:`~py_sc_fermi.defect_species.DefectSpecies.tl_profile` and
    :py:meth:`~py_sc_fermi.defect_species.DefectSpecies.get_transition_level_and_energy` all accept an optional
    ``temperature`` argument, which Boltzmann-averages over any metastable
    :py:class:`~py_sc_fermi.defect_charge_state.DefectChargeState`\\ s that share a formal charge.
