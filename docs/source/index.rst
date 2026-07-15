
py-sc-fermi
=======================================

.. toctree::
   :hidden:
   :maxdepth: 4
   :caption: contents

   installation
   migrating_to_v3
   theory
   usage_notes
   tutorials
   FAQs
   py_sc_fermi

:py:mod:`py-sc-fermi` is an open-source Python package for calculating
equilibrium point-defect and carrier concentrations in crystalline
semiconductors. From the cell volume, the bulk density of states, and the
formation energies and degeneracies of a set of point defects, it solves
self-consistently for the Fermi level that enforces charge neutrality, and
returns the equilibrium defect, electron, and hole concentrations, together
with the defect transition levels, at a given temperature.

.. image:: figures/outline.png
   :width: 800px
   :height: 400px
   :align: center

|

Quickstart
==========

.. code-block:: python

    from py_sc_fermi.dos import DOS
    from py_sc_fermi.defect_charge_state import DefectChargeState
    from py_sc_fermi.defect_species import DefectSpecies
    from py_sc_fermi.defect_system import DefectSystem

    # bulk density of states (the electron count and band gap are read from the vasprun)
    dos = DOS.from_vasprun("vasprun.xml")

    defects = [
        DefectSpecies(
            "V_O",
            nsites=2,
            charge_states=[
                DefectChargeState(charge=0, energy=4.7, degeneracy=1),
                DefectChargeState(charge=2, energy=1.1, degeneracy=1),
            ],
        ),
        # ... one DefectSpecies per defect
    ]

    system = DefectSystem(
        defect_species=defects,
        dos=dos,
        volume=250.0,      # cubic Angstroms
        temperature=1000,  # Kelvin
    )

    result = system.result
    print(result)           # self-consistent Fermi level, carriers, and defect concentrations
    result.fermi_energy     # Fermi level (eV, referenced to the VBM)
    result.concentrations   # {defect name: concentration in cm^-3}

The :doc:`tutorials` page works through this example and the features below in
full.

Features
========

- The self-consistent Fermi level and equilibrium defect, electron, and hole
  concentrations at a fixed temperature.
- Site-exclusion (Langmuir) statistics by default, capping each defect at its
  available sites; the dilute Boltzmann limit is recovered automatically.
- Temperature series and temperature-dependent band-edge shifts through
  ``DefectSystemFactory``.
- Frozen-defect and anneal-and-quench workflows through ``fixed_concentrations``.
- Constraints on the defect populations: a shared budget of sites across
  species (``site_pools``), and a fixed total content of an element solved
  through its chemical potential (``element_pools``).
- Metastable charge states, and Boltzmann-weighted effective formation energies.
- Defect transition-level diagrams.

The :doc:`theory` page sets out the model behind these features: the
site-exclusion statistics, the charge-neutrality solve, and the pool
constraints. Upgrading from v2? The :doc:`migrating_to_v3` guide covers the
breaking changes and the new default statistics.

py-sc-fermi's self-consistent Fermi-level algorithm was originally based on that
of the FORTRAN code `SC-Fermi <https://github.com/jbuckeridge/sc-fermi>`_
(`Buckeridge 2019 <https://www.sciencedirect.com/science/article/pii/S0010465519302048>`_).

Papers that use :py:mod:`py-sc-fermi`
======================================

- Haouari et al., Impact of Solution Chemistry on Growth and Structural Features of Mo-Substituted Spinel Iron Oxides, 2021, `10.1021/acs.inorgchem.1c00278 <https://pubs.acs.org/doi/abs/10.1021/acs.inorgchem.1c00278>`_
- Squires et al., Low Electronic Conductivity of Li\ :sub:`7`\ La\ :sub:`3`\ Zr\ :sub:`2`\ O\ :sub:`12`\  Solid Electrolytes from First Principles, 2022, `10.1103/PhysRevMaterials.6.085401 <https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.6.085401>`_
- Jackson, Parret and Willis et al., Computational Prediction and Experimental Realization of Earth-Abundant Transparent Conducting Oxide Ga-Doped ZnSb\ :sub:`2`\ O\ :sub:`6`\ , 2022, `10.1021/acsenergylett.2c01961 <https://doi.org/10.1021/acsenergylett.2c01961>`_
- Cen et al., Cation disorder dominates the defect chemistry of high-voltage LiMn\ :sub:`1.5`\ Ni \ :sub:`0.5`\ O \ :sub:`4`\ (LMNO) spinel cathodes, 2023, `10.1039/D3TA00532A <https://doi.org/10.1039/D3TA00532A>`_
- Nicolson et al., Cu\ :sub:`2`\SiSe\ :sub:`3`\  as a promising solar absorber: harnessing cation dissimilarity to avoid killer antisites, 2023, `10.1039/D3TA02429F <https://doi.org/10.1039/D3TA02429F>`_


Searching
=========

:ref:`genindex` | :ref:`modindex` | :ref:`search`
