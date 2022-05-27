
py-sc-fermi
=======================================

:py:mod:`py-sc-fermi` is an open-source Python package for calculating the concentration of point defects in (semiconducting) crystalline materials.
The required inputs are the volume, density of states of the bulk material, and the formation energies and degeneracies of the point defects. 
The outputs include the self consistent Fermi energy, defect transition levels, and concentrations of the point defects, electrons and holes at a given temperature.
:py:mod:`py-sc-fermi` uses a numerical method to solve for the self-consistent Fermi level in a material, necessary for accurately quantifing the populations of point defects in such materials. 

The approach used in this code was initially based off the algorithm used by the FORTRAN code `SC-Fermi <https://github.com/jbuckeridge/sc-fermi>`_, as described in this `Paper <https://www.sciencedirect.com/science/article/pii/S0010465519302048>`_.

.. image:: source/figures/outline.png
   :width: 800px
   :height: 400px
   :align: center

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Usage notes
-------------------------------------
- At different points in the documentation, a "unit cell" is referred to, this is the cell for which the density
  of states data was calculated, and any volumes, degeneracies and numbers of electrons provided should be 
  consitent with this structure for the reported defect concentrations to be accurate. In most cases, inconsistent
  specification of these values will lead to incorrect Fermi energies.
- All Fermi energies reported by the code, and provided in inputs are referenced to the VBM, which in turn
  must be referenced to 0 in the input density of states data.
- The code operates using the following units,
  and all user-input must be consistent
  - energy: electron volts
  - temperature: Kelvin
  - volume: Angstroms :superscript:`3`
- Concentrations are a special case, internally, the :py:mod:`py-sc-fermi` operates in the concentration of sites in the unit cell 
  which are defective, but will report in cm :superscript:`-3`. The documentation should always specify what kind of concentration
  is expected by a particular function, if it does not, this is a bug! Please report it. 
  

Papers that use :py:mod:`py-sc-fermi`
-------------------------------------
- `10.1021/acs.inorgchem.1c00278 <https://pubs.acs.org/doi/abs/10.1021/acs.inorgchem.1c00278>`_

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
