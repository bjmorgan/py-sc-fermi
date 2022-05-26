
py-sc-fermi
=======================================

:py:mod:`py-sc-fermi` is an open-source Python package for calculating the concentration of point defects in (semiconducting) crystalline materials.
The required inputs are the volume, density of states of the bulk material, and the formation energies and degeneracies of the point defects. 
The outputs include the self consistent Fermi energy, defect transition levels, and concentrations of the point defects, electrons and holes at a given temperature.
:py:mod:`py-sc-fermi` uses a numerical method to solve for the self-consistent Fermi level in a material, necessary for accurately quantifing the populations of point defects in such materials. 

The approach used in this code was initially based off the algorithm used by the FORTRAN code `SC-Fermi <https://github.com/jbuckeridge/sc-fermi>`_, as described in this `Paper <https://www.sciencedirect.com/science/article/pii/S0010465519302048>`_.

.. image:: figures/outline.jpg
   :width: 200px
   :height: 100px
   :scale: 50 %
   :alt: alternate text
   :align: right

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Papers that use :py:mod:`py-sc-fermi`
-------------------------------------
- `10.1021/acs.inorgchem.1c00278 <https://pubs.acs.org/doi/abs/10.1021/acs.inorgchem.1c00278>`_

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
