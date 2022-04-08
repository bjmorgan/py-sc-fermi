
py-sc-fermi
=======================================

:py:mod:`py-sc-fermi` is an open-source Python package focusessed on calculating the concentrations of point defects in (semiconducting) crystalline materials. Typically the inputs will be comprised of the outputs of DFT codes such as VASP, CASTEP etc. :py:mod:`py-sc-fermi` uses a numerical method to solve for the self-consistent Fermi level in a material, necessary for accurately quantifing the populations of point defects in such materials. 

The approach used in this code was initially based off the algorithm used by the `FORTRAN` code [SC-Fermi](https://github.com/jbuckeridge/sc-fermi), as described in this [paper](https://www.sciencedirect.com/science/article/pii/S0010465519302048). The original intention of the python reimplementation given here was to reimplement the functionalility of SC-Fermi for ease of incoporation into common python materials modelling workflows, and to provide a flexible python API for looping over input parameters and applying various custom constraints to the self-consitent Fermi energy solution.

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
