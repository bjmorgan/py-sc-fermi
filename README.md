# py-sc-fermi

[![Build Status](https://travis-ci.org/bjmorgan/py-sc-fermi.svg?branch=master)](https://travis-ci.org/bjmorgan/py-sc-fermi)
[![Test Coverage](https://api.codeclimate.com/v1/badges/e2ee22eaa4387f072ce7/test_coverage)](https://codeclimate.com/github/bjmorgan/py-sc-fermi/test_coverage)  

Originally an implementation of the FORTRAN code [SC-FERMI](https://github.com/jbuckeridge/sc-fermi), `py-sc-fermi` calculates self-consistent Fermi energies and concentrations given defect formation energies (in ionic crystals). For the theory, see https://doi.org/10.1016/j.cpc.2019.06.017.   

The inputs are (charged) defect formation energies, an (electronic) density of states and a structural information, such as a `VASP` POSCAR, a `.cif` file etc (TODO: implement generic structure parser). Having this data, a `DefectSystem` object can be inititalised, properties of which include the self consistent Fermi energy, defect concentrations and defect transition levels. Basic usage can be found in `examples/example_notebook.ipynb` (TODO: write up example notebooks).

