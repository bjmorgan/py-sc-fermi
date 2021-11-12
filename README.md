# py-sc-fermi

![Build Status](https://github.com/bjmorgan/py-sc-fermi/actions/workflows/build.yml/badge.svg)
[![Test Coverage](https://api.codeclimate.com/v1/badges/e2ee22eaa4387f072ce7/test_coverage)](https://codeclimate.com/github/bjmorgan/py-sc-fermi/test_coverage)  

Originally an implementation of the FORTRAN code [SC-FERMI](https://github.com/jbuckeridge/sc-fermi), `py-sc-fermi` calculates self-consistent Fermi energies and defect concentrations under thermodynamic equlibrium given defect formation energies (in ionic crystals). For the theory, see https://doi.org/10.1016/j.cpc.2019.06.017.   

The inputs are (charged) defect formation energies, an (electronic) density of states and a file which describes the bulk crystal structure, such as a `VASP` POSCAR, a `.cif` file etc. Having this data, a `DefectSystem` object can be inititalised, properties of which include the self consistent Fermi energy, defect concentrations and defect transition levels. Basic usage can be found in `examples/example_workflow.ipynb`.

### Citing

If you use py-sc-Fermi in your work, please consider citing 
- this repository (see `cite this repository` in the sidebar)
- the paper associated with the FORTRAN implementation which provides an excellent discussion of both the underlying theory and the self-consitent Fermi energy searching algorithm  

   > J. Buckeridge, Equilibrium point defect and charge carrier concentrations in a material determined through calculation of the self-consistent Fermi energy, Computer Physics      Communications, Volume 244, 2019, Pages 329-342, ISSN 0010-4655, https://doi.org/10.1016/j.cpc.2019.06.017.
