# py-sc-fermi

![Build Status](https://github.com/bjmorgan/py-sc-fermi/actions/workflows/build.yml/badge.svg)
[![Test Coverage](https://api.codeclimate.com/v1/badges/e2ee22eaa4387f072ce7/test_coverage)](https://codeclimate.com/github/bjmorgan/py-sc-fermi/test_coverage)
[![Documentation Status](https://readthedocs.org/projects/py-sc-fermi/badge/?version=latest)](https://py-sc-fermi.readthedocs.io/en/latest/?badge=latest)
      
`py-sc-fermi` is a materials modelling code for calculating self-consistent Fermi energies and defect concentrations under thermodynamic equlibrium given defect formation energies. For the theory, see [this paper](https://doi.org/10.1016/j.cpc.2019.06.017).   

The necessary inputs are (charged) defect formation energies, an (electronic) density of states and the volume of the unit cell. Having this data, a `DefectSystem` object can be inititalised, properties of which include the self consistent Fermi energy, defect concentrations and defect transition levels. 

Documentation and usage guides can be found [here](https://py-sc-fermi.readthedocs.io/en/latest/).

### Citing

If you use py-sc-Fermi in your work, please consider citing the following: 
- this repository (see `cite this repository` in the sidebar)
- the paper associated with the FORTRAN code [`SC-Fermi`](https://github.com/jbuckeridge/sc-fermi) on which this code was initially based, and provides an excellent discussion of both the underlying theory and the self-consistent Fermi-energy searching algorithm  

   > J. Buckeridge, Equilibrium point defect and charge carrier concentrations in a material determined through calculation of the self-consistent Fermi energy, Computer Physics      Communications, Volume 244, 2019, Pages 329-342, ISSN 0010-4655, https://doi.org/10.1016/j.cpc.2019.06.017.
