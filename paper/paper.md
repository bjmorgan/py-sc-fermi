---
title: 'py-sc-fermi: self-consitent Fermi energies in functional materials'
tags:
  - Python
  - materials modelling
  - materials physics
  - materials chemistry
  - thermodynamics
authors:
  - name: Alexander G. Squires 
    orcid: 0000-0001-6967-3690
    affiliation: "1, 4" # (Multiple affiliations must be quoted)
  - name: David O. Scanlon
    orcid: 0000-0001-9174-8601
    affiliation: "1, 3, 4"
  - name: Benjamin J. Morgan 
    orcid: 0000-0002-3056-8233
    affiliation: "2, 4"
affiliations:
 - name: Department of Chemistry, Universtiy College London, 20 Gordon Street, London WC1H 0AJ, UK
   index: 1
 - name: Deparment of Chemistry, University of Bath, Bath BA2 7AY, UK
   index: 2
 - name: Thomas Young Centre, University College London, Gower Street, London WC1E 6BT, UK
   index: 3
 - name: The Faraday Instiution, Didcot OX11 ORA, UK
   index: 4
date: 20 oct 20
bibliography: paper.bib

---

# Summary

`py-sc-fermi` is a Python package for calculating point defect concentrations in crystalline materials under thermodynamic equilibrium. 

Point defects are atomic-scale impefections that influence the properties of all functional materials: they participate in electronic and optical events critical to solar energy conversion [@TLC]; diffusion charged defects underpins modern electrochemical energy storage systems [@batteries] and they control the extent to which the conductivities of thermoelectrics and transparent conductors can be increased via doping [@thermoelectrics,@TCOs]. Attempts to quantify the concentrations of these defect species has become a common practice in materials modelling community in an attempt to desgin novel, highly efficient electronic materials and to rationlise and optimise the properties of known materials [@LLZO-elect].

`py-sc-fermi` provides a numerical approach for calculating defect concentrations under the condition that the removal and addition on ions throughout the material maintain overall net charge neutralility: a "self-consistent Fermi energy" approach. The inputs for the code are formation energies of all point defects in the system (typically obtained from a density funtional theory study) and the electronic structure (an electronic density of states, including the band gap and energy of the valence band maximum). The temperature under which the defects form can then be treated as a free parameter. 

`py-sc-fermi` also allows the user to implement constraints on the concentrations of defects, for example, if we wished to simulate a system in which defects formed at high temperature, but then the system was cooled such that kinetic barriers were too high for full requilibration of some or all defect concentrations to change, the system can then be resampled at lower temperatures allowing for estimating the room temperature conductivity of a material that was sysnthesied under much higher temperatures [@LLZO-elect]. In addition, `py-sc-fermi` is fully "agnositic" towards the software used to generate the input data, and as such, the user is not limitied in which materials modelling code is used to generate the input data.* `py-sc-fermi` strives to be flexible enough to incorporate into any existing point-defect modelling workflow, while remaining user and developer friendly to promote both a wide user base, and ease of extension.

 *There is some minor convenience functionaility for users of the Vienna Ab initio Software Package (VASP).

# Statement of need

`py-sc-fermi` is a Python package for determining the self-consistent Fermi energy of a material from knowledge of the point defect
energetics, allowing for quantification of point defect populations. While we aware of other scientific software that allows for these and related calculations [@Neilson2022-cj,Arrigoni2021-oc,@Ogawa2022-sn,@Buckeridge2019-fm], to our knowledge `py-sc-fermi` is unique amongst these in providing each of the following features:

  - built on a flexible Python API which allows for rapid prototyping and convergence testing with respect to the calculated self-consistent Fermi energy
  - is not part of a larger "point-defects workflow" package, making it as flexible as possible for the end user
  - agnostic towards the choice of simulation code used to generate the input data
  - possess both a command-line interface and user-friendly API
  - allows for constraints on the concentrations of any combination of defects and specific charge states of that defect, including the simulation of effective dopants
  - fully documented and tested.

# Acknowledgements

The authors are grateful for feature requests and user testing from Seán Kavanagh, Joe Willis, and Jiayi Cen.

# References
