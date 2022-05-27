---
title: 'py-sc-fermi: self-consitent Fermi energies in functional materials'
tags:
  - Python
  - materials modelling
  - materials physics
  - materials chemistry
  - thermodynamics
authors:
  - name: Alexander G. Squires # note this makes a footnote saying 'Co-first author'
    orcid: 0000-0000-0000-0000
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Benjamin J. Morgan # note this makes a footnote saying 'Co-first author'
    orcid: 0000-0000-0000-0000
    affiliation: 2
affiliations:
 - name: Department of Chemistry, Universtiy College London, UK
   index: 1
 - name: Deparment of Chemistry, University of Bath, UK
   index: 2
 - name: The Faraday Instiution, UK
date: 20 oct 20
bibliography: paper.bib

---

# Summary

Atomic scale impefections--known as point defects--influence the properties of all functional materials. Point defects participate in electronic and optical events critical to solar energy conversion [@TLC]. Diffusion charged defects underpins modern electrochemical energy storage systems [@batteries], and they control the extent to which the conductivities of thermoelectrics and transparent conductors can be increased via doping [@thermoelectrics,@TCOs]. Attempts to quantify the concentrations of these defect species has become a common 
practice in materials modelling community in an attempt to desgin new, highly efficient electronic materials and to rationlise and optimise the properties of known materials [cite Joe's preprint and all the electronic conductivty in solid electrolyte papers]. `py-sc-fermi` is a Python package which provides a numerical approach for calculating defect concentrations under the condition that the removal and addition on ions throughout the material maintain overall net charge neutralility. This is known as a "self-consistent Fermi energy" approach. The inputs for the code are formation energies of all point defects in the system (these are typically obtained from a density funtional theory study) and the electronic structure. The temperature under which the defects form can then be treated as a free parameter. 

# Statement of need

`py-sc-fermi` is a Python package for determining the self-consistent Fermi energy of a material from knowledge of the point defect
energetics, allowing for quantification of point defect populations. While we aware of other scientific software that allows for these and related calculations, we are not aware of any equivalent packages that are not

    - built in to other, larger "workflow" packages (reducing their overall flexibility)[Spinney,pycdt(?),pylada(?)]  
    - or are built on a flexible Python API which allows for rapid prototyping [SC-Fermi]

`py-sc-fermi` also allows the user to implement constraints on the concentrations of defects, for example, if we wished to simulate a system in which defects formed at high temperature, but then the system was cooled such that kinetic barriers were too high for full requilibration of some or all defect concentrations to change, the system can then be resampled at lower temperatures allowing for estimating the room temperature conductivity of a material that was sysnthesied under much higher temperatures. In addition, `py-sc-fermi` is fully "agnositic" towards the software used to generate the input data, and as such, the user is not limitied in which materials modelling code is used to generate the input data.* `py-sc-fermi` strives to be flexible enough to incorporate into any existing point-defect modelling workflow, while remaining user and developer friendly to promote both a wide user base, and ease of extension.

 *There is some minor convenience functionaility for users of VASP.

# Acknowledgements

The authors are grateful for feature requests and user testing from Seán Kavanagh, Joe Willis, and Jiayi Cen.

# References
