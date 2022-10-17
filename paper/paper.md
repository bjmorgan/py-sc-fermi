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

`py-sc-fermi` is a Python package for calculating point defect concentrations in crystalline materials under the constraint of net&mdash;charge&ndash;neutrality given the formation energies of all defect species in the system and an electronic density of states both obtained from a set of electronic structure calculations (e.g. density functional theory).

Point defects are atomic-scale imperfections in functional materials that influence energy conversion [@TLC], charge transport [@batteries] and the thermodynamics of their formation controls the limit to which we are able to tune materials properties via synthesis conditions and doping strategies [@thermoelectrics,@TCOs]. Attempts to quantify the concentrations of these defect species has become a common practice in materials modelling community in an attempt to desgin novel, highly efficient electronic materials and to rationlise and optimise the properties of known materials [@LLZO-elect].

The principle of charge conservation tells us that the total electric charge in an isolated system never changes; point defects can introduce a local (integer) charge but these must sum to zero across the full defective system. `py-sc-fermi` provides a numerical approach for calculating point defect populations in functional materials under the condition that the removal and addition of ions throughout the material maintains overall net charge neutrality.

The concentration $c$ of point defect $X$ carrying charge $q$ is given by a Boltzmann distribution,

$$
c[X^q] = n\exp\left(\frac{-E_\mathrm{f}[X^q]}{k_\mathrm{B}T}\right)
$$

where $n$ is the defect's degeneracy (a function of the number of symmetrically equivalent ways defect $X^q$ can form), $E_\mathrm{f}[X^q]$ is the formation energy of the defect and $k_\mathrm{B}$ and $T$ are the Boltzmann constant and temperature respectively. Crucially, $E_\mathrm{f}[X^q]$ is a function of the Fermi energy, and therefore as is $c[X^q]$. The Fermi energy is unknown, but can be solved for self-consistently by observing the condition of charge neutrality,

$$
0 = \sum_{X^𝑞} c[𝑋^𝑞] + n_0 − p_0,
$$

where $n_0$ and $p_0$ are the concentrations of free electrons and holes respectively. $n_0$ and $p_0$ are given by a Fermi-Dirac distribution, which is also a function of the Fermi energy. This gives us a means to solve for the concentrations of electronic charge carriers (holes and electrons) and point defects in which an initial Fermi energy is defined, and then iteratively updated until a self-consistent charge neutral solution is found.

The value of the Fermi energy itself can be used as a general desciptor for the electronic transport properties of the material [@SbTCOs]; the calculated concentrations of electronic charge carriers can be used–in combination with a method to solve for mobility [@amset]–to calculate electronic conductivity, a key figure of merit in many functional materials, and the concentration of the point defects can be used to make inferences about defect processes and the doping response of the material [@LLZO,@BiSI].

 <!-- In other words, the charge contributions of all the charged defects plus any positive holes and negative electrons must sum to zero. The concentration of free electrons and holes are determined by the Fermi-Dirac distribuition, which is a function of the Fermi level.

Likewise the concentration of a charged defect is a function of its formation energy, which is in turn a function of the Fermi level $E_\mathrm{Fermi}$ (the input defect formation energies for `py-sc-Fermi` are given for $E_\mathrm{Fermi} = 0$). It is possible to construct simultaneous equations using the Fermi-Dirac distribution and the formation energy of a defect, adjusting the Fermi energy until the charge neutrality condition is satisfied, also known as a "self consistent Fermi energy" approach. An excellent discussion of the theory is avaiable in the paper published alongside the 
FORTRAN code that formed the intitial inspiration for `py-sc-fermi`, `SC-FERMI` [@Buckeridge2019-fm]. -->

# Statement of need

`py-sc-fermi` is a Python package for determining the self-consistent Fermi energy of a material from knowledge of the point defect
energetics, allowing for quantification of point defect concentrations. While we aware of other scientific software that allows for these and related calculations [@Neilson2022-cj,Arrigoni2021-oc,@Ogawa2022-sn,@Buckeridge2019-fm], to our knowledge `py-sc-fermi` is unique amongst these in providing all of the following features:

  - built on a flexible Python API which allows for rapid prototyping and convergence testing with respect to the calculated self-consistent Fermi energy
  - is not part of a larger "point-defects workflow" package, making it as flexible as possible for the end user
  - agnostic towards the choice of simulation code used to generate the input data
  - possess both command-line functionality for those with little python experience and an object-oriented API that allows the user to construct objects representing individual defects or complete systems of defects directly, allowing easy use with input data generated from any electronic structure simulation code.* 
  - allows for constraints on the concentrations of any combination of defects and specific charge states of that defect, including the simulation of the influence of effective dopants.
  - fully documented and unit-tested.

The code has already been used in a number of studies including some focussed on rationalising the 
properties of known materials [@LLZO-elect,@Squires2021-je] and also predicting the properties of novel materials [@SbTCOs], in addition, the code can also assist with the visualisation of defect energetics in a flexible manner [@Fe2O3].

One feature in particular we would like to draw attention to is the ability to arbitrarily fix the concentrations of specific defects (or individual defect charge states). For example, if we wished to simulate a system in which defects formed at high temperature, but then the system was cooled such that kinetic barriers were too high for full re-equilibration of some or all defect concentrations to change, the system can then be resampled at lower temperatures allowing for estimating the room temperature conductivity of a material that was synthesised under much higher temperatures [@LLZO-elect]. 




 *There is some minor convenience functionality for users of the Vienna Ab initio Software Package (VASP).



# Acknowledgements

The authors are grateful for feature requests and user testing from Seán Kavanagh, Joe Willis, and Jiayi Cen.

# References
