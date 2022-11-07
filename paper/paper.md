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

`py-sc-fermi` is a Python package for calculating point defect concentrations in functional materials under the constraint of net&mdash;charge&ndash;neutrality given the formation energies of all defect species in the system and an electronic density of states both obtained from a set of electronic structure calculations (e.g. density functional theory).

Point defects are atomic-scale imperfections in functional materials that influence charge transport [@batteries], energy conversion [@TLC], and optical properties [@paper-on-color-centres] amongst many others. Studying the thermodynamics of point defect formation guides us in understanding the limit to which we are able to tune materials properties via synthesis conditions and doping strategies [@thermoelectrics,@TCOs,@LLZO-orig]. Attempts to quantify the prevalence of different point defects species has become a common practice in materials modelling community in an attempt to design novel, highly efficient electronic materials and to rationlise and optimise the properties of known materials [@LLZO-elect]. 

The main challenge in calculating point defect populations comes from the fact that point defects carry integer charge but the principle of charge conservation, which tells us that the total electric charge in an isolated system never changes, tells us these local integer charge must sum to zero over the full defective system. In other words, all charged defect concentrations are mutually dependent. Stated mathematically,

$$
\begin{equation}
0 = \sum_{X^𝑞} q[𝑋^𝑞] + n_0 − p_0
\end{equation}
$$

where the first term is the sum over all the charge contributions from each defect $X$ in its set of all possible charge states $q$ is the concentration and $n_0$ and $p_0$ are the concentrations of free electrons and holes respectively. All the variables in equation 1 are directly or indirectly functions of the electron chemical potential: the Fermi energy. Under a fixed set of growth conditions, the only unknown variable in the calculation of each term is the Fermi energy, and so the populations of all charged species in the system can be solved for self consistently.

`py-sc-fermi` provides a numerical approach for the self-consistent solution. An initial Fermi energy is guessed, and this is updated over multiple cycles until the value is found which satisfies charge neutrality (within a specified tolerance). The value of the Fermi energy itself can be used as a general descriptor for the electronic transport properties of the material [@SbTCOs]; the calculated concentrations of electronic charge carriers can be used&mdash;in combination with a method to solve for mobility [@amset]&mdash;to calculate electronic conductivity (a key figure of merit in many functional materials) and the concentration of the point defects can be used to make inferences about defect processes and the doping response of the material [@LLZO,@BiSI].

# Statement of need

`py-sc-fermi` is a Python package for determining the self-consistent Fermi energy of a material from knowledge of the point defect
energetics, allowing for quantification of point defect concentrations. While we aware of other scientific software that allows for these and related calculations [@Neilson2022-cj,Arrigoni2021-oc,@Ogawa2022-sn,@Buckeridge2019-fm], to our knowledge `py-sc-fermi` is unique amongst these in providing all of the following features:

  - built on a flexible Python API which allows for rapid prototyping and convergence testing of the solved-for Fermi energy.
  - is not part of a larger "point-defects workflow" package, making it as flexible as possible for the end user
  - agnostic towards the choice of simulation code used to generate the input data*
  - possess both command-line functionality for those with little python experience and an object-oriented API that allows the user to construct objects representing individual defects or complete systems of defects directly.* 
  - allows for constraints on the concentrations of any combination of defects and specific charge states of that defect, including the simulation of the influence of effective dopants.
  - fully documented with good unit-test coverage.

The code has already been used in a number of studies including some focussed on rationalising the 
properties of known materials [@LLZO-elect,@Squires2021-je] and also predicting the properties of novel materials [@SbTCOs], in addition, the code can also assist with the visualisation of defect energetics in a flexible manner [@Fe2O3].

One feature in particular we would like to draw attention to is the ability to arbitrarily fix the concentrations of specific defects (or individual defect charge states). For example, if we wished to simulate a system in which defects formed at high temperature, but then the system was cooled such that kinetic barriers were too high for full re-equilibration of some or all defect concentrations to change, the system can then be resampled at lower temperatures allowing for estimating the room temperature conductivity of a material that was synthesised under much higher temperatures [@LLZO-elect]. 




 *There is some minor convenience functionality for users of the Vienna Ab initio Software Package (VASP).



# Acknowledgements

The authors are grateful for feature requests and user testing from Seán Kavanagh, Joe Willis, and Jiayi Cen, Sabrine, and Lavan.

# References
