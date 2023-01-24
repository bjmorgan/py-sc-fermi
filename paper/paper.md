---
title: 'py-sc-fermi: self-consistent Fermi energies and defect concentrations from electronic structure calculations'
tags:
  - Python
  - materials modelling
  - materials physics
  - materials chemistry
  - thermodynamics
authors:
  - name: Alexander G. Squires 
    orcid: 0000-0001-6967-3690
    affiliation: "1, 4"
  - name: David O. Scanlon
    orcid: 0000-0001-9174-8601
    affiliation: "1, 3, 4"
  - name: Benjamin J. Morgan 
    orcid: 0000-0002-3056-8233
    affiliation: "2, 4"
affiliations:
 - name: Department of Chemistry, University College London, 20 Gordon Street, London WC1H 0AJ, UK
   index: 1
 - name: Department of Chemistry, University of Bath, Bath BA2 7AY, UK
   index: 2
 - name: Thomas Young Centre, University College London, Gower Street, London WC1E 6BT, UK
   index: 3
 - name: The Faraday Institution, Didcot OX11 ORA, UK
   index: 4
date: 20 oct 20
bibliography: paper.bib

---

# Summary

`py-sc-fermi` is a Python package for calculating point defect concentrations in crystalline materials, under the constraint of net&ndash;charge-neutrality and the assumption of thermodynamic (quasi-)equilibrium.
The required inputs are the formation energies of all point defects of interest and the electronic density of states.
These can be obtained from electronic structure calculations, e.g., density functional theory.

Point defects are atomic-scale imperfections in crystalline materials, which can strongly affect a range of material properties, including electronic structure [@Pastor2022-xa], charge transport [@batteries], and energy-conversion processes [@TLC].
Point defect concentrations are generally sensitive to synthesis conditions and to the inclusion of extrinsic dopants.
By understanding the thermodynamics of point-defect formation, it is possible to predict the response of specific material properties to changes in synthesis or doping protocols [@thermoelectric; @Squires2019-is; @Park2018-np].
The quantitative modelling of point-defect populations in functional materials has become common practice in the materials modelling community, often undertaken with the aim of informing the design of novel, highly performant materials [@Jackson2022-la; @Ganose2018-eu; @Toriyama2021-ch] or to rationlise and optimise the properties of known materials [@Kim2018-eg; @Shimoda2022-vs; @LLZO-elect]. 

While the formation energies of charge-neutral point defects can be calculated simply, e.g., from electronic structure calculations, the formation energies of charged defects depend on the electronic chemical potential (Fermi energy), which, in turn, depends on the net population of all charged defects in the system.
The equilibrium concentrations of charged point defects are therefore mutually dependent, and must be solved as a unified self-consistent system.
In practice, self-consistent defect concentrations are calculated by imposing the constraint of net&ndash;charge-neutrality [@Buckeridge2019-fm],
\begin{equation}
0 = \sum_{X,q} q \left[ X,q \right] - n_0 + p_0.
\end{equation}
Here, the first term is a sum over charge contributions from each defect $X$ in each of its accessible charge states $q$.
$n_0$ and $p_0$ are the concentrations of free electrons and electron holes respectively.

`py-sc-fermi` provides a Python implementation of an iterative numerical approach for calculating the self-consistent solution for an arbitrary set of point-defects.
The resulting self-consistent Fermi energy can be used as a general descriptor for the electronic-transport properties of a material [@Jackson2022-la]; e.g., the calculated electronic charge-carrier concentrations can be used to calculate electronic conductivities [@LLZO-elect]; and the point defect concentrations, and how these vary with synthesis conditions or doping protocol, can be used to model how the formation of competing defects affects material properties of interest [@Kavanagh2021-bj; @Ganose2018-eu] as well as quantifiying the doping response of a material [@Squires2019-is; @Squires2021-je].

# Statement of need

`py-sc-fermi` is a Python package for determining the self-consistent Fermi energy of a material from pre-calculated point defect data.
This code allows the quantification of point defect and electronic carrier concentrations.
Other software exists that is written to perform similar calculations [@Neilson2022-cj; @Arrigoni2021-oc; @Ogawa2022-sn; @Buckeridge2019-fm].
To the best of our knowledge `py-sc-fermi` is unique in providing all of the following features:

- It is built on a flexible Python API which allows for rapid prototyping and convergence testing of the solved-for Fermi energy.
- It is not part of a larger &ldquo;point-defects workflow&rdquo; package, making it as flexible as possible for the end user.
- It is agnostic towards the choice of simulation code used to generate the input data.
- It possess both command-line functionality for those with little Python experience and an object-oriented API that allows the user to easily construct objects representing individual defects or complete systems of defects.
- It allows for constraints on the concentrations of any combination of defects and specific charge states of that defect, including the simulation of the influence of effective dopants.
- It is fully documented with good unit-test coverage.

The code has already been used in a number of studies, including some focussed on rationalising the 
properties of known materials [@LLZO-elect; @Squires2021-je] and others focussed on predicting the properties of novel materials [@Jackson2022-la].
The code can also assist with visualising defect energetics [@Haouari2021-xz].

One feature we believe is notable is the ability to arbitrarily fix the concentrations of specific defects or individual defect charge states.
Defect and electronic carrier concentrations are functions of temperature.
For a material with initial equilibrium defect concentrations determined by high-temperature synthesis conditions, rapid cooling may result in some or all of these defect populations being &ldquo;frozen in&rdquo; (kinetically trapped) at their initial high-temperature values, while other defect populations remain free to vary, to establish a new pseudo-equilibrium [@Maier2003-pc]. 
`py-sc-fermi` allows this scenario to be modelled through a two-step calculation, where in the second step some defect concentrations are fixed to their high temperature &ldquo;synthesis&rdquo; values.
The new pseudo-equilibrium defect concentrations can be quickly computed as a function of quenching temperature, allowing defect concentrations to be calculated under a range of synthesis and quench scenarios [@LLZO-elect].  

There is some minor convenience functionality for users of the Vienna Ab initio Software Package (VASP).

# Acknowledgements

The authors are grateful for feature requests and user testing from Seán Kavanagh, Joe Willis, Jiayi Cen, Sabrine Hachmioune, and Lavan Ganeshkumar.
This work was supported by the Faraday Institution grant number FIRG017. D.O.S. acknowledges support from the European Research Council, ERC (Grant 758345).
B.J.M. acknowledges support from the Royal Society (UF13032 and URF\\R\\191006).

# References
