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

Point defects are determinants of the properties of functional materials. Point defects influnce the effiency of solar cell materials [@TLC], the capacities of
battery materials [@batteries], and performance thermoelectrics[@thermoelectric] to name a few. `py-sc-fermi` is a Python package which provides a numerical
approach for calculating a "self-consistent Fermi energy", which allows for the quantification of point defect poulations thereby allowing accurate calculation of figures of merit such as electronic conducitivity.

# Statement of need

`py-sc-fermi` is a Python package for determining the self-consistent Fermi energy of a material from knowledge of the point defect
energetics, allowing for quantification of point defect populations. While we aware of other scientific software that allows for these 
kinds of calculations, we are not aware of any equivalent packages that are not built in to other, larger "workflow" packages (reducing their overall flexibility)
or are built on a flexible python API which allows for rapid prototyping and exploration of different "what-if" scenarios, while striving to remain "code-agnositic", such that the program does not rely on the input data having being caluclated from a particular materials modelling code.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

The authors are grateful for feature requests and user testing from Seán Kavanagh, Joe Willis, and Jiayi Cen.

# References
