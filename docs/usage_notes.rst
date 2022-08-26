Usage notes
-------------------------------------
- At different points in the documentation, a "unit cell" is referred to, this is the cell for which the density
  of states data was calculated, and any volumes, degeneracies and numbers of electrons provided should be 
  consitent with this structure for the reported defect concentrations to be accurate. In most cases, inconsistent
  specification of these values will lead to incorrect Fermi energies.
- All Fermi energies reported by the code, and provided in inputs are referenced to the VBM, which in turn
  must be referenced to 0 in the input density of states data.
- The code operates using the following units,
  and all user-input must be consistent:
  
  - energy: electron volts. 
  
  - temperature: Kelvin.  
  
  - volume: Angstroms :superscript:`3`. 
  
- Concentrations are a special case, internally, the :py:mod:`py-sc-fermi` operates in the concentration of sites in the unit cell 
  which are defective, but will report in cm :superscript:`-3`. The documentation should always specify what kind of concentration
  is expected by a particular function, if it does not, this is a bug! Please report it. 
