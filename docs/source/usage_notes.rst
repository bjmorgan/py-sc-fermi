Usage notes
-------------------------------------
- At different points in the documentation, a "unit cell" is referred to, this is the cell for which the density
  of states data was calculated, and any volumes, degeneracies and numbers of electrons provided should be 
  consitent with this structure for the reported defect concentrations to be accurate. In most cases, inconsistent
  specification of these values will lead to incorrect Fermi energies.
- The reported Fermi energy and transition levels reported are referenced to 0 eV; the code expects that the input
  density of states data is zeroed on the valence band maximum.
- The code operates using the following units,
  and all user-input must be consistent:
  
  - energy: electron volts. 
  
  - temperature: Kelvin.  
  
  - volume: Angstroms :superscript:`3`. 
  
- Concentrations are a special case, internally, the :py:mod:`py-sc-fermi` operates in the concentration of sites in the unit cell 
  which are defective, but will report in cm :superscript:`-3`. The documentation should always specify what kind of concentration
  is expected by a particular function, if it does not, this is a bug! Please report it. 
