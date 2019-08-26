import numpy as np

from .constants import kboltz

class DefectChargeState(object):
    """Class for individual defect charge states"""
    
    def __init__(self, charge, energy, degeneracy):
        """Instantiate a DefectChargeState.

        Args:
            charge (int): Charge of this charge state.
            energy (float): Energy of this charge state when E_Fermi = E(VBM) (in eV),
                including any chemical potential contributions and energy
                corrections.
            degeneracy (int): Degeneracy of this charge state (e.g. spin degeneracy).

        Returns:
            None

        """ 
        self._charge = charge
        self._energy = energy
        self._degeneracy = degeneracy
        self._fixed_concentration = False
       
    @property
    def charge(self):
        return self._charge

    @property
    def energy(self):
        return self._energy

    @property
    def degeneracy(self):
        return self._degeneracy

    @property
    def concentration_is_fixed(self):
        return self._fixed_concentration
 
    def get_formation_energy(self, e_fermi):
        """Calculate the formation energy of this charge state at a 
        specified Fermi energy.

        Args:
            e_fermi (float): Position of the Fermi energy, relative to the VBM (in eV).

        Returns:
            (float): The charge state energy (in eV).

        """
        return self.energy + self.charge * e_fermi
    
    def get_concentration(self, e_fermi, temperature):
        """Calculate the concentration of this charge state at a
        specified Fermi energy and temperature, per site in the unit
        cell..

        Args:
            e_fermi (float): Position of the Fermi energy, relative to the VBM (in eV).
            temperature (float): Temperature (in K).
        
        Returns:
            (float): Per-site concentration of this charge state.

        """
        expfac = -self.get_formation_energy(e_fermi) / (kboltz * temperature)
        return self.degeneracy * np.exp(expfac)
        
    def __repr__(self):
        return f'q={self.charge:+2}, e={self.energy}, deg={self.degeneracy}'

class FrozenDefectChargeState(object):
    """Class for individual defect charge states with fixed concentrations"""

    def __init__(self, charge, concentration):
        """Instantiate a FrozenDefectChargeState.

        Args:
            charge (int): Charge of this charge state.
            concentration (float): Fixed concentration of this charge state per site.

        Returns:
            None

        """
        self._charge = charge
        self._concentration = concentration
        self._fixed_concentration = True

    @property
    def charge(self):
        return self._charge

    @property
    def concentration(self):
        return self._concentration

    @property
    def concentration_is_fixed(self):
        return self._fixed_concentration
 
    def get_concentration(self, e_fermi, temperature):
        """Convenience method to return the concentration of this charge state,
        per site in the unit cell.

        Args:
            e_fermi (float): Position of the Fermi energy, relative to the VBM (in eV).
            temperature (float): Temperature (in K).
        
        Returns:
            (float): Per-site concentration of this charge state.

        """
        return self.concentration

    def __repr__(self):
        return f'q={self.charge:+2}, [c]={self.concentration}'

  
