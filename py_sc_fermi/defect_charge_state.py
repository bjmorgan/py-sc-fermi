import numpy as np

from .constants import kboltz

class DefectChargeState(object):
    """Class for individual defect charge states"""
    
    def __init__(self, charge, energy, deg):
        """Instantiate a DefectChargeState.

        Args:
            charge (int): Charge of this charge state.
            energy (float): Energy of this charge state when E_Fermi = E(VBM) (in eV),
                including any chemical potential contributions and energy
                corrections.
            deg (int): Degeneracy of this charge state (e.g. spin degeneracy).

        Returns:
            None

        """ 
        self.charge = charge
        self.energy = energy
        self.deg = deg
        self.fixed_concentration = False
        
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
            (float): Mole fraction of this charge state.

        """
        expfac = -self.get_formation_energy(e_fermi) / (kboltz * temperature)
        return self.deg * np.exp(expfac)
        
    def __repr__(self):
        return f'q={self.charge:+2}, e={self.energy}, deg={self.deg}'

class FrozenDefectChargeState(object):

    def __init__(self, charge, concentration):
        self.charge = charge
        self.concentration = concentration
        self.fixed_concentration = True

    def get_concentration(self, e_fermi, temperature):
        return self.concentration

    def __repr__(self):
        return f'q={self.charge:+2}, [c]={self.concentration}'

  
