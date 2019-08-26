import numpy as np

from .constants import kboltz

class DefectChargeState(object):
    
    def __init__(self, charge, energy, deg):
        self.charge = charge
        self.energy = energy
        self.deg = deg
        self.fixed_concentration = False
        
    def get_formation_energy(self, e_fermi):
        return self.energy + self.charge * e_fermi
    
    def get_concentration(self, e_fermi, temperature, nsite):
        expfac = -self.get_formation_energy(e_fermi) / (kboltz * temperature)
        return nsite * self.deg * np.exp(expfac)
        
    def __repr__(self):
        return f'q={self.charge:+2}, e={self.energy}, deg={self.deg}'

class FrozenDefectChargeState(object):

    def __init__(self, charge, concentration):
        self.charge = charge
        self.concentration = concentration
        self.fixed_concentration = True

    def get_concentration(self, e_fermi, temperature, nsite):
        return self.concentration

    def __repr__(self):
        return f'q={self.charge:+2}, [c]={self.concentration}'

  
