import numpy as np

from .constants import kboltz

class DefectChargeState(object):
    
    def __init__(self, charge, energy, deg):
        self.charge = charge
        self.energy = energy
        self.deg = deg
        
    def get_formation_energy(self, e_fermi):
        return self.energy + self.charge * e_fermi
    
    def get_concentration(self, e_fermi, temperature, nsite):
        expfac = -self.get_formation_energy(e_fermi) / (kboltz * temperature)
        return nsite * self.deg * np.exp(expfac)
        
    def __repr__(self):
        return f'q={self.charge}, e={self.energy}, deg={self.deg}'
