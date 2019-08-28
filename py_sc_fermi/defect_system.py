import numpy as np
from scipy.optimize import minimize_scalar
from .constants import kboltz
from py_sc_fermi.dos import DOS
import multiprocessing

class DefectSystem(object):
    
    def __init__(self, defect_species, dos, volume, temperature):
        """Initialise a DefectSystem instance.

        Args:
            defect_species (list(DefectSpecies)): List of DefectSpecies objects.
            volume (float): Cell volume in A^3.
            dos (:obj:`DOS`): A DOS object.
            temperature (float): Temperature in K.

        Returns:
            None

        """ 
        self.defect_species = defect_species
        self.volume = volume
        self.dos = dos
        self.temperature = temperature
        self.kT = kboltz * temperature

    def __repr__(self):
        to_return = [f'DefectSystem\n',
                     f'  nelect: {self.dos.nelect} e\n',
                     f'  egap:   {self.dos.egap} eV\n',
                     f'  volume: {self.volume} A^3\n',
                     f'  temperature: {self.temperature} K\n',
                     f'\nContains defect species:\n']
        for ds in self.defect_species:
            to_return.append(str(ds))
        return ''.join(to_return)

    def defect_species_by_name(self, name):
        return [ ds for ds in self.defect_species if ds.name == name ][0]

    @property
    def defect_species_names(self):
        return [ ds.name for ds in self.defect_species ]
    
    def get_sc_fermi(self, conv=1e-16, emin=None, emax=None, verbose=True):
        if not emin:
            emin = self.dos.emin()
        if not emax:
            emax = self.dos.emax()
        phi_min = minimize_scalar(self.abs_q_tot, method='bounded', bounds=(emin, emax),
                        tol=conv, options={'disp': False} )
        if verbose:
            print(phi_min)
        return phi_min.x
    
    def report(self, emin=None, emax=None, conv=1e-16):
        if not emin:
            emin = self.dos.emin()
        if not emax:
            emax = self.dos.emax()
        e_fermi = self.get_sc_fermi(verbose=False, emin=emin, emax=emax, conv=conv)
        print(f'SC Fermi level :      {e_fermi}  (eV)\n')
        p0, n0 = self.dos.carrier_concentrations(e_fermi, self.kT)
        print( 'Concentrations:')
        print( f'n (electrons)  : {n0 * 1e24 / self.volume} cm^-3')
        print( f'p (holes)      : {p0 * 1e24 / self.volume} cm^-3')
        for ds in self.defect_species:
            concall = ds.get_concentration(e_fermi, self.temperature)
            print( f'{ds.name:9}      : {concall * 1e24 / self.volume} cm^-3')
        print()
        print('Breakdown of concentrations for each defect charge state:')
        for ds in self.defect_species:
            charge_state_concentrations = ds.charge_state_concentrations(e_fermi, self.temperature)
            concall = ds.get_concentration(e_fermi, self.temperature)
            print('---------------------------------------------------------')
            if concall == 0.0:
                print( f'{ds.name:11}: Zero total - cannot give breakdown')
                continue
            print(f'{ds.name:11}: Charge Concentration(cm^-3) Total')
            for q, conc in charge_state_concentrations.items():
                if ds.charge_states[q].concentration_is_fixed:
                    fix_str = ' [fixed]'
                else:
                    fix_str = ''
                print(f'           : {q: 1}  {conc * 1e24 / self.volume:5e}          {(conc * 100 / concall):.2f} {fix_str}')

    def defect_charge_contributions(self, e_fermi):
        contrib = np.array([ ds.defect_charge_contributions( e_fermi, self.temperature ) 
                             for ds in self.defect_species ])
        lhs = np.sum( contrib[:,0] )
        rhs = np.sum( contrib[:,1] )
        return lhs, rhs
    
    def q_tot(self, e_fermi):
        p0, n0 = self.dos.carrier_concentrations(e_fermi, self.kT)
        lhs_def, rhs_def = self.defect_charge_contributions(e_fermi)
        lhs = p0 + lhs_def
        rhs = n0 + rhs_def
        diff = rhs - lhs
        return diff

    def abs_q_tot(self, e_fermi):
        return abs( self.q_tot(e_fermi) )

    def check_concentrations(self):
        for ds in self.defect_species:
            if not ds.fixed_concentration:
                continue
            fixed_concentrations = [ cs.concentration for cs in ds.charge_states.values() 
                                         if cs.concentration_is_fixed ]
            if sum(fixed_concentrations) > ds.fixed_concentration:
                raise ValueError(f'ERROR: defect {ds.name} has a fixed'
                                 +'total concentration less than'
                                 +'the sum of its fixed concentration charge states.')
            if len(fixed_concentrations) == len(ds.charge_states):
                if sum(fixed_concentrations) != ds.fixed_concentration:
                    raise ValueError(f'ERROR: defect {ds.name} has fixed concentrations'
                                     +'for all charge states, but the sum of these concentrations'
                                     +'does not equal the fixed total concentration.')

