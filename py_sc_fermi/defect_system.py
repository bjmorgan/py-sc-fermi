import numpy as np
from scipy.optimize import minimize_scalar
from .constants import kboltz

class DefectSystem(object):
    
    def __init__(self, defect_species, volume, nelect, edos, dos, egap, temperature):
        self.defect_species = defect_species
        self.volume = volume
        self.nelect = nelect
        self.edos = edos
        self.dos = dos
        self.egap = egap
        self.temperature = temperature
        self.kT = kboltz * temperature
        self.normalise_dos()
        
    def normalise_dos(self):
        vbm_index = np.where(self.edos <= 0)[0][-1]
        sum1 = np.trapz(self.dos[:vbm_index+2], self.edos[:vbm_index+2]) # BJM: possible off-by-one error?
        # np_edos[vbm_index+2] has a positive energy.
        self.dos = self.dos / sum1 * self.nelect
            
    def get_sc_fermi(self, conv=1e-16, emin=None, emax=None, verbose=True):
        if not emin:
            emin = self.edos.min()
        if not emax:
            emax = self.edos.max()
        phi_min = minimize_scalar(self.abs_q_tot, method='bounded', bounds=(emin, emax),
                        tol=conv, options={'disp': False} )
        if verbose:
            print(phi_min)
        return phi_min.x
    
    def report(self):
        e_fermi = self.get_sc_fermi(verbose=False)
        print(f'SC Fermi level :      {e_fermi}  (eV)\n')
        p0, n0 = self.carrier_concentrations(e_fermi)
        print( 'Concentrations:')
        print( f'n (electrons)  : {n0* 1e24 / self.volume} cm^-3')
        print( f'p (holes)      : {p0 * 1e24 / self.volume} cm^-3')
        concall = []
        for ds in self.defect_species:
            concall.append( ds.get_concentration(e_fermi, self.temperature))
            print( f'{ds.name}      : {concall[-1] * 1e24 / self.volume} cm^-3')
        print()
        print('Breakdown of concentrations for each defect charge state:')
        for i, ds in enumerate(self.defect_species):
            print('---------------------------------------------------------')
            if concall[i] == 0:
                print( f'{ds.name}      : Zero total - cannot give breakdown')
                continue
            print(f'{ds.name}      : Charge Concentration(cm^-3) Total')
            for j, (q, cs) in enumerate(ds.charge_states.items()): 
                conc = cs.get_concentration(e_fermi,self.temperature, ds.nsite)
                print(f'           : {q} {conc * 1e24 / self.volume} {(conc*100 / concall[i]):.2f}')


    def carrier_concentrations(self, e_fermi):
        # get n0 and p0 using integrals (equations 28.9 in Ashcroft Mermin)
        p0_index = np.where(self.edos <= 0)[0][-1]
        n0_index = np.where(self.edos > self.egap)[0][0]
        p0 = np.trapz( p_func(e_fermi, self.dos[:p0_index+2], self.edos[:p0_index+2], self.kT ), 
                       self.edos[:p0_index+2]) # BJM: possible off-by-one error?
        n0 = np.trapz( n_func(e_fermi, self.dos[n0_index:], self.edos[n0_index:], self.kT ),
                       self.edos[n0_index:])
        return p0, n0

    def defect_charge_contributions(self, e_fermi):
        lhs = 0.0
        rhs = 0.0
        # get defect concentrations at E_F
        for ds in self.defect_species:
            for q, cs in ds.charge_states.items():
                if q != 0:
                    concd = cs.get_concentration( e_fermi, self.temperature, ds.nsite)
                    if q < 0:
                        rhs += concd * abs(cs.charge)
                    else:
                        lhs += concd * abs(cs.charge)
        return lhs, rhs
    
    def q_tot(self, e_fermi):
        p0, n0 = self.carrier_concentrations(e_fermi)
        lhs_def, rhs_def = self.defect_charge_contributions(e_fermi)
        lhs = p0 + lhs_def
        rhs = n0 + rhs_def
        diff = rhs - lhs
        return diff

    def abs_q_tot(self, e_fermi):
        return abs( self.q_tot(e_fermi) )

# Possibly better as class methods that integrate over a certain energy range.
def p_func(e_fermi, dos, edos, kT):
    return dos / (1.0 + np.exp((e_fermi - edos)/kT))

def n_func(e_fermi, dos, edos, kT):
    return dos / (1.0 + np.exp((edos - e_fermi)/kT))

