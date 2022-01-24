import numpy as np
from scipy.optimize import minimize_scalar, minimize # type: ignore
from py_sc_fermi.constants import kboltz
from py_sc_fermi.dos import DOS
from py_sc_fermi.defect_charge_state import  FrozenDefectChargeState
import multiprocessing
import pandas as pd # type: ignore

class DefectSystem(object):

    def __init__(self, defect_species, dos, volume, temperature, spin_pol):
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
        self.spin_pol = spin_pol

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

    def get_sc_fermi(self, conv=1e-20, emin=None, emax=None, verbose=True, n_trial_steps=1000):
        """
        Solve to find value of E_fermi for which the DefectSystem is
        charge neutral

        args:
            conv (float): convergence tolerance for total charge when
            assessing the condition of charge neutrality
            emin (float): minmum energy for E_fermi search, cannot be lower than 
            the lowest energy in the DOS input, defaults to None, which will use
            the lowest energy in the DOS
            emax (float): maximum energy for Fermi energy search, cannot exceed the
            highest energy in the DOS input, default to None, which will use the
            highest energy in the DOS

        returns:
            e_fermi (float): self consistent Fermi energy in electron volts
        
        """
        if not emin:
            emin = self.dos.emin()
        if not emax:
            emax = self.dos.emax()
        direction = +1.0
        e_fermi = (emin + emax)/2.0
        step = 1.0
        converged = False
        reached_e_min = False
        reached_e_max = False
      #  TODO: need to check whether emin and emax are reached.
        for i in range(n_trial_steps):
            q_tot = self.q_tot(e_fermi=e_fermi)
            if e_fermi > emax:
                if reached_e_min or reached_e_max:
                    raise RuntimeError(f'No solution found between {emin} and {emax}')
                reached_e_max = True
                direction = -1.0
            if e_fermi < emin:
                if reached_e_max or reached_e_min:
                    raise RuntimeError(f'No solution found between {emin} and {emax}')
                reached_e_min = True
                direction = +1.0
            if abs(q_tot) < conv:
                converged = True
                break
            if q_tot > 0.0:
                if direction == +1.0:
                    step *= 0.25
                    direction = -1.0
            else:
                if direction == -1.0:
                    step *= 0.25
                    direction = +1.0
            e_fermi += step * direction
        e_fermi_err = (self.q_tot(e_fermi=e_fermi+step) - self.q_tot(e_fermi=e_fermi-step))/2.0
        if verbose:
            print(f'converged: {converged}')
            print(f'residual: {abs(q_tot)}')
            print(f'e_fermi_err: {e_fermi_err}')
            print(f'e_fermi: {e_fermi}')
        return e_fermi

    def report(self, emin=None, emax=None, conv=1e-20):
        """
        print a report in the style of SC-FERMI which summarises key properties of
        the defect system.

        Args:
            emin: minimum energy considered in Fermi energy solver, the default is the bottom of the DOS
            emax: maximum energy considered in Fermi energy solver, the default is the top of the DOS
            conv: convergence tolerance for charge neutrality condition 

        Returns:
            None

        """
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
        """
        for a given Fermi energy, calculate the net charge as the difference between 
        all positive species (including holes) and all negative species (including)
        electrons.

        Args:
            e_fermi (float): Fermi energy
        
        Returns:
            diff (float): net charge
        """
        p0, n0 = self.dos.carrier_concentrations(e_fermi, self.kT)
        lhs_def, rhs_def = self.defect_charge_contributions(e_fermi)
        lhs = p0 + lhs_def
        rhs = n0 + rhs_def
        diff = rhs - lhs
        return diff

    def abs_q_tot(self, e_fermi):
        return abs( self.q_tot(e_fermi) )

    def check_concentrations(self):
        """
        Assess whether defects with fixed concentrations sum to reasonable values,
        i.e. does the sum of all fixed charge states of a defect equal the total fixed 
        concentration of the defect.

        Args: None

        Returns: None
        """
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

    def get_transition_levels(self):
        """
        method which returns the transition levels as a dictionary:

        returns:
         transition_levels (dict): {defect_species_label: [[x_values],[y_values]]}
        """
        transition_levels = {}
        for defect_species in self.defect_species_names:
            transition_level = self.defect_species_by_name(defect_species).tl_profile(ef_min=self.dos.emin(), ef_max=self.dos.emax())
            x = [[x_value][0][0] for x_value in transition_level]
            y = [[y_value][0][1] for y_value in transition_level]
            transition_levels.update({defect_species:[x,y]})
        return transition_levels

    def to_dict_per_volume(self, emin=None, emax=None, conv=1e-20, decomposed=False):
        """
        returns a dictionary of relevent properties of the DefectSystem
        concentrations are reported in cm^-3
        
        args:
            emin (float): minimum energy for the Fermi energy search
            emax (float): maximum energy for the Fermi energy search
            conv (float): convergence tolerance for the charge neutrality condition
            decomposed (bool): False (default) returns concentrations of defect species,
            true returns the concentraion of individual charge states.
            
        returns:
            {**run_stats, **concs} (dict): {'Fermi Energy': self consistent fermi energy value,
                                            'p0': concentration of holes in cm^-3
                                            'n0': concentration of electrons in cm^-3
                                             concs: {defect concentrations in cm^-3}}
        """
        if not emin:
            emin = self.dos.emin()
        if not emax:
            emax = self.dos.emax()
        e_fermi = self.get_sc_fermi(verbose=False, emin=emin, emax=emax, conv=conv)
        p0, n0 = self.dos.carrier_concentrations(e_fermi, self.kT)
        concs = {}
        for ds in self.defect_species:
            conc = ds.get_concentration(e_fermi, self.temperature)
            if decomposed == True:
                    chg_states = ds.charge_state_concentrations(e_fermi, self.temperature)
                    all_chg_states = {str(k) : v * 1e24/ self.volume for k,v in chg_states.items()}
                    concs.update({ds.name:all_chg_states})
            else:
                    concs.update({ds.name:conc* 1e24/ self.volume})
        run_stats = {'Fermi Energy': e_fermi, 'p0': p0 * 1e24/ self.volume, 'n0':n0 * 1e24/ self.volume}

        return {**run_stats, **concs}


    def to_dict(self, emin=None, emax=None, conv=1e-20, decomposed=False):
        """
        returns a dictionary of relevent properties of the DefectSystem
        concentrations are reported per unit cell
        
        args:
            emin (float): minimum energy for the Fermi energy search
            emax (float): maximum energy for the Fermi energy search
            conv (float): convergence tolerance for the charge neutrality condition
            decomposed (bool): False (default) returns concentrations of defect species,
            true returns the concentraion of individual charge states.
            
        returns:
            {**run_stats, **concs} (dict): {'Fermi Energy': self consistent fermi energy value,
                                            'p0': concentration of holes per unit cell
                                            'n0': concentration of electrons per unit cell
                                             concs: {defect concentrations in per unit cell}}
        """
        if not emin:
            emin = self.dos.emin()
        if not emax:
            emax = self.dos.emax()
        e_fermi = self.get_sc_fermi(verbose=False, emin=emin, emax=emax, conv=conv)
        p0, n0 = self.dos.carrier_concentrations(e_fermi, self.kT)
        concs = {}
        for ds in self.defect_species:
            conc = ds.get_concentration(e_fermi, self.temperature)
            if decomposed == True:
                    chg_states = ds.charge_state_concentrations(e_fermi, self.temperature)
                    all_chg_states = {str(k) : v for k,v in chg_states.items()}
                    concs.update({ds.name:all_chg_states})
            else:
                    concs.update({ds.name:conc})
        run_stats = {'Fermi Energy': e_fermi, 'p0': p0, 'n0':n0}

        return {**run_stats, **concs}

    def write_inputs( self, filename='input-fermi.dat' ):
        """
        writes an input file which is compatible with the FORTRAN code
        SC-FERMI on which py-sc-fermi was initially based. 

        args:
            filename (string): name of file to write. defaults to 'input-fermi.dat' as
            this is what SC-FERMI expects to read.
        
        returns:
            None
        """
        with open(filename, 'w') as f:
            f.write( str(self.spin_pol) + '\n' )
            f.write( str(self.dos._nelect) + '\n' )
            f.write( str(self.dos._egap) + '\n')
            f.write( str(self.temperature) + '\n')
            i = 0
            for d in self.defect_species:
                free_chg_states = []
                for c in d.charge_states:
                    if type(d.charge_states[c]) != FrozenDefectChargeState:
                        free_chg_states.append(c)
                if len(free_chg_states) > 0:
                    i=i+1
            #print(i)
            f.write(str(i) +'\n')
            frozen_defects = []
            frozen_charge_states = []
            free_defects_to_write = []
            for d in self.defect_species:
                free_chg_states = []
                for c in d.charge_states:
                    if type(d.charge_states[c]) != FrozenDefectChargeState:
                        free_chg_states.append(c)
                if len(free_chg_states) > 0:
                        f.write( '{} {} {}'.format( d.name, len(free_chg_states), d.nsites ) + '\n')
                if d._fixed_concentration is not None:
                        frozen_defects.append(d)
                for c in d.charge_states:
                    if type(d.charge_states[c]) != FrozenDefectChargeState:
                        f.write( '{} {} {}'.format( c, d.charge_states[c].energy, d.charge_states[c].degeneracy ) + '\n')
                    if d.charge_states[c]._fixed_concentration is not False:
                        frozen_charge_states.append((d.name, d.charge_states[c]))
            f.write( str(len(frozen_defects)) + '\n' )
            if frozen_defects is not []:
                for fd in frozen_defects:
                    f.write( '{} {}'.format( fd.name, fd.fixed_concentration * 1e24 / self.volume   ) + '\n')  #
            f.write( str(len(frozen_charge_states)) + '\n' )
            if frozen_charge_states is not []:
                for fc in frozen_charge_states:
                    #print(fc)
                    f.write( '{} {} {}'.format( fc[0], fc[1].charge, fc[1]._concentration * 1e24 / self.volume ) + '\n')  #
        f.close()
