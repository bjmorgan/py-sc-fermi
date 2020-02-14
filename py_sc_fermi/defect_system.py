import numpy as np
from scipy.optimize import minimize_scalar, minimize
from .constants import kboltz
from py_sc_fermi.dos import DOS
from py_sc_fermi.defect_charge_state import  FrozenDefectChargeState
import multiprocessing
import pandas as pd

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

    def get_constrained_sc_fermi(self, constraint, total, conv=1e-16, emin=None, emax=None, verbose=True):
        """Calculate the self-consistent Fermi energy, subject to constrained net
        defect concentrations. The constraint to be satisfied is defined as a set of
        coefficients for the relevant defect species concentrations, and their
        constrained sum, e.g. [V_Li] - [Li_i] = 0.
        Args:
            constraint (dict): Dictionary of non-zero concentration coefficients,
                               i.e. {'v_Li': +1, 'Li_i': -1}.
            total (float): Summed concentration of constrained defect species.
            conv (float, optional): TODO
            emin (float, optional): TODO
            emax (float, optional): TODO
            verbose (bool, optional): TODO
        Returns:
            (float): self-consistent Fermi energy.
        """
        if not emin:
            emin = self.dos.emin()
        if not emax:
            emax = self.dos.emax()
        def constraint_func(phi):
            summed_total = total
            for name, c in constraint.items():
                summed_total += self.defect_species_by_name(name).get_concentration(phi, self.temperature)
            summed_total *= 1e10
            return summed_total
        phi_min = minimize_scalar(self.abs_q_tot, method='bounded', bounds=(emin, emax),
                        tol=conv, options={'disp': False, 'xatol': 1e-4} )
        phi_min = minimize(self.abs_q_tot, x0=phi_min.x,
                           constraints={'fun': constraint_func, 'type': 'eq'})
        return phi_min.x[0]

    def get_sc_fermi(self, conv=1e-16, emin=None, emax=None, verbose=True):
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
        # TODO: need to check whether emin and emax are reached.
        for i in range(1000):
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

    def get_sc_fermi_new(self, conv=1e-16, emin=None, emax=None, verbose=True, niter=100):
        if not emin:
            emin = self.dos.emin()
        if not emax:
            emax = self.dos.emax()
        bounds = [emin, emax]
        q_tot_min = self.q_tot(e_fermi=bounds[0])
        q_tot_max = self.q_tot(e_fermi=bounds[1])
        if verbose:
            print(f'Initial E_Fermi bracketing values are {emin} {emax}')
            print(f'E_F = {emin} => Q_tot = {q_tot_min}')
            print(f'E_F = {emax} => Q_tot = {q_tot_max}')
        if (q_tot_min > 0) or (q_tot_max < 0):
            raise ValueError(f'emin and emax do not appear to bound a zero-charge solution: [{emin} {emax}]')
        for i in range(niter):
            e_fermi = np.mean(bounds)
            q_tot = self.q_tot(e_fermi=e_fermi)
            if verbose:
                print(f'iteration {i}')
                print(bounds[0], e_fermi, bounds[1])
                print(q_tot_min, q_tot, q_tot_max)
                print()
            if (q_tot > 0) and (q_tot <= q_tot_max):
                q_tot_max = q_tot
                bounds[1] = e_fermi
            elif (q_tot < 0) and (q_tot >= q_tot_min):
                q_tot_min = q_tot
                bounds[0] = e_fermi
            else:
                raise RuntimeError('Warning! You probably shouldn\'t reach here!')
            if (abs(q_tot_min) < conv) and (abs(q_tot_max) < conv):
                break
        return { 'e_fermi': e_fermi,
                 'n_iter': i,
                 'bracket': bounds }

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

    def get_transition_levels(self):
        tls = {}
        for ds in self.defect_species_names:
            tl = self.defect_species_by_name(ds).tl_profile(ef_min=self.dos.emin(), ef_max=self.dos.emax())
            x = [[j][0][0] for j in tl]
            y = [[j][0][1] for j in tl]
            tls.update({ds:[x,y]})
        return tls

    def to_dict(self, emin=None, emax=None, conv=1e-16, decomposed=False):
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


    def to_dict_per_site(self, emin=None, emax=None, conv=1e-16, decomposed=False):
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

            with open(filename, 'w') as f:

                f.write( str(self.spin_pol) + '\n' )
                f.write( str(self.dos._nelect) + '\n' )
                f.write( str(self.dos._egap) + '\n')
                f.write( str(self.temperature) + '\n')
                #f.write( str(len(self.defect_species_names)) + '\n' )
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
