import numpy as np

class DefectSpecies(object):
    
    def __init__(self, name, nsites, charge_states):
        self.name = name
        self.nsites = nsites
        self.charge_states = {cs.charge: cs for cs in charge_states}
        self.fixed_concentration = None
        
    def __repr__(self):
        to_return =  f'\n{self.name}, nsites={self.nsites}' 
        if self.fixed_concentration:
            to_return += f'\nfixed [c] = {self.fixed_concentration} cm^-3'
        to_return += '\n'+''.join([ f'  {cs.__repr__()}\n' 
            for cs in self.charge_states.values() ])
        return to_return
    
    def charge_states_by_energy( self, delphi ):
        return sorted( defect_species[0].charge_states, key=lambda x: x.energy(delphi) )
    
    def min_energy_charge_state( self, e_fermi ):
        return self.charge_states_by_energy( e_fermi )[0]   
    
    def get_formation_energy(self, charge, e_fermi):
        return self.charge_states[charge].get_formation_energy(e_fermi)
    
    def get_formation_energies(self, e_fermi):
        form_eners = {}
        for q, d in self.charge_states.items():
            form_eners[q] = self.get_formation_energy(
                q, e_fermi)
        return form_eners
    
    def charge_state_at_fermi_energy(self, e_fermi):
        form_e = self.get_formation_energies(e_fermi)
        return min([(q, d) for q, d in form_e.items()], key=lambda x: x[1])
    
    def tl_profile(self, ef_min,  ef_max):
        q1, d = self.charge_state_at_fermi_energy(
            ef_min)
        form_e = self.get_formation_energies(ef_min)
        points = [(ef_min, form_e[q1])]
        # TODO: is there a better way to do this?
        while q1 != min(self.charges):
            qlist = [q for q in self.charges if q < q1]
            nextp, nextq = min(((self.get_transition_level_and_energy(
                q1, q2),
                q2)
                for q2 in qlist), key=lambda p: p[0][0])
            if nextp[0] < ef_max:
                points.append(nextp)
                q1 = nextq
            else:
                break
        form_e = self.get_formation_energies(ef_max)
        points.append((ef_max, form_e[q1]))
        return np.array(points)
    
    def get_transition_level_and_energy(self, q1, q2):
        """
        Get transition level between two charge states
        Calculates the Fermi energy (relative to the host VBM) and formation
        energy for the transition level between charge states q1 and q2, as a
        function of elemental chemical potentials.
        Args:
            q1 (int): Charge of charge state 1.
            q2 (int): Charge of charge state 2.
            chem_pots (dict or None): Dictionary of elemental
                chemical potentials.  Keys are strings identifying
                elements. Values are the corresponding elemental chemical
                potentials. If None, set all chemical potentials to zero.
            correction_keys (list): List of keyword strings for
                the set of correction terms to include. Defaults to 'default'.
        Returns:
            2-tuple: Fermi energy (relative to the host VBM) and
            formation energy of the transition level as tuple of floats
            ``(trans_level, energy)``.
        """

        form_eners = self.get_formation_energies(0)
        trans_level = (form_eners[q2] - form_eners[q1]) / (q1 - q2)
        energy = q1 * trans_level + form_eners[q1]
        return (trans_level, energy)

    
    @property
    def charges(self):
        """
        Return a list of charges for the charge states of this defect.
        Args:
            None
        Returns:
            (list(int)): A list of integer charges.
        """
        return self.charge_states.keys()
    
    def get_concentration(self, e_fermi, temperature):
        if self.fixed_concentration:
            return self.fixed_concentration
        else:
            return sum( self.charge_state_concentrations( e_fermi, temperature ).values() )

    def fixed_conc_charge_states(self):
        return { q: cs for q, cs in self.charge_states.items() if cs.fixed_concentration }

    def variable_conc_charge_states(self):
        return { q: cs for q, cs in self.charge_states.items() if not cs.fixed_concentration }

    def charge_state_concentrations(self, e_fermi, temperature):
        cs_concentrations = { cs.charge: cs.get_concentration(e_fermi, temperature, self.nsites) 
                     for cs in self.charge_states.values() }
        if self.fixed_concentration:
            fixed_conc = sum( [ c for q, c in cs_concentrations.items() 
                                if q in self.fixed_conc_charge_states() ] )
            variable_conc = sum( [ c for q, c in cs_concentrations.items() 
                                   if q in self.variable_conc_charge_states() ] )
            constrained_conc = self.fixed_concentration - fixed_conc
            scaling = constrained_conc / variable_conc
            for q in cs_concentrations:
                if q in self.variable_conc_charge_states():
                    cs_concentrations[q] *= scaling
        return cs_concentrations    
