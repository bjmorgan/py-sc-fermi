import re
from vasppy.calculation import *
import subprocess
import pandas as pd
from scipy.stats import linregress
from scipy.constants import physical_constants
import numpy as np
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
import itertools
from math import exp, erfc
from pymatgen.io.vasp import Outcar
import yaml

with open('automator_config.yaml', 'r') as stream:
    try:
        config = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

# Define necessary physical constants

kb_e = physical_constants['Boltzmann constant in eV/K'][0]
j_to_ev = 1/physical_constants['electron volt-joule relationship'][0]
S_0 = 205 * j_to_ev / physical_constants['Avogadro constant'][0]
Cp = (7/2)*kb_e

def iden_defect(defect_label):
    """
    iden_defect is a conveneince function which utilises the decomposes title of a vasppy caclulation for utilisation in other funtions assuming naming conventions are followed
    args: defect label (string)
    returns: defect identifiers (list of strings)
    """
    return re.findall(r"[^]_[^_]+",defect_label)

def sc_fermi_vacancy_wrap( calc, lattice_site, stoich, limit, im_cor, charge ):

    """
    sc_fermi_vacancy_wrap determines the formation energy of a vacancy defect, and formats a 'ChargeState' object to be fed into sc_fermi
    args: calc = DFT calculation summary of the defective material (vasspy.calculation)
          lattice_site = The species in the lattice that is removed to form the vacancy (string)
          stoich = DFT calculation summary of the stoichiometric material that the vacancy is formed in (vasspy.caclulation)
          limit = chemical potential of the vacant species (float)
          im_cor = image charge correction for the relevant charge state (float)
          charge = charge associated with the defect (integer)
    """

    formation_energy = ( calc.energy - stoich.energy ) - ( -1 * (lattice_site.energy/sum(lattice_site.stoichiometry.values()) + limit)) + (charge * ( config['e_vbm']+ config['potals'][calc.title] )) + im_cor

    return ChargeState(  charge, formation_energy, int(2**config['ue'][calc.title]))

def sc_fermi_interstitial_wrap( calc, lattice_site, stoich, limit, im_cor, charge ):

    """
    sc_fermi_interstital_wrap determines the formation energy of a interstital defect, and formats a 'ChargeState' object to be fed into sc_fermi
    args: calc = DFT calculation summary of the defective material (vasspy.calculation)
          lattice_site = The species in the lattice that is added to form the interstital (string)
          stoich = DFT calculation summary of the stoichiometric material that the interstital is formed in (vasspy.caclulation)
          limit = chemical potential of the interstital species (float)
          im_cor = image charge correction for the relevant charge state (float)
          charge = charge associated with the defect (integer)
    """

    formation_energy = ( calc.energy - stoich.energy ) - ( 1 * (lattice_site.energy/sum(lattice_site.stoichiometry.values()) + limit)) + (charge * ( config['e_vbm'] + config['potals'][calc.title])) + im_cor

    return ChargeState(  charge, formation_energy, int(2**config['ue'][calc.title]))

def sc_fermi_sub_wrap( calc, non_native, native, stoich, non_native_limit, native_limit, im_cor, charge ):
    """
    sc_fermi_sub_wrap determines the formation energy of a substitutional defect, and formats a 'ChargeState' object to be fed into sc_fermi
    args: calc = DFT calculation summary of the defective material (vasspy.calculation)
          non_native = The species that is added to form the substitution (string)
          non_native = The species that is removed to form the substitution (string)
          stoich = DFT calculation summary of the stoichiometric material that the interstital is formed in (vasspy.caclulation)
          non_native_limit = chemical potential of the interstital species (float)
          native_limit = chemical potential of the vacant species (float)
          im_cor = image charge correction for the relevant charge state (float)
          charge = charge associated with the defect (integer)

          Returns ChargeState class with defect charge, E_F = 0 formation energy, and the spin degeneracy of the defect.
    """
    formation_energy = ( calc.energy - stoich.energy ) - ( 1 * (non_native.energy/sum(non_native.stoichiometry.values()) + non_native_limit) + -1 * (native.energy/sum(native.stoichiometry.values()) + native_limit)) + (charge * (config['e_vbm'] ) + config['potals'][calc.title]) + im_cor
    return ChargeState(  charge, formation_energy, int(2**config['ue'][calc.title]))

class ChargeState:
    """
    ChargeState class has three properties
        charge = charge associated with defect
        formation_energy = defect formation energy at E_F = 0
        degeneracy = spin degeneracy of defect
    """

    def __init__( self, charge, formation_energy, degeneracy ):
        self.charge = charge
        self.formation_energy = formation_energy
        self.degeneracy = degeneracy


class Defect:
    """
    Class describing a defect as understood by SCFermi, made up of charge states, each defect has a label, a charge_state, and a number of sites on which it can form
    """

    def __init__( self, label, charge_states, n_sites ):
        self.label = label
        self.charge_states = charge_states
        self.n_sites = n_sites

    @property
    def n_charge_states( self ):
        return len( self.charge_states )

class SCFermi:
    """
    SCFermi class containing all information needed to carry out an SC-Fermi calculationm including defects, number of eletrons, band gap magnitude, temperature, and whether the DFT calculation was spin polarised
    """

    def __init__( self, defects, nelect, e_gap, temperature, spin_polarised=False):
        self.defects = defects
        self.nelect = nelect
        self.e_gap = e_gap
        self.temperature = temperature
        self.spin_polarised = spin_polarised

    @property
    def n_defects( self ):
        return len( self.defects )

    def output( self ):

            with open('input-fermi.dat', 'w') as f:

                if self.spin_polarised:
                    f.write( '2' + '\n')
                else:
                    f.write( '1' + '\n' )
                f.write( str(self.nelect) + '\n' )
                f.write( str(self.e_gap) + '\n')
                f.write( str(self.temperature) + '\n')
                f.write( str(self.n_defects) + '\n' )
                for d in self.defects:
                    f.write( '{} {} {}'.format( d.label, d.n_charge_states, d.n_sites ) + '\n')
                    for c in d.charge_states:
                        f.write( '{} {} {}'.format( c.charge, c.formation_energy, c.degeneracy ) + '\n')

            f.close()

class frozen_SCFermi:
    """
    frozen_SCFermi class containing all information needed to carry out a frozen-SC-Fermi calculationm including defects, number of eletrons, band gap magnitude, temperature, the charge and concentration of the defect to be frozen and whether the DFT calculation was spin polarised.
    """

    def __init__( self, defects, nelect, e_gap, temperature, dopant_chg, dopant_conc, spin_polarised = False):
        self.defects = defects
        self.nelect = nelect
        self.e_gap = e_gap
        self.temperature = temperature
        self.spin_polarised = spin_polarised
        self.dopant_chg = float(dopant_chg)
        self.dopant_conc = float(dopant_conc)

    @property
    def n_defects( self ):
        return len( self.defects )

    def output( self ):

            with open('input-fermi-frozen.dat', 'w') as f:

                if self.spin_polarised:
                    f.write( '2' + '\n')
                else:
                    f.write( '1' + '\n' )
                f.write( str(self.nelect) + '\n' )
                f.write( str(self.e_gap) + '\n')
                f.write( str(self.temperature) + '\n')
                f.write( str(self.n_defects) + '\n' )
                for d in self.defects:
                    f.write( '{} {} {}'.format( d.label, d.n_charge_states, d.n_sites ) + '\n')
                    for c in d.charge_states:
                        f.write( '{} {} {}'.format( c.charge, c.formation_energy, c.degeneracy ) + '\n')
                f.write( '0' + '\n')
                f.write( '1' + '\n')
                f.write( 'dopant ' + str(self.dopant_chg) + ' ' + str(self.dopant_conc) + '\n')
            f.close()

def make_defect(defects,elements,stoich,delta_mu,corr=0,sites=1):
    """
    a function to create a defect from a vasppy.Calculation object given a little extra information:
        args:
        defects (list): a list of the the defective cell calculations that comprise this specific defect,
        elements (list): a list of the elemental reference calculations that comprise this system,
        stoich (vasppy.caluclation): the perfect reference calculation
        delta_mu (dict): a dictionary of elements, and the delta_mu values associated with them
        corr (float) = image charge correction
        sites (int) = number of sites in the perfect crystal this material can occupy

        returns: a defect object
    """
    chg_states = []
    for i in defects:
        defect = iden_defect(i.title)
        defect.append('0')
        if defect[0] is 'v':
            if re.search('\+', defect[-2]) is not None:
                formation_energy = sc_fermi_vacancy_wrap(i, elements['{}'.format(defect[1])], stoich, delta_mu['mu_'+'{}'.format(defect[1])], corr[len(defect[-2])], len(defect[-2]) )
            elif re.search('\-', defect[-2]) is not None:
                formation_energy = sc_fermi_vacancy_wrap(i, elements['{}'.format(defect[1])], stoich,  delta_mu['mu_'+'{}'.format(defect[1])], corr[len(defect[-2])], -1*len(defect[-2]) )
            else:
                formation_energy = sc_fermi_vacancy_wrap(i, elements['{}'.format(defect[1])], stoich,  delta_mu['mu_'+'{}'.format(defect[1])], 0, 0)
        elif defect[1] is 'i':
            if re.search('\+', defect[-2]) is not None:
                formation_energy = sc_fermi_interstitial_wrap(i, elements['{}'.format(defect[0])], stoich,  delta_mu['mu_'+'{}'.format(defect[0])], corr[len(defect[-2])], len(defect[-2]) )
            elif re.search('\-', defect[-2]) is not None:
                formation_energy = sc_fermi_interstitial_wrap(i, elements['{}'.format(defect[0])], stoich,  delta_mu['mu_'+'{}'.format(defect[0])], corr[len(defect[-2])], -1*len(defect[-2]) )
            else:
                formation_energy = sc_fermi_interstitial_wrap(i, elements['{}'.format(defect[0])], stoich,  delta_mu['mu_'+'{}'.format(defect[0])], 0, 0)
        else:
            if re.search('\+', defect[-2]) is not None:
                formation_energy = sc_fermi_sub_wrap(i, elements['{}'.format(defect[0])], elements['{}'.format(defect[1])], stoich, delta_mu['mu_'+'{}'.format(defect[0])], delta_mu['mu_'+'{}'.format(defect[1])], corr[len(defect[-2])], len(defect[-2]) )
            elif re.search('\-', defect[-2]) is not None:
                formation_energy = sc_fermi_sub_wrap(i, elements['{}'.format(defect[0])], elements['{}'.format(defect[1])], stoich, delta_mu['mu_'+'{}'.format(defect[0])], delta_mu['mu_'+'{}'.format(defect[1])], corr[len(defect[-2])], -1*len(defect[-2]) )
            else:
                formation_energy = sc_fermi_sub_wrap(i, elements['{}'.format(defect[0])], elements['{}'.format(defect[1])], stoich, delta_mu['mu_'+'{}'.format(defect[0])], delta_mu['mu_'+'{}'.format(defect[1])], 0, 0)
        chg_states.append(formation_energy)
        if len(defect) > 4:
            label = str (defect[0] + '_' + defect[1] + '_' + defect[2])
        else:
            label = str (defect[0] + '_' + defect[1])

    return Defect( label, chg_states, sites)