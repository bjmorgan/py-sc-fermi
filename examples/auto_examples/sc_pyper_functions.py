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

# Define necessary physical constants

kb_e = physical_constants['Boltzmann constant in eV/K'][0]
j_to_ev = 1/physical_constants['electron volt-joule relationship'][0]
S_0 = 205 * j_to_ev / physical_constants['Avogadro constant'][0]
Cp = (7/2)*kb_e


E_vbm=3.647 # Defining the valence band maximum from DFT calculation of the stoichiometric material

#### potential_alignment_corrections #####

# these were determined using the pot_al() funtion in this script 

potals = {'La_Zr' : 0,
          'Li_i' : 0,
          'La_Zr_-' : 0.10350315789471409,
          'Li_i_+' : 0.03406739130434033,
          'O_i' : 0,
          'O_i_-' : -0.004940000000004829,
          'v_La_---' : -0.3375894117647178,
          'v_Li' : 0,
          'v_Li_-' : -0.03239890109890098,
          'v_O' : 0,
          'v_O_+' : -0.043869148936181546,
          'v_O_++' : 0.03406739130434033,
          'v_Zr_----' : -0.2980481927711054,
          'Zr_La' : 0,
          'Zr_La_+' : -0.08226595744680765,
          'Zr_Li' : 0,
          'Zr_Li_+' : 0.1214978947368408,
          'Zr_Li_++' : 0.13771505376341509,
          'Zr_Li_+++' : 0.144820430107508,
          'Zr_Li_tet' : 0,
          'Zr_Li_tet_+' : 0.18679101123593966,
          'Zr_Li_tet_++' : 0.13769032258063163,
          'Zr_Li_tet_+++' : 0.18368409090908244,
          'Zr_i' : 0,
          'Zr_i_+' : 0.21727624999999762,
          'Zr_i_++' : 0.23355232558137118,
          'Zr_i_+++' : 0.23345833333332422,
          'Zr_i_++++' : 0.2488537500000021,
          'Li_Zr' : 0,
          'Li_Zr_-' : -0.004041111111121154,
          'Li_Zr_--' : -0.017765555555556034,
          'Li_Zr_---' : -0.026111111111113416,
          'Li_La' : 0,
          'Li_La_-' : -0.06163763440859782,
          'Li_La_--' : -0.06490326086959186,
          'La_Li' : 0,
          'La_Li_+' : 0.4787106060606092,
          'La_Li_++' : 0.4804250000000039}

##### unpaired electrons ####

# Unpaired electrons, for use in determining spin degeneracy of defects, were determined from magenetic moments in OUTCAR files found in the data set

ue = {'v_Zr_----' : -0.0,
      'v_La_---' : 0.0,
      'v_Li_-' : -0.0,
      'v_Li' : 1.0,
      'v_O' : -0.0,
      'v_O_+' : 1.0,
      'v_O_++' : 0.0,
      'Zr_La' : 1.0,
      'Zr_La_+' : -0.0,
      'La_Zr' : 1.0,
      'La_Zr_-' : 0.0,
      'Zr_Li' : 1.0,
      'Zr_Li_+' : 2.0,
      'Zr_Li_++' : 1.0,
      'Zr_Li_+++' : 0.0,
      'Zr_Li_tet' : 3.0,
      'Zr_Li_tet_+' : 2.0,
      'Zr_Li_tet_++' : 1.0,
      'Zr_Li_tet_+++' : -0.0,
      'Li_i' : 1.0,
      'Li_i_+' : 0.0,
      'O_i' : -0.0,
      'O_i_-' : 0.0,
      'Zr_i' : -0.0,
      'Zr_i_+' : -0.0,
      'Zr_i_++' : -0.0,
      'Zr_i_+++' : 1.0,
      'Zr_i_++++' : 0.0,
      'Li_Zr' : -0.0,
      'Li_Zr_-' : 2.0,
      'Li_Zr_--' : 1.0,
      'Li_Zr_---' : 0.0,
      'Li_La_-' : 0.0,
      'Li_La_--' : 0.0,
      'Li_La' : 0.0,
      'La_Li_+' : 0.0,
      'La_Li_++' : 0.0,
      'La_Li' : 0.0}



def dependance(P,T):
    """
    
    This function gives dependance of mu_O(T,P)
    args: P = pressure (float)
          T = temperature (float)
          
    returning:
          mu_O(T,P) (float)
          
    assuming oxygen is an ideal gas
    
    """
    
    chem_pot = 0.5 * ( (Cp * (T - 298))
                      - T * ( (S_0 + (Cp * np.log(T/298)) + (kb_e * np.log((0.21/P)) ) ) ))  
    return chem_pot                                                                      



def encircle(x,y, ax=None, **kw):
    if not ax: ax=plt.gca()
    p = np.c_[x,y]
    hull = ConvexHull(p)
    poly = plt.Polygon(p[hull.vertices,:], **kw)
    ax.add_patch(poly)
    
def con_hull(x,y, ax=None, **kw):
    p = np.c_[x,y]
    hull = ConvexHull(p)
    return hull

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
    
    formation_energy = ( calc.energy - stoich.energy ) - ( -1 * (lattice_site.energy/sum(lattice_site.stoichiometry.values()) + limit)) + (charge * ( E_vbm + potals[calc.title] )) + im_cor 
    
    return ChargeState(  charge, formation_energy, 2**ue[calc.title])

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
    
    formation_energy = ( calc.energy - stoich.energy ) - ( 1 * (lattice_site.energy/sum(lattice_site.stoichiometry.values()) + limit)) + (charge * ( E_vbm + potals[calc.title])) + im_cor 
    
    return ChargeState(  charge, formation_energy, 2**ue[calc.title])

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
    formation_energy = ( calc.energy - stoich.energy ) - ( 1 * (non_native.energy/sum(non_native.stoichiometry.values()) + non_native_limit) + -1 * (native.energy/sum(native.stoichiometry.values()) + native_limit)) + (charge * ( E_vbm ) + potals[calc.title]) + im_cor 
    return ChargeState(  charge, formation_energy, 2**ue[calc.title])

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

def run_some_fermi(defects,T,nelect, e_gap, spin_polarised):
    """
    a function to run an SCFermi calculation, and retun the concentrations of all the relevent defects, charge carriers, and the Fermi Level. 
    args:
        defects (list): a list of defects to be considered as part of this calculation
        T (float): the temperature (in Kelvin) the concentrations are to be calculated at
        e_gap (float): magnitude of bandgap
        spin_polarised (bool): whether the DFT data is spin polarised
    returns:
        flat (dict): a dictionary specifiying defect concentration, fermi level and charge carrier concentrations
        
        Note: function currently assumes executable sc-fermi in directory.
    """
    out = []
    scf = SCFermi( defects, nelect, e_gap, T, spin_polarised)
    scf.output()
    with open("out.txt", 'w') as f:
        sp = subprocess.run(["./sc-fermi"],stdout=f)
        text_file = open(("out.txt") , "r")
        lines =  text_file.readlines()
        for i in defects:
                for line in lines:
                    if re.search(str(i.label)+' '+r'.*?Charge',line) is not None:
                        for h in range(i.n_charge_states):
                            joop = lines[lines.index(line)+(h+1)]
                            coop = joop.split()
                            x = coop[1]
                            y = float(coop[-2])
                            a_dict = {(i.label)+'_'+str(x) : y}
                            out.append(a_dict)
                            flat = {k: v for d in out for k, v in d.items()}
        for line in lines:
            if re.search('\(electrons\)', line) is not None:
                flat.update({'c_e' : float(line.split()[3])})
            if re.search('\(holes\)', line) is not None:
                flat.update({'c_h' : float(line.split()[3])})
            if re.search('\(eV\)', line) is not None:
                flat.update({'Fermi_level' : float(line.split()[4])})
        return flat
    
def run_some_frozen_fermi(defects,T,nelect, e_gap, spin_polarised, dopant_chg, dopant_conc):
    """
    a function to run a frozen SCFermi calculation, and retun the concentrations of all the relevent defects. 
    args:
        defects (list): a list of defects to be considered as part of this calculation
        T (float): the temperature (in Kelvin) the concentrations are to be calculated at
        e_gap (float): magnitude of bandgap
        spin_polarised (bool): whether the DFT data is spin polarised
        dopant_chg (int): charge on the frozen defect
        dopant_conc (int): concentration of the frozen defect in n/cm-3
    returns:
        flat (dict): a dictionary specifiying defect concentrations.
        
        Note: function currently assumes executable frozen sc-fermi in directory.
    """
    out = []
    scf = frozen_SCFermi( defects, nelect, e_gap, T, dopant_chg, dopant_conc, spin_polarised)
    scf.output()
    with open("out.txt", 'w') as f:
        sp = subprocess.run(["./frozen-sc-fermi"],stdout=f)
        text_file = open(("out.txt") , "r")
        lines =  text_file.readlines()
        for i in defects:
                #print(i.label)
                for line in lines:
                    if re.search(str(i.label)+' '+r'.*?Charge',line) is not None:
                        for h in range(i.n_charge_states):
                            joop = lines[lines.index(line)+(h+1)]
                            coop = joop.split()
                            x = coop[1]
                            y = float(coop[-2])
                            a_dict = {(i.label)+'_'+str(x) : y}
                            out.append(a_dict)
                            #print(a_dict)
                            flat = {k: v for d in out for k, v in d.items()}
                for line in lines:
                    if re.search('\(electrons\)', line) is not None:
                        flat.update({'c_e' : float(line.split()[3])})
                    if re.search('\(holes\)', line) is not None:
                        flat.update({'c_h' : float(line.split()[3])})
                    if re.search('\(eV\)', line) is not None:
                        flat.update({'Fermi_level' : float(line.split()[4])})
        return flat





def get_image_charge_correction(lattice, dielectric_matrix, conv=0.3,
                                factor=30, motif=[0.0, 0.0, 0.0]):
    """Calculates the anisotropic image charge correction by Sam Murphy in eV.
    References:
        [1] S. T. Murphy and N. D. H. Hine, Phys. Rev. B 87, 094111 (2013).
    Args:
        lattice (list): The defect cell lattice as a 3x3 matrix.
        dielectric_matrix (list): The dielectric tensor as 3x3 matrix.
        conv (float): A value between 0.1 and 0.9 which adjusts how much real
                      space vs reciprocal space contribution there is.
        factor: The cut-off radius, defined as a mutliple of the longest cell
            parameter.
        motif: The defect motif (doesn't matter for single point defects, but
            included in case we include the extended code for defect clusters).
        verbose (bool): If True details of the correction will be printed.
    Returns:
        The image charge correction as {charge: correction}
    """
    inv_diel = np.linalg.inv(dielectric_matrix)
    det_diel = np.linalg.det(dielectric_matrix)
    latt = np.sqrt(np.sum(lattice**2, axis=1))

    # calc real space cutoff
    longest = max(latt)
    r_c = factor * longest

    # Estimate the number of boxes required in each direction to ensure 
    # r_c is contained (the tens are added to ensure the number of cells
    # contains r_c). This defines the size of the supercell in which
    # the real space section is performed, however only atoms within rc
    # will be conunted.
    axis = np.array([int(r_c/a + 10) for a in latt])

    # Calculate supercell parallelpiped and dimensions
    sup_latt = np.dot(np.diag(axis), lattice)

    # Determine which of the lattice parameters is the largest and determine
    # reciprocal space supercell
    recip_axis = np.array([int(x) for x in factor * max(latt)/latt])
    recip_volume = abs(np.dot(np.cross(lattice[0], lattice[1]), lattice[2]))

    # Calculatate the reciprocal lattice vectors (need factor of 2 pi)
    recip_latt = np.linalg.inv(lattice).T * 2 * np.pi

    real_space = _get_real_space(conv, inv_diel, det_diel, latt, longest,
                                 r_c, axis, sup_latt)
    reciprocal = _get_recip(conv, inv_diel, det_diel, latt, recip_axis,
                            recip_volume, recip_latt, dielectric_matrix)

    # calculate the other terms and the final Madelung potential
    third_term = -2*conv/np.sqrt(np.pi*det_diel)
    fourth_term = -3.141592654/(recip_volume*conv**2)
    madelung = -(real_space + reciprocal + third_term + fourth_term)

    # convert to atomic units
    conversion = 14.39942
    real_ev = real_space * conversion / 2
    recip_ev = reciprocal * conversion / 2
    third_ev = third_term * conversion / 2
    fourth_ev = fourth_term * conversion / 2
    madelung_ev = madelung * conversion / 2

    correction = {}
    for q in range(1, 8):
        makov = 0.5 * madelung * q**2 * conversion
        lany = 0.65 * makov
        correction[q] = lany
    return correction


def _get_real_space(conv, inv_diel, det_diel, latt, longest, r_c, axis,
                    sup_latt):
    # Calculate real space component
    real_space = 0.0
    axis_ranges = [range(-a, a) for a in axis]

    # Pre-compute square of cutoff distance for cheaper comparison than
    # separation < r_c
    r_c_sq = r_c**2

    def _real_loop_function(mno):
        # Calculate the defect's fractional position in extended supercell
        d_super = np.array(mno, dtype=float) / axis
        d_super_cart = np.dot(d_super, sup_latt)

        # Test if the new atom coordinates fall within r_c, then solve
        separation_sq = np.sum(np.square(d_super_cart))
        # Take all cases within r_c except m,n,o != 0,0,0
        if separation_sq < r_c_sq and any(mno):
            mod = np.dot(d_super_cart, inv_diel)
            dot_prod = np.dot(mod, d_super_cart)
            N = np.sqrt(dot_prod)
            contribution = 1/np.sqrt(det_diel) * erfc(conv * N)/N
            return contribution
        else:
            return 0.
    real_space = sum(_real_loop_function(mno) for mno in
                     itertools.product(*axis_ranges))
    return real_space


def _get_recip(conv, inv_diel, det_diel, latt, recip_axis, recip_volume,
               recip_latt, dielectric_matrix):
    # convert factional motif to reciprocal space and
    # calculate reciprocal space supercell parallelpiped
    recip_sup_latt = np.dot(np.diag(recip_axis), recip_latt)

    # Calculate reciprocal space component
    axis_ranges = [range(-a, a) for a in recip_axis]

    def _recip_loop_function(mno):
        # Calculate the defect's fractional position in extended supercell
        d_super = np.array(mno, dtype=float) / recip_axis
        d_super_cart = np.dot(d_super, recip_sup_latt)

        if any(mno):
            mod = np.dot(d_super_cart, dielectric_matrix)
            dot_prod = np.dot(mod, d_super_cart)
            contribution = (exp(-dot_prod / (4 * conv**2)) / dot_prod)
            return contribution
        else:
            return 0.
    reciprocal = sum(_recip_loop_function(mno) for mno in
                     itertools.product(*axis_ranges))
    scale_factor = 4 * np.pi / recip_volume
    return reciprocal * scale_factor

iccs = get_image_charge_correction(np.array([[11.1172443672,-0.0162675841,0.0000008494],[-4.0770844023,10.3426677250, -0.0000037259],[-3.5200794553,-5.1632033465,9.1947714849]]),np.array([[24.776440,7.924358,-0.002433],[7.924358,31.143939, 0.001167],[-0.002433,0.001167,19.424886]]))

def per_form_un(number):
    """
    works out, for the lattice parameters of t-LLZO, defects per cubic cm
    args:
        number (float): defect per formula unit
    returns:
        per_cubic_cm: defect per cubic cm
    """
    per_unit_cell = number * 8
    per_cubic_angstrom = per_unit_cell  / (13.003 * 13.003 * 12.498)
    per_cubic_cm = per_cubic_angstrom * 1e+24
    return per_cubic_cm

def off_stoich(conc):
    """
    works out, for the lattice parameters of t-LLZO, defects per formula unit from defects per cubic cm
    args:
        number (float): defect per cubic cm
    returns:
        per_cubic_cm: defect per formula unit
    """
    per_cubic_angstrom = conc / 1e+24
    per_unit_cell = per_cubic_angstrom * (13.003 * 13.003 * 12.498)
    return per_unit_cell

def pot_al(outcar, ref_outcar, tol=0.5):
    """
    given a reference outcar, and a outcar of interest, will calculate a potential alignment:
        args:
            outcar: outcar from defect calculation
            ref_outcar: outcar from stoichiometric calculation
            tol: tolerance factor for electrostatic potentials considered 'far' from the defect
    """
    avg_defect = []
    avg_ref = []
    for i,j in zip(Outcar.read_avg_core_poten(ref_outcar)[-1],Outcar.read_avg_core_poten(outcar)[-1]):
        diff = i - j
        if diff < tol and diff > -tol:
            avg_defect.append(j)
            avg_ref.append(i)
    plt.plot(avg_defect)
    plt.plot(avg_ref)
    pot_al =  np.mean(np.array(avg_defect)) - np.mean(np.array(avg_ref))
    return pot_al

def compile_defects(chemical_potentials,schottky=False):
        v_O = make_defect([defects['v_O'],defects['v_O_+'],defects['v_O_++']], elements, interest['LLZO'], delta_mu=chemical_potentials, corr=iccs, sites=3 )
        Li_i = make_defect([defects['Li_i'],defects['Li_i_+']], elements, interest['LLZO'], delta_mu=chemical_potentials , corr=iccs)
        v_Li = make_defect([defects['v_Li'],defects['v_Li_-']], elements, interest['LLZO'], delta_mu=chemical_potentials , corr=iccs, sites=3)
        v_La = make_defect([defects['v_La_---']], elements, interest['LLZO'], delta_mu=chemical_potentials, corr=iccs,sites=2 )
        v_Zr = make_defect([defects['v_Zr_----']], elements, interest['LLZO'], delta_mu=chemical_potentials , corr=iccs)
        Zr_La = make_defect([defects['Zr_La'],defects['Zr_La_+']], elements, interest['LLZO'], delta_mu=chemical_potentials, corr=iccs, sites=2 )
        La_Zr = make_defect([defects['La_Zr_-'],defects['La_Zr']], elements, interest['LLZO'], delta_mu=chemical_potentials, corr=iccs )
        Zr_Li_o = make_defect([defects['Zr_Li'],defects['Zr_Li_+'],defects['Zr_Li_++'],defects['Zr_Li_+++']], elements, interest['LLZO'], delta_mu=chemical_potentials, corr=iccs,sites=2 )
        O_i = make_defect([defects['O_i'],defects['O_i_-']], elements, interest['LLZO'], delta_mu=chemical_potentials, corr=iccs ) 
        Zr_Li_t = make_defect([defects['Zr_Li_tet'],defects['Zr_Li_tet_+'],defects['Zr_Li_tet_++'],defects['Zr_Li_tet_+++']], elements, interest['LLZO'], delta_mu=chemical_potentials, corr=iccs )
        Zr_i = make_defect([defects['Zr_i'],defects['Zr_i_+'],defects['Zr_i_++'],defects['Zr_i_+++'],defects['Zr_i_++++']], elements, interest['LLZO'], delta_mu=chemical_potentials, corr=iccs )
        Li_Zr = make_defect([defects['Li_Zr'],defects['Li_Zr_-'],defects['Li_Zr_--'],defects['Li_Zr_---']], elements, interest['LLZO'], delta_mu=chemical_potentials, corr=iccs )
        Li_La = make_defect([defects['Li_La'],defects['Li_La_-'],defects['Li_La_--']], elements, interest['LLZO'], delta_mu=chemical_potentials, corr=iccs )
        La_Li = make_defect([defects['La_Li'],defects['La_Li_+'],defects['La_Li_++']], elements, interest['LLZO'], delta_mu=chemical_potentials, corr=iccs )
        if schottky == True:
            return [v_Li,Li_i,v_O]
        else:
            return [v_O,v_Li,v_La,v_Zr,Li_i,Zr_i,O_i,La_Zr,Zr_La,Zr_Li_o,Zr_Li_t,Li_Zr,Li_La,La_Li]
        
def read_full_transition_levels():
    v_O_ts =  pd.read_csv('transition-levels.dat', delimiter='\s+', nrows=4)
    v_Li_ts = pd.read_csv('transition-levels.dat', delimiter='\s+', skiprows=5, nrows=3)
    v_La_ts = pd.read_csv('transition-levels.dat', delimiter='\s+', skiprows=10, nrows=2)
    v_Zr_ts = pd.read_csv('transition-levels.dat', delimiter='\s+', skiprows=14, nrows=2)
    Li_i_ts = pd.read_csv('transition-levels.dat', delimiter='\s+', skiprows=18, nrows=3)
    Zr_i_ts = pd.read_csv('transition-levels.dat', delimiter='\s+', skiprows=23, nrows=5)
    O_i_ts =  pd.read_csv('transition-levels.dat', delimiter='\s+', skiprows=30, nrows=3)
    La_Zr_ts =  pd.read_csv('transition-levels.dat', delimiter='\s+', skiprows=35, nrows=3)
    Zr_La_ts =  pd.read_csv('transition-levels.dat', delimiter='\s+', skiprows=40, nrows=3)
    Zr_Li_ts =  pd.read_csv('transition-levels.dat', delimiter='\s+', skiprows=45, nrows=5)
    Zr_Li_t_ts =  pd.read_csv('transition-levels.dat', delimiter='\s+', skiprows=52, nrows=4)
    Li_Zr_ts =  pd.read_csv('transition-levels.dat', delimiter='\s+', skiprows=59, nrows=5)
    Li_La_ts =  pd.read_csv('transition-levels.dat', delimiter='\s+', skiprows=66, nrows=4)
    La_Li_ts =  pd.read_csv('transition-levels.dat', delimiter='\s+', skiprows=71, nrows=4)
    transition_levels = v_O_ts,v_Li_ts,v_La_ts,v_Zr_ts,Li_i_ts,Zr_i_ts,O_i_ts,La_Zr_ts,Zr_La_ts,Zr_Li_ts,Zr_Li_t_ts,Li_Zr_ts,Li_La_ts,La_Li_ts
    return transition_levels

defects = import_calculations_from_file('defects.yaml') # DFT defect calculation information
elements = import_calculations_from_file('elements.yaml')     # DFT element calculation information
interest = import_calculations_from_file('interest.yaml')     # DFT LLZO calculation information
