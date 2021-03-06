import numpy as np
from .defect_species import DefectSpecies
from .defect_charge_state import DefectChargeState
from py_sc_fermi.dos import DOS

def read_unitcell_data(filename, verbose=True):
    with open(filename, 'r') as f:
        readin =  [ l for l in f.readlines() if l[0] != '#' ]
    factor = float(readin[0])
    lattvec = np.array([ [ float(s) for s in l.split() ] for l in readin[1:4] ], order='F')
    lattvec *= factor
    volume = ( lattvec[0,0] * (lattvec[1,1]*lattvec[2,2] - lattvec[1,2]*lattvec[2,1])
             + lattvec[0,1] * (lattvec[1,2]*lattvec[2,0] - lattvec[1,0]*lattvec[2,2])
             + lattvec[0,2] * (lattvec[1,0]*lattvec[2,1] - lattvec[1,1]*lattvec[2,0]) )
    if verbose:
        print(f"Volume of cell: {volume} A^3")
    return volume

def read_input_data(filename, verbose=True, frozen=False):
    with open(filename, 'r') as f:
        readin = f.readlines()
        pure_readin = [ l for l in readin if l[0] != '#']
    nspinpol = int(pure_readin.pop(0))
    if verbose:
        if nspinpol == 1:
            print( 'Found non-spin polarised system...' )
        elif nspinpol == 2:
            print( 'Found spin polarised system...' )
        else:
            raise ValueError('spin polarisation specification must be 1 or 2')
    nelect = int(pure_readin.pop(0))
    if verbose:
        print( f'Number of electrons in system: {nelect}')
    egap = float(pure_readin.pop(0))
    if verbose:
        print( f'Energy gap of system: {egap} eV')
        if(egap < 0.0):
            print("You have a negative gap - this is going to get weird...")
    temperature = float(pure_readin.pop(0))
    if verbose:
        print( f'Temperature: {temperature} K')
    ndefects = int(pure_readin.pop(0))
    if verbose:
        print( f'Number of defect species: {ndefects}')
        if ndefects < 0:
            print( '0 defects found...')
    # For each species read in name, number of sites in the unit cell, 
    # charge states, degeneracy and energies
    defect_species = []
    for i in range(ndefects):
        l = pure_readin.pop(0).split()
        name = l[0] 
        ncharge = int(l[1])
        nsites = int(l[2])
        if verbose:
            if ncharge <= 0:
                print(f"ERROR: defect {len[ndefects]+1} has idiotic number of charge states!!")
        charges = []
        energies = []
        degs = []
        for j in range(ncharge):
            l = pure_readin.pop(0).split()
            charges.append(float(l[0]))
            energies.append(float(l[1]))
            degs.append(int(l[2]))
        charge_states  = [ DefectChargeState(c, e, d) for c, e, d in zip( charges, energies, degs ) ]
        defect_species.append( DefectSpecies(name, nsites, charge_states) )
    return { 'defect_species': defect_species,
             'egap': egap,
             'temperature': temperature,
             'nspinpol': nspinpol,
             'nelect': nelect }

def read_dos_data(filename, egap, nelect):
    data = np.loadtxt(filename)
    if data.shape[1] == 2:
        print("Reading non-spin-polarised DOS")
    elif data.shape[1] == 3:
        print("Reading spin-polarised DOS")
    else:
        raise ValueError("DOS input file does not contain Nx2 or Nx3 elements")
    edos = data[:,0]
    dos = np.sum(data[:,1:], axis=1)
    if np.any( data[:,1:] < 0.0 ):
        print( "WARNING: Found negative value(s) of DOS!" )
        print("         Do you know what you are doing?")
        print("         These may cause serious problems...")
    return DOS(dos=dos, edos=edos, nelect=nelect, egap=egap)

def inputs_from_files( unitcell_filename, totdos_filename, input_fermi_filename ):
    inputs = {}
    inputs['volume'] = read_unitcell_data(unitcell_filename)
    inputs.update( read_input_data(input_fermi_filename) )
    inputs['dos'] = read_dos_data(totdos_filename, egap=inputs['egap'], nelect=inputs['nelect'])
    return inputs
