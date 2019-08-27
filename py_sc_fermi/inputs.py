import numpy as np
from .defect_species import DefectSpecies
from .defect_charge_state import DefectChargeState

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

def read_dos_data(filename, nspinpol, egap):
    with open( filename, 'r') as f:
        readin = f.readlines()
        pure_readin = [ l.split() for l in readin if l[0] != '#']
    edos = []
    dos = []
    dosu = []
    dosd = []
    countdos = 0
    if nspinpol == 1:
        for l in pure_readin:
            edos.append(float(l[0]))
            dos.append(float(l[1]))
            if dos[-1] < 0.0:
                countdos += 1
    elif nspinpol == 2:
        for l in pure_readin:
            edos.append(float(l[0]))
            dosu.append(float(l[1]))
            dosd.append(float(l[2]))
            if (dosu[-1] < 0.0) or (dosd[-1] < 0.0):
                countdos += 1
            dos.append(dosu[-1] + dosd[-1])
    if countdos > 0:
        print( "WARNING: Found negative value(s) of DOS!" )
        print("         Do you know what you are doing?")
        print("         These may cause serious problems...")
    numdos = len(edos)
    emin = edos[0]
    emax = edos[-1]
    if egap > emax:
        # BJM: Not sure why this check is here / implemented like this
        raise ValueError('ERROR: your conduction band is not present in energy range of DOS!!')
    dos = np.array(dos)
    edos = np.array(edos)
    return dos, edos

def inputs_from_files( unitcell_filename, totdos_filename, input_fermi_filename ):
    inputs = {}
    inputs['volume'] = read_unitcell_data(unitcell_filename)
    inputs.update( read_input_data(input_fermi_filename) )
    inputs['dos'], inputs['edos'] = read_dos_data(totdos_filename, inputs['nspinpol'], inputs['egap'])
    return inputs
