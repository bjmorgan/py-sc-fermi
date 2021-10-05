import numpy as np
from pymatgen.io.vasp import Vasprun # type: ignore
from pymatgen.electronic_structure.core import Spin

class DOS(object):
    """Class for handling density-of-states data and its integration"""

    def __init__(self, dos, edos, egap, nelect, spin_polarised=False, normalise=True):
        """Initialise a DOS instance.

        Args:
            dos (np.array): Density of states data.
            edos (np.array): Energies for the density of states (in eV).
            egap (float): Width of the band gap (in eV).
            nelect (int): Number of electrons.
            spin_polarise(`bool`, optional): is the calculated DOS spin polarisised? Default
                is `False`
            normalise (:obj:`bool`, optional): Normalise the DOS so that the integral
                from e_min --> 0.0 is equal to the `nelect`. Default is `True`.

        Returns:
            None

        Raises:
            ValueError: If `egap` > `max(edos)`.

        """
        if egap > edos[-1]:
            raise ValueError('ERROR: Your gap extends beyond your maximum energy')
        self._dos = dos
        self._edos = edos
        self._egap = egap
        self._nelect = nelect 
        if normalise:
            self.normalise_dos()

    @property
    def dos(self):
        """Get the dos array."""
        return self._dos

    @property
    def edos(self):
        """Get the edos array."""
        return self._edos

    @property
    def egap(self):
        """Get the bandgap."""
        return self._egap

    @property
    def nelect(self):
        """Get the number of electrons in the unit cell."""
        return self._nelect

    def sum_dos(self):
        """get integrated density of states"""
        vbm_index = np.where(self._edos <= 0)[0][-1]
        sum1 = np.trapz(self._dos[:vbm_index+1], self._edos[:vbm_index+1])
        return sum1
    
    def normalise_dos(self):
        """normalise the density of states w.r.t. number of electrons (self.nelect)""" 
        integrated_dos = self.sum_dos()  
        self._dos = self._dos / integrated_dos * self._nelect 

    def emin(self):
        """get lowest energy in the dos data"""
        return self._edos[0]

    def emax(self):
        """get highest energy in the dos data"""
        return self._edos[-1]

    def carrier_concentrations(self, e_fermi, kT):
        """get n0 and p0 using integrals (equations 28.9 in Ashcroft Mermin)
           for an abitrary value of kT and Fermi energy. Typically used internally 
           to calculate carrier concentrations at the self-cosistent Fermi energy
        Args:
            efermi (float): Fermi energy
            kT (float): k * temperature

        Returns:
            p0 (float): concentration of holes
            n0 (float): concentration of electrons
        """
        p0_index = np.where(self._edos <= 0)[0][-1]
        n0_index = np.where(self._edos > self.egap)[0][0]
        p0 = np.trapz( p_func(e_fermi, self._dos[:p0_index+1], self._edos[:p0_index+1], kT ),
                       self._edos[:p0_index+1])
        n0 = np.trapz( n_func(e_fermi, self._dos[n0_index:], self._edos[n0_index:], kT ),
                       self._edos[n0_index:])
        return p0, n0

# Possibly better as class methods that integrate over a certain energy range.
def p_func(e_fermi, dos, edos, kT):
    return dos / (1.0 + np.exp((e_fermi - edos)/kT))

def n_func(e_fermi, dos, edos, kT):
    return dos / (1.0 + np.exp((edos - e_fermi)/kT))

def DOS_from_vasprun(vasprun, nelect, bandgap = None):
    """
    generate DOS object (py_sc_fermi.dos) from a vasp 
    doscar given number of electrons in calculation
    
    Args:
        doscar (string): path to DOSCAR to parse
        nelect (float): number of electrons in vasp calculation
        bandgap (float): user supplied band-gap (defualts to `None`),
        in which case, the bandgap will be determined from the DOS provided
    
    Returns:
        dos (DOS): DOS object 
    """
    vr = Vasprun(vasprun)
    densities = vr.complete_dos.densities
    edos = vr.complete_dos.energies - vr.parameters['SIGMA']
    if len(densities) == 2:
        tdos_data = np.stack([edos,densities[Spin.up], densities[Spin.down]], axis=1)
    else:
        tdos_data = np.stack([edos,densities[Spin.up]], axis = 1)
    edos = tdos_data[:,0]
    dos = np.sum(tdos_data[:,1:], axis=1)
    if bandgap == None:
        bandgap = vr.tdos.get_gap()
    dos = DOS(dos, edos, bandgap, nelect=nelect)
    return dos
    

