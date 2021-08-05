import numpy as np

class DOS(object):
    """Class for handling density-of-states data and its integration"""

    def __init__(self, dos, edos, egap, nelect, normalise=True, replicate_sc_fermi=True):
        """Initialise a DOS instance.

        Args:
            dos (np.array): Density of states data.
            edos (np.array): Energies for the density of states (in eV).
            egap (float): Width of the band gap (in eV).
            nelect (int): Number of electrons.
            normalise (:obj:`bool`, optional): Normalise the DOS so that the integral
                from e_min --> 0.0 is equal to the `nelect`. Default is `True`.
            replicate_sc_fermi (:obj:`bool`, optional): If True will reproduce off-by-one
                errors in SC-Fermi when integrating the DOS. Default is `True`.

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
        self._replicate_sc_fermi = replicate_sc_fermi
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
        vbm_index = np.where(self._edos <= 0)[0][-1]
        if self._replicate_sc_fermi:
            sum1 = np.trapz(self._dos[:vbm_index+2], self._edos[:vbm_index+2]) # Off-by-one error
        else:
            sum1 = np.trapz(self._dos[:vbm_index+1], self._edos[:vbm_index+1]) # Correct integral
        return sum1
    
    def normalise_dos(self, verbose=False):
        integrated_dos = self.sum_dos()
        if verbose:
            print(f'Integration of DOS up to Fermi level: {integrated_dos}')
        self._dos = self._dos / integrated_dos * self._nelect
        if verbose:
            print(f'Renormalised integrated DOS        : {self.sum_dos()}')

    def emin(self):
        return self._edos[0]

    def emax(self):
        return self._edos[-1]

    def carrier_concentrations(self, e_fermi, kT):
        # get n0 and p0 using integrals (equations 28.9 in Ashcroft Mermin)
        p0_index = np.where(self._edos <= 0)[0][-1]
        n0_index = np.where(self._edos > self.egap)[0][0]
        if self._replicate_sc_fermi:
            p0 = np.trapz( p_func(e_fermi, self._dos[:p0_index+2], self._edos[:p0_index+2], kT ),
                           self._edos[:p0_index+2]) # Off-by-one error
        else:
            p0 = np.trapz( p_func(e_fermi, self._dos[:p0_index+2], self._edos[:p0_index+1], kT ),
                           self._edos[:p0_index+1]) # Off-by-one error
        n0 = np.trapz( n_func(e_fermi, self._dos[n0_index:], self._edos[n0_index:], kT ),
                       self._edos[n0_index:])
        return p0, n0

# Possibly better as class methods that integrate over a certain energy range.
def p_func(e_fermi, dos, edos, kT):
    return dos / (1.0 + np.exp((e_fermi - edos)/kT))

def n_func(e_fermi, dos, edos, kT):
    return dos / (1.0 + np.exp((edos - e_fermi)/kT))


