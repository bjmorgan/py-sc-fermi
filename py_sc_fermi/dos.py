import numpy as np

class DOS(object):

    def __init__(self, dos, edos, egap, nelect, normalise=True):
        if egap > edos[-1]:
            raise ValueError('ERROR: your conduction band is not present in energy range of DOS!!')
        self._dos = dos
        self._edos = edos
        self._egap = egap
        self._nelect = nelect
        if normalise:
            self.normalise_dos()

    @property
    def dos(self):
        return self._dos

    @property
    def edos(self):
        return self._edos

    @property
    def egap(self):
        return self._egap

    @property
    def nelect(self):
        return self._nelect

    def normalise_dos(self):
        vbm_index = np.where(self._edos <= 0)[0][-1]
        sum1 = np.trapz(self._dos[:vbm_index+2], self._edos[:vbm_index+2]) # BJM: possible off-by-one error?
        # np_edos[vbm_index+2] has a positive energy.
        self._dos = self._dos / sum1 * self._nelect

    def emin(self):
        return self._edos[0]

    def emax(self):
        return self._edos[-1]

    def carrier_concentrations(self, e_fermi, kT):
        # get n0 and p0 using integrals (equations 28.9 in Ashcroft Mermin)
        p0_index = np.where(self._edos <= 0)[0][-1]
        n0_index = np.where(self._edos > self.egap)[0][0]
        p0 = np.trapz( p_func(e_fermi, self._dos[:p0_index+2], self._edos[:p0_index+2], kT ),
                       self._edos[:p0_index+2]) # BJM: possible off-by-one error?
        n0 = np.trapz( n_func(e_fermi, self._dos[n0_index:], self._edos[n0_index:], kT ),
                       self._edos[n0_index:])
        return p0, n0

# Possibly better as class methods that integrate over a certain energy range.
def p_func(e_fermi, dos, edos, kT):
    return dos / (1.0 + np.exp((e_fermi - edos)/kT))

def n_func(e_fermi, dos, edos, kT):
    return dos / (1.0 + np.exp((edos - e_fermi)/kT))


