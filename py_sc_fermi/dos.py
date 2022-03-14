import numpy as np
from pymatgen.io.vasp import Vasprun  # type: ignore
from pymatgen.electronic_structure.core import Spin  # type: ignore
from py_sc_fermi.constants import kboltz


class DOS(object):
    """Class for handling density-of-states data and its integration"""

    def __init__(
        self,
        dos: np.array,
        edos: np.array,
        egap: float,
        nelect: int,
        spin_polarised=False,
        normalise=True,
    ):
        """Initialise a DOS instance.

        Args:
            dos (np.array): Density of states data.
            edos (np.array): Energies for the density of states (in eV).
            egap (float): Width of the band gap (in eV).
            nelect (int): Number of electrons.
            spin_polarised (`bool`, optional): is the calculated DOS spin polarisised? Default
                is `False`
            normalise (:obj:`bool`, optional): Normalise the DOS so that the integral
                from e_min --> 0.0 is equal to the `nelect`. Default is `True`.

        Returns:
            None

        Raises:
            ValueError: If `egap` > `max(edos)`.

        """
        if egap > edos[-1]:
            raise ValueError(
                "ERROR: Your conduction band minimum is not accounted for in your DOS \n, increase the number of bands in your DOS calculation)"
            )
        self._dos = dos
        self._edos = edos
        self._egap = egap
        self._nelect = nelect
        self._spin_polarised = spin_polarised
        if normalise:
            self.normalise_dos()

    @property
    def dos(self) -> np.array:
        """Get the dos array."""
        return self._dos

    @property
    def edos(self) -> np.array:
        """Get the edos array."""
        return self._edos

    @property
    def egap(self) -> float:
        """Get the bandgap."""
        return self._egap

    @property
    def spin_polarised(self) -> bool:
        """is the DOS calculation spin polarised?"""
        return self._spin_polarised

    @property
    def nelect(self) -> int:
        """Get the number of electrons in the unit cell."""
        return self._nelect

    @classmethod
    def from_vasprun(cls, path_to_vasprun: str, nelect: int, bandgap: float = None):
        """
        generate DOS object (py_sc_fermi.dos) from a vasp 
        doscar given number of electrons in calculation
        
        Args:
            vasprun (string): path to vasprun to parse
            nelect (int): number of electrons in vasp calculation
            bandgap (float): user supplied band-gap (defualts to `None`),
            in which case, the bandgap will be determined from the DOS provided
        
        Returns:
            dos (DOS): DOS object 
        """
        vr = Vasprun(path_to_vasprun, parse_potcar_file=False)
        densities = vr.complete_dos.densities
        cbm = vr.eigenvalue_band_properties[2]
        edos = vr.complete_dos.energies - vr.parameters["SIGMA"] - cbm
        if len(densities) == 2:
            tdos_data = np.stack(
                [edos, densities[Spin.up], densities[Spin.down]], axis=1
            )
        else:
            tdos_data = np.stack([edos, densities[Spin.up]], axis=1)
        edos = tdos_data[:, 0]
        dos = np.sum(tdos_data[:, 1:], axis=1)
        if bandgap == None:
            bandgap = vr.eigenvalue_band_properties[0]
        return cls(dos, edos, nelect, bandgap)

    def sum_dos(self) -> np.array:
        """get integrated density of states"""
        vbm_index = np.where(self._edos <= 0)[0][-1]
        sum1 = np.trapz(self._dos[: vbm_index + 1], self._edos[: vbm_index + 1])
        return sum1

    def normalise_dos(self) -> None:
        """normalise the density of states w.r.t. number of electrons (self.nelect)"""
        integrated_dos = self.sum_dos()
        self._dos = self._dos / integrated_dos * self._nelect

    def emin(self) -> float:
        """get lowest energy in the dos data"""
        return self._edos[0]

    def emax(self) -> float:
        """get highest energy in the dos data"""
        return self._edos[-1]

    def _p0_index(self) -> int:
        """get index of the vbm in self._edos"""
        return np.where(self._edos <= 0)[0][-1]

    def _n0_index(self) -> int:
        """get index of the cbm in self._edos"""
        return np.where(self._edos > self.egap)[0][0]

    def carrier_concentrations(self, e_fermi: float, temperature: float):
        """get n0 and p0 using integrals (equations 28.9 in Ashcroft Mermin)
           for an abitrary value of temperature and Fermi energy. Typically used internally 
           to calculate carrier concentrations at the self-cosistent Fermi energy
        Args:
            efermi (float): Fermi energy
            temperature (float): temperature for Fermi-Dirac distribution solution 

        Returns:
            p0 (float): concentration of holes
            n0 (float): concentration of electrons
        """
        p0 = np.trapz(
            self._p_func(e_fermi, temperature), self._edos[: self._p0_index() + 1]
        )
        n0 = np.trapz(
            self._n_func(e_fermi, temperature), self._edos[self._n0_index() :]
        )
        return p0, n0

    def _p_func(self, e_fermi: float, temperature: float):
        return self.dos[: self._p0_index() + 1] / (
            1.0
            + np.exp(
                (e_fermi - self.edos[: self._p0_index() + 1]) / (kboltz * temperature)
            )
        )

    def _n_func(self, e_fermi: float, temperature: float):
        return self.dos[self._n0_index() :] / (
            1.0
            + np.exp((self.edos[self._n0_index() :] - e_fermi) / (kboltz * temperature))
        )
