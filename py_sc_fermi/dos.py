import numpy as np
from typing import Tuple, Optional
from pymatgen.io.vasp import Vasprun  # type: ignore
from pymatgen.electronic_structure.core import Spin  # type: ignore
from scipy.constants import physical_constants  # type: ignore

kboltz = physical_constants["Boltzmann constant in eV/K"][0]


class DOS(object):
    """Class for handling density-of-states data and its integration.

    :param np.array dos: Density-of-states data.
    :param np.array edos: Energies for the density-of-states (in eV).
    :param float bandgap: Width of the band gap (in eV).
    :param int nelect: Number of electrons in the DOS calculation cell.
    :param bool spin_polarised: is the calculated DOS spin polarisised?
        (Default: ``False``)
    :raises ValueError: If `self.bandgap` > `max(self.edos)`.
    """

    def __init__(
        self,
        dos: np.ndarray,
        edos: np.ndarray,
        bandgap: float,
        nelect: int,
        spin_polarised=False,
    ):
        """Initialise a DOS instance."""
        self._dos = dos
        self._edos = edos
        self._bandgap = bandgap
        self._nelect = nelect
        self._spin_polarised = spin_polarised
        self.normalise_dos()

        if self.bandgap > self.emax():
            raise ValueError(
                """bandgap > max(self.edos). Please check your bandgap and
                 energy range (self.edos)."""
            )

    @property
    def dos(self) -> np.ndarray:
        """:return: Array representing the dos data."""
        return self._dos

    @property
    def edos(self) -> np.ndarray:
        """:return: Array representing the energy range data."""
        return self._edos

    @property
    def bandgap(self) -> float:
        """:return: band gap in eV."""
        return self._bandgap

    @property
    def spin_polarised(self) -> bool:
        """:return: true if dos data is spin polarised, else false"""
        return self._spin_polarised

    @property
    def nelect(self) -> int:
        """:return: number of electrons in the density-of-states calculation cell"""
        return self._nelect

    @classmethod
    def from_vasprun(
        cls, path_to_vasprun: str, nelect: int, bandgap: Optional[float] = None
    ):
        """
        generate ``py_sc_fermi.dos.DOS`` object from a ``VASP`` ``vasprun.xml``
        file. As this is parsed using pymatgen, the number of electrons is not
        contained in the vasprun data and must be passed in. On the other hand,
        If the bandgap is not passed in, it can be read from the vasprun file.

        :param str path_to_vasprun: path to vasprun.xml file
        :param int nelect: number of electrons in the calculation cell
        :param float bandgap: bandgap in eV. If not passed in, it will be read
            from the vasprun file.

        :return: DOS object
        :rtype: py_sc_fermi.dos.DOS
        """
        vr = Vasprun(path_to_vasprun, parse_potcar_file=False)
        densities = vr.complete_dos.densities
        cbm = vr.eigenvalue_band_properties[2]
        edos = vr.complete_dos.energies - cbm
        if len(densities) == 2:
            tdos_data = np.stack(
                [edos, np.abs(densities[Spin.up]), np.abs(densities[Spin.down])], axis=1
            )
            spin_pol = True
        else:
            tdos_data = np.stack([edos, np.abs(densities[Spin.up])], axis=1)
            spin_pol = False
        edos = tdos_data[:, 0]
        dos = np.sum(tdos_data[:, 1:], axis=1)
        if bandgap is None:
            bandgap = float(vr.eigenvalue_band_properties[0])
        return cls(
            dos=dos, edos=edos, nelect=nelect, bandgap=bandgap, spin_polarised=spin_pol 
        )

    @classmethod
    def from_dict(cls, dos_dict: dict):
        """
        return a DOS object from a dictionary containing the DOS data.
        If the density-of-states data is spin polarised, it should
        be stored as a list of two arrays, one for each spin.

        :param dict dos_dict: dictionary containing the DOS data.
        :return: :py:class:`DOS`
        :rtype: py_sc_fermi.dos.DOS
        """
        nelect = dos_dict["nelect"]
        bandgap = dos_dict["bandgap"]
        raw_dos = np.array(dos_dict["dos"])
        edos = np.array(dos_dict["edos"])
        shape = raw_dos.shape
        if len(shape) == 1:
            new_dos = raw_dos
            spin_pol = False
        elif shape[0] == 2:
            new_dos = np.sum(raw_dos, axis=0)
            spin_pol = True
        else:
             raise ValueError("dos_dict['dos'] is not in the correct format.")
        return cls(
            nelect=nelect,
            bandgap=bandgap,
            edos=edos,
            dos=new_dos,
            spin_polarised=spin_pol,
        )

    def sum_dos(self) -> np.ndarray:
        """
        :returns: integrated density-of-states up to the valence band maximum
        :rtype: np.array
        """
        vbm_index = np.where(self._edos <= 0)[0][-1]
        sum1 = np.trapz(self._dos[: vbm_index + 1], self._edos[: vbm_index + 1])
        return sum1

    def normalise_dos(self) -> None:
        """normalises the density of states w.r.t. number of electrons in the
        density-of-states calculation cell (self.nelect)"""
        integrated_dos = self.sum_dos()
        self._dos = self._dos / integrated_dos * self._nelect

    def emin(self) -> float:
        """:return: minimum energy in self.edos
        :rtype: float"""
        return self._edos[0]

    def emax(self) -> float:
        """:return: maximum energy in self.edos
        :rtype: float"""
        return self._edos[-1]

    def _p0_index(self) -> int:
        """:return: index of the valence band maximum in self._edos
        :rtype: int"""
        return np.where(self._edos <= 0)[0][-1]

    def _n0_index(self) -> int:
        """:return: index of the conduction band minimum in self._edos
        :rtype: int"""
        return np.where(self._edos > self.bandgap)[0][0]

    def carrier_concentrations(
        self, e_fermi: float, temperature: float
    ) -> Tuple[float, float]:
        """return electron and hole carrier concentrations from the Fermi-Dirac
        distribution multiplied by the density-of-states at a given Fermi energy
        and temperature.

        :param float e_fermi: Fermi energy in eV
        :param float temperature: temperature in K
        :return: electron and hole carrier concentrations
        :rtype: tuple[float, float]
        """
        p0 = np.trapz(
            self._p_func(e_fermi, temperature), self._edos[: self._p0_index() + 1]
        )
        n0 = np.trapz(
            self._n_func(e_fermi, temperature), self._edos[self._n0_index() :]
        )
        return p0, n0

    def _p_func(self, e_fermi: float, temperature: float) -> float:
        """Fermi Dirac distribution for holes."""
        return self.dos[: self._p0_index() + 1] / (
            1.0
            + np.exp(
                (e_fermi - self.edos[: self._p0_index() + 1]) / (kboltz * temperature)
            )
        )

    def _n_func(self, e_fermi: float, temperature: float) -> float:
        """Fermi Dirac distribution for electrons."""
        return self.dos[self._n0_index() :] / (
            1.0
            + np.exp((self.edos[self._n0_index() :] - e_fermi) / (kboltz * temperature))
        )
