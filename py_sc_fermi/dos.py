import numpy as np
from typing import Tuple, Optional
from scipy.constants import physical_constants  # type: ignore[import-untyped]
from scipy.integrate import trapezoid  # type: ignore[import-untyped]

from pymatgen.io.vasp import Vasprun  # type: ignore[import-untyped]
from pymatgen.electronic_structure.core import Spin  # type: ignore[import-untyped]

from py_sc_fermi.warnings import suppresses_numpy_overflow

kboltz = physical_constants["Boltzmann constant in eV/K"][0]


class DOS:
    """
    Class for handling density-of-states data and its integration.

    Args:
        dos (np.ndarray): density-of-states data.
        edos (np.ndarray): energies associated with density-of-states data.
        bandgap (float): band gap
        nelect (int): number of electrons in density-of-states calculation
        spin_polarised (bool): is the calculated density-of-states spin polarised?
    """

    def __init__(
        self,
        dos: np.ndarray,
        edos: np.ndarray,
        bandgap: float,
        nelect: int,
        spin_polarised: bool = False,
    ):
        self._edos = edos
        self._bandgap = bandgap
        self._nelect = nelect
        self._spin_polarised = spin_polarised

        if self.spin_polarised:
            new_dos = np.sum(dos, axis=0)
            self._dos = new_dos
        else:
            self._dos = dos

        # The band-edge indices are found relative to E=0 (the VBM), so edos must bracket zero:
        if not (self._edos[0] <= 0 <= self._edos[-1]):
            raise ValueError(
                f"DOS edos must bracket zero (i.e. extend above and below zero), as the valence band "
                f"maximum is taken as E=0. Got an energy range of [{self._edos[0]}, {self._edos[-1]}]."
            )

        # band-edge indices are determined as those closest to zero (for VBM) and bandgap (for CBM), rather
        # than previous first-index-below-or-equal-to-zero and first-index-above-or-equal-to-bandgap, to
        # avoid issues with accumulated noise in `edos` (e.g. ``edos`` of [-0.9999, 0.0001, 1.0001] should
        # give VBM (``_p0_idx``) of 0.0001 not -0.9999:
        self._p0_idx = int(np.argmin(np.abs(self._edos)))
        self._n0_idx = int(np.argmin(np.abs(self._edos - self._bandgap)))

        # Cache integration indices — depend only on edos and bandgap,
        # before normalise_dos(), which reads _p0_int_idx via sum_dos():
        mid_gap = (self._p0_idx + self._n0_idx) // 2
        self._p0_integration_idx = max(mid_gap, self._p0_idx + 1)
        self._n0_integration_idx = min(mid_gap, self._n0_idx - 1) + 1
        
        self.normalise_dos()

        if self.bandgap > self.emax():
            raise ValueError(
                """bandgap > max(self.edos). Please check your bandgap and
                 energy range (self.edos)."""
            )

    @property
    def dos(self) -> np.ndarray:
        """density-of-states array

        Returns:
            np.ndarray: density-of-states data
        """
        return self._dos

    @property
    def edos(self) -> np.ndarray:
        """energy associated with density-of-states data

        Returns:
            np.ndarray: energy associated with the density-of-states data
        """
        return self._edos

    @property
    def bandgap(self) -> float:
        """bandgap of the density of states data

        Returns:
            float: bandgap
        """
        return self._bandgap

    @property
    def spin_polarised(self) -> bool:
        """bool describing whether the density-of-states data is spin-polarised or not

        Returns:
            bool: True if the ``DOS`` is spin-polarised, else False
        """
        return self._spin_polarised

    @property
    def nelect(self) -> int:
        """number of electrons in density of states calculation with which to
        normalise the ``DOS`` with respect to.

        Returns:
            int: number of electrons
        """
        return self._nelect

    @classmethod
    def from_vasprun(
        cls,
        path_to_vasprun: str,
        nelect: Optional[int] = None,
        bandgap: Optional[float] = None,
    ) -> "DOS":
        """Generate DOS object from a VASP vasprun.xml
        file. As this is parsed using pymatgen, the number of electrons is not
        contained in the vasprun data and must be passed in. On the other hand,
        If the bandgap is not passed in, it can be read from the vasprun file.

        Args:
            path_to_vasprun (str): path to vasprun file
            nelect (int): number of electrons in vasp calculation associated with
              the vasprun
            bandgap (Optional[float], optional): bandgap. Defaults to None.
        """
        vr = Vasprun(
            path_to_vasprun,
            parse_potcar_file=False,
            separate_spins=False # This is the default, but it does not hurt to be explicit.
        )
        band_properties = vr.eigenvalue_band_properties
        if not (isinstance(band_properties, tuple) and len(band_properties) == 4 
            and all(isinstance(band_properties[i], float) for i in (0,1,2))
            and isinstance(band_properties[3], bool)):
            raise TypeError(
                "eigenvalue_band_properties from pymatgen has unexpected format. "
                "Expected tuple[float, float, float, bool]"
            )
        densities = vr.complete_dos.densities
        vbm = band_properties[2]
        edos = vr.complete_dos.energies - vbm
        if len(densities) == 2:
            dos = np.array([densities[Spin.up], densities[Spin.down]])
            spin_pol = True
        else:
            dos = np.array(densities[Spin.up])
            spin_pol = False

        if nelect is None:
            nelect = int(vr.parameters["NELECT"])
        if bandgap is None:
            bandgap = band_properties[0]

        return cls(
            dos=dos, edos=edos, nelect=nelect, bandgap=bandgap, spin_polarised=spin_pol
        )

    @classmethod
    def from_dict(cls, dos_dict: dict) -> "DOS":
        """return a ``DOS`` object from a dictionary containing the density-of-states
        data. If the density-of-states data is spin polarised, it should
        be stored as a list of two arrays, one for each spin. The order is not
        important.

        Args:
            dos_dict (dict): dictionary defining the density of states data
        """
        nelect = dos_dict["nelect"]
        bandgap = dos_dict["bandgap"]
        dos = np.array(dos_dict["dos"])
        edos = np.array(dos_dict["edos"])
        spin_pol = dos_dict["spin_pol"]

        return cls(
            nelect=nelect,
            bandgap=bandgap,
            edos=edos,
            dos=dos,
            spin_polarised=spin_pol,
        )

    def as_dict(self) -> dict:
        """Return a dictionary representation of the DOS object

        Returns:
            dict: DOS as dictionary

        Note:
            The defect dictionary will always report the DOS data is not spin
            polarised, even if the input data was. This is an artefact related
            to maintaining the ability of `py-sc-fermi` to read files formatted
            for the FORTRAN SC-Fermi code. Future versions will consider how
            the code parses these files such that this is no longer an issue.
        """

        return dict(
            nelect=int(self.nelect),
            bandgap=float(self.bandgap),
            edos=list(self.edos),
            dos=list(self.dos),
            spin_pol=False,
        )

    def sum_dos(self) -> np.ndarray:
        """
        Returns:
            np.ndarray: integrated density-of-states up to mid-gap
        """
        return trapezoid(self._dos[: self._p0_integration_idx], self._edos[: self._p0_integration_idx])

    def normalise_dos(self) -> None:
        """normalises the density of states w.r.t. number of electrons in the
        density-of-states unit cell (``self.nelect``)
        """
        integrated_dos = self.sum_dos()
        self._dos = self._dos / integrated_dos * self._nelect

    def emin(self) -> float:
        """minimum energy in ``self.edos``

        Returns:
            float: minimum energy in ``self.edos``
        """
        return self._edos[0]

    def emax(self) -> float:
        """maximum energy in ``self.edos``

        Returns:
            float: maximum energy in ``self.edos``
        """
        return self._edos[-1]

    @suppresses_numpy_overflow
    def carrier_concentrations(
        self, e_fermi: float, temperature: float
    ) -> Tuple[float, float]:
        """return electron and hole carrier concentrations from the Fermi-Dirac
        distribution multiplied by the density-of-states at a given Fermi energy
        and temperature.

        Args:
            e_fermi (float): fermi energy
            temperature (float): temperature

        Returns:
            Tuple[float, float]: concentration of holes, concentration of electrons
        """
        p0 = trapezoid(self._p_func(e_fermi, temperature), self._edos[: self._p0_integration_idx])
        n0 = trapezoid(self._n_func(e_fermi, temperature), self._edos[self._n0_integration_idx :])
        return p0, n0

    def _p_func(self, e_fermi: float, temperature: float) -> np.ndarray:
        """Fermi Dirac distribution for holes."""
        return self.dos[: self._p0_integration_idx] / (
            1.0
            + np.exp((e_fermi - self.edos[: self._p0_integration_idx]) / (kboltz * temperature))
        )

    def _n_func(self, e_fermi: float, temperature: float) -> np.ndarray:
        """Fermi Dirac distribution for electrons."""
        return self.dos[self._n0_integration_idx :] / (
            1.0
            + np.exp((self.edos[self._n0_integration_idx :] - e_fermi) / (kboltz * temperature))
        )
