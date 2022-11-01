from collections import namedtuple
from dataclasses import dataclass
import numpy as np
from py_sc_fermi.defect_species import DefectSpecies
from py_sc_fermi.defect_charge_state import DefectChargeState
from py_sc_fermi.dos import DOS
from pymatgen.core import Structure
from typing import List
import yaml
import os

InputFermiData = namedtuple(
    "InputFermiData",
    "spin_pol nelect bandgap temperature defect_species",
)


@dataclass
class InputSet:
    dos: DOS
    volume: float
    defect_species: List[DefectSpecies]
    temperature: float
    convergence_tolerance: float = 1e-18
    n_trial_steps: int = 1500

    @classmethod
    def from_yaml(cls, input_file: str, structure_file: str = "", dos_file: str = ""):
        """
        Generate an InputSet object from a given yaml file

        Args:
            input_file (str): path to yaml file to read
            structure_file (str): path to structure file to read
            dos_file (str): path to dos file to read

        Returns:
            InputSet: full set of inputs for ``py-sc-fermi.DefectSystem``.

        Note:
            Only the ``.yaml`` file is required. If the structure_file and dos_file
            are not specified, the ``.yaml`` file must contain the volume and
            the density-of-states data.
        """
        with open(input_file, "r") as f:
            input_dict = yaml.safe_load(f)

        if dos_file != "":
            if dos_file.endswith(".dat"):
                dos_data = read_dos_data(
                    input_dict["bandgap"], input_dict["nelect"], dos_file
                )
                dos = DOS(
                    dos=dos_data.dos,
                    edos=dos_data.edos,
                    nelect=input_dict["nelect"],
                    bandgap=input_dict["bandgap"],
                )

            # or if DOS file is an .xml, try and read it as a vasprun
            elif dos_file.endswith(".xml"):
                dos = DOS.from_vasprun(
                    dos_file, input_dict["nelect"], input_dict["bandgap"]
                )

        elif "edos" in input_dict.keys() and "dos" in input_dict.keys():
            dos = DOS.from_dict(input_dict)

        # or if there is a `totdos.dat` in the current folder
        elif "totdos.dat" in os.listdir("."):
            dos = read_dos_data(
                filename="totdos.dat",
                bandgap=input_dict["bandgap"],
                nelect=input_dict["nelect"],
            )
        # or if there is a vasprun in the current folder
        elif "vasprun.xml" in os.listdir("."):
            dos = DOS.from_vasprun("vasprun.xml", nelect=input_dict["nelect"])

        # if all else fails, raise an Error

        else:
            raise ValueError(
                """No DOS file specified, or dos and edos entry found in .yaml file, 
                and the dos could not be read from any other file  in 
                the current directory."""
            )

        # read volume
        if structure_file != "":
            volume = read_volume_from_structure_file(structure_file)

        elif "volume" not in input_dict.keys():
            if "unitcell.dat" in os.listdir("."):
                volume = volume_from_unitcell("unitcell.dat")
            elif "POSCAR" in os.listdir("."):
                volume = volume_from_structure("POSCAR")
            else:
                raise ValueError(
                    """No volume found in input file and no file defining the 
                    structure detected in this directory. We recommend specifying
                    the volume of the cell in the input .yaml file."""
                )
        else:
            volume = input_dict["volume"]

        # if the solver parameters are not in the .yaml file, set them
        if "convergence_tol" not in list(input_dict.keys()):
            input_dict["convergence_tol"] = 1e-18
        if "n_trial_steps" not in list(input_dict.keys()):
            input_dict["n_trial_steps"] = 1500

        defect_species = [
            DefectSpecies.from_dict(d, volume) for d in input_dict["defect_species"]
        ]

        return cls(
            dos=dos,
            volume=volume,
            defect_species=defect_species,
            temperature=input_dict["temperature"],
            convergence_tolerance=input_dict["convergence_tolerance"],
            n_trial_steps=input_dict["n_trial_steps"],
        )

    @classmethod
    def from_sc_fermi_inputs(
        cls,
        input_file: str,
        structure_file: str,
        dos_file: str,
        n_trial_steps: int = 1000,
        convergence_tolerance: float = 1e-18,
        frozen: bool = False,
    ) -> "InputSet":
        """Generate an InputSet object from a
        `SC-Fermi <https://github.com/jbuckeridge/sc-fermi>`_ -formatted input file.

        Args:
            input_file (str): path to file to read
            structure_file (str): path to structure file to read
            dos_file (str): path to dos file to read
            n_trial_steps (int, optional): number of trial steps for py-sc-fermi
              solver. Defaults to 1000.
            convergence_tolerance (float, optional): convergence tolerance for
              py-sc-fermi solver. Defaults to 1e-18.
            frozen (bool, optional): True if any defects or defect charge states in
              in the input file have fixed concentrations. Defaults to False.

        Returns:
            InputSet: full set of inputs for ``py-sc-fermi.DefectSystem``.
        """

        volume = read_volume_from_structure_file(structure_file)
        input_data = read_input_fermi(input_file, volume, frozen)
        dos = read_dos_data(input_data.bandgap, input_data.nelect, dos_file)
        return cls(
            dos=dos,
            volume=volume,
            defect_species=input_data.defect_species,
            temperature=input_data.temperature,
            n_trial_steps=n_trial_steps,
            convergence_tolerance=convergence_tolerance,
        )


def is_yaml(filename: str) -> bool:
    """True if file is readable as a yaml file

    Args:
        filename (str): path to file to check

    Returns:
        bool: ``True`` if file is readable yaml, else ``False``.
    """
    try:
        with open(filename, "r") as f:
            yaml.safe_load(f)
    except yaml.YAMLError:
        return False
    return True


def volume_from_unitcell(filename: str) -> float:
    """Get volume in A^3 from ``unitcell.dat`` file-type used in
    `SC-Fermi <https://github.com/jbuckeridge/sc-fermi>`_

    Args:
        filename (str): path to ``unitcell.dat`` file

    Returns:
        float: volume in A^3
    """

    with open(filename, "r") as f:
        readin = [l for l in f.readlines() if l[0] != "#"]
    factor = float(readin[0])
    lattvec = np.array([[float(s) for s in l.split()] for l in readin[1:4]], order="F")
    lattvec *= factor
    volume = np.linalg.det(lattvec)
    return volume


def read_input_fermi(
    filename: str, volume: float = None, frozen: bool = False
) -> InputFermiData:
    """Return all information from a input file correctly formatted to work
    for `SC-Fermi <https://github.com/jbuckeridge/sc-fermi>`_.

    Args:
        filename (str): path to ``SC-Fermi`` -formatted input file.
        volume (float, optional): unit cell volume. Only required if there are
          fixed-concentration defects in input file. Defaults to None.
        frozen (bool, optional): whether there are fixed-concentration defects
          in the input file. Defaults to False.

    Raises:
        ValueError: if the ``volume`` is not specified, but ``frozen == True``
        ValueError: if fixed-concentration defect does not have any charge-states
          defined.

    Returns:
        InputFermiData: input for generating a ``DefectSystem``.
    """

    if volume == None and frozen == True:
        raise ValueError("Volume must be specified if input contains 'frozen' defects.")

    with open(filename, "r") as f:
        readin = f.readlines()
        pure_readin = [line.strip() for line in readin if line[0] != "#"]

    # read in general defect system information
    spin_pol = int(pure_readin.pop(0))
    nelect = int(pure_readin.pop(0))
    bandgap = float(pure_readin.pop(0))
    temperature = int(pure_readin.pop(0))

    # read in defect species information
    defect_species = []
    ndefects = int(pure_readin.pop(0))
    for i in range(ndefects):
        ds = DefectSpecies._from_list_of_strings(pure_readin)
        defect_species.append(ds)

    # read in frozen defect concentrations
    if frozen is True and volume is not None:
        # fix defect concentrations
        n_frozen_defects = int(pure_readin.pop(0))
        for n in range(n_frozen_defects):
            l = pure_readin.pop(0).split()
            name, concentration = l[0], float(l[1])
            try:
                defect = [ds for ds in defect_species if ds.name == name][0]
                defect.fix_concentration(concentration / 1e24 * volume)
            except:
                raise ValueError(
                    f"Frozen defect {name} not found in defect species list"
                )

        # read fixed concentration charge states
        n_frozen_charge_states = int(pure_readin.pop(0))
        defect_names = [ds.name for ds in defect_species]
        for n in range(n_frozen_charge_states):
            full_string = str(pure_readin.pop(0))
            defect_info = full_string.split()
            name, charge_state, concentration = (
                defect_info[0],
                int(defect_info[1]),
                float(defect_info[2]),
            )
            if name in defect_names:
                defect = [ds for ds in defect_species if ds.name == name][0]
                defect.charge_states[charge_state].fix_concentration(
                    concentration / 1e24 * volume
                )
            else:
                defect_charge_state = DefectChargeState.from_string(
                    full_string, volume, frozen=True
                )
                defect_species.append(
                    DefectSpecies(name, 1, {charge_state: defect_charge_state})
                )

    return InputFermiData(spin_pol, nelect, bandgap, temperature, defect_species)


def read_dos_data(
    bandgap: float,
    nelect: int,
    filename: str = "totdos.dat",
) -> DOS:
    """read density of states data from an `SC-Fermi <https://github.com/jbuckeridge/sc-fermi>`_
    formatted ``totdos.dat`` file.

    Args:
        bandgap (float): bandgap of density-of-states data.
        nelect (int): number of electrons in the density-of-states data
        filename (str, optional): path to ``todos.dat`` file. Defaults to "totdos.dat".

    Returns:
        DOS: py-sc-Fermi ``DOS`` object
    """
    data = np.loadtxt(filename)
    edos = data[:, 0]
    dos = np.sum(np.abs(data[:, 1:]), axis=1)
    dos = DOS(dos=dos, edos=edos, nelect=nelect, bandgap=bandgap)
    return dos


def volume_from_structure(structure_file: str) -> float:
    """get volume of any structure file readable by ``pymatgen``

    Args:
        structure_file (str): path to file defining structure

    Returns:
        float: volume of structure
    """
    return Structure.from_file(structure_file).volume


def read_volume_from_structure_file(structure_file: str) -> float:
    """read the volume of a structure defined in a given file.

    Args:
        structure_file (str): path to the structure file.

    Raises:
        ValueError: if the structure file is not readable by ``py-sc-fermi``

    Returns:
        float: volume of structure
    """
    # if the structure is specified in the SC-Fermi format, calculate
    # the volume from that
    if structure_file.endswith(".dat"):
        volume = volume_from_unitcell(structure_file)
    # if all else fails, use pymatgen to try and get the
    # volume from the structure file
    else:
        try:
            volume = volume_from_structure(structure_file)
        except:
            raise ValueError(
                "Volume could not be read from structure file. Please check the file format."
            )
    return volume
