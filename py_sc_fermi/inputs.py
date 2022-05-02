from collections import namedtuple
from dataclasses import dataclass
import numpy as np
from py_sc_fermi.defect_species import DefectSpecies
from py_sc_fermi.defect_charge_state import DefectChargeState
from py_sc_fermi.dos import DOS
from pymatgen.core import Structure
from typing import List
import yaml

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
    def from_yaml(
        cls, input_file: str, structure_file: str = '', dos_file: str = ''
    ):
        """
        Generate a ``py_sc_fermi.inputs.InputSet`` object from a yaml file.

        :param str input_file: path to yaml file
        :param str structure_file: path to structure file
            (Default: ``None``)
        :param str dos_file: path to dos file
            (Default: ``None``)
        :return: InputSet object
        :rtype: py_sc_fermi.inputs.InputSet

        .. note::
            Only the ``.yaml`` file is required. If the structure_file and dos_file
            are not specified, the ``.yaml`` file must contain the volume and the density-of-states data.
        """
        with open(input_file, "r") as f:
            input_dict = yaml.safe_load(f)

        # if edos and dos are specified in the yaml file, generate
        # the dos object
        if "edos" in input_dict.keys() and "dos" in input_dict.keys():
            dos = DOS.from_dict(input_dict)
        # if the dos is specified in the SC-Fermi format, generate
        # from that
        elif dos_file.endswith(".dat"):
            dos_data = read_dos_data(input_dict["bandgap"], input_dict["nelect"], dos_file)
            dos = DOS(
                dos=dos_data.dos,
                edos=dos_data.dos,
                nelect=input_dict["nelect"],
                bandgap=input_dict["bandgap"],
            )
        # or if DOS file is an .xml, try and read it as a vaspun
        elif dos_file.endswith(".xml"):
            dos = DOS.from_vasprun(dos_file, input_dict["nelect"], input_dict["bandgap"])

        # if volume is specified in the yaml file, read it
        if "volume" in input_dict.keys():
            volume = input_dict["volume"]
        # otherwise, read it from the structure file
        else:
            volume = read_volume_from_structure_file(structure_file)

        return cls(
            dos=dos,
            volume=volume,
            defect_species=input_dict["defect_species"],
            temperature=input_dict["temperature"],
            convergence_tolerance=input_dict["convergence_tolerance"],
            n_trial_steps=input_dict["n_trial_steps"],
        )

    @classmethod
    def from_sc_fermi_inputs(
        cls, input_file: str, structure_file: str, dos_file: str, frozen: bool = False
    ):
        """
        Generate an InputSet object from a SC-Fermi input file.

        :param str input_file: path to SC-Fermi input file
        :param str structure_file: path to structure file
        :param str dos_file: path to totdos file
        :param bool frozen: whether the input file contains
            any fixed concentration defects (Default: ``False``)
        :return: InputSet object
        :rtype: InputSet
        """

        volume = read_volume_from_structure_file(structure_file)
        input_data = read_input_fermi(input_file, volume, frozen)
        dos = read_dos_data(input_data.bandgap, input_data.nelect, dos_file)
        return cls(
            dos=dos,
            volume=volume,
            defect_species=input_data.defect_species,
            temperature=input_data.temperature
        )


def is_yaml(filename: str) -> bool:
    """
    Check if a file is a yaml file.

    :param str filename: path to file
    :return: True if file is a yaml file, False otherwise
    :rtype: bool
    """
    try:
        with open(filename, "r") as f:
            yaml.safe_load(f)
    except yaml.YAMLError:
        return False
    return True

def volume_from_unitcell(filename: str) -> float:
    """
    Get volume in A^3 from `unitcell.dat` file-type used in SC-Fermi.

    :param str filename: path to `unitcell.dat` file
    :return: volume in A^3
    :rtype: float
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
    """
    Return all information from a input file correctly formated to work
    for ``SC-Fermi``.

    :param str filename: path to ``SC-Fermi` input file.
    :param float volume: volume of unit cell in A^3.
    :param bool frozen: whether the input file contains any fixed concentration
        ``DefectSpecies`` or ``DefectChargeState``s. I.e. whether this is an
        input file for ``SC-Fermi`` or ``Frozen-SC-Fermi``.
    :return: inputs (dict): set of input data for py_sc_fermi.DefectSystem
    :rtype: Dict[str, Union[List[:py:class:`DefectSpecies`], float, int]]
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
            except ValueError:
                raise ValueError(
                    f"Frozen defect {name} not found in defect species list"
                )

        # read fixed concentration charge states
        n_frozen_charge_states = int(pure_readin.pop(0))
        defect_names = [ds.name for ds in defect_species]
        print(defect_names)
        for n in range(n_frozen_charge_states):
            l = list(pure_readin.pop(0))
            defect_info = str(l).split(" ")
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
                print(l)
                defect_charge_state = DefectChargeState.from_string(
                    str(l), volume, frozen=True
                )
                defect_species.append(DefectSpecies(name, 1, {charge_state: defect_charge_state}))

    return InputFermiData(spin_pol, nelect, bandgap, temperature, defect_species)


def read_dos_data(
    bandgap: float,
    nelect: int,
    filename: str = "totdos.dat",
) -> DOS:
    """
    Read in total DOS from `totdos.dat` file used in SC-Fermi.

    :param float bandgap: bandgap of defect system in eV.
    :param int nelect: number of electrons in defect system.
    :param str filename: path to `totdos.dat` file.
    :return: :py:class:`DOS`
    :rtype: py_sc_fermi.dos.DOS
    """
    data = np.loadtxt(filename)
    edos = data[:, 0]
    dos = np.sum(np.abs(data[:, 1:]), axis=1)
    # if np.any(data[:, 1:] < 0.0):
    #     raise ValueError("Negative DOS values found")
    dos = DOS(dos=dos, edos=edos, nelect=nelect, bandgap=bandgap)
    return dos


def volume_from_structure(structure_file: str) -> float:
    """
    Calculate the volume of a structure from any file readable by
    :py:mod:pymatgen.

    :param str structure_file: path to file defining the structure.
    :return: volume of structure in Angstrom^3
    :rtype: float
    """
    return Structure.from_file(structure_file).volume


def read_volume_from_structure_file(structure_file: str) -> float:
    """
    Read the volume of a structure from a file.

    :param str file: path to file defining the structure.
    :return: volume of structure in Angstrom^3
    :rtype: float
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
