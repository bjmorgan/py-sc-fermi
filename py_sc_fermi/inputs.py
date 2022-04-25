import numpy as np
from .defect_system import DefectSystem
from .defect_species import DefectSpecies
from .defect_charge_state import DefectChargeState
from py_sc_fermi.dos import DOS
from pymatgen.core import Structure
from typing import Union, Dict, List
import yaml, os


def read_unitcell_data(filename: str) -> float:
    """
    Get volume in A^3 from `unitcell.dat` file used in SC-Fermi.
    
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


def read_input_data(
    filename: str, volume: float = None, frozen: bool = False
) -> Dict[str, Union[List["py_sc_fermi.defect_species.DefectSpecies"], float, int]]:
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
        raise ValueError("Volume must be specified if frozen is True")

    with open(filename, "r") as f:
        readin = f.readlines()
        pure_readin = [line.strip() for line in readin if line[0] != "#"]

    input_data = {}

    # read in general defect system information
    input_data["spinpol"] = int(pure_readin.pop(0))
    input_data["nelect"] = int(pure_readin.pop(0))
    input_data["bandgap"] = float(pure_readin.pop(0))
    input_data["temperature"] = int(pure_readin.pop(0))

    # read in defect species information
    defect_species = []
    ndefects = int(pure_readin.pop(0))
    for i in range(ndefects):
        defect_description = pure_readin.pop(0).split()
        defect_name = defect_description[0]
        n_charge_states = int(defect_description[1])
        nsites = int(defect_description[2])
        charge_states = []
        for charge_state in range(n_charge_states):
            charge_state_description = pure_readin.pop(0).split()
            charge_state = int(charge_state_description[0])
            charge_state_formation_energy = float(charge_state_description[1])
            charge_state_degeneracy = int(charge_state_description[2])
            charge_states.append(
                DefectChargeState(
                    charge=charge_state,
                    energy=charge_state_formation_energy,
                    degeneracy=charge_state_degeneracy,
                )
            )
        defect_species.append(
            DefectSpecies(name=defect_name, charge_states=charge_states, nsites=nsites)
        )

    # read in frozen defect concentrations
    if frozen == True:

        # fix defect concentrations
        n_frozen_defects = int(pure_readin.pop(0))
        for n in range(n_frozen_defects):
            l = pure_readin.pop(0).split()
            name, concentration = l[0], float(l[1])
            for defect in defect_species:
                if name == defect.name:
                    defect.fix_concentration(concentration / 1e24 * volume)

        # read fixed concentration charge states
        n_frozen_charge_states = int(pure_readin.pop(0))
        for n in range(n_frozen_charge_states):
            l = pure_readin.pop(0).split()
            name, charge_state, concentration = l[0], int(l[1]), float(l[2])
            for defect in defect_species:
                if name == defect.name:
                    defect.charge_states[charge_state].fix_concentration(
                        concentration / 1e24 * volume
                    )
            else:
                defect_charge_state = DefectChargeState(
                    charge=charge_state,
                    fixed_concentration=concentration / 1e24 * volume,
                    degeneracy=1,
                )
                defect_species.append(
                    DefectSpecies(
                        name=name, charge_states=[defect_charge_state], nsites=1
                    )
                )

    input_data.update({"defect_species": defect_species})

    return input_data


def read_dos_data(
    bandgap: float, nelect: int, filename: str = "totdos.dat",
) -> "py_sc_fermi.dos.DOS":
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
    dos = np.sum(data[:, 1:], axis=1)
    if np.any(data[:, 1:] < 0.0):
        raise ValueError("Negative DOS values found")
    dos = DOS(dos=dos, edos=edos, nelect=nelect, bandgap=bandgap)
    return dos


def inputs_from_files(
    structure_filename: str = "unitcell.dat",
    totdos_filename: str = "totdos.dat",
    input_fermi_filename: str = "input-fermi.dat",
    frozen: bool = False,
) -> dict[str, Union[list["py_sc_fermi.defect_species.DefectSpecies"], float, int]]:
    """
    Read a set of inputs descrbing a full defect system from SC-Fermi input files
    and return them as a dictionary.

    :param str structure_filename: path to file defining the structure.
    :param str totdos_filename: path to ``totdos.dat`` file.
    :param str input_fermi_filename: path to ``SC-Fermi`` input file.
    :param bool frozen: whether the input file contains any fixed concentration
        ``DefectSpecies`` or ``DefectChargeState``s. I.e. whether this is an
        input file for ``SC-Fermi`` or ``Frozen-SC-Fermi``.
    :return: inputs (dict): set of input data for py_sc_fermi.DefectSystem
    :rtype: Dict[str, Union[List[:py:class:`DefectSpecies`], float, int]]

    """
    inputs = {}
    try:
        inputs["volume"] = read_unitcell_data(structure_filename)
    except:
        inputs["volume"] = volume_from_structure(structure_filename)
    inputs.update(
        read_input_data(input_fermi_filename, frozen=frozen, volume=inputs["volume"])
    )
    inputs["dos"] = read_dos_data(
        filename=totdos_filename, bandgap=inputs["bandgap"], nelect=inputs["nelect"]
    )
    return inputs


def volume_from_structure(structure_file: str) -> float:
    """
    Calculate the volume of a structure from any file readable by 
    :py:mod:pymatgen.

    :param str structure_file: path to file defining the structure.
    :return: volume of structure in Angstrom^3
    :rtype: float
    """
    return Structure.from_file(structure_file).volume


def dos_from_dict(dos_dict: dict) -> "py_sc_fermi.dos.DOS":
    """
    Return a DOS object from a dictionary containing the DOS data.

    :param dict dos_dict: dictionary containing the DOS data.
    :return: :py:class:`DOS`
    :rtype: py_sc_fermi.dos.DOS
    """
    nelect = dos_dict["nelect"]
    bandgap = dos_dict["bandgap"]
    dos = dos_dict["dos"]
    edos = dos_dict["edos"]
    if type(dos) == dict:
        dos = np.array(dos["up"]) + np.array(dos["down"])
        spin_pol = True
    else:
        dos = dos
        spin_pol = False
    return DOS(
        nelect=nelect,
        bandgap=bandgap,
        edos=np.array(edos),
        dos=np.array(dos),
        spin_polarised=spin_pol,
    )


def defect_species_from_dict(
    defect_species_dict: dict, volume: float
) -> list["py_sc_fermi.defect_species.DefectSpecies"]:
    """
    return a DefectSpecies object from a dictionary containing the defect
    species data.

    :param dict defect_species_dict: dictionary containing the defect species
        data.
    :param float volume: volume of the defect system in Angstrom^3.
    :return: :py:class:`DefectSpecies`
    :rtype: py_sc_fermi.defect_species.DefectSpecies
    """
    charge_states = []
    name = list(defect_species_dict.keys())[0]
    for n, c in defect_species_dict[name]["charge_states"].items():
        if "fixed_concentration" not in list(c.keys()):
            fixed_concentration = None
        else:
            fixed_concentration = float(c["fixed_concentration"]) / 1e24 * volume
        if "formation_energy" not in list(c.keys()):
            formation_energy = None
        else:
            formation_energy = float(c["formation_energy"])
        if formation_energy == None and fixed_concentration == None:
            raise ValueError(
                f"{name, n} must have one or both fixed concentration or formation energy"
            )
        charge_state = DefectChargeState(
            charge=n,
            energy=formation_energy,
            degeneracy=c["degeneracy"],
            fixed_concentration=fixed_concentration,
        )
        charge_states.append(charge_state)

    if "fixed_concentration" in defect_species_dict[name].keys():
        fixed_concentration = (
            float(defect_species_dict[name]["fixed_concentration"]) / 1e24 * volume
        )
        return DefectSpecies(
            name,
            defect_species_dict[name]["nsites"],
            charge_states=charge_states,
            fixed_concentration=fixed_concentration,
        )
    else:
        return DefectSpecies(
            name, defect_species_dict[name]["nsites"], charge_states=charge_states,
        )


def defect_system_from_yaml(filename: str) -> "py_sc_fermi.defect_system.DefectSystem":
    """
    Return a DefectSystem object from a yaml file

    :param str filename: path to yaml file containing the defect system data.
    :return: :py:class:`DefectSystem`
    :rtype: py_sc_fermi.defect_system.DefectSystem

    """
    with open(filename, "r") as f:
        data = yaml.load(f, Loader=yaml.SafeLoader)

    if "volume" not in data.keys():
        if "unitcell.dat" in os.listdir("."):
            volume = read_unitcell_data("unitcell.dat")
        elif "POSCAR" in os.listdir("."):
            volume = volume_from_structure("POSCAR")
        else:
            raise ValueError(
                "No volume found in input file and no file defining the structure detected in this directory"
            )
    else:
        volume = data["volume"]

    if "edos" in data.keys() and "dos" in data.keys():
        dos = dos_from_dict(data)
    elif "totdos.dat" in os.listdir("."):
        dos = read_dos_data(bandgap=data["bandgap"], nelect=data["nelect"])
    elif "vasprun.xml" in os.listdir("."):
        dos = DOS.from_vasprun(nelect=data["nelect"])
    else:
        raise ValueError("No DOS data found")

    defect_species = [
        defect_species_from_dict(d, volume) for d in data["defect_species"]
    ]

    if "convergence_tol" not in list(data.keys()):
        data["convergence_tol"] = 1e-18
    if "n_trial_steps" not in list(data.keys()):
        data["n_trial_steps"] = 1500

    return DefectSystem(
        dos=dos,
        volume=volume,
        defect_species=defect_species,
        temperature=data["temperature"],
        convergence_tolerance=float(data["convergence_tol"]),
        n_trial_steps=int(data["n_trial_steps"]),
    )
