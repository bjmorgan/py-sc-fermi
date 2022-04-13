import numpy as np
from .defect_system import DefectSystem
from .defect_species import DefectSpecies
from .defect_charge_state import DefectChargeState
from py_sc_fermi.dos import DOS
from pymatgen.core import Structure
from typing import Union
import yaml, os


def read_unitcell_data(filename: str) -> float:
    """
    get volume in A^3 from `unitcell.dat` file used in SC-Fermi

    Args:
        filename (str): `unitcell.dat` file to parse
    Returns:
        volume (float): cell volume in A^3
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
) -> dict[str, Union[list["py_sc_fermi.defect_species.DefectSpecies"], float, int]]:
    """
    return all information from a input file correctly formated to work
    for `SC-Fermi`.

    Args:
        filename (str): path to input file to parse
        volume (float): volume of structure in A^3
        frozen (bool): if the file to be read contains frozen defect concentrations
    Returns:
        input_data (dict): a dictionary of data required to initialise a defect system
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
    return dos information from a `totdos.dat` file.

    Args:
        filename (str): path to `totdos.dat` to parse
        bandgap (float): bandgap of host material
        nelect (int): number of the electrons in host material calculation
    Returns:
        dos (DOS): py_sc_fermi.DOS object
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
    return a set of inputs descrbing a full defect system from py_sc_fermi input files

    Args:
        structure_filename (str): path to `unitcell.dat` to parse
        totdos_filename (str): path to `totdos.dat` to parse
        input_fermi_filename (str): path to `SC-Fermi` to parse
    Returns:
        inputs (dict): set of input data for py_sc_fermi.DefectSystem
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
    return volume in A^3 for a given structure file. Relies on pymatgen structure parser
    so accepts a wide range of formats inc. POSCAR, .cif, vasp output files (OUTCAR, vasprun.xml) etc.

    Args:
        structure_file (str): path to structure file to parse

    Returns:
        volume (float): volume of structure in A^3
    """
    structure = Structure.from_file(structure_file)
    volume = structure.volume
    return float(volume)


def dos_from_dict(dos_dict: dict) -> "py_sc_fermi.dos.DOS":
    """
    return a DOS object from a dictionary
    """
    nelect = dos_dict["nelect"]
    bandgap = dos_dict["bandgap"]
    edos = dos_dict["edos"]
    if type(dos) == dict:
        dos = dos["up"] + dos["down"]
    else:
        dos = dos
    return DOS(nelect=nelect, bandgap=bandgap, edos=edos, dos=dos)


def defect_species_from_dict(
    defect_species_dict: dict, volume: float
) -> list["py_sc_fermi.defect_species.DefectSpecies"]:
    """
    return a DefectSpecies object from a dictionary
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
    return a DefectSystem object from a yaml file
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

    # if "convergence_tolerance" not in list(data.keys()):
    #     data["convergence_tolerance"] = 1e-19
    # if "max_iterations" not in list(data.keys()):
    #     data["max_iterations"] = 1500

    return DefectSystem(
        dos=dos,
        volume=volume,
        defect_species=defect_species,
        temperature=data["temperature"],
        convergence_tolerance=float(data["convergence_tolerance"]),
        n_trial_steps=int(data["n_trial_steps"]),
    )
