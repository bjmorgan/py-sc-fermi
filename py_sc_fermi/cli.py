from ast import Pass
from py_sc_fermi.inputs import (
    inputs_from_files,
    read_input_data,
    read_unitcell_data,
    volume_from_structure,
    read_dos_data,
)
from py_sc_fermi.defect_system import DefectSystem
from py_sc_fermi.dos import DOS
import os, yaml, argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i", "--input_file", type=str, help="input file defining the defect system"
)
parser.add_argument(
    "-s", "--structure_file", help="Structure file giving the volume of a defect system"
)
parser.add_argument("-n", "--nelect", help="number of electrons in DOS calculation")
parser.add_argument(
    "-c", "--convergence_tolerance", help="convergence tolerance for sc-Fermi search", default= 1e-20, type = float
)
parser.add_argument(
    "-f", "--frozen_defects", help="frozen defects present in the defect system", type= bool
)
parser.add_argument("-b", "--band_gap", help="band gap of bulk system")
args = parser.parse_args()


def read_structure_file(
    files: list[str], structure_file: str = args.structure_file
) -> float:
    """
    tries to read the structure file as a unitcell.dat, 
    then will try to read via pymatgen
    """
    volume = None
    if "unitcell.dat" in files:
        volume = read_unitcell_data("unitcell.dat")
    if "POSCAR" in files:
        volume = volume_from_structure("POSCAR")
    else:
        try:
            volume = volume_from_structure(structure_file)
        except:
            Pass
        try:
            volume = read_unitcell_data(structure_file)
        except:
            Pass
    if volume is None:
        raise ValueError(
            "Could not read structure file from current directory. Please check that the file exists and is formatted correctly (see documentation)"
        )
    return volume


def read_DOS_file(
    files: list[str], nelect: int = args.nelect, bandgap: float = args.band_gap
) -> "py_sc_fermi.dos.DOS":
    """
    tries to the DOS as totdos.dat, then vasprun.xml
    """
    if "totdos.dat" in files:
        dos = read_dos_data("totdos.dat")
    elif "vasprun.xml" in files:
        dos = DOS().from_vasprun("vasprun.xml", nelect=nelect, bandgap=bandgap)
    else:
        print(
            "Could not read DOS information from current directory. Please define DOS in an alternate format (see documentation)"
        )
    return dos


def main():
    """
    read in input files for a defect system and return the defect system dumped
    to a yaml file.
    """
    input = args.input_file
    files = os.listdir(".")

    try:
        input_data = inputs_from_files(frozen=args.frozen_defects)
    except:
        if input.endswith(".yaml"):
            None
        else:
            volume = read_structure_file(files)
            dos = read_DOS_file(files)
            input_data = read_input_data(input, volume, args.frozen_defects)
            input_data["dos"] = dos

    defect_system = DefectSystem(
        defect_species=input_data["defect_species"],
        dos=input_data["dos"],
        volume=input_data["volume"],
        temperature=input_data["temperature"],
    )
    defect_system.report(conv=args.convergence_tolerance)


if __name__ == "__main__":
    main()
