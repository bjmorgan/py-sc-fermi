from ast import Pass
from py_sc_fermi.inputs import (
    read_input_data,
    read_unitcell_data,
    volume_from_structure,
    read_dos_data,
    inputs_from_files,
    defect_system_from_yaml,
)
from py_sc_fermi.defect_system import DefectSystem
from py_sc_fermi.dos import DOS
import os, yaml, argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i", "--input_file", type=str, help="input file defining the defect system"
)
parser.add_argument(
    "-s",
    "--structure_file",
    help="Structure file giving the volume of a defect system",
    default="unitcell.dat",
)
parser.add_argument(
    "-d",
    "--dos_file",
    help="File specifying the totdos of the system",
    default="totdos.dat",
)
parser.add_argument(
    "-c",
    "--convergence_tolerance",
    help="convergence tolerance for sc-Fermi search",
    default=1e-20,
    type=float,
)
parser.add_argument(
    "-f",
    "--frozen_defects",
    help="frozen defects present in the defect system",
    type=bool,
)
parser.add_argument("-b", "--band_gap", help="band gap of bulk system")
args = parser.parse_args()


def main():
    """
    read in input files for a defect system and return the defect system dumped
    to a yaml file.
    """
    input = args.input_file
    structure = args.structure_file
    totdos = args.dos_file
    frozen = args.frozen_defects

    if input.endswith(".yaml"):
        defect_system = defect_system_from_yaml(input)
    else:
        input_data = inputs_from_files(
            input_fermi_filename=input,
            structure_filename=structure,
            totdos_filename=totdos,
            frozen=frozen,
        )

    defect_system = DefectSystem(
        defect_species=input_data["defect_species"],
        dos=input_data["dos"],
        volume=input_data["volume"],
        temperature=input_data["temperature"],
    )
    defect_system.report(conv=args.convergence_tolerance)


if __name__ == "__main__":
    main()
