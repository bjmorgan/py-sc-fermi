from py_sc_fermi.inputs import InputSet
from py_sc_fermi.defect_system import DefectSystem
import argparse
import yaml

parser = argparse.ArgumentParser()
parser.add_argument(
    "input_file", type=str, help="Path to input file defining the defect system"
)
parser.add_argument(
    "-s",
    "--structure_file",
    help="Path to structure file giving the volume of a defect system",
    default="unitcell.dat",
)
parser.add_argument(
    "-d",
    "--dos_file",
    help="Path to file specifying the totdos of the system",
    default="totdos.dat",
)
parser.add_argument(
    "-f",
    "--frozen_defects",
    help="frozen defects present in the defect system",
    action="store_true",
)
parser.add_argument(
    "-c", "--convergence_tol", help="convergence tolerance", type=float, default=1e-19
)
parser.add_argument(
    "-n", "--n_trial", help="maximum number of trial steps", type=int, default=1500
)
args = parser.parse_args()
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
    conv = args.convergence_tol
    n_trial = args.n_trial

    if input.endswith(".yaml"):
        defect_system = DefectSystem.from_yaml(input)
    else:
        input_data = InputSet.from_sc_fermi_inputs(
            input_file=input,
            structure_file=structure,
            dos_file=totdos,
            frozen=frozen,
        )
        defect_system = DefectSystem.from_input_set(input_data)
    defect_system.report()

    dump_dict = defect_system.as_dict(decomposed=True)
    dump_dict["temperature"] = defect_system.temperature
    with open("py_sc_fermi_out.yaml", "w") as f:
        yaml.dump(dump_dict, f)


if __name__ == "__main__":
    main()
