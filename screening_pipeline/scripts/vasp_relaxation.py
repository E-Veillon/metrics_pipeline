#!/usr/bin/python
"""
A script that writes and performs VASP relaxation on structure data read from CIF files,
using pymatgen as a setting interface between raw data and VASP.
The relaxation results then may be used in other scripts for material properties analysis.
"""

########################################
# SYSTEM I/O MODULES

import os
from datetime import datetime
from argparse import ArgumentParser, Namespace, RawTextHelpFormatter
from monty.os import cd


########################################
# PYTHON MATERIAL GENOMICS PACKAGE

from pymatgen.core.structure import SiteCollection

########################################
# LOCAL MODULES

from screening_pipeline.utils import (
    _yaml_loader, PMGRelaxSet, read_cif, vasp_relaxation_settings
)

########################################
# LOCAL FUNCTIONS

def assert_args(args: Namespace) -> None:
    """Asserting input arguments validity."""

    assert os.path.exists(args.input_file), (
    f"{args.input_file}: path to input file not found."
    )
    assert os.path.isfile(args.input_file), (
    f"{args.input_file} found but it is not a file."
    )
    assert args.input_file.endswith(".cif"), (
    "Input structure data must be in CIF format."
    )
    assert args.executable_path.startswith("vasp") or os.path.exists(args.executable_path), (
    f"{args.executable_path}: executable file not found."
    )
    assert os.path.exists(args.output), (
    f"{args.output}: path to output directory not found."
    )
    assert os.path.isdir(args.output), (
    f"{args.output} found but it is not a directory."
    )
    assert args.preset in PMGRelaxSet, (
    "Provided relaxation preset must be one of the following:\n"
    f"{PMGRelaxSet}"
    )
    assert os.path.isfile(args.user_settings), (
    f"{args.user_settings}: file not found."
    )
    assert args.user_settings.endswith(".yaml"), (
    "user settings file must be of .yaml format."
    )
    assert args.workers >= 1, (
    "'workers' argument value must be strictly positive."
    )

    assert args.task_index >= 0 or args.task_index is None, (
    "'task_index' argument value must be positive or zero."
    )

########################################
# MAIN FUNCTION

def main():
    """Main function."""
    start = datetime.now()

    # ARGUMENTS PARSING BLOCK

    prog_name = "vasp_relaxation.py"
    prog_description = """
    A script that writes and performs VASP relaxation on structure data read from CIF files,
    using pymatgen as a setting interface between raw data and VASP.
    The relaxation results then may be used in other scripts for material properties analysis.
    """

    helper_format = RawTextHelpFormatter

    parser = ArgumentParser(
        prog=prog_name,
        description=prog_description,
        formatter_class=helper_format
    )

    parser.add_argument(
        "input_file",
        type=str,
        help="Path to the CIF file containing structure data to read.",
        metavar="input_file.cif"
    )
    parser.add_argument(
        "executable_path",
        type=str,
        help="Path to the VASP executable.",
        metavar="PATH"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="./",
        help=(
            "Path to the output directory where VASP files will be written. "
            "A subdirectory will be created in output directory "
            "for each structure found in input_file."
        ),
        metavar="outdir"
    )
    parser.add_argument(
        "-p",
        "--preset",
        type=str,
        default="MPRelaxSet",
        help=(
            "The pymatgen preset to use for VASP relaxation. "
            "More info on possible presets in pymatgen documentation:\n"
            "https://pymatgen.org/pymatgen.io.vasp.html#pymatgen.io.vasp.sets."
        ),
        metavar="RelaxSet"
    )
    parser.add_argument(
        "-u",
        "--user-settings",
        type=str,
        default="user_settings.yaml",
        help="Path to the .yaml file containing tags overrides to put over the PMG preset.",
        metavar="file.yaml",
        dest="user_settings"
    )
    parser.add_argument(
        "-t",
        "--task_index",
        type=int,
        default=0,
        help=(
            "If a job array is used, provide here the structure index "
            "to treat according to task IDs (e.g. if task ID 0 treats "
            "structure 0 and so on, just provide the task ID).\n"
            "If not given, it will default to the first index possible, i.e. index 0."
        ),
    )

    args: Namespace = parser.parse_args()

    assert_args(args)

    user_settings = _yaml_loader(args.user_settings)
    struct_idx    = args.task_index


    # MAIN BLOCK

    # Convert CIF data into Structure objects
    structures, *_ = read_cif(
        filename=args.input_file,
        keep_rare_gases=True, # Avoid calling rare gaz screening function
        keep_rare_earths=True # Avoid calling rare earth screening function
    )

    structure: SiteCollection = structures[struct_idx]
    dir_name = f"{struct_idx}_{structure.composition.reduced_formula}"

    vasp_input = vasp_relaxation_settings(
        structure=structure,
        preset=args.preset,
        user_corrections=user_settings
    )

    #vasp_launcher(
    #   vasp_exe=exe_path,
    #   vasp_input=vasp_input,
    #   path=os.path.join(outdir, dir_name)
    #)
    run_dir = os.path.join(args.output, dir_name)
    os.makedirs(run_dir, exist_ok=True)
    vasp_input.write_input(output_dir=run_dir)
    with cd(run_dir):
        os.system(f"{args.executable_path}")

    stop = datetime.now()
    print(f"Elapsed time: {stop-start}")


if __name__=="__main__":
    main()
