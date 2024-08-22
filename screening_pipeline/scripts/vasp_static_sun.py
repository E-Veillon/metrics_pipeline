#!/usr/bin/python
"""
A script that writes and performs VASP static run on structure data read from CIF files,
using pymatgen as a setting interface between raw data and VASP.
Energy results may be used for phase diagrams computations.
"""

import os
from datetime import datetime
import argparse

# PYTHON MATERIAL GENOMICS
from pymatgen.core.structure import SiteCollection

# LOCAL IMPORTS
from screening_pipeline.utils import (
    CONFIGPATH, check_file_or_dir, add_new_dir, yaml_loader,
    PMGStaticSet, read_cif, vasp_static_settings, write_and_run_vasp
)


########################################


def _assert_args(args: argparse.Namespace) -> None:
    """Asserting input arguments validity."""

    check_file_or_dir(args.input_file, "file", allowed_formats="cif")

    assert (
        args.executable_path.startswith("vasp")
        or os.path.exists(args.executable_path)
    ), f"{args.executable_path}: executable file not found."

    check_file_or_dir(args.output, "dir")

    assert args.preset in PMGStaticSet, (
    "Provided relaxation preset must be one of the following:\n"
    f"{PMGStaticSet.values}"
    )

    settings_path = os.path.join(CONFIGPATH, args.user_settings)
    check_file_or_dir(settings_path, "file", allowed_formats="yaml")

    assert args.workers >= 1, (
    "'workers' argument value must be strictly positive."
    )

    if args.task_index is not None:
        assert args.task_index >= 0, (
        "'task_index' argument value must be positive or zero."
        )


########################################


def main():
    """Main function."""
    start = datetime.now()

    # ARGUMENTS PARSING BLOCK

    prog_name = "vasp_static_sun.py"
    prog_description = (
        "A script that writes and performs VASP static runs on structure data "
        "read from CIF files, using pymatgen as a setting interface between "
        "raw data and VASP.\n"
        "Energy results may be used for phase diagrams computations."
    )

    parser = argparse.ArgumentParser(
        prog=prog_name,
        description=prog_description
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
        "-o", "--output",
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
        "-p", "--preset",
        type=PMGStaticSet,
        default=PMGStaticSet.MPSTATICSET,
        help=(
            "The pymatgen preset to use for VASP static run. "
            "More info on possible presets in pymatgen documentation:\n"
            "https://pymatgen.org/pymatgen.io.vasp.html#pymatgen.io.vasp.sets."
        ),
        metavar="StaticSet"
    )
    parser.add_argument(
        "-u",
        "--user-settings",
        type=str,
        default="default_settings.yaml",
        help=(
            "Name of the .yaml file containing tags overrides to put over the PMG preset.\n"
            "The file must be at location screening_pipeline/config to be found."
        ),
        metavar="FILENAME",
        dest="user_settings"
    )
    parser.add_argument(
        "-w", "--workers",
        type=int,
        default=1,
        help="Number of parallel processes to spawn for parallelized steps."
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

    args: argparse.Namespace = parser.parse_args()

    _assert_args(args)

    settings = os.path.join(CONFIGPATH, args.user_settings)
    user_settings = yaml_loader(settings)
    struct_idx    = args.task_index


    # MAIN BLOCK

    # Convert CIF data into Structure objects
    structures, *_ = read_cif(
        filename=args.input_file,
        workers=args.workers,
        keep_rare_gases=True, # Avoid calling rare gaz screening function
        keep_rare_earths=True # Avoid calling rare earth screening function
    )

    structure: SiteCollection = structures[struct_idx]
    dir_name = f"{struct_idx}_{structure.composition.reduced_formula}"

    vasp_input = vasp_static_settings(
        structure=structure,
        preset=args.preset,
        user_corrections=user_settings
    )
    run_dir = add_new_dir(args.output, dir_name)

    write_and_run_vasp(vasp_input, run_dir, args.executable_path)

    stop = datetime.now()
    print(f"Elapsed time: {stop-start}")


if __name__=="__main__":
    main()
