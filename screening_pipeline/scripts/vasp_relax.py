#!/usr/bin/python
"""
A script that writes and performs VASP relaxation on structure data read from CIF files,
using pymatgen as a setting interface between raw data and VASP.
The relaxation results then may be used in other scripts for material properties analysis.
"""

import os
import json
from datetime import datetime
from argparse import ArgumentParser, Namespace

# PYTHON MATERIAL GENOMICS
from pymatgen.core import SiteCollection

# LOCAL IMPORTS
from screening_pipeline.utils import (
    CONFIGPATH, check_file_or_dir, add_new_dir, yaml_loader,
    PMGRelaxSet, read_cif, vasp_relaxation_settings, write_and_run_vasp,
    get_struct_from_vasp
)


########################################


def _assert_args(args: Namespace) -> None:
    """Asserting input arguments validity."""

    check_file_or_dir(args.input_file, "file", allowed_formats=("cif", "json"))

    assert (
        args.executable_path.startswith("vasp")
        or os.path.exists(args.executable_path)
    ), f"{args.executable_path}: executable file not found."

    if args.output is not None:
        check_file_or_dir(args.output, "dir")

    assert args.preset in PMGRelaxSet, (
    "Provided relaxation preset must be one of the following:\n"
    f"{PMGRelaxSet.values}"
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

    parser = ArgumentParser(
        prog=prog_name,
        description=prog_description
    )

    parser.add_argument(
        "input_file",
        help=(
            "Path to the file containing structure data to read. "
            "It can be a CIF file to read structure data directly, or "
            "a JSON summary file from a previous pipeline step to extract "
            "structure data from a previous VASP run."
        )
    )
    parser.add_argument(
        "executable_path",
        help="Path to the VASP executable.",
        metavar="PATH"
    )
    parser.add_argument(
        "-o", "--output",
        help=(
            "Path to the output directory where VASP files will be written. "
            "A subdirectory will be created in output directory "
            "for each structure found in input_file."
        ),
        metavar="outdir"
    )
    parser.add_argument(
        "-p", "--preset",
        type=PMGRelaxSet,
        default=PMGRelaxSet.MPRELAXSET,
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
        default="default_settings.yaml",
        help="Path to the .yaml file containing tags overrides to put over the PMG preset.",
        metavar="file.yaml",
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

    args: Namespace = parser.parse_args()

    _assert_args(args)

    if args.output is None:
        outdir = os.path.join(os.path.dirname(args.input_file), "Relaxations")
    else:
        outdir = args.output

    user_settings = yaml_loader(os.path.join(CONFIGPATH, args.user_settings))
    struct_idx    = args.task_index


    # MAIN BLOCK

    if args.input_file.endswith(".cif"):
        # Convert CIF data into Structure objects
        structures, *_ = read_cif(
            filename=args.input_file,
            workers=args.workers,
            keep_rare_gases=True, # Avoid calling rare gaz screening function
            keep_rare_earths=True # Avoid calling rare earth screening function
        )
        # Get structure of interest
        structure: SiteCollection = structures[struct_idx]
        dir_name = f"{struct_idx}_{structure.composition.reduced_formula}"

    elif args.input_file.endswith(".json"):

        with open(args.input_file, "rt", encoding="utf-8") as fp:
            data = json.load(fp)

        paths = [d["path"] for d in data]

        try:
            struct_dir = next(filter(
                lambda path: os.path.basename(path).startswith(f"{struct_idx}_"),
                paths
            ))
        except StopIteration as exc:
            raise ValueError(
                f"Given 'task_index' value ({struct_idx}) do not match any data "
                f"in file '{args.input_file}'."
            ) from exc

        structure = get_struct_from_vasp(struct_dir, try_xdatcar=False)
        dir_name = os.path.basename(struct_dir)

    vasp_input = vasp_relaxation_settings(
        structure=structure,
        preset=args.preset,
        user_corrections=user_settings
    )
    run_dir = add_new_dir(outdir, dir_name)

    write_and_run_vasp(vasp_input, run_dir, args.executable_path)

    stop = datetime.now()
    print(f"Elapsed time: {stop-start}")


if __name__=="__main__":
    main()
