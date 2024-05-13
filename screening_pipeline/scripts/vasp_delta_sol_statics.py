#!/usr/bin/python
"""
Parses previous VASP data and launches Δ-Sol static calculations for structures not already rejected.
"""

########################################
# SYSTEM I/O MODULES

import os
import sys
from monty.os import cd
from datetime import datetime
from argparse import ArgumentParser, Namespace, RawTextHelpFormatter

########################################
# PYTHON MATERIALS GENOMICS PACKAGE


########################################
# LOCAL MODULES

from screening_pipeline.utils.utils import _yaml_loader
from screening_pipeline.utils.custom_types import PMGStaticSet
from screening_pipeline.utils.vasp_io import (
    extract_vasp_data_for_delta_sol_init, delta_sol_calculation_init, vasp_launcher
)
########################################
# LOCAL FUNCTIONS

def assert_args(args: Namespace) -> None:

    assert os.path.isdir(args.input_dir), (
        f"{args.input_dir}: No directory found."
    )
    assert args.executable_path.startswith("vasp") or os.path.exists(args.executable_path), (
        f"{args.executable_path}: Executable file not found."
    )
    assert args.task_id >= 0, (
        "task-id argument must be positive or zero."
    )
    assert os.path.isdir(args.output), (
        f"{args.output}: No directory found."
    )
    if args.prev_summary is not None:
        assert os.path.isfile(args.prev_summary), (
            f"{args.previous_results}: No such file found."
        )

    #PMGStaticSet.add("GenMatStatic54Set")

    assert args.preset in PMGStaticSet, (
    "Provided static preset must be one of the following:\n"
    f"{PMGStaticSet}"
    )
    assert os.path.isfile(args.user_settings) or args.user_settings is None, (
    f"{args.user_settings}: file not found."
    )

    assert args.user_settings.endswith(".yaml"), (
    "user settings file must be of .yaml format."
    )

def calc_idx_to_dir_name(calc_index: int) -> str:
    assert isinstance(calc_index, int), f"Expected 'int' type, got '{type(calc_index)}' instead."
    if calc_index == 0: return "_neutral"
    if calc_index == 1: return "_best_plus"
    if calc_index == 2: return "_best_minus"
    if calc_index == 3: return "_min_plus"
    if calc_index == 4: return "_min_minus"
    if calc_index == 5: return "_max_plus"
    if calc_index == 6: return "_max_minus"
    else: raise ValueError(f"Only int from 0 to 6 supported, got {calc_index}")

########################################
# MAIN FUNCTION

def main():
    start = datetime.now()

    # ARGUMENTS PARSING BLOCK
    
    prog_name = "vasp_delta_sol_statics.py"
    prog_desc = """
        Parses previous VASP data and launches Δ-Sol static calculations for structures not already rejected.
        Uses job arrays properties to maximize the parallelization efficiency.

        Reference for Δ-Sol method:
            M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010)
            (reference 32 in screening_pipeline/Bibliography)
        """
    prog_missing_steps = """
        Missing steps to complete this script:
            - Mandatory job array id
            - First verify that corresponding structure is not rejected
            - If rejected, stop the sub-job with a simple message in output about it
            - Else, extract the structure, prepare corresponding calculation and launch it
        """
    helper_format = RawTextHelpFormatter

    parser = ArgumentParser(
        prog=prog_name, 
        description=prog_desc, 
        epilog=prog_missing_steps, 
        formatter_class=helper_format
    )

    parser.add_argument(
        "input_dir",
        type=str,
        help="Base directory containing structure directories.", 
    )
    parser.add_argument(
        "executable_path", 
        type=str,
        help="Path to the VASP executable.", 
    )
    parser.add_argument( 
        "task_id", 
        type=int, 
        help=(
            "Provide here the job array task ID that will treat one calculation for one structure.\n"
            "The total number of jobs should be the number of structures multiplied by the number of\n"
            "calculations for one structure (3 for a direct estimation only, 7 with uncertainties)."
        )
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default=None,
        help=(
            "Path to the output directory where VASP files will be written.\n"
            "A subdirectory will be created in this directory for each Δ-Sol calculation.\n"
            "If not provided, a 'Band_gaps' directory is created at the same path as the directory "
            "provided in the 'input_dir' argument."
        ), 
        metavar="outdir"
    )
    parser.add_argument(
        "-R", "--read-previous-summary",
        type=str,
        default=None,
        help=(
            "Path to a JSON summary file produced by a previous screening step.\n"
            "If given, the file will be checked to filter structures that are already rejected."
        ),
        metavar="/path/to/summary.json",
        dest="prev_summary"
    )
    parser.add_argument(
        "-p", "--preset", 
        type=str, 
        default="MPStaticSet", 
        help=(
            "The pymatgen preset to use for VASP static calculations.\n"
            f"Supported presets: {PMGStaticSet}"
            "More info on possible presets in pymatgen documentation:\n"
            "https://pymatgen.org/pymatgen.io.vasp.html#pymatgen.io.vasp.sets."
        )
    )
    parser.add_argument(
        "-u", "--user-settings", 
        type=str, 
        default=None, 
        help="Path to the .yaml file containing user defined VASP tags that will override those of the preset.", 
        metavar="file.yaml",
    )
    parser.add_argument(
        "--with-uncertainties", 
        action="store_true", 
        help=(
            "Pass this flag to enable computation of minimal and maximal Δ-Sol band gaps.\n"
            "This will need two more VASP static total energy computation for each limit."
        )
    )

    args: Namespace = parser.parse_args()

    assert_args(args)

    # Positional args
    input_dir = args.input_dir
    exe_path  = args.executable_path
    task_id   = args.task_id

    # Optional args
    if args.output is None:
        outdir = os.path.join(os.path.dirname(input_dir), "Band_gaps")
    else:
        outdir = args.output

    prev_summary = args.prev_summary or None
    preset = args.preset
    user_settings = _yaml_loader(args.user_settings)

    # Variables coming from args
    tasks_per_struct = 7 if args.with_uncertainties else 3
    struct_idx = task_id // tasks_per_struct
    calc_idx = task_id % tasks_per_struct

    try:
        struct_dir = next(
            filter(
            lambda dirname: dirname.startswith(f"{struct_idx}_"), 
            os.listdir(input_dir)
            )
        )
    except StopIteration:
        print(
            "No structure directory found with index corresponding to given 'task-id' argument.\n"
            f"Searched directory: {input_dir}\n"
            f"Given task-id argument: {args.task_id}\n"
            f"Corresponding structure index: {struct_idx}\n"
            f"Corresponding calculation ID: {calc_idx}\n"
            "(0 = E(N0), 1-2 = E(N0 +/- n(best)), 3-4 = E(N0 +/- n(min)), 5-6 = E(N0 +/- n(max)))."
        )
        sys.exit(0)

    struct_path = os.path.join(input_dir, struct_dir)


    # MAIN BLOCK

    struct_data = extract_vasp_data_for_delta_sol_init(
        struct_dir=struct_path, path_to_summary=prev_summary
    )
    if not struct_data:
        print(
            "The structure data corresponding to given 'task-id' argument "
            "was not found or is already rejected in the summary file from previous step.\n"
            f"Given task-id argument: {args.task_id}\n"
            f"Corresponding structure index: {struct_idx}\n"
            f"Corresponding calculation ID: {calc_idx}\n"
            "(0 = E(N0), 1-2 = E(N0 +/- n(best)), 3-4 = E(N0 +/- n(min)), 5-6 = E(N0 +/- n(max)))."
        )
        sys.exit(0)

    dir_name   = f"{struct_idx}_{struct_data[1]['structure'].composition.reduced_formula}"
    calc_name  = calc_idx_to_dir_name(calc_idx)
    calc_dir   = os.path.join(outdir, dir_name, ''.join((dir_name, calc_name)))

    input_data = delta_sol_calculation_init(
        structure=struct_data[1]["structure"], 
        calc_index=calc_idx, 
        preset=preset, 
        user_corrections=user_settings
    )

    os.makedirs(calc_dir, exist_ok=True)

    #vasp_launcher(vasp_exe=exe_path, path=calc_dir, vasp_input=input_data)
    input_data.write_input(output_dir=calc_dir)
    with cd(calc_dir):
        os.system(f"{exe_path}")

    stop = datetime.now()
    print(f"elapsed time: {stop-start}")


if __name__ == "__main__":
    main()
