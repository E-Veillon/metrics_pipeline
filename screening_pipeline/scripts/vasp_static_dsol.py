#!/usr/bin/python
"""
Parses previous VASP data and launches Δ-Sol static calculations 
for structures not already rejected.
"""

import os
from datetime import datetime
import argparse

# LOCAL IMPORTS
from screening_pipeline.utils import (
    CONFIGPATH, check_file_or_dir, add_new_dir, yaml_loader,
    PMGStaticSet, extract_vasp_data_for_delta_sol_init,
    get_dsol_struct_dir, calc_idx_to_dir_name, dsol_calc_init,
    write_and_run_vasp
)


########################################


def _assert_args(args: argparse.Namespace) -> None:
    """Asserting input arguments validity."""

    check_file_or_dir(args.input_dir, "dir")

    assert (
        args.executable_path.startswith("vasp")
        or os.path.exists(args.executable_path)
    ), f"{args.executable_path}: Executable file not found."

    assert args.task_id >= 0, (
        "task-id argument must be positive or zero."
    )

    check_file_or_dir(args.output, "dir")

    if args.prev_summary is not None:
        check_file_or_dir(args.prev_summary, "file", allowed_formats="json")

    assert args.preset in PMGStaticSet.values or args.preset == "DsolStaticSet", (
    "Provided static preset must be one of the following:\n"
    f"{PMGStaticSet.values} or 'DsolStaticSet'."
    )

    if args.user_settings is not None:
        settings_path = os.path.join(CONFIGPATH, args.user_settings)
        check_file_or_dir(settings_path, "file", allowed_formats="yaml")


########################################


def main() -> None:
    """Main function."""
    start = datetime.now()

    # ARGUMENTS PARSING BLOCK

    prog_name = "vasp_delta_sol_statics.py"
    prog_desc = """
        Parses previous VASP data and launches Δ-Sol static calculations for structures not already rejected.
        Uses job arrays properties to maximize the parallelization efficiency.

        Reference for Δ-Sol method:
            M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010)
        """
    prog_missing_steps = """
        Missing steps to complete this script:
            - Mandatory job array id
            - First verify that corresponding structure is not rejected
            - If rejected, stop the sub-job with a simple message in output about it
            - Else, extract the structure, prepare corresponding calculation and launch it
        """

    parser = argparse.ArgumentParser(
        prog=prog_name,
        description=prog_desc,
        epilog=prog_missing_steps
    )

    parser.add_argument(
        "input_dir",
        help="Base directory containing structure directories.",
    )
    parser.add_argument(
        "executable_path",
        help="Path to the VASP executable.",
    )
    parser.add_argument(
        "task_id",
        type=int,
        help=(
            "Provide here the job array task ID that will treat "
            "one calculation for one structure.\n"
            "The total number of jobs should be the number of structures "
            "multiplied by the number of calculations for one structure "
            "(3 for a direct estimation only, 7 with uncertainties)."
        )
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
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
        help=(
            "Path to a JSON summary file produced by a previous screening step.\n"
            "If given, the file will be checked to filter structures that are already rejected."
        ),
        metavar="/path/to/summary.json",
        dest="prev_summary"
    )
    parser.add_argument(
        "-k", "--key-to-check",
        help=(
            "The dict key associated to the bool used to verify eligibility "
            "in the previous summary file.\n"
            "If --read-previous-summary is given, it must be given too."
        )
    )
    parser.add_argument(
        "-p", "--preset",
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
        help=(
            "Path to the .yaml file containing user defined VASP tags "
            "that will override those of the preset."
        ),
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

    args: argparse.Namespace = parser.parse_args()

    _assert_args(args)

    # Positional args
    input_dir = args.input_dir
    exe_path  = args.executable_path
    task_id   = args.task_id

    if preset != "DsolStaticSet":
        preset = PMGStaticSet(preset)

    # Optional args
    if args.output is None:
        outdir = os.path.join(os.path.dirname(input_dir), "Band_gaps")
    else:
        outdir = args.output

    prev_summary = args.prev_summary or None
    preset = args.preset
    user_settings = yaml_loader(args.user_settings)

    # Variables deduced from args
    struct_path, struct_idx, calc_idx = get_dsol_struct_dir(
        input_dir, task_id, with_uncertainties=args.with_uncertainties
    )


    # MAIN BLOCK

    struct_data = extract_vasp_data_for_delta_sol_init(
        struct_dir=struct_path, path_to_summary=prev_summary
    )
    if not struct_data:
        raise ValueError(
            "The structure data corresponding to given 'task-id' argument "
            "was not found or is already rejected in the summary file from previous step.\n"
            f"Given task-id argument: {args.task_id}\n"
            f"Corresponding structure index: {struct_idx}\n"
            f"Corresponding calculation ID: {calc_idx}\n"
            "(0 = E(N0), 1-2 = E(N0 +/- n(best)), 3-4 = E(N0 +/- n(min)), 5-6 = E(N0 +/- n(max)))."
        )

    input_data = dsol_calc_init(
        structure=struct_data[1]["structure"],
        calc_index=calc_idx,
        preset=preset,
        user_corrections=user_settings
    )

    dir_name   = f"{struct_idx}_{struct_data[1]['structure'].composition.reduced_formula}"
    calc_name  = calc_idx_to_dir_name(dir_name, calc_idx)
    calc_dir   = add_new_dir(outdir, dir_name, calc_name)

    write_and_run_vasp(input_data, calc_dir, exe_path)

    stop = datetime.now()
    print(f"elapsed time: {stop-start}")


if __name__ == "__main__":
    main()
