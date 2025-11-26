#!/usr/bin/python
"""
Parses previous VASP data and launches Δ-Sol static calculations 
for structures not already rejected.
"""

import os
import typing as tp
from datetime import datetime
import argparse as argp

from pymatgen.core import Structure

# LOCAL IMPORTS
from . import CONFIGPATH, _parse_input_args
from .utils import (
    check_type, check_num_value, check_file_or_dir, add_new_dir, yaml_loader,
    PMGStaticSet, extract_vasp_data_for_delta_sol_init,
    get_dsol_struct_dir, calc_idx_to_dir_name, dsol_calc_init,
    write_and_run_vasp
)


########################################
# ARGUMENTS HANDLING

def _get_command_line_args() -> argp.Namespace:
    parser = argp.ArgumentParser(prog=os.path.basename(__file__), description=__doc__)
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
        default=PMGStaticSet.MPSTATICSET.value,
        help=(
            "The pymatgen preset to use for VASP static calculations.\n"
            f"Supported presets: {PMGStaticSet}."
            "More info on possible presets in pymatgen documentation:\n"
            "https://pymatgen.org/pymatgen.io.vasp.html#pymatgen.io.vasp.sets."
        )
    )
    parser.add_argument(
        "-u", "--user-settings",
        default="DSolStaticSet.yaml",
        help=(
            "Name of the YAML file containing user defined VASP tags that will override "
            f"those of the preset. Given file must be located in {CONFIGPATH} to be found. "
            "Defaults to %(default)s."
        ),
        metavar="<filename>",
    )
    parser.add_argument(
        "--with-uncertainties",
        action="store_true",
        help=(
            "Pass this flag to enable computation of minimal and maximal Δ-Sol band gaps.\n"
            "This will need two more VASP static total energy computation for each limit."
        )
    )

    args: argp.Namespace = parser.parse_args()
    return args


def _process_input_args(args_dict: dict[str, tp.Any]) -> dict[str, tp.Any]:
    """Handle input arguments assertions and processing."""
    if args_dict is None:
        raise ValueError(f"No arguments found at '{os.path.basename(__file__)}' script call.")
    
    check_type(args_dict, "args_dict", (dict,))

    # Set default values for unset optional arguments
    args_dict = {k: v for k, v in args_dict.items() if v is not None}
    args_dict.setdefault("executable_path", "vasp")
    # WARNING:
    # Usual VASP shortcut, but may activate wrong VASP version if several are installed.
    # Prefer giving a true VASP executable path for unambiguous computation.
    default_output = os.path.join(os.path.dirname(args_dict.get("input_dir", "")), "Band_gaps")
    args_dict.setdefault("output", default_output)
    args_dict.setdefault("preset", PMGStaticSet.MPSTATICSET.value)
    args_dict.setdefault("user_settings", "DSolStaticSet.yaml")
    args_dict.setdefault("with_uncertainties", False)

    # Assert set arguments conformity
    check_file_or_dir(args_dict.get("input_dir"), "dir")

    if not str(args_dict.get("executable_path")).startswith("vasp"):
        check_file_or_dir(args_dict.get("executable_path"), "file")
    
    check_type(args_dict.get("task_id"), "task_id", (int,))
    check_num_value(args_dict.get("task_id"), "task_id", ">=", 0)
    check_file_or_dir(args_dict.get("output"), "dir")

    if args_dict.get("prev_summary") is not None:
        check_file_or_dir(args_dict.get("prev_summary"), "file", allowed_formats="json")

        if args_dict.get("key_to_check") is None:
            raise ValueError(
                f"'prev_summary' argument was provided ({args_dict.get('prev_summary')}), "
                "therefore 'key-to-check' argument has to be given as well."
            )

    if args_dict.get("key_to_check") is not None:
        check_type(args_dict.get("key_to_check"), "key_to_check", (str,))
    
    assert args_dict["preset"] in PMGStaticSet.values or args_dict.get("preset") == "DSolStaticSet", (
    "Provided static preset must be one of the following:\n"
    f"{PMGStaticSet.values} or 'DsolStaticSet'."
    )
    config_path = os.path.join(CONFIGPATH, args_dict.get("user_settings", "DSolStaticSet.yaml"))
    check_file_or_dir(config_path, "file", allowed_formats=("yml", "yaml"))
    check_type(args_dict["with_uncertainties"], "with-uncertainties", (bool,))

    # Additional arguments processing
    os.makedirs(args_dict["output"], exist_ok=True)

    if args_dict.get("preset") != "DsolStaticSet":
        args_dict["preset"] = PMGStaticSet(args_dict.get("preset"))

    args_dict["settings"] = yaml_loader(config_path)

    return args_dict


########################################


def main(standalone: bool = True, **kwargs) -> None:
    """
    Parses previous VASP data and launches Δ-Sol static calculations for structures not
    already rejected.

    Args:
        standalone (bool):          Whether parsed script is used directly through
                                    command-line (stand-alone script) or in an external
                                    pipeline script.

        input_dir (str|Path):       Base directory containing structure directories.

        executable_path (str|Path): Path to the VASP executable.

        task_id (int):              Provide here the job array task ID that will treat one
                                    calculation for one structure. The total number of jobs
                                    should be the number of structures multiplied by the number
                                    of calculations for one structure (3 for a direct estimation
                                    only, 7 with uncertainties).

        output (str|Path):          Path to the output directory where VASP files will be written.
                                    A subdirectory will be created in this directory for each Δ-Sol
                                    calculation. If not provided, a 'Band_gaps' directory is
                                    created at the same path as the directory provided in the
                                    'input_dir' argument.

        prev_summary (str|Path):    Path to a JSON summary file produced by a previous screening
                                    step. If given, the file will be checked to filter structures
                                    that are already rejected.

        key_to_check (str):         The dict key associated to the bool used to verify eligibility
                                    in the previous summary file. If prev_summary is given, it must
                                    be given too.

        preset (str):               The pymatgen preset to use for VASP static calculations.
                                    List of currently supported presets can be checked in command
                                    line within the --help documentation of this argument for this
                                    script. More info on possible presets in pymatgen documentation:
                                    https://pymatgen.org/pymatgen.io.vasp.html#pymatgen.io.vasp.sets.

        user_settings (str):        Name of the YAML file containing user defined VASP tags that will
                                    override those of the preset. Given file must be located in
                                    'metrics_pipeline/config' to be found.
                                    Defaults to 'DSolStaticSet.yaml'.

        with_uncertainties (bool):  Pass this flag to enable computation of minimal and maximal Δ-Sol
                                    band gaps. This will need two more VASP static total energy
                                    computation for each limit.
    """
    start = datetime.now()
    args = _parse_input_args(_get_command_line_args, _process_input_args, standalone, **kwargs)

    # Variables deduced from args
    struct_path, struct_idx, calc_idx = get_dsol_struct_dir(
        path=args["input_dir"],
        task_id=args["task_id"],
        with_uncertainties=args["with_uncertainties"]
    )

    struct_data = extract_vasp_data_for_delta_sol_init(
        struct_dir=struct_path, path_to_summary=args.get("prev_summary")
    )
    if not struct_data:
        raise ValueError(
            "The structure data corresponding to given 'task-id' argument "
            "was not found or is already rejected in the summary file from previous step.\n"
            f"Given task-id argument: {args['task_id']}\n"
            f"Corresponding structure index: {struct_idx}\n"
            f"Corresponding calculation ID: {calc_idx}\n"
            "(0 = E(N0), 1-2 = E(N0 +/- n(best)), 3-4 = E(N0 +/- n(min)), 5-6 = E(N0 +/- n(max)))."
        )

    input_data = dsol_calc_init(
        structure=struct_data[1]["structure"],
        calc_index=calc_idx,
        preset=args["preset"],
        user_corrections=args.get("settings")
    )
    structure: Structure = struct_data[1]["structure"]
    dir_name   = f"{struct_idx}_{structure.composition.reduced_formula}" # type: ignore
    calc_name  = calc_idx_to_dir_name(dir_name, calc_idx)
    calc_dir   = add_new_dir(args["output"], dir_name, calc_name)

    write_and_run_vasp(input_data, calc_dir, args["executable_path"])

    stop = datetime.now()
    print(f"elapsed time: {stop-start}")


if __name__ == "__main__":
    main()
