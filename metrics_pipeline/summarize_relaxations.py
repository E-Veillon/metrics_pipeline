#!/usr/bin/python
"""
Make a JSON convergence summary file out of VASP relax run directories.
"""

import typing as tp
import argparse as ap
from pathlib import Path

from src.utils import parse_input_args, check_type, check_num_value, VisualIterator
from src.io import (
    check_file_or_dir, check_file_format,
    VaspParser, JsonWriter
)


def _get_cmd_line_args() -> ap.Namespace:
    """Handle Command Line Interface (CLI) arguments."""
    parser = ap.ArgumentParser(prog=Path(__file__).name, description=__doc__)
    parser.add_argument(
        "base_dir",
        help="Base directory containing VASP run directories."
    )
    parser.add_argument(
        "-p", "--prev-summary", metavar="<path>",
        help=(
            "Path to a JSON summary file produced by a previous screening step. "
            "If given, the file will be checked to filter out structures that are "
            "already rejected."
        )
    )
    parser.add_argument(
        "-k", "--summary-key",
        help=(
            "The dict key associated to the bool used to verify eligibility in previous "
            "summary file. If --prev-summary is given, it must be given too."
        )
    )
    parser.add_argument(
        "-s", "--summary", default="relax_summary.json",
        help=(
            "Name of output JSON summary file indicating calculation results for this step. "
            "Also usable by further steps to filter out structures rejected in this step. "
            "The file will be automatically written in the 'base_dir' directory."
        )
    )
    parser.add_argument(
        "-w", "--workers", type=int, metavar="int",
        help=(
            "Number of processes to use in parallel. If not given, will use default of "
            "`tqdm.contrib.concurrent.process_map()`. Pass 0 to disable `process_map()` "
            "and execute sequentially."
        )
    )
    args = parser.parse_args()
    return args


def _process_input_args(args_dict: dict[str, tp.Any]) -> dict[str, tp.Any]:
    """Handle input arguments assertions and processing."""
    if args_dict is None:
        raise ValueError(f"No arguments found at '{Path(__file__).name}' script call.")

    check_type(args_dict, "args_dict", (dict,))

    # Set default values for unset optional arguments
    args_dict = {k: v for k, v in args_dict.items() if v is not None}
    args_dict.setdefault("summary", "relax_summary.json")

    # Assert set arguments conformity
    check_file_or_dir(args_dict.get("base_dir"), "dir")

    if args_dict.get("prev_summary") is not None:
        check_file_or_dir(args_dict.get("prev_summary"), "file", allowed_formats="json")

        if args_dict.get("summary_key") is None:
            raise ValueError(
                f"'prev_summary' argument was provided ({args_dict.get('prev_summary')}), "
                "therefore 'key-to-check' argument has to be given as well."
            )

    if args_dict.get("summary_key") is not None:
        check_type(args_dict.get("summary_key"), "summary_key", (str,))

    check_file_format(args_dict.get("summary"), allowed_formats="json")

    if args_dict.get("workers") is not None:
        check_type(args_dict.get("workers"), "workers", (int,))
        check_num_value(args_dict.get("workers"), "workers", ">", 0)

    # Additional arguments processing
    args_dict["summarypath"] = str(Path(args_dict["base_dir"], args_dict["summary"]).resolve())

    return args_dict


def main(standalone: bool = True, **kwargs) -> None:
    """
    Make a JSON convergence summary file out of VASP relax run directories.

    Parameters
    ----------
    standalone: bool
        Whether parsed script is used directly through command-line (stand-alone script)
        or in an external pipeline script.

    base_dir: str | Path
        Base directory containing structure directories with VASP runs inside.

    prev_summary: str | Path, optional
        Path to a JSON summary file produced by a previous screening step. If given, the file
        will be checked to filter out structures that are already rejected.

    summary_key: str, optional
        The dict key associated to the bool used to verify eligibility in previous summary file.
        If --prev-summary is given, it must be given too.

    summary: str | Path, optional
        Name of output JSON summary file indicating calculation results for this step. Also usable
        by further steps to filter out structures rejected in this step. The file will be
        automatically written in the 'base_dir' directory.

    workers: int, optional
        Number of parallel processes to spawn for parallelized steps. If not given, default value
        is the 'max_workers' default value from `tqdm.contrib.concurrent.process_map()` function.
        Pass 0 to disable the use of `process_map()` and execute sequentially.
    """
    args = parse_input_args(_get_cmd_line_args, _process_input_args, standalone, **kwargs)
    all_runs = VaspParser(args["base_dir"], workers=args.get("workers")).all_struct_dirs
    converged_runs_set = {
        run for run in all_runs if VaspParser.get_safe_vasprun(
            run, converged=True, parse_dos=False, parse_eigen=False, parse_potcar_file=False
        ) is not None
    }
    results = {}
    for run in VisualIterator(all_runs, desc="Writing relax summary", percent=True):
        results[run.name] = {
            "path": str(run.resolve()),
            "name": run.name,
            "converged": run in converged_runs_set
        }
    JsonWriter(args["summarypath"], data=results, indent=4).write_as_dict()


if __name__ == "__main__":
    main()

