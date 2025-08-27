#!/usr/bin/python
"""
A script using previous VASP static computations to compute the phase stability
of given structures by comparing their energy with a convex hull of
energies of reference structures.

Algorithm:

    1 - Group structures by dimension and composition.

    2 - For each group, search in the dataset all corresponding structures.

    3 - Build the phase diagram corresponding to dataset structures.

    4 - Compute ΔH for all generated structures of the group.

    5 - Reject structures too much above reference convex hull.
"""


import os
import json
import argparse as argp
import typing as tp
from datetime import datetime

# LOCAL IMPORTS
from . import _parse_input_args
from .utils import (
    check_type, check_file_format, check_file_or_dir, check_num_value,
    get_max_dim, init_entries_from_dict, filter_database_entries,
    get_elements_from_entries, get_lacking_elts_entries,
    group_by_composition, batch_compute_e_above_hull,
    batch_extract_vasp_data, load_phase_diagram_entries,
)


########################################
# ARGUMENTS HANDLING

def _get_command_line_args() -> argp.Namespace:
    """Command Line Interface (CLI)."""
    parser = argp.ArgumentParser(prog=os.path.basename(__file__), description=__doc__)
    parser.add_argument(
        "run_dir",
        help="Base directory containing structure directories with VASP runs inside.",
    )
    parser.add_argument(
        "-r", "--reference",
        help=(
            "Dataset of already known structure to construct a reference convex hull "
            "and compare generated structure against it (json format).\n"
            "If not given, the default reference hull will be defined only with "
            "elemental entries of energy 0.0 eV/atom."
        ),
    )
    parser.add_argument(
        "-R", "--read-previous-summary",
        help=(
            "Path to a JSON summary file produced by a previous screening step.\n"
            "If given, the file will be checked to filter out structures that are "
            "already rejected."
        ),
        metavar="<path>",
        dest="prev_summary"
    )
    parser.add_argument(
        "-k", "--key-to-check",
        help=(
            "The dict key associated to the bool used to verify eligibility "
            "in the previous summary file.\n"
            "If --read-previous-summary is given, it must be given too."
        ),
        metavar="<str>"
    )
    parser.add_argument(
        "-s", "--summary",
        default="summary.json",
        help=(
            "Output file indicating calculation results for this step (json format).\n"
            "Also usable by further steps to filter out structures rejected in this step."
            "This argument only affects the name of the file, it is automatically written "
            "at the location given by the 'run_dir' argument."
        ),
        metavar="<str>"
    )
    parser.add_argument(
        "-l", "--limit",
        type=float,
        default=0.1,
        help=(
            "Maximum value of ΔH (in eV/atom) above which structures "
            "are considered too unstable and rejected.\n"
            "Defaults to 0.1 eV/atom, as it is commonly assumed to be sufficient.\n"
        ),
        metavar="<float>",
    )
    parser.add_argument(
        "-w", "--workers",
        type=int,
        help=(
            "Number of parallel processes to spawn for parallelized steps. "
            "If not given, If not given, default value is the 'max_workers' "
            "default value from tqdm.contrib.concurrent.process_map function."
        ),
        metavar="<int>",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Whether to print each reference entry used when building a phase diagram."
    )
    parser.add_argument(
        "--pause-after-init",
        action="store_true",
        help="Pauses the program after finishing data preparations. Press Enter to unpause."
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
    args_dict.setdefault("summary", "summary.json")
    args_dict.setdefault("limit", 0.1)
    args_dict.setdefault("verbose", False)
    args_dict.setdefault("pause_after_init", False)

    # Assert set arguments conformity
    check_file_or_dir(args_dict.get("run_dir"), "dir")
    if args_dict.get("reference") is not None:
        check_file_or_dir(args_dict.get("reference"), "file", allowed_formats="json")

    if args_dict.get("prev_summary") is not None:
        check_file_or_dir(args_dict.get("prev_summary"), "file", allowed_formats="json")

        if args_dict.get("key_to_check") is None:
            raise ValueError(
                f"'prev_summary' argument was provided ({args_dict.get('prev_summary')}), "
                "therefore 'key-to-check' argument has to be given as well."
            )

    if args_dict.get("key_to_check") is not None:
        check_type(args_dict.get("key_to_check"), "key_to_check", (str,))

    check_file_format(args_dict.get("summary"), allowed_formats="json")
    check_type(args_dict.get("limit"), "limit", (float,))
    check_num_value(args_dict.get("limit"), "limit", ">=", 0.0)

    if args_dict.get("workers") is not None:
        check_type(args_dict.get("workers"), "workers", (int,))
        check_num_value(args_dict.get("workers"), "workers", ">", 0)

    # Additional arguments processing
    args_dict["limit"] = round(args_dict["limit"], 8)
    args_dict["summarypath"] = os.path.join(args_dict["run_dir"], args_dict["summary"])

    return args_dict


########################################


def main(standalone: bool = True, **kwargs):
    """
    A script using previous VASP static computations to compute the phase stability
    of given structures by comparing their energy with a convex hull of
    energies of reference structures.

    Algorithm:

        1 - Group structures by dimension and composition.
        2 - For each group, search in the dataset all corresponding structures.
        3 - Build the phase diagram corresponding to dataset structures.
        4 - Compute ΔH for all generated structures of the group.
        5 - Reject structures too much above reference convex hull.

    Args:
        standalone (bool):          Whether parsed script is used directly through
                                    command-line (stand-alone script) or in an external
                                    pipeline script.

        run_dir (str|Path):         Base directory containing structure directories with VASP runs
                                    inside.

        reference (str|Path):       Dataset of already known structure to construct a reference
                                    convex hull and compare generated structure against it
                                    (json format). If not given, the default reference hull will
                                    be defined only with elemental entries of energy 0.0 eV/atom.

        prev_summary (str|Path):    Path to a JSON summary file produced by a previous screening
                                    step. If given, the file will be checked to filter out
                                    structures that are already rejected.

        key_to_check (str):         The dict key associated to the bool used to verify eligibility
                                    in the previous summary file. If prev_summary is given, it must
                                    be given too.

        summary (str|Path):         Output file indicating calculation results for this step
                                    (json format). Also usable by further steps to filter out
                                    structures rejected in this step. This argument only affects
                                    the name of the file, it is automatically written at the
                                    location given by the 'run_dir' argument.

        limit (float):              Maximum value of ΔH (in eV/atom) above which structures are
                                    considered too unstable and rejected. Defaults to 0.1 eV/atom,
                                    as it is commonly assumed to be sufficient.

        workers (int):              Number of parallel processes to spawn for parallelized steps.
                                    If not given, default value is the 'max_workers' default value
                                    from tqdm.contrib.concurrent.process_map function.

        verbose (bool):             Whether to print each reference entry used when building a
                                    phase diagram.

        pause_after_init (bool):    Pauses the program after finishing data preparations.
                                    Press Enter to unpause.
"""
    start = datetime.now()
    args = _parse_input_args(_get_command_line_args, _process_input_args, standalone, **kwargs)

    # Extract generated data
    structs_data = batch_extract_vasp_data(
        method="convex_hull",
        base_dir=args["run_dir"],
        path_to_summary=args.get("prev_summary"),
        key_to_check=args.get("key_to_check"),
        workers=args.get("workers"),
    )
    generated_entries = init_entries_from_dict(
        entries_dict=structs_data, attribute="generated"
    )
    # Eliminate high dimension structures (> 10) from computations to avoid softlock
    gen_entries = list(
        filter(
            lambda entry: len(entry.composition) < 11,
            generated_entries
        )
    )

    max_dim_generated = get_max_dim(gen_entries)
    used_elts = get_elements_from_entries(gen_entries)

    # Extract reference dataset
    if args.get("reference") is not None:
        ref_data = load_phase_diagram_entries(args["reference"])
        ref_entries = init_entries_from_dict(
            entries_dict=ref_data, attribute="ref_structs"
        )
        ref_entries = filter_database_entries(
            entries=ref_entries,
            max_dim=max_dim_generated,
            ref_elts=used_elts
        )
    else:
        ref_entries = []

    if args.get("pause_after_init"):
        input("Tap Enter to continue:")

    # Generate and add default elemental references if they are not in the reference dataset
    auto_elts_entries = get_lacking_elts_entries(ref_entries, ref_elts=used_elts)
    ref_entries += auto_elts_entries

    # Group generated entries by composition
    grouped_entries = group_by_composition(comps=gen_entries)

    # Compute energy above hulls in each group
    screening_results = batch_compute_e_above_hull(
        entries=grouped_entries,
        ref_entries=ref_entries,
        stable_limit=args["limit"],
        workers=args.get("workers"),
        verbose=args.get("verbose", False)
    )

    for dct in screening_results:
        path = os.path.join(args["run_dir"], dct["name"]) # type: ignore
        dct.update({"path": os.path.abspath(path)})

    # Too high dimension structures are added as not stable in results
    high_dim_entries = list(
        filter(
            lambda entry: len(entry.composition) >= 11,
            generated_entries
        )
    )
    for entry in high_dim_entries:
        msg = f"Uncomputable due to its too high dimension ({len(entry.composition)} > 10)."
        screening_results.append(
            {
                "name": entry.name,
                "path": os.path.abspath(os.path.join(args["run_dir"], entry.name)),
                "e_above_hull": None, # type: ignore
                "stable": False,
                "comment": msg
            }
        )

    def sort_by_path(dct: dict) -> str:
        return dct["path"]

    screening_results = sorted(screening_results, key=sort_by_path)

    with open(args["summarypath"], mode="wt", encoding="utf-8") as fp:
        json.dump(screening_results, fp, indent=4)

    stop = datetime.now()
    print(f"elapsed time: {stop-start}")


if __name__ == "__main__":
    main()
