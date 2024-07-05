#!/usr/bin/python
"""
A script using previous VASP relaxations to compute the phase stability
of given structures by comparing their energy with a convex hull of
energies of refererence structures.

Algorithm:

    1 - Group structures by dimension and composition.

    2 - For each group, search in the dataset all corresponding structures.

    3 - Build the phase diagram corresponding to dataset structures.

    4 - Compute ΔH for all generated structures of the group.

    5 - Reject structures too much above reference convex hull.
"""


import os
import json
import argparse
from datetime import datetime

# LOCAL IMPORTS
from screening_pipeline.utils import (
    check_file_format, check_file_or_dir, check_num_value,
    get_max_dim, init_entries_from_dict, filter_database_entries,
    get_elements_from_entries, get_lacking_elts_entries,
    group_by_composition, batch_compute_e_above_hull,
    batch_extract_vasp_data, load_phase_diagram_entries,
)


########################################


def _assert_args(args: argparse.Namespace) -> None:
    """Asserting input arguments validity."""
    check_file_or_dir(args.run_dir, "dir")

    if args.reference is not None:
        check_file_or_dir(args.reference, "file", format="json")

    if args.prev_summary is not None:
        check_file_or_dir(args.prev_summary, "file", format="json")

    check_file_format(args.summary, format="json")
    check_num_value(args.limit, "--limit", ">=", float(1e-6))
    check_num_value(args.workers, "--workers", ">", 0)


########################################


def main():
    """Main function."""
    start = datetime.now()

    # ARGUMENTS PARSING BLOCK

    prog_name = "stability_screening.py"
    prog_description = (
        "A script using previous VASP relaxations to compute the phase "
        "stability of given structures by comparing their energy with "
        "a convex hull of energies of refererence structures."
    )

    parser = argparse.ArgumentParser(
        prog=prog_name,
        description=prog_description
    )

    parser.add_argument(
        "run_dir",
        type=str,
        help="Base directory containing structure directories.",
    )
    parser.add_argument(
        "-r", "--reference",
        help=(
            "dataset of already known structure to construct a reference convex hull "
            "and compare generated structure against it (json format).\n"
            "If not given, the default reference hull will be defined only with "
            "elemental entries of energy 0.0 eV/atom."
        ),
    )
    parser.add_argument(
        "-R", "--read-previous-summary",
        type=str,
        default=None,
        help=(
            "Path to a JSON summary file produced by a previous screening step.\n"
            "If given, the file will be checked to filter out structures that are "
            "already rejected."
        ),
        metavar="/path/to/summary.json",
        dest="prev_summary"
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
        metavar="float",
    )
    parser.add_argument(
        "-w", "--workers",
        type=int,
        default=1,
        help="Number of parallel processes to spawn for parallelized steps.",
        metavar="int",
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

    args: argparse.Namespace = parser.parse_args()

    _assert_args(args)

    prev_summary = args.prev_summary or None
    delta_H_limit = round(args.limit, 8)
    summary = os.path.join(args.run_dir, args.summary)


    # MAIN BLOCK

    # Extract generated data
    structs_data = batch_extract_vasp_data(
        method="convex_hull",
        base_dir=args.run_dir,
        path_to_summary=prev_summary,
        workers=args.workers,
    )
    generated_entries = init_entries_from_dict(
        entries_dict=structs_data, attribute="generated"
    )
    max_dim_generated = get_max_dim(generated_entries)
    used_elts = get_elements_from_entries(generated_entries)

    # Extract reference dataset
    if args.reference is not None:
        ref_data = load_phase_diagram_entries(args.reference)
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

    if args.pause_after_init:
        input("Tap Enter to continue:")

    # Generate and add default elemental references if they are not in the reference dataset
    auto_elts_entries = get_lacking_elts_entries(ref_entries, ref_elts=used_elts)
    ref_entries += auto_elts_entries

    # Group generated entries by composition
    grouped_entries = group_by_composition(comps=generated_entries)

    # Compute energy above hulls in each group
    screening_results = batch_compute_e_above_hull(
        entries=grouped_entries,
        ref_entries=ref_entries,
        stable_limit=delta_H_limit,
        workers=args.workers,
        verbose=args.verbose
    )

    for dct in screening_results:
        path = os.path.join(args.run_dir, dct["name"])
        dct.update({"path": os.path.abspath(path)})

    def sort_by_path(dct: dict) -> str:
        return dct.get("path")

    screening_results = sorted(screening_results, key=sort_by_path)

    with open(summary, mode="wt", encoding="utf-8") as fp:
        json.dump(screening_results, fp, indent=4)

    stop = datetime.now()
    print(f"elapsed time: {stop-start}")


if __name__ == "__main__":
    main()
