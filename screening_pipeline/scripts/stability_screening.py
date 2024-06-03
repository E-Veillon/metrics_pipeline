#!/usr/bin/python
"""
A script using previous VASP relaxations to compute the relative stability of given structures.

Algorithm:

    1/  Group structures by dimension and composition. - OK
        
    2/  Get reference structures from a dataset and add them to the pool
        with a specific attribute. - Not done yet

    3/  Build one convex hull per structure group 
        (initialize ref elements from group composition) - OK
        
    4/  Inside each convex hull, compute ΔH for each generated entry. - OK

    5/  Reject structures that have a ΔH above given threshold value. - OK

Possible alternatives:

    Comparing to known references:
        1 - For each group, search in the dataset all corresponding structures.
        2 - Build the phase diagram corresponding to dataset structures.
        3 - Compute ΔH for all generated structures of the group.
        4 - Reject structures too much above reference convex hull.

    Opti:   Store some common phase spaces in a local util file.

    Opti:   Match and ignore structures equivalents to references with StructureMatcher.
"""


import os
import json
import argparse
from datetime import datetime


from screening_pipeline.utils import (
    get_max_dim,
    init_entries_from_dict,
    filter_database_entries,
    get_elements_from_entries,
    get_lacking_elts_entries,
    group_by_composition,
    batch_compute_e_above_hull,
    batch_extract_vasp_data,
    load_phase_diagram_entries,
)


########################################
# LOCAL FUNCTIONS


def assert_args(args: argparse.Namespace) -> None:
    """Asserting input arguments validity."""
    import os

    assert os.path.isdir(args.run_dir), f"{args.run_dir}: No such directory found."

    assert os.path.isfile(args.reference), f"{args.reference}: No such file found."

    if args.prev_summary is not None:
        assert os.path.isfile(args.prev_summary), (
            f"{args.prev_summary}: No such file found."
        )

    assert args.limit >= 1e-8, (
        "Instability energy elimination criterion must be strictly positive.\n"
        "Moreover, any value below 1.10-8 eV/atom is too low and not supported."
    )

    assert args.workers >= 1, "The number of workers must be stricly positive."

    assert args.summary.endswith(".json"), (
        f"'{args.summary}' is not a valid JSON file, verify the .json extension."
    )


########################################
# MAIN FUNCTION


def main():
    """Main function."""
    start = datetime.now()

    # ARGUMENTS PARSING BLOCK

    prog_name = "stability_screening.py"
    prog_description = "A script using previous VASP relaxations to compute the relative \
                          stability of given structures."
    prog_missing_steps = """
        Missing steps to complete the script:
            - Determine the critical formation reaction of the structure
            - (Relax reference structures for calculations consistency)
            - Compare structure energy with the sum of reference structure energies
        """
    helper_format = argparse.RawTextHelpFormatter

    parser = argparse.ArgumentParser(
        prog=prog_name,
        description=prog_description,
        epilog=prog_missing_steps,
        formatter_class=helper_format,
    )

    parser.add_argument(
        "run_dir",
        type=str,
        help="Base directory containing structure directories.",
    )
    parser.add_argument(
        "-r",
        "--reference",
        help=(
            "dataset of already known structure to construct a reference convex hull "
            "and compare generated structure against it (json format).\n"
            "If not given, the default reference hull will be defined only with elemental "
            "entries of energy 0.0 eV/atom."
        ),
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
        "-s",
        "--summary",
        default="summary.json",
        help=(
            "Output file indicating calculation results for this step (json format).\n"
            "Also usable by further steps to filter out structures rejected in this step."
            "This argument only affects the name of the file, it is automatically written "
            "at the location given by the 'run_dir' argument."
        ),
    )
    parser.add_argument(
        "-l",
        "--limit",
        type=float,
        default=0.1,
        help=(
            "Maximum value of ΔH (in eV/atom) above which structures "
            "are considered too unstable and rejected.\n"
            "Defaults to 0.1 eV/atom, as it is commonly assumed to be sufficient."
        ),
        metavar="float",
    )
    parser.add_argument(
        "-w",
        "--workers",
        type=int,
        default=1,
        help="Number of parallel processes to spawn for parallelized steps.",
        metavar="int",
    )

    args: argparse.Namespace = parser.parse_args()

    assert_args(args)

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
        entries_dict=structs_data, attribute="generated", workers=args.workers
    )
    max_dim_generated = get_max_dim(generated_entries)
    used_elts = get_elements_from_entries(generated_entries)

    # Extract reference dataset
    if args.reference is not None:
        ref_data = load_phase_diagram_entries(args.reference)
        ref_entries = init_entries_from_dict(
            entries_dict=ref_data, attribute="ref_structs", workers=args.workers
        )
        ref_entries = filter_database_entries(
            entries=ref_entries,
            max_dim=max_dim_generated,
            ref_elts=used_elts
        )
    else:
        ref_entries = []

    # Generate and add default elemental references if they are not in the reference dataset
    auto_elts_entries = get_lacking_elts_entries(ref_entries, ref_elts=used_elts)
    ref_entries += auto_elts_entries

    # Group generated entries by composition
    grouped_entries = group_by_composition(comps=generated_entries)

    # Compute energy above hulls in each group
    screening_results = batch_compute_e_above_hull(
        entries_to_compute=grouped_entries, ref_entries=ref_entries,
        stable_limit=delta_H_limit, workers=args.workers
    )

    for dct in screening_results:
        path = os.path.join(args.run_dir, dct["name"])
        dct.update({"path": path})

    def sort_by_path(dct: dict) -> str:
        return dct.get("path")

    screening_results = sorted(screening_results, key=sort_by_path)

    with open(summary, mode="wt", encoding="utf-8") as fp:
        json.dump(screening_results, fp, indent=4)

    stop = datetime.now()
    print(f"elapsed time: {stop-start}")


if __name__ == "__main__":
    main()
