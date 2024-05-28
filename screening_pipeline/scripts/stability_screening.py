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

########################################
# SYSTEM I/O MODULES

from datetime import datetime
from argparse import ArgumentParser, Namespace, RawTextHelpFormatter
from typing import Dict
from pathlib import Path
import os
import json

########################################
# OPTIMIZATION MODULES

########################################
# PYTON MATERIALS GENOMICS PACKAGE

########################################
# LOCAL MODULES

from screening_pipeline.utils import (
    batch_extract_vasp_data,
    batch_calculate_instability_energies,
    load_phase_diagram_entries,
)


########################################
# LOCAL FUNCTIONS


def assert_args(args: Namespace) -> None:
    """Asserting input arguments validity."""

    assert os.path.isdir(args.run_dir), f"{args.run_dir}: No such directory found."

    if args.prev_summary is not None:
        assert os.path.isfile(args.prev_summary), (
            f"{args.prev_summary}: No such file found."
        )

    assert args.limit >= 1e-8, (
        "Instability energy elimination criterion must be strictly positive.\n"
        "Moreover, any value below 1.10-8 eV/atom is too low and not supported."
    )

    assert args.workers >= 1, "The number of workers cannot be negative or zero."


########################################
# MAIN FUNCTION


def main():
    """Main function."""
    start = datetime.now()

    # ARGUMENTS PARSING BLOCK

    prog_name = "stability_screening"
    prog_description = "A script using previous VASP relaxations to compute the relative \
                          stability of given structures."
    prog_missing_steps = """
        Missing steps to complete the script:
            - Determine the critical formation reaction of the structure
            - (Relax reference structures for calculations consistency)
            - Compare structure energy with the sum of reference structure energies
        """
    helper_format = RawTextHelpFormatter

    parser = ArgumentParser(
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
            "If the user wants to calculate the energy above the hull from an existing dataset "
            "used as a reference (json format)."
        ),
    )
    parser.add_argument(
        "-a",
        "--accept",
        type=str,
        default=None,
        help=(
            "Defines an optional file whose presence in a structure directory means it passed "
            "this screening step successfully and can be kept for further calculations.\n"
            "This file contains a small text giving result of the step "
            "for the corresponding structure."
        ),
        metavar="filename",
    )
    parser.add_argument(
        "-i",
        "--ignore",
        type=str,
        default=None,
        help=(
            "Defines a file whose presence in a structure directory means it did not pass "
            "previous screening steps and should not be used in this calculation.\n" 
            "This file will also be written in structure directories that did not pass this step.\n"
            "WARNING:\n" 
            "If this feature is used, care must be taken that the same file name is used throughout"
            "every screening steps to make sure rejected structures don't go further."
        ),
        metavar="ignore_file.txt",
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
        "-l",
        "--limit",
        type=float,
        default=0.1,
        help=(
            "Maximum value of ΔH (in eV/atom) above which structures "
            "are considered too unstable and rejected.\n"
            "Defaults to 36 meV/atom, as used in the following paper, and seems fairly strict:\n" 
            "Y. Wu, P. Lazic, G. Hautier, K. Persson, and G. Ceder,\n"
            "First principles high throughput screening of oxynitrides "
            "for water-splitting photocatalysts,\n" 
            "Energy & Environmental Science 6, no. 1 (2012) 157."
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
    parser.add_argument(
        "-s",
        "--summary",
        default="summary.json",
        help=(
            "Output file indicating calculation results for this step (json format).\n"
            "Also usable by further steps to filter out structures that were rejected in this step."
        ),
    )
    args: Namespace = parser.parse_args()

    assert_args(args)

    run_dir = Path(args.run_dir)
    accept_file = args.accept or None
    ignore_file = args.ignore or None
    prev_summary = args.prev_summary or None
    delta_H_limit = round(args.limit, 8)
    workers = args.workers
    summary = os.path.join(str(run_dir), args.summary)

    # MAIN BLOCK

    structs_data = batch_extract_vasp_data(
        method="convex_hull",
        base_dir=run_dir,
        ignore_file=ignore_file,
        path_to_summary=prev_summary,
        workers=workers,
    )

    if args.reference is not None:
        structs_reference = load_phase_diagram_entries(args.reference)
    else:
        structs_reference = None

    structs_data = batch_calculate_instability_energies(
        structs_data=structs_data, structs_ref=structs_reference, workers=workers
    )

    stable_structs = filter(
            lambda tup: tup[1]["delta_H"] <= delta_H_limit, list(structs_data.items())
    )

    unstable_structs = filter(
            lambda tup: tup[1]["delta_H"] > delta_H_limit, list(structs_data.items())
    )

    screening_results = []

    for struct in stable_structs:
        name = struct[0]
        data = struct[1]

        if accept_file is not None:
            accept_msg = (
                "STABILITY TEST PASSED\n"
                f"Instability energy for this structure is estimated at {data['delta_H']} eV/atom,\n"
                f"which is below or equal to the fixed instability limit of {delta_H_limit} eV/atom.\n"
                "Therefore, it is considered suitable for wanted application, "
                "and should be considered for further screening steps.\n"
            )
            accept_file_path = os.path.join(str(run_dir), name, accept_file)
            with open(accept_file_path, mode="wt", encoding="utf-8") as fp:
                fp.write(accept_msg)

        screening_results.append(
            {
                "path": os.path.join(str(run_dir), name),
                "energy_above_hull": data["delta_H"],
                "stable": True,
            }
        )

    for struct in unstable_structs:
        name = struct[0]
        data = struct[1]

        if ignore_file is not None:
            reject_msg = (
                "STABILITY REJECTION\n"
                f"Instability energy for this structure is estimated at {data['delta_H']} eV/atom,\n"
                f"which is above the fixed instability limit of {delta_H_limit} eV/atom.\n"
                "Therefore, it is considered not suitable for wanted application, "
                "and should not be considered in further screening steps.\n"
            )
            reject_file_path = os.path.join(str(run_dir), name, ignore_file)
            with open(reject_file_path, mode="wt", encoding="utf-8") as fp:
                fp.write(reject_msg)

        screening_results.append(
            {
                "path": os.path.join(str(run_dir), name),
                "energy_above_hull": data["delta_H"],
                "stable": False,
            }
        )

    def sort_by_path(dct: Dict) -> str:
        return dct.get("path")

    screening_results = sorted(screening_results, key=sort_by_path)

    with open(summary, mode="w", encoding="utf-8") as fp:
        json.dump(screening_results, fp, indent=4)

    stop = datetime.now()
    print(f"elapsed time: {stop-start}")


if __name__ == "__main__":
    main()
