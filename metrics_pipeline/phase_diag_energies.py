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
import argparse as argp
import typing as tp
from datetime import datetime

from pymatgen.core import Element

from . import _parse_input_args
from utils.utils import check_type, check_num_value
from utils.io import (
    VaspParser, VaspExtractor, PDDataset, JsonWriter,
    check_file_or_dir, check_file_format
)
from utils.metrics import Stability
from utils.computations.local import get_lacking_elts_entries


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
        "--compact", action="store_true",
        help=(
            "Only used if 'process-dataset is given. "
            "If passed, tells the parser that given JSON is organized by lists of "
            "attributes, e.g. {'id': [id1, id2, ...], 'composition': [comp1, comp2, ...], ...} "
            "instead of being organized by individual objects (default), "
            "e.g. {'data1': {'id': id1, ...}, 'data2': {'id': id2, ...}, ...}."
        )
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
    args_dict.setdefault("compact", False)
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
    check_type(args_dict.get("compact"), "compact", (bool,))

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

    Algorithm
    ---------
    1 - Group structures by dimension and composition.
    2 - For each group, search in the dataset all corresponding structures.
    3 - Build the phase diagram corresponding to dataset structures.
    4 - Compute ΔH for all generated structures of the group.
    5 - Reject structures too much above reference convex hull.

    Notes
    -----
    This script assumes that loaded reference dataset is already properly formatted for
    phase diagram computations. Make sure to pass a compatible dataset file.
    JSON dataset files can be made compatible by using the 'dataset_loader.py' script.

    Parameters
    ----------
    standalone: bool
        Whether parsed script is used directly through command-line (stand-alone script)
        or in an external pipeline script.

    run_dir: str | Path
        Base directory containing structure directories with VASP runs inside.

    reference: str | Path
        Dataset of already known structure to construct a reference convex hull and compare
        generated structure against it (json format). If not given, the default reference hull
        will be defined only with elemental entries of energy 0.0 eV/atom.

    prev_summary: str | Path
        Path to a JSON summary file produced by a previous screening step. If given, the file
        will be checked to filter out structures that are already rejected.

    key_to_check: str
        The dict key associated to the bool used to verify eligibility in the previous summary
        file. If prev_summary is given, it must be given too.

    summary: str | Path
        Output file indicating calculation results for this step (json format). Also usable by
        further steps to filter out structures rejected in this step. This argument only affects
        the name of the file, it is automatically written at the location given by the 'run_dir'
        argument.

    limit: float
        Maximum value of ΔH (in eV/atom) above which structures are considered too unstable and
        rejected. Defaults to 0.1 eV/atom, as it is commonly assumed to be sufficient.

    compact: bool
        Only used if 'process-dataset is given. If passed, tells the parser that given JSON is
        organized by lists of attributes, e.g.
        {'id': [id1, id2, ...], 'composition': [comp1, comp2, ...], ...} instead of being organized
        by individual objects (default), e.g.
        {'data1': {'id': id1, ...}, 'data2': {'id': id2, ...}, ...}.
        Defaults to False.

    workers: int
        Number of parallel processes to spawn for parallelized steps. If not given, default value
        is the 'max_workers' default value from `tqdm.contrib.concurrent.process_map()` function.
        Pass 0 to disable the use of `process_map()` and execute sequentially.

    verbose: bool
        Whether to print each reference entry used when building a phase diagram.
        Defaults to False.

    pause_after_init: bool
        Pauses the program after finishing data preparations. Press Enter to unpause.
        Defaults to False.
"""
    start = datetime.now()
    args = _parse_input_args(_get_command_line_args, _process_input_args, standalone, **kwargs)

    # Extract generated data
    gen_data = VaspExtractor(
        vasp_parser=VaspParser(args["run_dir"]),
        method="convex_hull",
        summary_name=args.get("prev_summary"),
        summary_key=args.get("key_to_check"),
        workers=args.get("workers")
    ).get_data()

    generated_dataset = PDDataset(
        gen_data,
        composition_key="composition",
        energy_key="final_energy",
        attribute="generated",
        compact=False
    )

    # Eliminate high dimension structures (> 10) from computations to avoid softlock
    true_max_dim = generated_dataset.max_dim
    if true_max_dim > 10:
        max_dim_generated = generated_dataset.max_computable_dim
        gen_entries = generated_dataset.get_computable_entries()
    else: # Save filtering computation
        max_dim_generated = true_max_dim
        gen_entries = generated_dataset.get_all_entries()

    used_elts = generated_dataset.computable_elements

    # Extract reference dataset
    if args.get("reference") is not None:
        ref_dataset = PDDataset.from_file(
            filepath=args["reference"],
            composition_key="composition",
            energy_key="final_energy",
            attribute="reference",
            compact=args["compact"]
        )
        ref_entries = ref_dataset.get_filtered_entries(
            elts=used_elts,
            dims=set(list(range(max_dim_generated)))
        )
    else:
        ref_entries = {}

    if args.get("pause_after_init"):
        input("Tap Enter to continue:")

    # Generate and add default elemental references if they are not in the reference dataset
    used_elts_pmg = set(Element(elt) for elt in used_elts)
    auto_elts_entries = get_lacking_elts_entries(
        list(ref_entries.values()), ref_elts=used_elts_pmg
    )
    ref_entries.update({f"auto_{entry.elements[0].symbol}": entry for entry in auto_elts_entries})

    # Compute Stability
    stability = Stability.from_entries(
        entries=list(gen_entries.values()),
        ref_entries=list(ref_entries.values()),
        stable_tol=args["stable_tol"],
        workers=args.get("workers")
    )
    results = {}
    stable_ids = set(id(entry) for entry in stability.stable_entries)
    for name, entry in gen_entries.items():
        dct = {
            "name": name,
            "path": os.path.join(args["run_dir"], name),
            "e_above_hull": entry.attribute[Stability._delta_e_attr], # type: ignore
            "is_stable": id(entry) in stable_ids
        }
        results[name] = dct

    # Too high dimension structures are added as not stable in results
    high_dim_entries = generated_dataset.get_uncomputable_entries()
    for name, entry in high_dim_entries.items():
        msg = f"Uncomputable due to its too high dimension ({len(entry.composition)} > 10)."
        results[name] = {
                "name": entry.name,
                "path": os.path.abspath(os.path.join(args["run_dir"], entry.name)),
                "e_above_hull": None, # type: ignore
                "stable": False,
                "comment": msg
            }

    sorted_results = dict(sorted(results.items(), key=lambda tup: tup[1]["path"]))

    JsonWriter(args["summarypath"], data=sorted_results, indent=4).write_as_dict()

    stop = datetime.now()
    print(f"elapsed time: {stop-start}")


if __name__ == "__main__":
    main()
