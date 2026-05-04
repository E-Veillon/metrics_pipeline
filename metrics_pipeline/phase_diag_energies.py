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
import argparse as ap
import typing as tp
from datetime import datetime

from pymatgen.core import Element

from src.utils import parse_input_args, check_type, check_num_value
from src.io import (
    VaspParser, VaspExtractor, ExtractMethod, PDDataset, GenMatPDDataset, JsonWriter, GenMatFile,
    check_file_or_dir, check_file_format
)
from src.metrics import Stability
from src.computations.local import get_lacking_elts_entries
from src.computations.models import vectors_from_alignn


def _get_command_line_args() -> ap.Namespace:
    """Command Line Interface (CLI)."""
    parser = ap.ArgumentParser(prog=os.path.basename(__file__), description=__doc__)
    parser.add_argument(
        "base_dir",
        help="Base directory containing VASP run directories."
    )
    parser.add_argument(
        "-r", "--reference", metavar="<path>",
        help=(
            "JSON dataset file of already known structure to construct a reference "
            "convex hull and compare generated structure against it. If not given, "
            "the default reference hull will be defined only with elemental entries "
            "of energy 0.0 eV/atom."
        ),
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
        "-s", "--summary", default="stability_summary.json",
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
    parser.add_argument(
        "-l", "--limit", type=float, default=0.1,
        help=(
            "Maximum value of ΔH (in eV/atom) above which structures are considered too "
            "unstable and rejected. Defaults to 0.1 eV/atom, as it is commonly assumed "
            "to be sufficient."
        )
    )
    parser.add_argument(
        "--alignn", action=ap.BooleanOptionalAction, default=False,
        help=(
            "Whether to predict generated data energies with ALIGNN instead of extracting them "
            "from VASP runs. Defaults to %(default)s. If enabled, 'base_dir' must be a path "
            "to a preprocessed JSON file instead of a VASP base directory."
        )
    )
    parser.add_argument(
        "--compact", action="store_true",
        help=(
            "Flag to pass if reference dataset is organized by data type instead of by structure. "
            "By default, dataset is assumed to be organized by structure, i.e. "
            "{'struct_name_1': struct_dict_1, 'struct_name_2': struct_dict_2, ...}. "
            "If the flag is passed, dataset is instead assumed to be organized by data key, i.e. "
            "{'entry_id': [struct_name_1, struct_name_2, ...], "
            "'composition': [struct_comp_1, struct_comp_2, ...], ...}."
        )
    )
    parser.add_argument(
        "--match-all", action=ap.BooleanOptionalAction, default=True,
        help=(
            "Whether to raise an error if some run indices do not match any structure directory "
            "in base directory. Defaults to %(default)s. Deactivate it if some structures "
            "could not be written as VASP input when using 'vasp_rundir_writer.py' script."
        )
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Whether to print each reference entry used when building a phase diagram."
    )
    parser.add_argument(
        "--pause-after-init", action="store_true",
        help="Pauses the program after finishing data preparations. Press Enter to unpause."
    )
    args: ap.Namespace = parser.parse_args()
    return args


def _process_input_args(args_dict: dict[str, tp.Any]) -> dict[str, tp.Any]:
    """Handle input arguments assertions and processing."""
    if args_dict is None:
        raise ValueError(f"No arguments found at '{os.path.basename(__file__)}' script call.")
    
    check_type(args_dict, "args_dict", (dict,))

    # Set default values for unset optional arguments
    args_dict = {k: v for k, v in args_dict.items() if v is not None}
    args_dict.setdefault("summary", "stability_summary.json")
    args_dict.setdefault("limit", 0.1)
    args_dict.setdefault("alignn", False)
    args_dict.setdefault("compact", False)
    args_dict.setdefault("match_all", False)
    args_dict.setdefault("verbose", False)
    args_dict.setdefault("pause_after_init", False)

    # Assert set arguments conformity
    check_file_or_dir(args_dict.get("base_dir"), "dir")
    if args_dict.get("reference") is not None:
        check_file_or_dir(args_dict.get("reference"), "file", allowed_formats="json")

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
    check_type(args_dict["limit"], "limit", (float,))
    check_num_value(args_dict["limit"], "limit", ">=", 0.0)
    check_type(args_dict.get("compact"), "compact", (bool,))

    if args_dict.get("workers") is not None:
        check_type(args_dict["workers"], "workers", (int,))
        check_num_value(args_dict["workers"], "workers", ">=", 0)

    # Additional arguments processing
    args_dict["limit"] = round(args_dict["limit"], 8)
    args_dict["summarypath"] = os.path.join(args_dict["base_dir"], args_dict["summary"])

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

    base_dir: str | Path
        Base directory containing structure directories with VASP runs inside.

    reference: str | Path, optional
        JSON dataset file of already known structure to construct a reference  convex hull
        and compare generated structure against it. If not given, the default reference hull
        will be defined only with elemental entries of energy 0.0 eV/atom.

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

    limit: float, optional
        Maximum value of ΔH (in eV/atom) above which structures are considered too unstable and
        rejected. Defaults to 0.1 eV/atom, as it is commonly assumed to be sufficient.

    alignn: bool
        Whether to predict generated data energies with ALIGNN instead of extracting them
        from VASP runs. Defaults to False. If enabled, 'base_dir' must be a path
        to a preprocessed JSON file instead of a VASP base directory.

    compact: bool
        Whether reference dataset is organized by data type instead of by structure.
        By default, dataset is assumed to be organized by structure, i.e.
        {'struct_name_1': struct_dict_1, 'struct_name_2': struct_dict_2, ...}.
        If the flag is passed, dataset is instead assumed to be organized by data key, i.e.
        {
            'entry_id': [struct_name_1, struct_name_2, ...],
            'composition': [struct_comp_1, struct_comp_2, ...],
            ...
        }. Defaults to False.

    match_all: bool
        Whether to raise an error if some run indices do not match any structure directory
        in base directory. Defaults to %(default)s. Deactivate it if some structures
        could not be written as VASP input when using 'vasp_rundir_writer.py' script.

    verbose: bool
        Whether to print each reference entry used when building a phase diagram.
        Defaults to False.

    pause_after_init: bool
        Pauses the program after finishing data preparations. Press Enter to unpause.
        Defaults to False.
"""
    start = datetime.now()
    args = parse_input_args(_get_command_line_args, _process_input_args, standalone, **kwargs)

    if args["alignn"]:
        # Predict generated data energies with ALIGNN
        gfile = GenMatFile.from_file(args["input_file"])
        alignn_results = vectors_from_alignn(gfile.parse_pmg_structures(), output="energy")
        structures = {
            struct.name: {
                "entry_id": struct.name,
                "composition": struct.composition,
                "final_energy": energy.item()
            } for struct, energy in zip(gfile, alignn_results)
        }
    else:
        # Extract generated data from VASP runs
        structures = VaspExtractor(
            vasp_parser=VaspParser(args["base_dir"], match_all=args["match_all"]),
            method=ExtractMethod.CONVEX_HULL,
            summary_name=args.get("prev_summary"),
            summary_key=args.get("summary_key"),
            workers=args.get("workers")
        ).get_data()

    generated_dataset = GenMatPDDataset(
        structures,
        composition_key="composition",
        energy_key="final_energy",
        attribute="generated"
    )

    # Eliminate high dimension structures (> 10) from computations to avoid softlock
    true_max_dim = generated_dataset.max_dim
    if true_max_dim > 10:
        max_dim_generated = generated_dataset.max_computable_dim
        gen_entries = generated_dataset.computable_entries
    else: # Save filtering computation
        max_dim_generated = true_max_dim
        gen_entries = generated_dataset.all_entries

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
            dims=set(range(max_dim_generated + 1))
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
        stable_tol=args["limit"],
        workers=args.get("workers"),
        verbose=args["verbose"]
    )
    results = {}
    for data in stability.computed_data:
        dct = {
            "name": data.entry.name,
            "path": os.path.join(args["base_dir"], data.entry.name),
            "e_above_hull": data.entry.energy_above_hull,
            "stable": data.is_stable,
            "metastable": data.is_metastable
        }
        results[data.entry.name] = dct

    # Too high dimension structures are added as not stable in results
    for name, entry in generated_dataset.uncomputable_entries.items():
        msg = f"Uncomputable due to its too high dimension ({len(entry.composition)} > 10)."
        results[name] = {
                "name": entry.name,
                "path": os.path.abspath(os.path.join(args["base_dir"], entry.name)),
                "e_above_hull": None,
                "stable": False,
                "metastable": False,
                "comment": msg
            }

    sorted_results = dict(sorted(results.items(), key=lambda tup: tup[1]["path"]))

    JsonWriter(args["summarypath"], data=sorted_results, indent=4).write_as_dict()

    stop = datetime.now()
    print(f"elapsed time: {stop-start}")


if __name__ == "__main__":
    main()
