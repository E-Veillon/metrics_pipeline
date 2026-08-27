#!/usr/bin/python
"""
Download and process the Materials Project phase diagram database directly from their API,
or process a local JSON database file to make it compatible with GenMat phase diagram computations.
"""

import os
import argparse as argp
import typing as typ

from metrics_pipeline.core.utils import parse_input_args, check_type
from metrics_pipeline.core.genmat_io import check_file_or_dir, PDDataset, MPDatasetDownloader


def _get_command_line_args() -> argp.Namespace:
    """Command Line Interface (CLI)."""
    parser = argp.ArgumentParser(prog=os.path.basename(__file__), description=__doc__)
    parser.add_argument(
        "--from-mp-api",
        help=(
            "Path to create a JSON file to store the Materials Project data."
            "Fetches entries with thermodynamic data via the Materials Project's REST API, "
            "and store them in a json file usable for the phase_diag_energies.py script as its "
            "'--reference' argument."
        ),
        metavar="<path>"
    )
    parser.add_argument(
        "-k", "--mp-api-key",
        help=(
            "If the --from-mp-api arg is used, you can either enter manually a valid MP "
            "API key here or set 'PMG_MAPI_KEY' in .pmgrc.yaml for the program to be able "
            "to fetch the data from the Materials Project API."
        )
    )
    parser.add_argument(
        "--process-dataset",
        help=(
            "Path to JSON dataset file with structure data suitable for phase diagrams, "
            "i.e. an id, a composition or a formula and a number of atoms, and an energy "
            "or energy per atom. The data will be processed to fit data formatting needed "
            "by phase stability script. The processed data is written to a file with the "
            "same name as the one given but with a '_genmat' suffix, at the same location "
            "as the one given."
        ),
        metavar="<path>"
    )
    parser.add_argument(
        "--id_key",
        help=(
            "Only used and mandatory if 'process-dataset' is given. "
            "Give the name of the key storing structure IDs."
        )
    )
    parser.add_argument(
        "--composition_key",
        help=(
            "Only used if 'process-dataset' is given. "
            "Give the name of the key storing structure compositions. "
            "If not given, 'forumla_key' and 'natoms_key' must be given."
        )
    )
    parser.add_argument(
        "--formula_key",
        help=(
            "Only used if 'process-dataset' is given. "
            "Give the name of the key storing structure formulas."
            "If not given, 'composition_key' must be given."
        )
    )
    parser.add_argument(
        "--natoms_key",
        help=(
            "Only used if 'process-dataset' is given. "
            "Give the name of the key storing structure number of atoms."
            "If not given, 'composition_key' must be given."
        )
    )
    parser.add_argument(
        "--energy_key",
        help=(
            "Only used if 'process-dataset' is given. "
            "Give the name of the key storing structure total energy."
            "If not given, 'energy_per_atom_key' must be given."
        )
    )
    parser.add_argument(
        "--energy_per_atom_key",
        help=(
            "Only used if 'process-dataset' is given. "
            "Give the name of the key storing structure energy per atom."
            "If not given, 'energy_key' must be given."
        )
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
    args = parser.parse_args()
    return args


def _process_input_args(args_dict: dict[str, typ.Any]) -> dict[str, typ.Any]:
    """Handle input arguments assertions and processing."""
    if args_dict is None:
        raise ValueError(f"No arguments found at '{os.path.basename(__file__)}' script call.")

    check_type(args_dict, "args_dict", (dict,))

    # Set default values for unset optional arguments
    args_dict = {k: v for k, v in args_dict.items() if v is not None}
    args_dict.setdefault("compact", False)

    # Assert set arguments conformity
    if args_dict.get("from_mp_api") is not None:
        check_file_or_dir(os.path.dirname(args_dict["from_mp_api"]), "dir")

    if args_dict.get("process_dataset") is not None:
        check_file_or_dir(args_dict.get("process_dataset"), "file", allowed_formats="json")
        args_dict.setdefault("outfile", args_dict["process_dataset"].replace(".json", "_genmat.json"))
        check_type(args_dict.get("id_key"), "id_key", (str,))
        assert (
            args_dict.get("composition_key") or
            (args_dict.get("formula_key") and args_dict.get("natoms_key"))
        ), ValueError(f"Either 'composition_key' or 'formula_key' and 'natoms_key' must be given.")
        assert args_dict.get("energy_key") or args_dict.get("energy_per_atom_key"), ValueError(
            f"Either 'energy_key' or 'energy_per_atom_key' must be given."
        )
        check_type(args_dict.get("compact"), "compact", (bool,))

    return args_dict


def main(standalone: bool = True, **kwargs) -> None:
    """
    Get MP and/or process OQMD database entries for the pipeline.
    
    Parameters
    ----------
    standalone: bool                  
        Whether parsed script is used directly through command-line (stand-alone script) or
        in an external pipeline script.

    from_mp_api: str | Path, optional   
        Path to create a JSON file to store the Materials Project data. Fetches entries with
        thermodynamic data via the Materials Project's REST API, and store them in a json file
        usablefor the phase_diag_energies.py script as its'--reference' argument.

    mp_api_key: str, optional
        If `from_mp_api` is used, you can either enter manually a valid MP API key here or set
        `PMG_MAPI_KEY` in .pmgrc.yaml for the program to be able to fetch the data from the
        Materials Project API.

    process_oqmd: str | Path, optional
        Path to JSON dataset file with structure data suitable for phase diagrams, i.e. an id,
        a composition or a formula and a number of atoms, and an energy or energy per atom.
        The data will be processed to fit data formatting needed by phase stability script.
        The processed data is written to a file with the same name as the one given but with a
        '_genmat' suffix, at the same location as the one given."
    """
    args = parse_input_args(_get_command_line_args, _process_input_args, standalone, **kwargs)

    if args.get("from_mp_api") is not None:
        MPDatasetDownloader(args["from_mp_api"], api_key=args.get("mp_api_key"))

    if args.get("process_dataset") is not None:
        processed_dataset = PDDataset.from_file(
            filepath=args["process_dataset"],
            id_key=args["id_key"],
            composition_key=args.get("composition_key"),
            formula_key=args.get("formula_key"),
            natoms_key=args.get("natoms_key"),
            energy_key=args.get("energy_key"),
            energy_per_atom_key=args.get("energy_per_atom_key"),
            compact=args["compact"]
        )
        processed_dataset.write_dataset(args["outfile"], compact=args["compact"])


if __name__ == "__main__":
    main()
