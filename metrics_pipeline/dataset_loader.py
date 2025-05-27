#!/usr/bin/python
"""Get MP and/or process OQMD database entries for the pipeline."""

import os
import argparse as argp
import typing as typ

# LOCAL IMPORTS
from . import _parse_input_args
from .utils import (
    check_type, check_file_or_dir, mp_api_download, process_oqmd_json_file
)


########################################


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
        "--process-oqmd",
        help=(
            "Path to JSON file in OQMD RESTful API format with at least the 'entry_id', "
            "'composition', 'natoms' and 'delta_e' fields. The data will be processed to "
            "fit data formatting needed by phase stability computing script. "
            "The processed data is written to a file with the same name as the one given "
            "but with a '_proc' suffix, at the same location as the one given."
        ),
        metavar="<path>"
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

    # Assert set arguments conformity
    if args_dict.get("from_mp_api") is not None:
        check_file_or_dir(os.path.dirname(args_dict["from_mp_api"]), "dir")
    check_file_or_dir(args_dict.get("process_oqmd"), "file", allowed_formats="json")

    return args_dict


########################################


def main(standalone: bool = True, **kwargs) -> None:
    """
    Get MP and/or process OQMD database entries for the pipeline.
    
    Args:
        standalone (bool):                  Whether parsed script is used directly through
                                            command-line (stand-alone script) or in an external
                                            pipeline script.

        from_mp_api (str|Path, optional):   Path to create a JSON file to store the
                                            Materials Project data. Fetches entries with
                                            thermodynamic data via the Materials Project's
                                            REST API, and store them in a json file usable
                                            for the phase_diag_energies.py script as its
                                            '--reference' argument.

        mp_api_key (str):                   If from_mp_api is used, you can either enter manually
                                            a valid MP API key here or set 'PMG_MAPI_KEY' in
                                            .pmgrc.yaml for the program to be able to fetch the
                                            data from the Materials Project API.

        process_oqmd (str|Path, optional):  Path to JSON file in OQMD RESTful API format with at
                                            least the 'entry_id', 'composition', 'natoms' and
                                            'delta_e' fields. The data will be processed to fit
                                            data formatting needed by phase stability computing
                                            script. Processed data is written to a file with the
                                            same name as the one given but with a '_proc' suffix,
                                            at the same location as the one given.
    """
    args = _parse_input_args(_get_command_line_args, _process_input_args, standalone, **kwargs)

    if args.get("from_mp_api") is not None:
        dir_path = os.path.dirname(args["from_mp_api"])
        check_file_or_dir(dir_path, "dir")
        mp_api_download(args["from_mp_api"], api_key=args.get("mp_api_key"))

    if args.get("process_oqmd") is not None:
        check_file_or_dir(args.get("process_oqmd"), "file", allowed_formats="json")
        process_oqmd_json_file(path=args["process_oqmd"])


if __name__ == "__main__":
    main()
