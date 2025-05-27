#!/usr/bin/python
"""A script to generate a JSON summary from manually selected structures."""

import os
import json
import argparse as argp
import typing as typ

# LOCAL IMPORTS
from . import _parse_input_args
from .utils import (
    check_type, check_file_format, check_file_or_dir, match_struct_dirs
)


########################################
# ARGUMENTS HANDLING

def _get_command_line_args() -> argp.Namespace:
    """Command Line Interface (CLI)."""
    parser = argp.ArgumentParser(prog=os.path.basename(__file__), description=__doc__)

    parser.add_argument(
        "input_dir",
        help="Directory containing structure directories."
    )
    parser.add_argument(
        "-i", "--indices",
        nargs="*",
        type=int,
        help=(
            "The indices of the structures of interest "
            "(the ID number before each structure directory name). "
            "By default, they are the structures to keep "
            "(i.e. having a 'selected = True' key). "
            "All other structures detected in input_dir will get the opposite value."
            "If not given, all structures will be given the same bool value "
            "(i.e. all set to 'true' by default, or all set to 'false' with --reject flag)."
        )
    )
    parser.add_argument(
        "-o", "--output",
        default="manual_summary.json",
        help=(
            "Name of the output file. Must be a JSON format. "
            "Only affect the name of the file, its location is the path given as input_dir. "
            "Defaults to 'manual_summary.json'."
            )
    )
    parser.add_argument(
        "-n", "--key-name",
        default="selected",
        help="Name for the dict key containing the selection boolean. Defaults to 'selected'."
    )
    parser.add_argument(
        "--reject",
        action="store_true",
        help=(
            "Pass this flag to reject (i.e. having a 'selected = False' key) "
            "given structures instead of saving them (i.e. having a 'selected = True' key)."
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
    args_dict.setdefault("output", "manual_summary.json")
    args_dict.setdefault("key_name", "selected")
    args_dict.setdefault("reject", False)

    # Assert set arguments conformity
    check_file_or_dir(args_dict.get("input_dir"), "dir")
    check_file_format(args_dict.get("output"), allowed_formats="json")
    match_struct_dirs(args_dict["input_dir"], args_dict.get("indices"), no_return=True)

    # Additional arguments processing
    if args_dict.get("indices") is not None:
        args_dict["indices"] = sorted(args_dict["indices"])

    args_dict["outfile"] = os.path.join(args_dict["input_dir"], args_dict["output"])

    return args_dict


########################################


def main(standalone: bool = True, **kwargs) -> None:
    """
    A script to generate a JSON summary from manually selected structures.
    
    Args:
        standalone (bool):      Whether parsed script is used directly through
                                command-line (stand-alone script) or in an external
                                pipeline script.

        input_dir (str|Path):   Directory containing structure directories.

        indices ([int]):        List of the indices of the structures of interest
                                (the ID number before each structure directory name).
                                By default, they are the structures to keep
                                (i.e. having a 'selected = True' key).
                                All other structures detected in input_dir will get
                                the opposite value. If not given, all structures will be given
                                the same bool value (i.e. all set to 'true' by default, or all
                                set to 'false' with reject = True).

        output (str|Path):      Name of the output file. Must be a JSON format.
                                Only affect the name of the file, its location is the path
                                given as input_dir. Defaults to 'manual_summary.json'.

        key_name (str):         Name for the dict key containing the selection boolean.
                                Defaults to 'selected'.

        reject (bool):          If True, reject (i.e. having a 'selected = False' key)
                                given structures instead of saving them
                                (i.e. having a 'selected = True' key).
    """
    args = _parse_input_args(_get_command_line_args, _process_input_args, standalone, **kwargs)

    all_struct_dirs = match_struct_dirs(args["input_dir"])
    wanted_struct_dirs = match_struct_dirs(args["input_dir"], args.get("indices"))
    summary_result = []

    for struct_dir in sorted(all_struct_dirs):
        summary_result.append(
            {
                "path": struct_dir,
                "name": os.path.basename(struct_dir),
                f"{args.get('key_name')}": (struct_dir in wanted_struct_dirs) ^ args["reject"]
            }
        )

    with open(args["outfile"], "wt", encoding="utf-8") as fp:
        json.dump(summary_result, fp, indent=4)


if __name__ == "__main__":
    main()