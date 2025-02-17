#!/usr/bin/python
"""A script to generate a JSON summary from manually rejected structures."""

import os
import json
import argparse as argp
import typing as typ

from .utils import (
    check_type, check_file_format, check_file_or_dir, match_struct_dirs
)


########################################
# ARGUMENTS HANDLING

def _get_command_line_args() -> argp.Namespace:
    """Command-line arguments UI."""
    parser = argp.ArgumentParser(prog="create_summary.py", description=__doc__)

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
        raise ValueError("No arguments found at 'create_summary.py' script call.")
    
    check_type(args_dict, "args_dict", (dict))

    args_dict.setdefault("output", "manual_summary.json")
    args_dict.setdefault("key_name", "selected")
    args_dict.setdefault("reject", False)

    check_file_or_dir(args_dict.get("input_dir"), "dir")
    check_file_format(args_dict.get("output"), allowed_formats="json")
    match_struct_dirs(args_dict.get("input_dir"), args_dict.get("indices"), no_return=True)

    if args_dict.get("indices") is not None:
        args_dict["indices"] = sorted(args_dict.get("indices"))

    args_dict["outfile"] = os.path.join(args_dict.get("input_dir"), args_dict.get("output"))

    return args_dict


def _parse_input_args(**kwargs) -> dict[str, typ.Any]:
    """Gather and process input arguments, either from command-line or external script."""
    if __name__ == "__main__": # Direct use of the script through command line
        args = _get_command_line_args()
        args = args.__dict__
    else: # Indirect use of the script as part of a pipeline script in another file
        args = kwargs

    args: dict[str, typ.Any] = _process_input_args(args)

    return args


########################################


def main(**kwargs) -> None:
    """
    create_summary.py: A script to generate a JSON summary from manually rejected structures.
    
    Possible kwargs:
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
    args = _parse_input_args(**kwargs)

    all_struct_dirs = match_struct_dirs(args.get("input_dir"))
    wanted_struct_dirs = match_struct_dirs(args.get("input_dir"), args.get("indices"))
    summary_result = []

    for struct_dir in sorted(all_struct_dirs):
        summary_result.append(
            {
                "path": struct_dir,
                "name": os.path.basename(struct_dir),
                f"{args.get('key_name')}": (struct_dir in wanted_struct_dirs) ^ args.get("reject")
            }
        )

    with open(args.get("outfile"), "wt", encoding="utf-8") as fp:
        json.dump(summary_result, fp, indent=4)


if __name__ == "__main__":
    main()