#!/usr/bin/python
"""
Update paths keys in a summary json file when the corresponding structure directories are moved.
"""

import os
import typing as typ
import argparse as argp

from src.utils import parse_input_args, check_type
from src.io import check_file_or_dir, JsonLoader, JsonWriter


def _get_command_line_args() -> argp.Namespace:
    """Command Line Interface (CLI)."""
    parser = argp.ArgumentParser(prog=os.path.basename(__file__), description=__doc__)
    parser.add_argument(
        "input_file", help="Path to the JSON file to update paths in."
    )
    parser.add_argument(
        "new_path", help="Path to the directory where structure directories are actually stored."
    )
    parser.add_argument(
        "-a", "--absolute",
        action="store_true",
        help="The new path is absolutized before replacing the old one."
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
    args_dict.setdefault("absolute", False)

    # Assert set arguments conformity
    check_file_or_dir(args_dict.get("input_file"), "file", allowed_formats="json")
    check_file_or_dir(args_dict.get("new_path"), "dir")
    check_type(args_dict.get("absolute"), "absolute", (bool,))

    return args_dict


def main(standalone: bool = True, **kwargs) -> None:
    """
    Update paths keys in a summary json file when the corresponding structure directories
    are moved.

    Parameters
    ----------
    standalone: bool
        Whether parsed script is used directly through command-line (stand-alone script)
        or in an external pipeline script.

    input_file: str | Path
        Path to the JSON file to update paths in.

    new_path: str | Path
        Path to the directory where structure directories are actually stored.

    absolute: bool
        If set to True, the new path is absolutized before replacing the old one.
        Defaults to False.
    """
    args = parse_input_args(_get_command_line_args, _process_input_args, standalone, **kwargs)

    data = JsonLoader(args["input_file"]).load_as_list()

    assert all(isinstance(struct, dict) for struct in data), (
        "The data inside the JSON file must be a list of structure dicts."
    )

    if args.get("absolute"):
        new_path = os.path.realpath(args["new_path"])
    else:
        new_path = args["new_path"]

    for struct in data:
        old_path = struct["path"]
        struct_name = os.path.basename(old_path)
        struct["path"] = os.path.join(new_path, struct_name)

    JsonWriter(args["input_file"], data, indent=4).write_as_list()

    file = args["input_file"]
    print(f"the file '{os.path.basename(file)}' was successfully modified.")


if __name__ == "__main__":
    main()
