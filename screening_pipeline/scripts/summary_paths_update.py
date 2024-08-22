#!/usr/bin/python
"""
Update paths keys in a summary json file when the corresponding structure directories are moved.
"""

import os
import json
from typing import List, Dict
import argparse


########################################


def _assert_args(args: argparse.Namespace) -> None:
    """Asserting input arguments validity."""
    if not os.path.isfile(args.input_file):
        raise FileNotFoundError(
            f"{args.input_file}: No such file found."
        )
    if not args.input_file.endswith(".json"):
        raise ValueError(
            f"{args.input_file}: allowed file format is 'json', "
            f"got '{args.input_file.split(sep='.')[-1]}' format instead."
        )
    if not os.path.isdir(args.new_path):
        raise FileNotFoundError(
            f"{args.input_file}: No such directory found."
        )


########################################


def main() -> None:
    """Main function."""

    # ARGUMENTS PARSING BLOCK

    parser = argparse.ArgumentParser(
        description="A command-line tool to modify path designation in JSON summary files."
    )
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

    _assert_args(args)


    # MAIN BLOCK

    with open(args.input_file, "rt", encoding="utf-8") as fp:
        data = json.load(fp)

    assert (
        isinstance(data, List)
        and all(isinstance(struct, Dict) for struct in data)
    ), "The data inside the JSON file must be a list of structure dicts."

    if args.absolute:
        new_path = os.path.abspath(args.new_path)
    else:
        new_path = args.new_path

    for struct in data:
        old_path = struct["path"]
        struct_name = os.path.basename(old_path)
        struct["path"] = os.path.join(new_path, struct_name)

    with open(args.input_file, "wt", encoding="utf-8") as fp:
        json.dump(data, fp, indent=4)

    print(f"the file '{os.path.basename(args.input_file)}' was successfully modified.")


if __name__ == "__main__":
    main()
