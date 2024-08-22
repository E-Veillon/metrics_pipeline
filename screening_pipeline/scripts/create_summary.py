#!/usr/bin/python
"""A script to generate a JSON summary from manually rejected structures."""

import os
import json
import argparse

# LOCAL IMPORTS
from screening_pipeline.utils import (
    check_file_format, check_file_or_dir, match_struct_dirs
)


########################################


def _assert_args(args: argparse.Namespace) -> None:
    """Check input args validity."""
    check_file_or_dir(args.input_dir, "dir")
    check_file_format(args.output, allowed_formats="json")
    match_struct_dirs(args.input_dir, args.indices, no_return=True)


########################################


def main() -> None:
    """Main function."""

    # ARGUMENTS PARSING BLOCK

    parser = argparse.ArgumentParser(
        prog="create_summary",
        description="A script to generate a JSON summary from manually rejected structures."
    )

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
            "Name of the output file. Must be a JSON format (i.e. with '.json' suffix). "
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

    _assert_args(args)

    if args.indices is None:
        indices = args.indices
    else:
        indices = sorted(args.indices)

    outfile = os.path.join(args.input_dir, args.output)


    # MAIN BLOCK

    all_struct_dirs = match_struct_dirs(args.input_dir)
    wanted_struct_dirs = match_struct_dirs(args.input_dir, indices)
    summary_result = []

    for struct_dir in sorted(all_struct_dirs):
        summary_result.append(
            {
                "path": struct_dir,
                "name": os.path.basename(struct_dir),
                f"{args.key_name}": (struct_dir in wanted_struct_dirs) ^ args.reject
            }
        )

    with open(outfile, "wt", encoding="utf-8") as fp:
        json.dump(summary_result, fp, indent=4)


if __name__ == "__main__":
    main()
