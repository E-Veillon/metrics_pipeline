#!/usr/bin/python
"""Counts the number of structures in given CIF files."""

import os
import argparse as argp
import typing as typ

# LOCAL IMPORTS
from . import _parse_input_args
from .utils import check_type


########################################
# ARGUMENTS HANDLING

def _get_command_line_args() -> argp.Namespace:
    """Command Line Interface (CLI)."""
    parser = argp.ArgumentParser(prog=os.path.basename(__file__), description=__doc__)

    parser.add_argument(
        "filenames", 
        nargs="+",
        help="files to count structures in."
    )
    parser.add_argument(
        "--cut", "-c",
        nargs="*",
        type=int,
        help=(
            "Create a truncated copy of the tested file containing only some of the "
            "structures data from original file. Takes 1 or more integer values depending "
            "on chosen cutting algorithm. See '--how-to-cut' for more infos."
        ),
        metavar="int"
    )
    parser.add_argument(
        "--how-to-cut",
        default="from_start",
        help=(
            "Which cutting algorithm to use. "
            "Supports 'from_start', keeping the first <cut value> data (default); "
            "'from_end', keeping the last <cut value> data; "
            "'interval', keeping data between <cut value 1> and <cut value 2> indices (limits included); "
            "'select', keeping only data at given <list of cut values> indices. "
            "Indexation of data is zero-based, and in the file order."
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
    args_dict.setdefault("how_to_cut", "from_start")

    # Assert set arguments conformity
    check_type(args_dict.get("filenames"), "filenames", (list, tuple))

    if args_dict.get("cut") is not None:
        for val in args_dict["cut"]:
            check_type(val, "--cut", (int,))

    valid_cut_algos = ("from_start", "from_end", "interval", "select")
    if args_dict.get("how_to_cut") not in valid_cut_algos:
        raise ValueError(
            f"'how-to-cut' arg must be one of {', '.join(valid_cut_algos)}, "
            f"got {args_dict.get('how_to_cut')} instead."
        )
    if args_dict.get("cut") is not None and (
        args_dict.get("how_to_cut") in ("from_start", "from_end") and len(args_dict["cut"]) != 1
        or args_dict.get("how_to_cut") == "interval" and len(args_dict["cut"]) != 2
    ):
        raise ValueError(
            "The number of values given to the 'cut' argument does not match with the "
            "algorithm given to the 'how-to-cut' argument. See --help for more infos "
            "on how to set those values properly."
        )

    return args_dict


########################################


def write_cut_file(
        filename: str,
        file_lines: list[str],
        data_breakpoints: list[tuple[int, str]],
        cut_idx: list[int],
        cut_algo: str,
    ) -> None:
    """
    Create a truncated copy of given CIF file according to given algorithm.

    Parameters:
        filename (str):                     The name of the original file.
    
        file_lines ([str]):                 The original file data as a list of its lines.

        data_breakpoints ([(int, str)]):    List of the lines starting with "data_" denoting
                                            the beginning of the data of each structure.

        cut_idx ([int]):                    Indices indicating which data to keep.

        cut_algo (str):                     Cutting algorithm to use for the truncation.
    """
    valid_cut_idx = list(filter(lambda idx: 0 <= idx < len(data_breakpoints), cut_idx))
    skipped_idx = list(filter(lambda idx: idx < 0 or idx >= len(data_breakpoints), cut_idx))
    skipped_idx = list(map(str, skipped_idx))

    if not valid_cut_idx:
        print(
            "All given '--cut' argument indices are greater or equal "
            f"to the total number of structures in '{filename}'.\n"
            "Cutting algorithm skipped."
        )
        return

    if skipped_idx:
        print(
            "WARNING: Some of given '--cut' argument indices are out of the range "
            f"of the total number of structures in '{filename}'.\n"
            f"These indices are ignored (skipped indices: {', '.join(skipped_idx)})."
        )

    if (
        cut_algo in ("from_start", "from_end") and len(valid_cut_idx) != 1
        or cut_algo == "interval" and len(valid_cut_idx) != 2
    ):
        print(
            "After removing out-of-range indices, the number of left indices"
            "does not match with the cutting algorithm anymore. "
            "Cutting algorithm skipped."
        )
        return

    print("Cutting file...")

    match cut_algo:
        case "from_start":
            cut_line_idx = data_breakpoints[valid_cut_idx[0]][0]
            lines_left = file_lines[:cut_line_idx]
            cut_text = "".join(lines_left)
            suffix = f"_first_{valid_cut_idx[0]}"

        case "from_end":
            cut_line_idx = data_breakpoints[-valid_cut_idx[0]][0]
            lines_left = file_lines[cut_line_idx:]
            cut_text = "".join(lines_left)
            suffix = f"_last_{valid_cut_idx[0]}"

        case "interval":
            first_line_idx = data_breakpoints[min(valid_cut_idx)][0]
            last_line_idx = data_breakpoints[max(valid_cut_idx) + 1][0]
            lines_left = file_lines[first_line_idx:last_line_idx]
            cut_text = "".join(lines_left)
            suffix = f"_{min(valid_cut_idx)}_to_{max(valid_cut_idx)}"

        case "select":
            valid_cut_idx = sorted(valid_cut_idx)
            kept_cifs = []
            for idx in valid_cut_idx:
                data_lines_idx = (data_breakpoints[idx][0], data_breakpoints[idx + 1][0])
                selected_data = file_lines[data_lines_idx[0]:data_lines_idx[1]]
                kept_cifs.append("".join(selected_data))
            cut_text = "".join(kept_cifs)
            suffix = f"_select_{'-'.join(list(map(str, valid_cut_idx)))}"

        case _:
            valid_cut_algos = ("from_start", "from_end", "interval", "select")
            raise ValueError(
                f"'how-to-cut' arg must be one of {', '.join(valid_cut_algos)}, "
                f"got {cut_algo} instead."
            )

    cut_file = filename.replace(".cif", f"{suffix}.cif")
    with open(cut_file, mode="wt", encoding="utf-8") as outfile:
        outfile.write(cut_text)

    print(f"The file {cut_file} truncated from {filename} was successfully created.")


########################################


def main(standalone: bool = True, **kwargs) -> None:
    """
    Counts the number of structures in given CIF files.
    
    Args:
        standalone (bool):      Whether parsed script is used directly through
                                command-line (stand-alone script) or in an external
                                pipeline script.

        filenames ([str|Path]): Files to count structures in.

        cut ([int]):            Create a truncated copy of the tested file containing only some
                                of the structures data from original file. Takes 1 or more integer
                                values depending on chosen cutting algorithm.
                                See 'how-to-cut' for more infos.

        how_to_cut (str):       Which cutting algorithm to use.
                                Supports 'from_start', keeping the first <cut value> data (default);
                                'from_end', keeping the last <cut value> data;
                                'interval', keeping data between <cut value 1> and <cut value 2>
                                indices (limits included);
                                'select', keeping only data at given <list of cut values> indices.
                                Indexation of data is zero-based, and in the file order.
    """
    args = _parse_input_args(_get_command_line_args, _process_input_args, standalone, **kwargs)

    for file in args["filenames"]:

        if not os.path.isfile(file):
            print(f"FileNotFoundError: {file}: No such file found, skipped.")
            continue

        if not str(file).endswith(".cif"):
            print(f"ValueError: {file} is not a valid CIF file, skipped.")
            continue

        with open(file, mode="r", encoding="utf-8") as data:
            lines = data.readlines()

        data_breakpoints = list(filter(lambda idxline: idxline[1].startswith("data_"), enumerate(lines)))
        nbr_structs = len(data_breakpoints)
        print(f"{nbr_structs} structures found in file '{file}'.")

        if args.get("cut") is not None:
            write_cut_file(file, lines, data_breakpoints, args["cut"], args["how_to_cut"])


if __name__ == "__main__":
    main()
