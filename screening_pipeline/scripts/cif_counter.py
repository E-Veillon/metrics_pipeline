#!/usr/bin/python
"""Counts the number of structures in given CIF files."""

import os
from argparse import ArgumentParser
from typing import List


########################################


def has_str(string: str, pattern: str) -> bool:
    """
    Returns True if pattern is in string, otherwise returns False.

    Parameters:
        string (str): string to search in.
        pattern (str): which string to search.
    
    Returns:
        Whether the searched string was found or not.
    """
    return string.find(pattern) > -1

def write_cut_file(
        filename: str,
        file_lines: List[str],
        data_breakpoints: List[str],
        cut_nbr: int,
    ) -> None:
    """
    Searches in the file the right line where it should be truncated, 
    then writes a truncated copy with a suffix denoting the number of 
    structures kept in the file name.

    Parameters:
        filename (str):             The name of the original file.
    
        file_lines ([str]):         The original file data as a list of its lines.

        data_breakpoints ([str]):   List of the lines starting with "data_" denoting
                                    the beginning of the data of each structure.

        cut_nbr (int):              The number of structures to keep in the truncated file.
    """
    wrong_cut = False
    idx_changed_lines = []
    cut_line = data_breakpoints[cut_nbr]

    # Case where an identical line appears before the one wanted
    if data_breakpoints.index(cut_line) < cut_nbr:
        nbr_iter = 0

        while True:
            idx_same_line = data_breakpoints.index(cut_line)

            if idx_same_line == cut_nbr:
                break

            if idx_same_line > cut_nbr:
                print(
                    "Unexpected behaviour happened while trying to cut the file.\n"
                    "The cutting option will stop here without doing anything."
                )
                wrong_cut = True
                break

            if idx_same_line < cut_nbr:
                idx_changed_lines.append(idx_same_line)
                data_breakpoints[idx_same_line] = f"data_xx{nbr_iter}xx"
                file_lines[file_lines.index(cut_line)] = data_breakpoints[idx_same_line]
                nbr_iter += 1

    if not wrong_cut:
        cut_lines = file_lines[:file_lines.index(cut_line)]

        if cut_lines[-1] == "# generated using pymatgen":
            cut_lines = cut_lines[:-1]

        for idx in idx_changed_lines:
            line_to_restore = data_breakpoints[idx]
            idx_to_restore = file_lines.index(line_to_restore)
            cut_lines[idx_to_restore] = cut_line

        cut_text = "\n".join(cut_lines)
        cut_file = filename.replace(".cif", f"_{cut_nbr}.cif")

        with open(cut_file, mode="wt", encoding="utf-8") as outfile:
            outfile.write(cut_text)

        print(
            f"The file '{cut_file}' containing the first {cut_nbr} structure(s) "
            f"from '{filename}' was successfully created."
        )


########################################


def main() -> None:
    """Main function."""

    # ARGUMENTS PARSING BLOCK

    parser = ArgumentParser(
        prog="cif_counter.py",
        description="Counts the number of structures in given CIF files."
    )

    parser.add_argument(
        "filenames", 
        nargs="+",
        help="files to count structures in."
    )
    parser.add_argument(
        "--cut",
        type=int,
        help=(
            "create a truncated copy of the tested file containing only given "
            "number of structures, going from the first one. If given value is not "
            "integer compatible, it will be ignored."
        ),
        metavar="int"
    )

    args = parser.parse_args()
    cut_nbr = args.cut if isinstance(args.cut, int) and (args.cut >= 1) else None


    for file in args.filenames:

        if not os.path.isfile(file):
            raise FileNotFoundError(f"{file}: No such file found.")
        if not file.endswith(".cif"):
            raise ValueError(f"{file} is not a valid CIF file.")

        with open(file, mode="r", encoding="utf-8") as data:
            lines = data.read().splitlines()

        data_breakpoints = list(filter(lambda line: has_str(line, "data_"), lines))
        nbr_structs = len(data_breakpoints)

        if cut_nbr is not None and cut_nbr >= nbr_structs:
            print(
                f"Provided 'cut' arg ({cut_nbr}) is larger or equal "
                f"to the total number of structures in '{file}'.\n"
                "Therefore, as it will not change anything, "
                "the cut option will now be deactivated."
            )

        elif cut_nbr is not None:
            write_cut_file(file, lines, data_breakpoints, cut_nbr)

        print(f"{nbr_structs} structures found in file '{file}'.")


if __name__ == "__main__":
    main()
