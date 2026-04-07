"""
Compute space group symmetry of crystals and produce a report of distribution at each symmetry level.
"""

import os
import typing as tp
import argparse as ap

from src.utils import check_type, check_num_value, parse_input_args
from src.io import check_file_or_dir, check_file_format, CIFFile, PoscarFile
from src.metrics import SymmetryClassifier


def _get_cmd_line_args() -> ap.Namespace:
    """Command line handler."""
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input_file",
        help="Path to the concatenated CIF or POSCAR file to extrtact crystals from."
    )
    parser.add_argument(
        "-o", "--output",
        help="path to write symmetry distribution report (plain '.txt' file)."
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
        "--symprec", type=float, default=0.01, metavar="float",
        help="Fractional coordinates tolerance for symmetry finding (Default: %(default)s).",
    )
    parser.add_argument(
        "--angleprec", type=float, default=5.0, metavar="float",
        help="Angle tolerance for symmetry finding in degrees (Default: %(default)s degrees).",
    )
    parser.add_argument(
        "--on-error", type=str, default="warn", choices=["raise", "warn", "ignore"], metavar="str",
        help=(
            "What to do in case the symmetry of a structure cannot be found. Defaults to %(default)s. "
            "If not set to 'raise', the structure will be added to a special 'uncomputable' category."
        )
    )
    parser.add_argument(
        "-v", "--verbose", action=ap.BooleanOptionalAction, default=False,
        help=(
            "Whether to add the list of headers in each symmetry class in report. "
            "Defaults to %(default)s."
        )
    )
    args = parser.parse_args()
    return args


def _process_input_args(args_dict: dict[str, tp.Any]) -> dict[str, tp.Any]:
    """Handle input arguments assertions and processing."""
    if args_dict is None:
        raise ValueError(f"No arguments found at '{os.path.basename(__file__)}' script call.")

    check_type(args_dict, "args_dict", (dict,))

    # Set default values for unset optional arguments
    args_dict = {k: v for k, v in args_dict.items() if v is not None}
    default_output = os.path.join(
        os.path.dirname(args_dict["input_file"]),
        "symmetry_distribution.txt"
    )
    args_dict.setdefault("output", default_output)
    args_dict.setdefault("symprec", 0.01)
    args_dict.setdefault("angleprec", 5.0)
    args_dict.setdefault("on_error", "warn")
    args_dict.setdefault("verbose", False)

    # Assert set arguments conformity
    check_file_or_dir(args_dict["input_file"], "file", allowed_formats="cif")
    check_file_or_dir(args_dict["reference"], "file", allowed_formats="cif")
    check_file_format(args_dict["output"], allowed_formats="txt")

    if args_dict.get("workers") is not None:
        check_type(args_dict["workers"], "workers", (int,))
        check_num_value(args_dict["workers"], "workers", ">=", 0)

    check_type(args_dict["symprec"], "symprec", (float,))
    check_num_value(args_dict["symprec"], "symprec", ">=", 0.0)
    check_num_value(args_dict["symprec"], "symprec", "<=", 1.0)
    check_type(args_dict["angleprec"], "angleprec", (float,))
    check_num_value(args_dict["angleprec"], "angleprec", ">=", 0.0)
    check_num_value(args_dict["angleprec"], "angleprec", "<=", 90.0)
    check_type(args_dict["on_error"], "on_error", (str,))
    if args_dict["on_error"] not in {"raise", "warn", "ignore"}:
        raise ValueError(
            f"'on_error' only supports 'raise', 'warn' or 'ignore', got {args_dict['on_error']!r}."
        )
    return args_dict


def main(standalone: bool = True, **kwargs) -> None:
    """
    Compute space group symmetry of crystals and produce a report of distribution at each
    symmetry level.

    Parameters
    ----------
    standalone: bool
        Whether parsed script is used directly through command-line (stand-alone script)
        or in an external pipeline script.

    input_file: str | Path
        Path to the concatenated CIF or POSCAR file to extrtact crystals from.

    output: str | Path, optional
        path to write symmetry distribution report (plain text file).

    workers: int, optional
        Number of processes to use in parallel. If not given, will use default of
        `tqdm.contrib.concurrent.process_map()`. Pass 0 to disable `process_map()`
        and execute sequentially.

    symprec: float
        Fractional coordinates tolerance for symmetry finding. Defaults to 0.01.

    angleprec: float
        Angle tolerance for symmetry finding in degrees. Defaults to 5.0 degrees.

    on_error: 'raise' | 'warn' | 'ignore'
        What to do in case the symmetry of a structure cannot be found. Defaults to `warn`.
        If not set to 'raise', the structure will be added to a special 'uncomputable' category.

    verbose: bool
        Whether to add the list of headers in each symmetry class in report.
        Defaults to False.
    """
    args = parse_input_args(_get_cmd_line_args, _process_input_args, standalone, **kwargs)

    # Read input file
    match os.path.splitext(args["input_file"]):
        case ".cif":
            cfile = CIFFile.from_file(
                args["input_file"],
                special_keys=["header"],
                workers=args.get("workers")
            )
            structures, _ = cfile.parse_structures()
        case ".poscar":
            pfile = PoscarFile(
                args["input_file"],
                workers=args.get("workers")
            )
            structures = pfile.parse_structures()

    # Compute symmetries
    classifier = SymmetryClassifier(
        structures,
        symprec=args["symprec"],
        angleprec=args["angleprec"],
        workers=args.get("workers"),
        on_error=args["on_error"]
    )

    # Write results to file
    classifier.write_result(args["output"], verbose=args["verbose"])


if __name__ == "__main__":
    main()