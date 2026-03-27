"""
Compare the distribution of spacegroup symmetries between generated and reference datasets.
"""

import typing as tp
import argparse as ap
import os

from src.utils import parse_input_args, check_type, check_num_value
from src.io import check_file_or_dir, check_file_format
from src.metrics import SymmetryClassifier


def _get_cmd_line_args() -> ap.Namespace:
    parser = ap.ArgumentParser()
    parser.add_argument(
        "generated", help="Path to generated structures file."
    )
    parser.add_argument(
        "reference", help="Path to reference structures file."
    )
    parser.add_argument(
        "-o", "--output", default="symmetry_distribution_comparison.json",
        help=(
            "Path to output file to write results to. "
            "If not given, defaults to %(default)s in 'generated' file directory."
        )
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
        os.path.dirname(args_dict["generated"]),
        "symmetry_distribution_comparison.json"
    )
    args_dict.setdefault("output", default_output)
    args_dict.setdefault("symprec", 0.01)
    args_dict.setdefault("angleprec", 5.0)
    args_dict.setdefault("on_error", "warn")

    # Assert set arguments conformity
    check_file_or_dir(args_dict["generated"], "file", allowed_formats="cif")
    check_file_or_dir(args_dict["reference"], "file", allowed_formats="cif")
    check_file_format(args_dict["output"], allowed_formats="json")

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
    Compare the distribution of spacegroup symmetries between generated and reference datasets.

    Parameters
    ----------
    standalone: bool
        Whether to run the function in standalone mode, i.e. parsing command line arguments.
        Defaults to True.

    **kwargs:
        Keyword arguments to be passed to `SymmetryClassifier` metric class.
    """
    raise NotImplementedError("This script is not yet implemented.")
    args = parse_input_args(_get_cmd_line_args, _process_input_args, standalone, **kwargs)

    gen_classifier = SymmetryClassifier(
        args["generated"],
        symprec=args["symprec"],
        angleprec=args["angleprec"],
        workers=args.get("workers"),
        on_error=args["on_error"]
    )
    ref_classifier = SymmetryClassifier(
        args["reference"],
        symprec=args["symprec"],
        angleprec=args["angleprec"],
        workers=args.get("workers"),
        on_error=args["on_error"]
    )
