#!/usr/bin/python
"""
A script to determine material fundamental band gap from VASP energies and Δ-Sol method.

Reference for Δ-Sol method:
    M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010).
"""

import os
import json
import argparse as argp
import typing as tp
from datetime import datetime


# LOCAL IMPORTS
from . import _parse_input_args
from .utils import (
    check_type, check_num_value, check_file_format, check_file_or_dir,
    batch_extract_vasp_data, batch_get_dsol_band_gaps
)


########################################


def _get_command_line_args() -> argp.Namespace:
    """Command Line Interface (CLI)."""
    parser = argp.ArgumentParser(prog=os.path.basename(__file__), description=__doc__)
    parser.add_argument(
        "input_dir",
        type=str,
        help=(
            "Base directory containing structures directories, "
            "themself containing static calculations directories."
        ),
    )
    parser.add_argument(
        "-f", "--functional",
        type=str,
        default="PBE",
        help=(
            "DFT functional used for static calculations. "
            "Can be chosen between 'LDA', 'PBE', and 'AM05' "
            "(default: '%(default)s')."
        ),
    )
    parser.add_argument(
        "-v", "--valid-interval",
        nargs=2,
        type=float,
        help=(
            "Valid band gaps interval in eV (default: [1.3 ; 3.6] eV).\n"
            "If another interval is given, the min AND max values must be given, "
            "even if one of them matches the default values."
        ),
        metavar="float"
    )
    parser.add_argument(
        "-w", "--workers",
        type=int,
        help=(
            "Number of parallel processes to spawn for parallelized steps. "
            "If not given, If not given, default value is the 'max_workers' "
            "default value from tqdm.contrib.concurrent.process_map function."
        ),
        metavar="int",
    )
    parser.add_argument(
        "-s", "--summary",
        type=str,
        default="summary.json",
        help=(
            "Output file indicating calculation results for this step (json format).\n"
            "Also usable by further steps to filter out structures "
            "that were rejected in this step.\n"
            "This arg only changes the file name, "
            "its path is automatically set in the step directory."
        ),
        metavar="<new_file_name>"
    )
    parser.add_argument(
        "--with-uncertainties", 
        action="store_true",
        help=(
            "Enables computation of minimal and maximal Δ-Sol band gaps.\n"
            "If enabled, it will search for uncertainty calculations results."
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
    args_dict.setdefault("functional", "PBE")
    args_dict.setdefault("valid_interval", (1.3, 3.6))
    args_dict.setdefault("summary", "summary.json")
    args_dict.setdefault("with_uncertainties", False)

    # Assert set arguments conformity
    check_file_or_dir(args_dict.get("input_dir"), "dir")

    assert args_dict.get("functional") in {"LDA", "PBE", "AM05"}, (
        f"{args_dict.get('functional')} is not supported by delta-Sol. "
        "See --help for valid functional argument values."
    )
    check_type(args_dict.get("valid_interval"), "valid_interval", (tuple, list))
    assert len(args_dict["valid_interval"]) == 2, (
        "An interval should be between 2 values, "
        f"got {len(args_dict['valid_interval'])} instead."
    )
    for idx, elt in enumerate(args_dict["valid_interval"]):
        check_type(elt, f"valid_interval[{idx}]", (float, int))

    assert all(value >= 0.0 for value in args_dict["valid_interval"]), (
        "Acceptable band gap values must be positive or zero."
    )
    assert args_dict["valid_interval"][0] != args_dict["valid_interval"][1], (
        "Acceptable band gap values cannot have the exact same value."
    )
    check_file_format(args_dict.get("summary"), allowed_formats="json")

    if args_dict.get("workers") is not None:
        check_num_value(args_dict["workers"], "workers", ">", 0)

    # Additional arguments processing
    args_dict["valid_interval"] = sorted(args_dict["valid_interval"])
    args_dict["summary"] = os.path.basename(args_dict["summary"])
    args_dict["summaryfile"] = os.path.join(args_dict["input_dir"], args_dict["summary"])

    return args_dict


########################################


def main(standalone: bool = True, **kwargs) -> None:
    """
    A script to determine material fundamental band gap from VASP energies and Δ-Sol method.

    Reference for Δ-Sol method:
        M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010).
    
    Args:
        standalone (bool):              Whether parsed script is used directly through
                                        command-line (stand-alone script) or in an external
                                        pipeline script.

        input_dir (str|Path):           Base directory containing structures directories,
                                        themself containing static calculations directories.

        functional (str):               DFT functional used for static calculations.
                                        Can be chosen between 'LDA', 'PBE', and 'AM05'.
                                        Defaults to 'PBE'.

        valid_interval ((int, int)):    Valid band gaps interval in eV (default: [1.3 ; 3.6] eV).
                                        If another interval is given, the min AND max values must
                                        be given, even if one of them matches the default values.

        workers (int):                  Number of parallel processes to spawn for parallelized
                                        steps. If not given, default value is the 'max_workers'
                                        default value from tqdm.contrib.concurrent.process_map
                                        function.

        summmary (str):                 Output file indicating calculation results for this step
                                        (json format). Also usable by further steps to filter out
                                        structures rejected in this step. This arg only changes
                                        the file name, its path is automatically set in the step
                                        directory.

        with_uncertainties (bool):      Enables computation of minimal and maximal Δ-Sol band gaps.
                                        If enabled, it will search for uncertainty calculations
                                        results.
    """
    start = datetime.now()
    args = _parse_input_args(_get_command_line_args, _process_input_args, standalone, **kwargs)

    # Extract VASP static calculations results
    bg_data = batch_extract_vasp_data(
        method="delta_sol_calc",
        base_dir=args["input_dir"],
        workers=args.get("workers")
    )

    e_band_gaps = batch_get_dsol_band_gaps(
        bg_data, args["functional"], args["with_uncertainties"], args.get("workers")
    )

    good_bg_structs = list(filter(
        lambda tup: min(args["valid_interval"]) <= tup[1]["E_band_gap"] <= max(args["valid_interval"]),
        list(e_band_gaps.items())
    ))

    bad_bg_structs = list(filter(
        lambda tup: not min(args["valid_interval"]) <= tup[1]["E_band_gap"] <= max(args["valid_interval"]),
        list(e_band_gaps.items())
    ))

    screening_results: list[dict[str, tp.Any]] = []

    # Keep good structures
    for struct in good_bg_structs:
        name       = struct[0]
        bgdict     = struct[1]
        e_band_gap = round(bgdict["E_band_gap"], 6)
        e_band_gap_rectified = max(e_band_gap, 0.0)
        true_neg_bg = f" (true measurement: {e_band_gap})" if e_band_gap_rectified == 0.0 else ""

        struct_dict: dict[str, tp.Any] = {
                "path": os.path.join(str(args.get("input_dir")), name),
                "bandgap (eV)": f"{e_band_gap_rectified}{true_neg_bg}",
                "valid_gap": True
        }

        if args.get("with_uncertainties"):
            e_band_gap_min = round(bgdict["E_band_gap_min"], 6)
            e_band_gap_max = round(bgdict["E_band_gap_max"], 6)
            e_band_gap_min_rectified = max(e_band_gap_min, 0.0)
            e_band_gap_max_rectified = max(e_band_gap_max, 0.0)
            e_min = min(e_band_gap_min_rectified, e_band_gap_max_rectified)
            e_max = max(e_band_gap_min_rectified, e_band_gap_max_rectified)
            true_neg_bg_min = (
                f" (true measurement: {min(e_band_gap_min, e_band_gap_max)})" 
                if e_min == 0.0 else ""
            )
            true_neg_bg_max = (
                f" (true measurement: {max(e_band_gap_min, e_band_gap_max)})" 
                if e_max == 0.0 else ""
            )
            struct_dict.update(
                {
                    "bandgap_min (eV)": f"{e_min}{true_neg_bg_min}",
                    "bandgap_max (eV)": f"{e_max}{true_neg_bg_max}"
                }
            )

        screening_results.append(struct_dict)

    # Reject unsuitable structures
    for struct in bad_bg_structs:
        name       = struct[0]
        bgdict     = struct[1]
        e_band_gap = round(bgdict["E_band_gap"], 6)
        e_band_gap_rectified = max(e_band_gap, 0.0)
        true_neg_bg = f" (true measurement: {e_band_gap})" if e_band_gap_rectified == 0.0 else ""

        struct_dict = {
                "path": os.path.join(str(args.get("input_dir")), name),
                "bandgap (eV)": f"{e_band_gap_rectified}{true_neg_bg}",
                "valid_gap": False
        }

        if args.get("with_uncertainties"):
            e_band_gap_min = round(bgdict["E_band_gap_min"], 6)
            e_band_gap_max = round(bgdict["E_band_gap_max"], 6)
            e_band_gap_min_rectified = max(e_band_gap_min, 0.0)
            e_band_gap_max_rectified = max(e_band_gap_max, 0.0)
            e_min = min(e_band_gap_min_rectified, e_band_gap_max_rectified)
            e_max = max(e_band_gap_min_rectified, e_band_gap_max_rectified)
            true_neg_bg_min = (
                f" (true measurement: {min(e_band_gap_min, e_band_gap_max)})"
                if e_min == 0.0 else ""
            )
            true_neg_bg_max = (
                f" (true measurement: {max(e_band_gap_min, e_band_gap_max)})"
                if e_max == 0.0 else ""
            )
            struct_dict.update(
                {
                    "bandgap_min (eV)": f"{e_min}{true_neg_bg_min}",
                    "bandgap_max (eV)": f"{e_max}{true_neg_bg_max}"
                }
            )

        screening_results.append(struct_dict)

    def sort_by_path(dct: dict) -> str:
        return dct["path"]

    screening_results = sorted(screening_results, key=sort_by_path)

    with open(args["summaryfile"], "wt", encoding="utf-8") as fp:
        json.dump(screening_results, fp, indent=4)

    stop = datetime.now()
    print(f"elapsed time: {stop-start}")


if __name__ == "__main__":
    main()
