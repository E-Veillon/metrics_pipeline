#!/usr/bin/python
"""
A script to determine material fundamental band gap from VASP energies and Δ-Sol method.

Reference for Δ-Sol method:
    M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010).
"""

import os
import argparse as argp
import typing as tp
from datetime import datetime

from metrics_pipeline.core.utils import parse_input_args, check_type, check_num_value
from metrics_pipeline.core.genmat_io import (
    check_file_format, check_file_or_dir, VaspParser, VaspExtractor, ExtractMethod, JsonWriter
)
from metrics_pipeline.core.computations.local import DSolStructure, batch_get_dsol_band_gaps


########################################


def _get_command_line_args() -> argp.Namespace:
    """Command Line Interface (CLI)."""
    parser = argp.ArgumentParser(prog=os.path.basename(__file__), description=__doc__)
    parser.add_argument(
        "input_dir", type=str,
        help=(
            "Base directory containing structures directories, "
            "themself containing static calculations directories."
        ),
    )
    parser.add_argument(
        "-f", "--functional",
        type=str, default="PBE",
        help=(
            "DFT functional used for static calculations. "
            "Can be chosen between 'LDA', 'PBE', and 'AM05' "
            "(default: '%(default)s')."
        ),
    )
    parser.add_argument(
        "-v", "--valid-interval", nargs=2, type=float, metavar="float",
        help=(
            "Valid band gaps interval in eV (default: [1.3 ; 3.6] eV).\n"
            "If another interval is given, the min AND max values must be given, "
            "even if one of them matches the default values."
        )
    )
    parser.add_argument(
        "--workers", "-w", type=int, metavar="int",
        help=(
            "Number of processes to use in parallel. If not given, will use default of "
            "`tqdm.contrib.concurrent.process_map()`. Pass 0 to disable `process_map()` "
            "and execute sequentially."
        )
    )
    parser.add_argument(
        "-s", "--summary", type=str, default="summary.json",
        help=(
            "Output file indicating calculation results for this step (json format).\n"
            "Also usable by further steps to filter out structures "
            "that were rejected in this step.\n"
            "This arg only changes the file name, "
            "its path is automatically set in the step directory."
        ),
        metavar="<new_file_name>"
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
    
    Parameters
    ----------
    standalone: bool
        Whether parsed script is used directly through command-line (stand-alone script)
        or in an external pipeline script.

    input_dir: str | Path
        Base directory containing structures directories, themself containing static
        calculations directories.

    functional: str
         DFT functional used for static calculations. Can be chosen between 'LDA', 'PBE',
         and 'AM05'. Defaults to 'PBE'.

    valid_interval: tuple[int, int]
        Valid band gaps interval in eV (default: [1.3 ; 3.6] eV). If another interval is given,
        the min AND max values must be given, even if one of them matches the default values.

    workers: int, optional
        Number of processes to use in parallel. If not given, will use default of
        `tqdm.contrib.concurrent.process_map()`. Pass 0 to disable `process_map()`
        and execute sequentially.

    summmary: str
        Output file indicating calculation results for this step (json format). Also usable
        by further steps to filter out structures rejected in this step. This arg only changes
        the file name, its path is automatically set in the step directory.
    """
    start = datetime.now()
    args = parse_input_args(_get_command_line_args, _process_input_args, standalone, **kwargs)

    # Extract VASP static calculations results
    bg_data = VaspExtractor(
        vasp_parser=VaspParser(args["input_dir"]),
        method=ExtractMethod.DSOL_BANDGAP,
        workers=args.get("workers")
    ).get_data()

    dsol_structs = [
        DSolStructure(
            data["structure"],
            name,
            args["functional"],
            data[name + "_neutral"],
            data[name + "_best_plus"],
            data[name + "_best_minus"],
            data.get(name + "_min_plus"),
            data.get(name + "_min_minus"),
            data.get(name + "_max_plus"),
            data.get(name + "_max_minus"),
        ) for name, data in bg_data.items()
    ]

    e_band_gaps = batch_get_dsol_band_gaps(dsol_structs, args.get("workers"))

    def is_good_bg(e_band_gap: float) -> bool:
        return min(args["valid_interval"]) <= e_band_gap <= max(args["valid_interval"])

    screening_results: list[dict[str, tp.Any]] = []

    for name, bgdict in e_band_gaps.items():
        e_band_gap = round(bgdict["E_band_gap"], 6)
        e_band_gap_rectified = max(e_band_gap, 0.0)
        true_neg_bg = f" (true measurement: {e_band_gap})" if e_band_gap_rectified == 0.0 else ""

        struct_dict: dict[str, tp.Any] = {
                "path": os.path.join(str(args.get("input_dir")), name),
                "bandgap (eV)": f"{e_band_gap_rectified}{true_neg_bg}",
                "valid_gap": is_good_bg(bgdict["E_band_gap"])
        }

        if len(bgdict) > 1:
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

    JsonWriter(args["summaryfile"], screening_results, indent=4).write_as_list()

    stop = datetime.now()
    print(f"elapsed time: {stop-start}")


if __name__ == "__main__":
    main()
