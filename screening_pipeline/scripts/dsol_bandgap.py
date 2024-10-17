#!/usr/bin/python
"""
A script to determine material fundamental band gap from VASP energies and Δ-Sol method.

Reference for Δ-Sol method:
    M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010).
"""

import os
import json
from typing import Dict
from datetime import datetime
from argparse import ArgumentParser, Namespace

# LOCAL IMPORTS
from screening_pipeline.utils import (
    check_file_format, check_file_or_dir,
    batch_extract_vasp_data, batch_get_dsol_band_gaps
)


########################################


def _assert_args(args: Namespace) -> None:

    check_file_or_dir(args.input_dir, "dir")

    assert args.functional in {"LDA", "PBE", "AM05"}, (
        f"{args.functional} is not supported by delta-Sol. "
        "See --help for valid functional argument values."
    )

    if args.valid_interval is not None:
        assert all(value >= 0.0 for value in args.valid_interval), (
        "Acceptable band gap values must be positive or zero."
        )

        assert args.valid_interval[0] != args.valid_interval[1], (
        "Acceptable band gap values cannot have the same value."
        )

    check_file_format(args.summary, allowed_formats="json")

    assert args.workers >= 1, (
        "The number of workers cannot be negative or zero."
    )


########################################


def main() -> None:
    """Main function."""
    start = datetime.now()

    # ARGUMENTS PARSING BLOCK

    prog_name = "dsol_bandgap"
    prog_desc = """
        A script to determine material fundamental band gap from VASP energies and Δ-Sol method.

        Reference for Δ-Sol method:
            M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010)
        """

    parser = ArgumentParser(
        prog=prog_name,
        description=prog_desc,
    )

    parser.add_argument(
        "input_dir",
        type=str,
        default=None,
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
        default=None,
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
        default=1,
        help="Number of parallel processes to spawn for parallelized steps.",
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
        metavar="summary.json"
    )
    parser.add_argument(
        "--with-uncertainties", 
        action="store_true",
        help=(
            "Enables computation of minimal and maximal Δ-Sol band gaps.\n"
            "If enabled, it will search for uncertainty calculations results."
        )
    )

    args: Namespace = parser.parse_args()

    _assert_args(args)

    input_dir = args.input_dir
    dft_func  = args.functional

    if args.valid_interval is None:
        valid_interval = (1.3, 3.6)
    else:
        valid_interval = sorted(args.valid_interval)

    workers = args.workers
    summary = os.path.join(input_dir, args.summary)


    # MAIN BLOCK

    # Extract VASP static calculations results
    bg_data = batch_extract_vasp_data(
        method="delta_sol_calc",
        base_dir=input_dir,
        workers=workers
    )

    e_band_gaps = batch_get_dsol_band_gaps(
        bg_data, dft_func, args.with_uncertainties, workers
    )

    good_bg_structs = list(filter(
        lambda tup: min(valid_interval) <= tup[1]["E_band_gap"] <= max(valid_interval),
        list(e_band_gaps.items())
    ))

    bad_bg_structs = list(filter(
        lambda tup: not min(valid_interval) <= tup[1]["E_band_gap"] <= max(valid_interval),
        list(e_band_gaps.items())
    ))

    screening_results = []

    # Keep good structures
    for struct in good_bg_structs:
        name       = struct[0]
        bgdict     = struct[1]
        e_band_gap = round(bgdict["E_band_gap"], 6)
        e_band_gap_rectified = max(e_band_gap, 0.0)
        true_neg_bg = f" (true measurement: {e_band_gap})" if e_band_gap_rectified == 0.0 else ""

        struct_dict = {
                "path": os.path.join(str(input_dir), name),
                "bandgap (eV)": f"{e_band_gap_rectified}{true_neg_bg}",
                "valid_gap": True
        }

        if args.with_uncertainties:
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
                "path": os.path.join(str(input_dir), name),
                "bandgap (eV)": f"{e_band_gap_rectified}{true_neg_bg}",
                "valid_gap": False
        }

        if args.with_uncertainties:
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

    def sort_by_path(dct: Dict) -> str:
        return dct.get("path")

    screening_results = sorted(screening_results, key=sort_by_path)

    with open(summary, "w", encoding="utf-8") as fp:
        json.dump(screening_results, fp, indent=4)

    stop = datetime.now()
    print(f"elapsed time: {stop-start}")


if __name__ == "__main__":
    main()
