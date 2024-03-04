#!/usr/bin/python
"""
A script to determine material fundamental band gap from VASP energies and Δ-Sol method.

Reference for Δ-Sol method:
    M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010)
    (reference 32 in screening_pipeline/Bibliography)
"""

########################################
# SYSTEM I/O MODULES

import os
import json
from typing import Dict
from pathlib import Path
from datetime import datetime
from argparse import ArgumentParser, Namespace, RawTextHelpFormatter

########################################
# PYTHON MATERIALS GENOMICS PACKAGE


########################################
# LOCAL MODULES

from screening_pipeline.utils.utils import _yaml_loader
from screening_pipeline.utils.typing import PMGStaticSet
from screening_pipeline.utils.vasp_io import vasp_batch_launch, vasp_launcher, batch_extract_vasp_data, \
                                             delta_sol_inputs_init
from screening_pipeline.utils.data_process import batch_calculate_delta_sol_band_gaps

########################################
# LOCAL FUNCTIONS

def assert_args(args: Namespace) -> None:

    assert os.path.isdir(args.input_dir), f"{args.input_dir}: No such directory found."

    if args.mini_maxi is not None:
        assert all([value >= 0.0 for value in args.mini_maxi]), (
        f"Acceptable band gap values must be positive or zero."
        )

        assert args.mini_maxi[0] != args.mini_maxi[1], (
        f"Acceptable band gap values cannot have the same value."
        )

    if args.accept is not None:
        assert args.accept.endswith(".txt"), (
        "Structure acceptance file must be a plain text file type (.txt)."
        )

    if args.ignore is not None:
        assert args.ignore.endswith(".txt"), (
        "Structure ignoring file must be a plain text file type (.txt)."
        )

    assert args.summary.endswith(".json"), (
        "'summary' argument value must be a JSON format file name (.json)."
    )

    assert args.workers >= 1, \
    "The number of workers cannot be negative or zero."

########################################
# MAIN FUNCTION

def main():
    start = datetime.now()

    # ARGUMENTS PARSING BLOCK
    
    prog_name = "band_gap_screening"
    prog_desc = """
        A script to determine material fundamental band gap from VASP energies and Δ-Sol method.

        Reference for Δ-Sol method:
            M.K.Y. Chan and G. Ceder, Phys. Rev. Lett., 105, 196403 (2010)
            (reference 32 in screening_pipeline/Bibliography)
        """
    prog_missing_steps = """
        Missing steps to complete this script: None
        """
    helper_format = RawTextHelpFormatter

    parser = ArgumentParser(
        prog=prog_name, 
        description=prog_desc, 
        #epilog=prog_missing_steps, 
        formatter_class=helper_format
    )

    parser.add_argument(
        "input_dir",
        type=str,
        default=None,
        help="Base directory containing static calculations results directories.",
    )
    parser.add_argument(
        "-v", "--valid-interval",
        narg=2,
        type=float,
        default=None,
        help=(
            "Valid band gaps interval in eV (default: [1.3 ; 3.6] eV)."
            "If another interval is given, the min AND max values must be given, even if one "
            "of them matches the default values."
        ), 
        metavar="float"
    )
    parser.add_argument(
        "-a", "--accept", 
        type=str, 
        default=None, 
        help="""Defines a file whose presence in a structure directory means it passed this
        screening step successfully and can be kept for further calculations.""",  
        metavar="accept_file.txt"
    )
    parser.add_argument(
        "-i", "--ignore", 
        type=str, 
        default=None, 
        help="""Defines a file whose presence in a structure directory means it did not pass
        previous screening steps and should not be used in this calculation. 
        This file will also be written in structure directories that did not pass this step.
        WARNING: 
        If the name of this file is overwritten, care must be taken that it is the same file
        throughout every used screening steps to make sure rejected structures don"t go further.""", 
        metavar="ignore_file.txt"
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
            "Also usable by further steps to filter out structures that were rejected in this step.\n"
            "This arg only changes the file name, its path is automatically in the step directory."
        ),
    )
    parser.add_argument(
        "--with-uncertainties", 
        action="store_true", 
        help="Pass this flag to enable computation of minimal and maximal Δ-Sol band gaps.\n\
              If enabled, it will search for uncertainty calculations results."
    )

    args: Namespace = parser.parse_args()

    assert_args(args)

    input_dir = args.input_dir

    if args.mini_maxi is None:
        valid_interval = (1.3, 3.6)
    else:
        valid_interval = sorted(args.mini_maxi)

    accept_file    = args.accept or None
    ignore_file    = args.ignore or None
    workers        = args.workers
    summary        = os.path.join(input_dir, args.summary)


    # MAIN BLOCK

    # Extract VASP static calculations results
    bg_data = batch_extract_vasp_data(
        method="delta_sol", 
        base_dir=input_dir, 
        workers=workers
    )

    E_band_gaps = batch_calculate_delta_sol_band_gaps(
        bg_data, args.with_uncertainties, workers
    )

    good_bg_structs = list(filter(
        lambda tup: min(valid_interval) <= tup[1]["E_band_gap"] <= max(valid_interval), 
        list(E_band_gaps.items())
    ))

    bad_bg_structs = list(filter(
        lambda tup: tup[1]["E_band_gap"] < min(valid_interval) or tup[1]["E_band_gap"] > max(valid_interval), 
        list(E_band_gaps.items())
    ))

    screening_results = []

    # Keep good structures
    for struct in good_bg_structs:
        name       = struct[0]
        bgdict     = struct[1]
        E_band_gap = max(round(bgdict["E_band_gap"], 6), 0.0)

        if accept_file is not None:
            accept_msg = [
                "BAND GAP TEST PASSED", 
                f"Δ-Sol band gap was estimated to {E_band_gap} eV, which is inside the interval [{min(valid_interval)}, {max(valid_interval)}].", 
                "Therefore, it is suitable for wanted application, and should be considered for further screening steps."
            ]

            if args.with_uncertainties:
                E_band_gap_min = max(round(bgdict["E_band_gap_min"], 6), 0.0)
                E_band_gap_max = max(round(bgdict["E_band_gap_max"], 6), 0.0)
                accept_msg += [
                    "\nUncertainty interval (does not affect acception or rejection):", 
                    f"Band Gap minimum = {E_band_gap_max} eV", 
                    f"Band Gap maximum = {E_band_gap_min} eV"
                ]

            accept_msg = "\n".join(accept_msg)

            for calc_dir in filter(lambda path: path.name.startswith(name), Path(input_dir).iterdir()):
                accept_path = os.path.join(calc_dir, accept_file)
                with open(accept_path, "wt") as fp:
                    fp.write(accept_msg)
        
        struct_dict = {
                "path": os.path.join(str(input_dir), name),
                "bandgap": E_band_gap,
                "valid_gap": True
        }

        if args.with_uncertainties:
            E_band_gap_min = max(round(bgdict["E_band_gap_min"], 6), 0.0)
            E_band_gap_max = max(round(bgdict["E_band_gap_max"], 6), 0.0)
            struct_dict.update(
                {
                    "bandgap_min": E_band_gap_max,
                    "bandgap_max": E_band_gap_min
                }
            )

        screening_results.append(struct_dict)

    # Reject unsuitable structures
    for struct in bad_bg_structs:
        name       = struct[0]
        bgdict     = struct[1]
        E_band_gap = max(round(bgdict["E_band_gap"], 6), 0.0)

        if ignore_file is not None:
            reject_msg = [
                "BAND GAP REJECTION", 
                f"Δ-Sol band gap was estimated to {E_band_gap} eV, which is not inside the interval [{min(valid_interval)}, {max(valid_interval)}].",
                "Therefore, it is not suitable for wanted application, and should not be considered in further screening steps."
            ]

            if args.with_uncertainties:
                E_band_gap_min = max(round(bgdict["E_band_gap_min"], 6), 0.0)
                E_band_gap_max = max(round(bgdict["E_band_gap_max"], 6), 0.0)
                reject_msg += [
                    "\nUncertainty interval (does not affect acception or rejection):", 
                    f"Band Gap minimum = {E_band_gap_max} eV", 
                    f"Band Gap maximum = {E_band_gap_min} eV"
                ]

            reject_msg = "\n".join(reject_msg)
            for calc_dir in filter(lambda path: path.name.startswith(name), Path(input_dir).iterdir()):
                reject_path = os.path.join(calc_dir, ignore_file)
                with open(reject_path, "wt") as fp:
                    fp.write(reject_msg)

        struct_dict = {
                "path": os.path.join(str(input_dir), name),
                "bandgap": E_band_gap,
                "valid_gap": False
        }

        if args.with_uncertainties:
            E_band_gap_min = max(round(bgdict["E_band_gap_min"], 6), 0.0)
            E_band_gap_max = max(round(bgdict["E_band_gap_max"], 6), 0.0)
            struct_dict.update(
                {
                    "bandgap_min": E_band_gap_max,
                    "bandgap_max": E_band_gap_min
                }
            )

        screening_results.append(struct_dict)

    def sort_by_path(dct: Dict) -> str:
        return dct.get("path")

    screening_results = sorted(screening_results, key=sort_by_path)

    with open(summary, "w") as fp:
        json.dump(screening_results, fp, indent=4)

    stop = datetime.now()
    print(f"elapsed time: {stop-start}")


if __name__ == "__main__":
    main()
