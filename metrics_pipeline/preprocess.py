#!/usr/bin/python
"""
Preprocess generated structures data for further pipeline steps:
- Convert POSCAR and CIF input data into standardized CIF for format compatibility with the pipeline;
- Generate a unique indexed name for each structure for data tracking throughout computations;
- (optional) Eliminate data containing specific unwanted elements;
- (optional) Compute spacegroup symmetry and add symmetry infos into standardized output data.
"""

import os
from pathlib import Path
import typing as tp
from time import time
import argparse as ap

from src.utils import (
    parse_input_args, check_type, check_num_value,
    ALL_ELTS_CATEGORIES,
    get_elts_from_symbol_or_z, get_elts_in_categories,
    filter_by_elements,
    generate_genmat_names
)
from src.io import check_file_or_dir, check_file_format, CIFFile, PoscarFile

if tp.TYPE_CHECKING:
    from pymatgen.core import Structure


def _get_command_line_args() -> ap.Namespace:
    """Command Line Interface (CLI)."""
    parser = ap.ArgumentParser(
        prog=os.path.basename(__file__), description=__doc__,
        formatter_class=ap.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "input_file",
        help=(
            "Path to the CIF or POSCAR formatted file containing structure data to process. "
            "If POSCAR format, it must have the extension '.poscar', and each structure must "
            "have a header comment line."
        )
    )
    parser.add_argument(
        "-o", "--output",
        help=(
            "Path to wanted CIF output file where processed structures "
            "will be stored. If not given, output file is written in input file "
            "directory with the same name but a '_preproc' suffix is added."
        )
    )
    parser.add_argument(
        "--special-keys", "-sk", nargs="*",
        help=(
            "Each structure will be attributed a unique name in output file for easy data "
            "tracking throughout the pipeline, composed of a unique index followed by the reduced "
            "formula of the structure. Pass 'header' here to keep original header line content "
            "(comment line for POSCARs) in '_original_header' label in output CIFs. "
            "For a CIF formatted file, also pass any non-looping CIF Label to keep in output data."
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
        "-s", "--symmetrize", action="store_true",
        help="Compute and add spacegroup symmetry data in output file."
    )
    parser.add_argument(
        "--symprec", type=float, default=0.1, metavar="float",
        help="Fractional coordinates tolerance for symmetry finding (Default: %(default)s).",
    )
    parser.add_argument(
        "--angleprec", type=float, default=5.0, metavar="float",
        help="Angle tolerance for symmetry finding in degrees (Default: %(default)s degrees).",
    )
    parser.add_argument(
        "-r", "--remove-elts", nargs="*",
        help=(
            "Structures containing given elements will not be saved in output file. "
            "Elements to track can be specified by their symbol, atomic number, or a mix of both."
        )
    )
    parser.add_argument(
        "-R", "--remove-elt-categories", nargs="*",
        help=(
            "Pass valid element categories to eliminate structures containing any element from "
            f"these categories. Supported categories are: {', '.join(sorted(ALL_ELTS_CATEGORIES))}."
        )
    )
    args: ap.Namespace = parser.parse_args()
    return args


def _process_input_args(args_dict: dict[str, tp.Any]) -> dict[str, tp.Any]:
    """Handle input arguments assertions and processing."""
    if args_dict is None:
        raise ValueError(f"No arguments found at '{os.path.basename(__file__)}' script call.")
    
    check_type(args_dict, "args_dict", (dict,))

    # Set default values for unset optional arguments
    args_dict = {k: v for k, v in args_dict.items() if v is not None}

    check_file_or_dir(args_dict.get("input_file"), "file", allowed_formats=("cif", "poscar"))
    args_dict["input_file"] = str(Path(args_dict["input_file"]).resolve(strict=True))

    default_output = os.path.splitext(args_dict["input_file"])[0] + "_preproc.cif"
    args_dict.setdefault("output", default_output)
    args_dict["output"] = str(Path(args_dict["output"]).resolve())

    args_dict.setdefault("special_keys", [])
    args_dict.setdefault("symmetrize", False)
    args_dict.setdefault("symprec", 0.1)
    args_dict.setdefault("angleprec", 5.0)
    args_dict.setdefault("remove_elts", [])
    args_dict.setdefault("remove_elt_categories", [])

    # Assert set arguments conformity
    check_file_format(args_dict["output"], allowed_formats="cif")

    if args_dict.get("special_keys") is not None:
        check_type(args_dict["special_keys"], "special_keys", (list,))
        for idx, elt in enumerate(args_dict["special_keys"]):
            check_type(elt, f"special_keys[{idx}]", (str,))

    if args_dict.get("workers") is not None:
        check_type(args_dict["workers"], "workers", (int,))
        check_num_value(args_dict["workers"], "workers", ">=", 0)

    check_type(args_dict.get("symmetrize", False), "symmetrize", (bool,))
    check_type(args_dict["symprec"], "symprec", (float,))
    check_num_value(args_dict["symprec"], "symprec", ">=", 0.0)
    check_num_value(args_dict["symprec"], "symprec", "<=", 1.0)
    check_type(args_dict["angleprec"], "angleprec", (float,))
    check_num_value(args_dict["angleprec"], "angleprec", ">=", 0.0)
    check_num_value(args_dict["angleprec"], "angleprec", "<=", 90.0)
    check_type(args_dict["remove_elts"], "remove_elts", (list,))
    for idx, elt in enumerate(args_dict["remove_elts"]):
        check_type(elt, f"remove_elts[{idx}]", (str,))

    check_type(args_dict["remove_elt_categories"], "remove_elt_categories", (list,))
    for idx, elt in enumerate(args_dict["remove_elt_categories"]):
        check_type(elt, f"remove_elt_categories[{idx}]", (str,))
        if not elt in ALL_ELTS_CATEGORIES:
            raise ValueError(f"The category {elt!r} is not supported.")

    # Parse forbidden elements
    forbidden_elts = get_elts_in_categories(args_dict["remove_elt_categories"])
    forbidden_elts.update(get_elts_from_symbol_or_z(args_dict["remove_elts"]))
    args_dict["forbidden_elts"] = forbidden_elts

    return args_dict


def _print_config(args_dict: dict) -> None:
    print(" ")
    print("------------------------------")
    print(" ")
    print(" - I/O ARGUMENTS - ")
    print(f"INPUT FILE: {args_dict.get('input_file')}")
    print(f"OUTPUT FILE: {args_dict.get('output')}")
    print(f"SPECIAL KEYS: {'None' if args_dict.get('special_keys') is None else ''}")
    if args_dict.get("special_keys") is not None:
        keys_lst = '\n'.join(args_dict["special_keys"])
        print(f"{keys_lst}")
    print(" ")
    print("------------------------------")
    print(" ")
    print(" - PROCESSING ARGUMENTS - ")
    print(
        "NUMBER OF WORKERS: "
        f"{'auto' if args_dict.get('workers') is None else args_dict.get('workers')} "
        f"{'(sequential)' if args_dict.get('workers') == 0 else ''}"
    )
    print(f"SYMMETRIZE: {args_dict.get('symmetrize', False)}")
    print(
        f"* Fractional coordinates tolerance: {args_dict.get('symprec')} "
        f"{'(ignored)' if not args_dict.get('symmetrize') else ''}"
    )
    print(
        f"* Angles tolerance: {args_dict.get('angleprec')} degrees "
        f"{'(ignored)' if not args_dict.get('symmetrize') else ''}"
    )
    print(
        "REMOVED ELEMENTS: "
        f"{', '.join(args_dict['remove_elts']) if args_dict['remove_elts'] != [] else None}."
    )
    print(
        "REMOVED CATEGORIES: "
        f"{', '.join(args_dict['remove_elt_categories']) if args_dict['remove_elt_categories'] != [] else None}."
    )
    print(
        "ALL REMOVED ELEMENTS (parsed categories): "
        f"{', '.join(args_dict['forbidden_elts'].keys()) if args_dict['forbidden_elts'] != {} else None}."
    )
    print(" ")
    print("------------------------------")
    print(" ")


def _get_structures_from_cif(args_dict: dict[str, tp.Any]) -> tuple[list[Structure], int, int, int]:
    cif_file = CIFFile.from_file(
        args_dict["input_file"],
        special_keys=args_dict.get("special_keys"),
        workers=args_dict.get("workers")
    )
    cifs = cif_file.get_cifs()
    nbr_loaded_data = len(cifs)
    assert nbr_loaded_data > 0, "No structure could be parsed from given data"
    print(f"{nbr_loaded_data} structures detected in total")

    # Remove structures containing unwanted elements
    if len(args_dict["forbidden_elts"]) > 0:
        cifs, nbr_discarded = filter_by_elements(
            cifs, list(args_dict["forbidden_elts"].keys()), format="cif"
        )
        print(f"{nbr_discarded} CIFs containing forbidden elements were discarded.")

    cif_file.clear()
    cif_file.add_cifs(cifs)
    structures, invalid_indices = cif_file.parse_structures()
    nbr_invalid = len(invalid_indices)
    print(f"{nbr_invalid} CIFs could not be loaded and will not be saved in output.")
    return structures, nbr_loaded_data, nbr_discarded, nbr_invalid


def _get_structures_from_poscar(args_dict: dict[str, tp.Any]) -> tuple[list[Structure], int, int]:
    pfile = PoscarFile(args_dict["input_file"], workers=args_dict.get("workers"))
    nbr_loaded_data = len(pfile)
    assert nbr_loaded_data > 0, "No structure could be parsed from given data"
    print(f"{nbr_loaded_data} structures detected in total")

    # Remove structures containing unwanted elements
    if len(args_dict["forbidden_elts"]) > 0:
        poscars, nbr_discarded = filter_by_elements(
            pfile.data, list(args_dict["forbidden_elts"].keys()), format="poscar"
        )
        pfile = PoscarFile.from_str("\n".join(poscars), strict=False, workers=args_dict.get("workers"))
        print(f"{nbr_discarded} POSCARs containing forbidden elements were discarded.")

    structures = pfile.parse_structures()
    if args_dict["special_keys"] and "header" in args_dict["special_keys"]:
        for header, structure in zip(pfile.headers, structures):
            structure.properties["header"] = header
    return structures, nbr_loaded_data, nbr_discarded


def _print_computation_summary(args_dict: dict[str, tp.Any], counters: dict[str, int | float]) -> None:
    nbr_loaded_data = counters.pop("nbr_loaded_data", 0)
    nbr_written_data = counters.pop("nbr_written_data", 0)
    nbr_discarded = counters.pop("nbr_discarded", 0)
    nbr_invalid = counters.pop("nbr_invalid", 0)
    computation_time = counters.pop("computation_time", 0)

    print("\n------------------------------")
    print("\nSUMMARY OF THE PREPROCESSING")
    print(f"{nbr_loaded_data} structures processed in total, including:")
    print(f"- {nbr_written_data} structure(s) written in output")

    if len(args_dict["forbidden_elts"]) > 0:
        print(f"- {nbr_discarded} structure(s) containing forbidden elements")
    
    if args_dict["input_file"].endswith("cif") and nbr_invalid > 0:
        print(f"- {nbr_invalid} structure(s) that could not be parsed")

    print(f"\nOutput results written in {args_dict['output']}")
    print(f"elapsed time: {computation_time}")


def main(standalone: bool = True, **kwargs) -> None:
    """
    Preprocess generated structures data for further pipeline steps:
    - Convert POSCAR and CIF input data into standardized CIF for format compatibility with the pipeline;
    - Generate a unique indexed name for each structure for data tracking throughout computations;
    - (optional) Eliminate data containing specific unwanted elements;
    - (optional) Compute spacegroup symmetry and add symmetry infos into standardized output data.

    Parameters
    ----------
    standalone: bool
        Whether parsed script is used directly through command-line (stand-alone script)
        or in an external pipeline script.

    input_file: str | Path
        Path to the CIF or POSCAR file containing structure data to process.
        If POSCAR format, it must have extension '.poscar'.

    output: str | Path
        Path to wanted CIF output file where processed structures will be stored. If not given,
        output file is written in input file directory with the same name with a '_out' suffix
        added.

    special_keys: list[str], optional
        Each structure will be attributed a unique name in output file for easy data
        tracking throughout the pipeline, composed of a unique index followed by the reduced
        formula of the structure. Pass 'header' here to keep original header line content
        (comment line for POSCARs) in '_original_header' label in output CIFs.
        For a CIF formatted file, also pass any non-looping CIF Label to keep in output data.

    workers: int, optional
        Number of processes to use in parallel. If not given, will use default of
        `tqdm.contrib.concurrent.process_map()`. Pass 0 to disable `process_map()`
        and execute sequentially.

    symmetrize: bool
        Compute and add spacegroup symmetry data in output file.

    symprec: float
        Fractional coordinates tolerance for symmetry finding. Defaults to 0.01.

    angleprec: float
        Angle tolerance for symmetry finding in degrees. Defaults to 5.0 degrees.

    remove_elts: list[str], optional
        Structures containing given elements will not be saved in output file.
        Elements to track can be specified by their symbol, atomic number, or a mix of both.

    remove_elt_categories: list[str], optional
        Pass valid element categories to eliminate structures containing any element from
        these categories.
    """
    start = time()
    args = parse_input_args(_get_command_line_args, _process_input_args, standalone, **kwargs)
    _print_config(args)

    # Data extraction and conversion to structures
    if args["input_file"].endswith(".cif"):
        structures, nbr_loaded_data, nbr_discarded, nbr_invalid = _get_structures_from_cif(args)

    elif args["input_file"].endswith(".poscar"):
        structures, nbr_loaded_data, nbr_discarded = _get_structures_from_poscar(args)
        nbr_invalid = 0

    print(f"{len(structures)} structures with valid composition were successfully parsed.")

    # Generate unique indexed names as headers
    structures = generate_genmat_names(structures)

    for key in ("header", "_original_header"):
        if key not in args["special_keys"]:
            args["special_keys"].append(key)

    # Symmetrize and write standardized CIF file
    if args["symmetrize"]:
        print("Symmetrization is activated, now symmetrizing and writing CIF file...")
    else:
        print("Symmetrization is deactivated, now writing CIF file...")

    cif_file = CIFFile(special_keys=args["special_keys"], workers=args.get("workers"))
    cif_file.add_structures(
        structures,
        symmetrize=args["symmetrize"],
        symprec=args["symprec"],
        angleprec=args["angleprec"]
    )
    cif_file.write_file(args["output"])

    # Time of the preprocessing
    stop = time()
    counters: dict[str, int | float] = {
        "nbr_loaded_data": nbr_loaded_data,
        "nbr_written_data": len(structures),
        "nbr_discarded": nbr_discarded,
        "nbr_invalid": nbr_invalid,
        "computation_time": stop-start,
    }
    _print_computation_summary(args, counters)


if __name__ == "__main__":
    main()
