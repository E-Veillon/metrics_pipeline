"""Convert CIF formatted structures to POSCAR format."""

import os
import typing as tp
import argparse as argp
from pymatgen.io.vasp import Poscar

from . import _parse_input_args
from .utils import (
    check_type, check_file_or_dir, check_file_format, check_num_value,
    read_cif, VisualIterator
)

########################################
# ARGUMENTS HANDLING

def _get_command_line_args() -> argp.Namespace:
    """Command Line Interface (CLI)."""
    parser = argp.ArgumentParser(prog=os.path.basename(__file__), description=__doc__)
    parser.add_argument(
        "input_file", help="Path to the CIF file containing structure data to process."
    )
    parser.add_argument(
        "--output", "-o",
        help=(
            "Path to the file to write POSCAR files in. If not given, "
            "output file is written in input file directory with the "
            "same name but with a '.poscar' extension instead of '.cif'."
        )
    )
    parser.add_argument(
        "--save-header", "-sh",
        action="store_true",
        help=(
            "Save the header of CIF data into the comment line of POSCAR format. "
            "If not passed, the pymatgen default is applied, i.e. writing full "
            "composition of the crystal in the comment line."
        )
    )
    parser.add_argument(
        "--significant-figures", "-sf",
        type=int,
        help="Number of significant digits to output all quantities. Defaults to 16."
    )
    parser.add_argument(
        "--workers", "-w",
        type=int,
        help=(
            "Number of parallel processes to spawn for parallelized steps. "
            "If not given, If not given, default value is the 'max_workers' "
            "default value from tqdm.contrib.concurrent.process_map function."
        ),
        metavar="int",
    )
    parser.add_argument(
        "--sequential", "-S",
        action="store_true",
        help=(
            "Pass this flag to deactivate multiprocessing handling and switch to "
            "sequential computation. Generally, multiprocess is faster, but sometimes "
            "(e.g. when data is very big) it can softlock into resource distribution. "
            "Switch to more stable sequential computing if such problem  were to arise."
        )
    )
    args: argp.Namespace = parser.parse_args()
    return args


def _process_input_args(args_dict: dict[str, tp.Any]) -> dict[str, tp.Any]:
    """Handle input arguments assertions and processing."""
    if args_dict is None:
        raise ValueError(f"No arguments found at '{os.path.basename(__file__)}' script call.")
    
    check_type(args_dict, "args_dict", (dict,))

    # Set default values for unset optional arguments
    args_dict = {k: v for k, v in args_dict.items() if v is not None}
    default_output = str(args_dict.get("input_file")).replace(".cif", ".poscar")
    args_dict.setdefault("output", default_output)
    args_dict.setdefault("save_header", False)
    args_dict.setdefault("significant_figures", 16)
    args_dict.setdefault("sequential", False)

    # Assert set arguments conformity
    check_file_or_dir(args_dict.get("input_file"), "file", allowed_formats="cif")
    check_file_format(args_dict.get("output"), allowed_formats="poscar")
    check_type(args_dict.get("save_header"), "save_header", (bool,))
    check_type(args_dict.get("significant_figures"), "significant_figures", (int,))
    check_num_value(args_dict.get("significant_figures"), "significant_figures", ">=", 4)
    check_num_value(args_dict.get("significant_figures"), "significant_figures", "<=", 20)

    args_dict["input_file"] = os.path.abspath(os.path.realpath(args_dict["input_file"]))
    args_dict["output"] = os.path.abspath(os.path.realpath(args_dict["output"]))

    if args_dict.get("workers") is not None and not args_dict.get("sequential", False):
        check_type(args_dict.get("workers"), "workers", (int,))
        check_num_value(args_dict.get("workers"), "workers", ">", 0)

    # Print final configuration
    print(" ")
    print("------------------------------")
    print(" ")
    print(" - I/O ARGUMENTS - ")
    print(f"INPUT FILE: {args_dict.get('input_file')}")
    print(f"OUTPUT FILE: {args_dict.get('output')}")
    print(" ")
    print("------------------------------")
    print(" ")
    print(" - ACTIVATED FEATURES - ")
    print(f"SAVE HEADER: {args_dict.get('save_header')}")
    print(f"SIGNIFICANT FIGURES: {args_dict.get('significant_figures')}")
    print(f"IS SEQUENTIAL: {args_dict.get('sequential')}")
    print(
        "NUMBER OF WORKERS: "
        f"{'auto' if args_dict.get('workers') is None else args_dict.get('workers')} "
        f"{'(ignored)' if args_dict.get('sequential') else ''}"
    )
    print(" ")
    print("------------------------------")
    print(" ")

    return args_dict


########################################


def main(standalone: bool = True, **kwargs) -> None:
    """
    Convert CIF formatted structures to POSCAR format.

    Args:
        standalone (bool):          Whether parsed script is used directly through
                                    command-line (stand-alone script) or in an external
                                    pipeline script.

        input_file (str|Path):      Path to the CIF file containing structure data to process.

        output (str|Path):          Path to the file to write POSCAR files in. If not given,
                                    output file is written in input file directory with the
                                    same name but with a '.poscar' extension instead of '.cif'.

        save_header (bool):         Save the header of CIF data into the comment line of POSCAR
                                    format. If not passed, the pymatgen default is applied, i.e.
                                    writing full composition of the crystal in the comment line.

        significant_figures (int):  Number of significant digits to output all quantities.
                                    Defaults to 16.

        workers (int):              Number of parallel processes to spawn for parallelized steps.
                                    If not given, If not given, default value is the 'max_workers'
                                    default value from tqdm.contrib.concurrent.process_map function.

        sequential (bool):          Pass this flag to deactivate multiprocessing handling and
                                    switch to sequential computation. Generally, multiprocess
                                    is faster, but sometimes (e.g. when data is very big) it
                                    can softlock into resource distribution. Switch to more
                                    stable sequential computing if such problem  were to arise.
    """
    args = _parse_input_args(_get_command_line_args, _process_input_args, standalone, **kwargs)

    structures, *_ = read_cif(
        filename=args["input_file"],
        keep_rare_gases=True,
        keep_rare_earths=True,
        special_keys=(["header"] if args.get("save_header") else None),
        workers=args.get("workers"),
        sequential=args["sequential"]
    )
    structures = VisualIterator(structures, desc="Converting to Poscar")
    outfile_content = [
        Poscar(
            structure=struct,
            comment=f"# {struct.properties['header']}" if args.get("save_header") else None
        ).get_str(significant_figures=args["significant_figures"]) for struct in structures
    ]
    with open(args["output"], "wt", encoding="utf-8") as fp:
        fp.write("".join(outfile_content)) # '\n' separating POSCAR data are given by get_str().


if __name__ == "__main__":
    main()