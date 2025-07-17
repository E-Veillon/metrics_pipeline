#!/usr/bin/python
"""
A script that writes and performs VASP static run on structure data read from CIF files,
using pymatgen as a setting interface between raw data and VASP.
Energy results may be used for phase diagrams computations.
"""

import os
import typing as typ
from datetime import datetime
import argparse as argp

# PYTHON MATERIAL GENOMICS
from pymatgen.core.structure import SiteCollection

# LOCAL IMPORTS
from . import CONFIGPATH, _parse_input_args
from .utils import (
    check_type, check_num_value, check_file_or_dir, add_new_dir, yaml_loader,
    PMGStaticSet, read_cif, vasp_static_settings, write_and_run_vasp
)


########################################
# ARGUMENTS HANDLING

def _get_command_line_args() -> argp.Namespace:
    """Command Line Interface (CLI)."""
    parser = argp.ArgumentParser(prog=os.path.basename(__file__), description=__doc__)
    parser.add_argument(
        "input_file",
        type=str,
        help="Path to the CIF file containing structure data to read.",
        metavar="input_file.cif"
    )
    parser.add_argument(
        "executable_path",
        type=str,
        help="Path to the VASP executable.",
        metavar="PATH"
    )
    parser.add_argument(
        "-o", "--output",
        help=(
            "Path to the output directory where VASP files will be written. "
            "A subdirectory will be created in output directory for each structure "
            "found in input_file. By default, it will create a 'Statics' directory "
            "in the input_file directory and write structures runs in it."
        ),
        metavar="outdir"
    )
    parser.add_argument(
        "-p", "--preset",
        type=PMGStaticSet,
        default=PMGStaticSet.MPSTATICSET.value,
        help=(
            "The pymatgen preset to use for VASP static run. "
            "More info on possible presets in pymatgen documentation:\n"
            "https://pymatgen.org/pymatgen.io.vasp.html#pymatgen.io.vasp.sets."
        ),
        metavar="StaticSet"
    )
    parser.add_argument(
        "-u",
        "--user-settings",
        type=str,
        default="default_settings.yaml",
        help=(
            "Name of the .yaml file containing tags overrides to put over the PMG preset.\n"
            "The file must be at location metrics_pipeline/config to be found."
        ),
        metavar="FILENAME",
        dest="user_settings"
    )
    parser.add_argument(
        "-w", "--workers",
        type=int,
        help=(
            "Number of parallel processes to spawn for parallelized steps. "
            "If not given, default value is the 'max_workers' default value "
            "from tqdm.contrib.concurrent.process_map function."
        )
    )
    parser.add_argument(
        "-t",
        "--task_index",
        type=int,
        default=0,
        help=(
            "If a job array is used, provide here the structure index "
            "to treat according to task IDs (e.g. if task ID 0 treats "
            "structure 0 and so on, just provide the task ID).\n"
            "If not given, it will default to the first possible index, i.e. index 0."
        ),
    )
    args: argp.Namespace = parser.parse_args()
    return args


def _process_input_args(args_dict: dict[str, typ.Any]) -> dict[str, typ.Any]:
    """Handle input arguments assertions and processing."""
    if args_dict is None:
        raise ValueError(f"No arguments found at '{os.path.basename(__file__)}' script call.")
    
    check_type(args_dict, "args_dict", (dict,))

    # Set default values for unset optional arguments
    args_dict = {k: v for k, v in args_dict.items() if v is not None}
    args_dict.setdefault("executable_path", "vasp")
    # WARNING:
    # Usual VASP shortcut, but may activate wrong VASP version if several are installed.
    # Prefer giving a true VASP executable path for unambiguous computation.
    default_output = os.path.join(os.path.dirname(args_dict.get("input_file", "")), "Statics")
    args_dict.setdefault("output", default_output)
    args_dict.setdefault("preset", PMGStaticSet.MPSTATICSET.value)
    args_dict.setdefault("user_settings", "default_settings.yaml")
    args_dict.setdefault("task_index", 0)

    # Assert set arguments conformity
    check_file_or_dir(args_dict.get("input_file"), "file", allowed_formats="cif")

    if not str(args_dict.get("executable_path")).startswith("vasp"):
        check_file_or_dir(args_dict.get("executable_path"), "file")

    assert args_dict.get("preset",PMGStaticSet.MPSTATICSET.value) in PMGStaticSet.values, (
    "Provided static preset must be one of the following:\n"
    f"{PMGStaticSet.values}."
    )
    config_path = os.path.join(CONFIGPATH, args_dict["user_settings"])
    check_file_or_dir(config_path, "file", allowed_formats=("yml", "yaml"))
    
    if args_dict.get("workers") is not None:
        check_type(args_dict.get("workers"), "workers", (int,))
        check_num_value(args_dict.get("workers"), "workers", ">", 0)
    
    check_type(args_dict.get("task_index"), "task_index", (int,))
    check_num_value(args_dict.get("task_index"), "task_index", ">=", 0)

    # Additional arguments processing
    os.makedirs(args_dict["output"], exist_ok=True)
    args_dict["preset"] = PMGStaticSet(args_dict.get("preset"))
    args_dict["settings"] = yaml_loader(config_path)

    return args_dict


########################################


def main(standalone: bool = True, **kwargs):
    """
    A script that writes and performs VASP static run on structure data read from CIF files,
    using pymatgen as a setting interface between raw data and VASP. Energy results may be used
    for phase diagrams computations.
    
    Args:
        standalone (bool):          Whether parsed script is used directly through
                                    command-line (stand-alone script) or in an external
                                    pipeline script.

        input_file (str|Path):      Path to the CIF file containing structure data to read.

        executable_path (str|Path): Path to the VASP executable.

        output (str|Path):          Path to the output directory where VASP files will be written.
                                    A subdirectory will be created in output directory for each
                                    structure found in input_file. By default, it will create a
                                    'Statics' directory in the input_file directory and write
                                    structures runs in it.

        preset (str):               The pymatgen preset to use for VASP static calculations.
                                    List of currently supported presets can be checked in command
                                    line within the --help documentation of this argument for this
                                    script. More info on possible presets in pymatgen documentation:
                                    https://pymatgen.org/pymatgen.io.vasp.html#pymatgen.io.vasp.sets.

        user_settings (str):        Name of the .yaml file containing tags overrides to put over
                                    the PMG preset. The file must be at location
                                    'metrics_pipeline/config' to be found.

        workers (int):              Number of parallel processes to spawn for parallelized steps.
                                    If not given, default value is the 'max_workers' default value
                                    from tqdm.contrib.concurrent.process_map function.

        task_index (int):           If a job array is used, provide here the structure index to
                                    treat according to task IDs (e.g. if task ID 0 treats structure
                                    0 and so on, just provide the task ID). If not given, it will
                                    default to the first possible index, i.e. index 0.
        
        standalone (bool):      Whether parsed script is used directly through
                                command-line (stand-alone script) or in an external
                                pipeline script.

        input_file (str|Path):  Path to the CIF file containing structure data to read.
    """
    start = datetime.now()
    args = _parse_input_args(_get_command_line_args, _process_input_args, standalone, **kwargs)

    # Convert CIF data into Structure objects
    structures, *_ = read_cif(
        filename=args["input_file"],
        workers=args.get("workers"),
        keep_rare_gases=True, # Avoid calling rare gaz screening function
        keep_rare_earths=True # Avoid calling rare earth screening function
    )
    task_idx: int = args["task_index"]
    structure: SiteCollection = structures[task_idx]
    dir_name = f"{args.get('task_index')}_{structure.composition.reduced_formula}"

    vasp_input = vasp_static_settings(
        structure=structure,
        preset=args["preset"],
        user_corrections=args.get("settings")
    )
    run_dir = add_new_dir(args["output"], dir_name)

    write_and_run_vasp(vasp_input, run_dir, args["executable_path"])

    stop = datetime.now()
    print(f"Elapsed time: {stop-start}")


if __name__=="__main__":
    main()
