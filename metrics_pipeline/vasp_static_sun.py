#!/usr/bin/python
"""
Write VASP static run directories for structure data read from CIF file
or summary JSON file using pymatgen as a setting interface between raw data and VASP.
"""

import os
import typing as typ
from datetime import datetime
import argparse as argp

from pymatgen.core.structure import SiteCollection

from . import CONFIGPATH, _parse_input_args
from utils.utils import check_type, check_num_value
from utils.io import (
    check_file_or_dir, read_cif, load_yaml_as_dict, VaspParser, VaspWriter, JsonLoader
)
from .utils.computations.vasp import vasp_static_settings, PMGStaticSet


def _get_command_line_args() -> argp.Namespace:
    """Command Line Interface (CLI)."""
    parser = argp.ArgumentParser(prog=os.path.basename(__file__), description=__doc__)
    parser.add_argument(
        "input_file", type=str,
        help="Path to the CIF file containing structure data to read.",
    )
    parser.add_argument(
        "-o", "--output",
        help=(
            "Path to the output directory where VASP files will be written. "
            "A subdirectory will be created in output directory for each structure "
            "found in input_file. By default, it will create a 'Statics' directory "
            "in the input_file directory and write structures runs in it."
        )
    )
    parser.add_argument(
        "-p", "--preset", type=PMGStaticSet, default=PMGStaticSet.MPSTATICSET.value,
        help=(
            "The pymatgen preset to use for VASP static run. "
            "More info on possible presets in pymatgen documentation:\n"
            "https://pymatgen.org/pymatgen.io.vasp.html#pymatgen.io.vasp.sets."
        )
    )
    parser.add_argument(
        "-u", "--user-settings", type=str, default="default_settings.yaml",
        help=(
            "Name of the .yaml file containing tags overrides to put over the PMG preset.\n"
            "The file must be at location metrics_pipeline/config to be found."
        ),
        dest="user_settings"
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
        "-t", "--task_index", type=int, default=0,
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
    default_output = os.path.join(os.path.dirname(args_dict.get("input_file", "")), "Statics")
    args_dict.setdefault("output", default_output)
    args_dict.setdefault("preset", PMGStaticSet.MPSTATICSET.value)
    args_dict.setdefault("user_settings", "default_settings.yaml")
    args_dict.setdefault("task_index", 0)

    # Assert set arguments conformity
    check_file_or_dir(args_dict.get("input_file"), "file", allowed_formats="cif")

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
    args_dict["settings"] = load_yaml_as_dict(config_path)

    return args_dict

# TODO: simplify by asking an input CIF or base dir and output base dir more explicitly
def main(standalone: bool = True, **kwargs):
    """
    Write VASP static run directories for structure data read from CIF file
    or summary JSON file using pymatgen as a setting interface between raw data and VASP.
    
    Parameters
    ----------
    standalone: bool
        Whether parsed script is used directly through command-line (stand-alone script)
        or in an external pipeline script.

    input_file: str | Path
        Path to the file containing structure data to read. It can be a CIF file to read
        structure data directly, or a JSON summary file from a previous pipeline step to
        extract structure data from a previous VASP run.

    output: str | Path
        Path to the output directory where VASP files will be written. A subdirectory will
        be created in output directory for each structure found in input_file.

    preset: str, optional
        The pymatgen preset to use for VASP static run. More info on possible presets in
        pymatgen documentation: https://pymatgen.org/pymatgen.io.vasp.html#pymatgen.io.vasp.sets.

    user_settings: str, optional
        Name of the YAML file containing tags to override the PMG preset. Given filename
        must be located in 'metrics_pipeline/config' to be found.

    workers: int, optional
        Number of processes to use in parallel. If not given, will use default of
        `tqdm.contrib.concurrent.process_map()`. Pass 0 to disable `process_map()`
        and execute sequentially.

    task_index: int, optional
        If a job array is used, provide here the structure index to treat according to task IDs
        (e.g. if task ID 0 treats structure 0 and so on, just provide the task ID). If not given,
        it will default to the first possible index, i.e. index 0.
    """
    start = datetime.now()
    args = _parse_input_args(_get_command_line_args, _process_input_args, standalone, **kwargs)

    if str(args.get("input_file")).endswith(".cif"):
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

    elif str(args["input_file"]).endswith(".json"):
        data = JsonLoader(args["input_file"]).load_as_list()
        paths = [d["path"] for d in data]
        task_idx: int = args["task_index"]

        try:
            struct_dir = next(filter(
                lambda path: os.path.basename(path).startswith(f"{task_idx}_"),
                paths
            ))
        except StopIteration as exc:
            raise ValueError(
                f"Given 'task_index' value ({task_idx}) do not match any data "
                f"in file {args.get('input_file')}."
            ) from exc

        structure = VaspParser(
            base_dir=os.path.dirname(args["input_file"]),
            indices=[task_idx]
        ).parse_structures()[struct_dir]

        dir_name = os.path.basename(struct_dir)

    vasp_input = vasp_static_settings(
        structure=structure,
        preset=args["preset"],
        user_corrections=args.get("settings")
    )
    VaspWriter(args["output"], {dir_name: vasp_input})

    stop = datetime.now()
    print(f"Elapsed time: {stop-start}")


if __name__=="__main__":
    main()
