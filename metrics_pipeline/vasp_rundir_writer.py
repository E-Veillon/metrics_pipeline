#!/usr/bin/python
"""
Write VASP run directories for structures read from CIF file or summary JSON file
using pymatgen as a setup interface between raw data and VASP. Can setup both static
and relaxation runs.
"""

import os
import typing as tp
from datetime import datetime
import argparse as ap

from src.utils import parse_input_args, check_type, check_num_value
from src.io import (
    check_file_or_dir, load_yaml_as_dict, CIFFile,
    VaspWriter, VaspParser, VaspExtractor, CONFIGPATH, EmptyDirectoryError
)
from src.computations.vasp import init_vasp_settings, PMGRelaxSet, PMGStaticSet


def _get_command_line_args() -> ap.Namespace:
    """Command Line Interface (CLI)."""
    parser = ap.ArgumentParser(prog=os.path.basename(__file__), description=__doc__)
    parser.add_argument(
        "input_file",
        help=(
            "Path to the file containing structure data to read. "
            "It can be a CIF file to read structure data directly, or "
            "a JSON summary file from a previous pipeline step to extract "
            "structure data from a previous VASP run."
        )
    )
    parser.add_argument(
        "--output", "-o", metavar="outdir",
        help=(
            "Path to the output directory where VASP files will be written. "
            "A subdirectory will be created in output directory "
            "for each structure found in input_file."
        )
    )
    parser.add_argument(
        "--indices", "-i", nargs="*", type=int, metavar="int",
        help=(
            "If only a few specific structures need to be parsed, pass here their respective index "
            "(0-based position in CIF file order or name index in JSON summary). If not given, all "
            "structures in input_file are parsed normally."
        ),
    )
    parser.add_argument(
        "--summary-key", "-k",
        help=(
            "If a JSON summary file is given as input_file, pass here the name of the key "
            "containing a boolean value telling if said structure passed previous step filter."
        )
    )
    parser.add_argument(
        "--preset", "-p",
        help=(
            "Name of the pymatgen preset to use as a base to write the VASP input files "
            "(case insensitive). More info on possible presets in pymatgen documentation: "
            "https://pymatgen.org/pymatgen.io.vasp.html#pymatgen.io.vasp.sets."
        )
    )
    parser.add_argument(
        "--user-settings", "-u", default="default_settings.yaml",
        help=(
            "Name of the YAML file containing tags to override the PMG preset. "
            f"Given filename must be located in {CONFIGPATH} to be found."
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
    args: ap.Namespace = parser.parse_args()
    return args


def _process_input_args(args_dict: dict[str, tp.Any]) -> dict[str, tp.Any]:
    """Handle input arguments assertions and processing."""
    if args_dict is None:
        raise ValueError(f"No arguments found at '{os.path.basename(__file__)}' script call.")
    check_type(args_dict, "args_dict", (dict,))

    # Set default values for unset optional arguments
    args_dict = {k: v for k, v in args_dict.items() if v is not None}
    default_output = os.path.join(os.path.dirname(args_dict.get("input_file", "")), "Relaxations")
    args_dict.setdefault("output", default_output)
    args_dict.setdefault("user_settings", "default_settings.yaml")

    # Assert set arguments conformity
    check_file_or_dir(args_dict.get("input_file"), "file", allowed_formats=("cif", "json"))
    check_file_or_dir(args_dict.get("output"), "dir")

    if args_dict.get("indices") is not None:
        for idx, struct_idx in enumerate(args_dict["indices"]):
            check_num_value(struct_idx, f"indices[{idx}]", ">=", 0)

    if str(args_dict["input_file"]).endswith(".json") and not args_dict.get("summary_key"):
        raise ValueError("Given input file is a JSON summary, but the summary key was not given.")

    allowed_presets = list(
        map(str.lower, PMGRelaxSet.names() + PMGStaticSet.names() + ["DsolStaticSet"])
    )
    assert str(args_dict.get("preset", "")).lower() in allowed_presets, (
        "Provided preset must be one of the following (case insensitive): "
        f"{', '.join(sorted(allowed_presets))}."
    )
    config_path = os.path.join(CONFIGPATH, args_dict["user_settings"])
    check_file_or_dir(config_path, "file", allowed_formats=("yml", "yaml"))
    
    if args_dict.get("workers") is not None:
        check_type(args_dict.get("workers"), "workers", (int,))
        check_num_value(args_dict.get("workers"), "workers", ">", 0)

    # Additional arguments processing
    os.makedirs(args_dict["output"], exist_ok=True)
    args_dict["settings"] = load_yaml_as_dict(config_path)

    return args_dict


def main(standalone: bool = True, **kwargs):
    """
    Write VASP run directories for structures read from CIF file or summary JSON file
    using pymatgen as a setup interface between raw data and VASP. Can setup both static
    and relaxation runs.

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

    indices: list[int], optional
        If only a few specific structures need to be parsed, pass here their respective index
        (0-based position in CIF file order or name index in JSON summary). If not given, all
        structures in input_file are parsed normally.

    summary_key: str, optional
        If a JSON summary file is given as input_file, pass here the name of the key
        containing a boolean value telling if said structure passed previous step filter.

    preset: str, optional
        Name of the pymatgen preset to use as a base to write the VASP input files
        (case insensitive). More info on possible presets in pymatgen documentation:
        https://pymatgen.org/pymatgen.io.vasp.html#pymatgen.io.vasp.sets.

    user_settings: str, optional
        Name of the YAML file containing tags to override the PMG preset. Given filename
        must be located in 'metrics_pipeline/config' to be found.

    workers: int, optional
        Number of processes to use in parallel. If not given, will use default of
        `tqdm.contrib.concurrent.process_map()`. Pass 0 to disable `process_map()`
        and execute sequentially.
    """
    start = datetime.now()
    args = parse_input_args(_get_command_line_args, _process_input_args, standalone, **kwargs)
    input_file: str = args["input_file"]
    indices: list[int] | None = args.get("indices")

    if input_file.endswith(".cif"):
        # Convert CIF data into Structure objects
        structures, _ = CIFFile.from_file(
            input_file,
            workers=args.get("workers")
        ).parse_structures()

        # Get structures of interest
        if indices is None:
            structs_data = [
                (
                    f"{struct_idx}_{structure.composition.reduced_formula}",
                    structure
                ) for struct_idx, structure in enumerate(structures)
            ]
        else:
            structs_data = [
                (
                    f"{struct_idx}_{structures[struct_idx].composition.reduced_formula}",
                    structures[struct_idx]
                ) for struct_idx in indices
            ]

    elif input_file.endswith(".json"):
        base_dir, summary_name = os.path.split(input_file)
        vparser = VaspParser(
            base_dir=base_dir,
            indices=indices
        )
        if not vparser.struct_dirs:
            if indices is None:
                raise EmptyDirectoryError(
                    f"{base_dir}: No run directory found in this directory."
                )
            else:
                raise ValueError(
                    f"None of given indices ({', '.join(sorted(map(str, indices)))}) correspond "
                    f"to a run directory inside {base_dir}."
                )

        data_dict = VaspExtractor(
            vparser, method="final_state",
            summary_name=summary_name,
            summary_key=args["summary_key"],
            workers=args.get("workers")
        ).get_data()

        if not data_dict: # No Structure passed last step
            print(
                "None of parsed structures passed previous step filter. "
                "Therefore, no VASP run directory is being written."
            )
            exit(0)

        structs_data = [
            (os.path.basename(dir_path), structure_data["structure"])
            for dir_path, structure_data in data_dict.items()
        ]

    vasp_inputs = {}
    for dir_name, structure in structs_data:
        vasp_input = init_vasp_settings(
            structure=structure,
            preset=args["preset"],
            user_corrections=args.get("settings")
        )
        vasp_inputs[dir_name] = vasp_input

    VaspWriter(base_dir=args["output"], vasp_inputs=vasp_inputs)
    print(f"{len(vasp_inputs)} run directories were successfully written in {args['output']}.")

    stop = datetime.now()
    print(f"Elapsed time: {stop-start}")


if __name__=="__main__":
    main()
