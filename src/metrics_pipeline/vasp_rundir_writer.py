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

from pymatgen.io.vasp import VaspInput

from metrics_pipeline.core.utils import parse_input_args, check_type, check_num_value
from metrics_pipeline.core.genmat_io import (
    check_file_or_dir, load_yaml_as_dict,
    VaspWriter, VaspParser, VaspExtractor, ExtractMethod, CONFIGPATH, GenMatFile
)
from metrics_pipeline.core.computations.local import DSolCalcType
from metrics_pipeline.core.computations.vasp import (
    init_vasp_settings, dsol_calc_init, ALL_PRESETS_NAMES, ALL_PRESETS_NAMES_LOWER,
    STATIC_PRESETS_NAMES_NO_DSOL_LOWER, GenMatSet
)


def _get_command_line_args() -> ap.Namespace:
    """Command Line Interface (CLI)."""
    parser = ap.ArgumentParser(prog=os.path.basename(__file__), description=__doc__)
    parser.add_argument(
        "input_file",
        help=(
            "Path to the JSON file containing structure data to read. "
            "It can be either a preprocessed structure data file to read directly, "
            "or a summary file from a previous pipeline step to extract "
            "structure data from a previous VASP run. It will be assumed to be a "
            "summary file only if the --summary-key argument is passed."
        )
    )
    parser.add_argument(
        "-o", "--output",
        help=(
            "Path to the output directory where VASP files will be written. "
            "A subdirectory will be created in output directory "
            "for each structure found in input_file."
        )
    )
    parser.add_argument(
        "--output-names",
        help="Path to a file to create and store all created run directories paths."
    )
    parser.add_argument(
        "-i", "--indices", nargs="*", type=int, metavar="int",
        help=(
            "If only a few specific structures need to be parsed, pass here their respective "
            "preprocessing indices. If not given, all structures in input_file are parsed."
        ),
    )
    parser.add_argument(
        "-k", "--summary-key",
        help=(
            "If a JSON summary file is given as input_file, pass here the name of the key "
            "containing a boolean value telling if said structure passed previous step filter."
        )
    )
    parser.add_argument(
        "-p", "--preset",
        help=(
            "Name of the pymatgen preset to use as a base to write the VASP input files "
            "(case insensitive). More info on possible presets in pymatgen documentation: "
            "https://pymatgen.org/pymatgen.io.vasp.html#pymatgen.io.vasp.sets."
        )
    )
    parser.add_argument(
        "-u", "--user-settings",
        help=(
            "Name of the optional YAML file containing tags to override in the preset. "
            f"Given filename must be located in {CONFIGPATH} to be found."
        )
    )
    parser.add_argument(
        "--delta-sol", action="store_true",
        help=(
            "Pass this flag to enable delta-sol method sub-runs generation. Inside each written "
            "structure directory, 3 VASP run subdirectories will be generated to compute energy "
            "for neutral and ionized versions of the structure in order to compute the band gap. "
            "Note that only structures that are relaxed may give accurate band gap values."
        )
    )
    parser.add_argument(
        "--dsol-uncertainty", action="store_true",
        help=(
            "If 'delta-sol' flag is passed, pass this flag to enable computations of delta-sol "
            "uncertainties on the band gap value, generating 7 sub-runs instead of 3. "
            "Ignored if 'delta-sol' flag is not passed. WARNING: uncertainty computations were "
            "not rigorously tested yet and might give weird results in some cases. To interpret "
            "with caution."
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
        "--match-all", action=ap.BooleanOptionalAction, default=True,
        help=(
            "Whether to raise an error if some run indices do not match any structure directory "
            "in base directory. Defaults to %(default)s. Deactivate it if some structures "
            "could not be written as VASP input when using this script before."
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
    args_dict.setdefault("delta_sol", False)
    args_dict.setdefault("dsol_uncertainty", False)
    args_dict.setdefault("match_all", True)

    # Assert set arguments conformity
    check_file_or_dir(args_dict.get("input_file"), "file", allowed_formats="json")

    if args_dict.get("indices") is not None:
        for idx, struct_idx in enumerate(args_dict["indices"]):
            check_num_value(struct_idx, f"indices[{idx}]", ">=", 0)

    assert str(args_dict.get("preset", "")).lower() in ALL_PRESETS_NAMES_LOWER, (
        "Provided preset must be one of the following (case insensitive): "
        f"{', '.join(ALL_PRESETS_NAMES)}, got {str(args_dict.get('preset', '')).lower()!r}."
    )
    if args_dict["preset"].lower() == GenMatSet.DSOLSTATICSET.name.lower():
        default_outdir = "Bandgaps"
    elif args_dict["preset"].lower() in STATIC_PRESETS_NAMES_NO_DSOL_LOWER:
        default_outdir = "Statics"
    else:
        default_outdir = "Relaxations"
    default_output = os.path.join(os.path.dirname(args_dict.get("input_file", "")), default_outdir)
    args_dict.setdefault("output", default_output)

    if args_dict.get("user_settings"):
        config_path = os.path.join(CONFIGPATH, args_dict["user_settings"])
        check_file_or_dir(config_path, "file", allowed_formats=("yml", "yaml"))
        args_dict["user_settings"] = load_yaml_as_dict(config_path)

    if args_dict.get("workers") is not None:
        check_type(args_dict["workers"], "workers", (int,))
        check_num_value(args_dict["workers"], "workers", ">", 0)

    # Additional arguments processing
    os.makedirs(args_dict["output"], exist_ok=True)
    if args_dict.get("output_names") is not None:
        os.makedirs(os.path.dirname(args_dict["output_names"]), exist_ok=True)

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
        Path to the JSON file containing structure data to read.
        It can be either a preprocessed structure data file to read directly, or a summary
        file from a previous pipeline step to extract structure data from a previous VASP run.
        It will be assumed to be a summary file only if summary_key is passed.

    output: str | Path
        Path to the output directory where VASP files will be written. A subdirectory will
        be created in output directory for each structure found in input_file.

    output_names: str | Path
        Path to a file to create and store all created run directories paths.

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

    delta_sol: bool
        Whether to enable delta-sol method sub-runs generation. Inside each written structure
        directory, 3 VASP run subdirectories will be generated to compute energy for neutral
        and ionized versions of the structure in order to compute the band gap. Note that only
        structures that are relaxed may give accurate band gap values.

    dsol_uncertainty: bool
        Only used if `delta_sol` is set to True. Enable computations of delta-sol uncertainties
        on the band gap value, generating 7 sub-runs instead of 3. WARNING: uncertainty
        computations were not rigorously tested yet and might give weird results in some cases.
        To interpret with caution.

    workers: int, optional
        Number of processes to use in parallel. If not given, will use default of
        `tqdm.contrib.concurrent.process_map()`. Pass 0 to disable `process_map()`
        and execute sequentially.

    match_all: bool
        Whether to raise an error if some run indices do not match any structure directory
        in base directory. Defaults to %(default)s. Deactivate it if some structures
        could not be written as VASP input when using 'vasp_rundir_writer.py' script.
    """
    start = datetime.now()
    args = parse_input_args(_get_command_line_args, _process_input_args, standalone, **kwargs)
    input_file: str = args["input_file"]
    indices: list[int] | None = args.get("indices")

    if args.get("summary_key") is None:
        # Load GenMatStructure data file
        gfile = GenMatFile.from_file(input_file)

        # Get structures of interest
        if indices is None:
            structs_data = [(structure.name, structure) for structure in gfile]
        else:
            structs_data = [
                (structure.name, structure) for structure in gfile
                if structure.name_index in indices
            ]
    else:
        base_dir, summary_name = os.path.split(input_file)
        vparser = VaspParser(
            base_dir=base_dir,
            indices=indices,
            match_all=args["match_all"]
        )
        data_dict = VaspExtractor(
            vparser, method=ExtractMethod.FINAL_STATE,
            summary_name=summary_name,
            summary_key=args["summary_key"],
            workers=args.get("workers")
        ).get_data()

        if not data_dict: # No Structure passed last step
            raise ValueError(
                "None of parsed structures passed previous step filter. "
                "Therefore, no VASP run directory is being written."
            )

        structs_data = [
            (os.path.basename(dir_path), structure_data["structure"])
            for dir_path, structure_data in data_dict.items()
        ]

    vasp_inputs: dict[str, VaspInput] = {}
    if args["delta_sol"]:
        for dir_name, structure in structs_data:
            for calc_idx in range(7 if args["dsol_uncertainty"] else 3):
                vasp_input = dsol_calc_init(
                    structure=structure,
                    calc_index=calc_idx,
                    preset=args["preset"],
                    user_corrections=args.get("user_settings")
                )
                run_name = "_".join((dir_name, DSolCalcType(calc_idx).name.lower()))
                vasp_inputs[os.path.join(dir_name, run_name)] = vasp_input

    else:
        for dir_name, structure in structs_data:
            vasp_input = init_vasp_settings(
                structure=structure,
                preset=args["preset"],
                user_corrections=args.get("user_settings")
            )
            vasp_inputs[dir_name] = vasp_input

    VaspWriter(base_dir=args["output"], vasp_inputs=vasp_inputs).write_run_dirs()

    if args.get("output_names"):
        dir_list = [os.path.join(args["output"], run_dir) for run_dir in vasp_inputs.keys()]
        with open(args["output_names"], "wt", encoding="utf-8") as fp:
            fp.write("\n".join(dir_list))

    print(f"{len(vasp_inputs)} run directories were successfully written in {args['output']}.")

    stop = datetime.now()
    print(f"Elapsed time: {stop-start}")


if __name__=="__main__":
    main()
