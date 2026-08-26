"""
Compute structures energies using ALIGNN to build a phase diagram reference dataset without using DFT.
"""

import argparse as ap
import typing as tp
from pathlib import Path

from src.utils import parse_input_args, check_type, check_num_value
from src.genmat_io import check_file_or_dir, check_file_format, CIFFile, JsonWriter
from src.computations.models import vectors_from_alignn


def _get_cmd_line_args() -> ap.Namespace:
    """Command Line Interface (CLI)."""
    parser = ap.ArgumentParser(prog=Path(__file__).name, description=__doc__)
    parser.add_argument(
        "input_file",
        help="CIF file containing structures to compute ALIGNN energies from."
    )
    parser.add_argument(
        "-o", "--output",
        help=(
            "Path to the JSON file to write dataset in. "
            "Defaults to same path as 'input_file' but with a '.json' extension."
        )
    )
    parser.add_argument(
        "-s", "--save-header", action="store_true",
        help=(
            "Pass this flag to save CIF headers as structures 'entry_id'. "
            "If not passed, 'entry_id' will default to '<index>_<reduced_formula>'."
        )
    )
    parser.add_argument(
        "-e", "--on-error", choices=["raise", "warn", "ignore"], default="warn",
        help=(
            "What to do in case an error occurs during structure parsing. "
            "Can be either 'raise', 'warn' or 'ignore'. Defaults to '%(default)s'. "
            "If not 'raise', unparsed structures will be excluded from ALIGNN computations "
            "and final dataset."
        )
    )
    parser.add_argument(
        "-c", "--compact", action="store_true",
        help=(
            "Flag to pass if JSON dataset should be organized by data type instead of by structure. "
            "By default, dataset is assumed to be organized by structure, i.e. "
            "{'struct_name_1': struct_dict_1, 'struct_name_2': struct_dict_2, ...}. "
            "If the flag is passed, dataset is instead assumed to be organized by data key, i.e. "
            "{'entry_id': [struct_name_1, struct_name_2, ...], "
            "'composition': [struct_comp_1, struct_comp_2, ...], ...}."
        )
    )
    parser.add_argument(
        "-l", "--load-bar-style", default="tqdm", choices=["tqdm", "local", "quiet"],
        help=(
            "Choose a style for the ALIGNN process bar, between 'tqdm', 'local' and 'none'. "
            "'tqdm' uses the default tqdm.tqdm() bar, 'local' uses GenMat VisualIterator(), "
            "and 'quiet' deactivates progress bar completely. Defaults to %(default)s."
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
    args = parser.parse_args()
    return args


def _process_input_args(args_dict: dict[str, tp.Any]) -> dict[str, tp.Any]:
    """Handle input arguments assertions and processing."""
    if args_dict is None:
        raise ValueError(f"No arguments found at '{os.path.basename(__file__)}' script call.")
    
    check_type(args_dict, "args_dict", (dict,))

    # Set default values for unset optional arguments
    args_dict = {k: v for k, v in args_dict.items() if v is not None}
    default_output = args_dict["input_file"].replace(".cif", ".json")
    args_dict.setdefault("output", default_output)
    args_dict.setdefault("on_error", "warn")
    args_dict.setdefault("load_bar_style", "tqdm")
    
    # Assert set arguments conformity
    check_file_or_dir(args_dict["input_file"], allowed_formats="cif")
    check_file_format(args_dict["output"], allowed_formats="json")

    if args_dict.get("workers") is not None:
        check_type(args_dict["workers"], "workers", (int,))
        check_num_value(args_dict["workers"], "workers", ">=", 0)
    
    check_type(args_dict["on_error"], "on_error", (str,))
    assert args_dict["on_error"].lower() in {"raise", "warn", "ignore"}, ValueError(
        f"'on_error' must be either 'raise', 'warn' or 'ignore', got {args_dict['on_error']!r}."
    )
    assert args_dict["load_bar_style"] in {"tqdm", "local", "quiet"}, ValueError(
        "'load_bar_style' must be either 'tqdm', 'local' or 'quiet', "
        f"got {args_dict['load_bar_style']}"
    )

    return args_dict


def main(standalone: bool = True, **kwargs) -> None:
    """
    Compute structures energies using ALIGNN to build a phase diagram reference dataset
    without using DFT.

    Parameters
    ----------
    standalone (bool
        Whether parsed script is used directly through command-line (stand-alone script)
        or in an external pipeline script.

    input_file: str | Path
        CIF file containing structures to compute ALIGNN energies from.

    output: str | Path, optional
        Path to the JSON file to write dataset in. Defaults to same path as 'input_file'
        but with a '.json' extension.

    save_header: bool
        Whether to save CIF headers as structures 'entry_id'. If set to `False`,
        'entry_id' will be '<index>_<reduced_formula>'. Defaults to `False`.

    on_error: 'raise', 'warn', 'ignore'
        What to do in case an error occurs during structure parsing.
        Can be either 'raise', 'warn' or 'ignore'. Defaults to 'warn'.
        If not 'raise', unparsed structures will be excluded from ALIGNN computations
        and final dataset.

    compact: bool
        Whather reference dataset is organized by data type instead of by structure.
        By default, dataset is assumed to be organized by structure, i.e.
        {'struct_name_1': struct_dict_1, 'struct_name_2': struct_dict_2, ...}.
        If the flag is passed, dataset is instead assumed to be organized by data key, i.e.
        {
            'entry_id': [struct_name_1, struct_name_2, ...],
            'composition': [struct_comp_1, struct_comp_2, ...],
            ...
        }. Defaults to False.

    load_bar_style: str
        Choose a style for the ALIGNN process bar, between 'tqdm', 'local' and 'none'.
        'tqdm' uses the default tqdm.tqdm() bar, 'local' uses GenMat VisualIterator(),
        and 'none' deactivates progress bar completely. Defaults to 'tqdm'.

    workers: int, optional
        Number of parallel processes to spawn for parallelized steps. If not given, default value
        is the 'max_workers' default value from `tqdm.contrib.concurrent.process_map()` function.
        Pass 0 to disable the use of `process_map()` and execute sequentially.
    """
    args = parse_input_args(_get_cmd_line_args, _process_input_args, standalone, **kwargs)

    cfile = CIFFile.from_file(
        args["input_file"],
        special_keys=(['header'] if args["save_header"] else None),
        workers=args["workers"]
    )
    structures, _ = cfile.parse_structures(on_error=args["on_error"])
    energies = vectors_from_alignn(structures, output="energy", load_bar=args["load_bar_style"])

    if args["compact"]:
        results = {}
        results["entry_id"] = [
            structure.properties.get("header", f"{idx}_{structure.reduced_formula}")
            for idx, structure in enumerate(structures)
        ]
        results["composition"] = [structure.composition.as_dict() for structure in structures]
        results["final_energy"] = energies.tolist()

    else:
        results = {}
        for idx, (structure, energy) in enumerate(zip(structures, energies)):
            name = structure.properties.get("header", f"{idx}_{structure.reduced_formula}")
            final_energy = float(energy.item()) if hasattr(energy, "item") else float(energy)
            results[name] = {
                "entry_id": name,
                "composition": structure.composition.as_dict(),
                "final_energy": final_energy
            }

    JsonWriter(args["output"], results).write_as_dict()


if __name__ == "__main__":
    main()
